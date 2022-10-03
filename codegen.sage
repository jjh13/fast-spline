#!/usr/bin/env python

"""
Main class that implements part I of the paper. Basically all of the analysis 
happens here, anything that needs a one-time computation is intended to be in 
part I.

Copyright © 2016-2022 Joshua Horacsek

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

__author__ = "Joshua Horacsek"
from llvmlite import ir
from types import MethodType
import os

"""
We need to monkey patch this into the VectorType package --- llvmlite 0.38.0 has bad support for vector operations in 
general
"""
def vector_gep(self, i):
    """
    Resolve the type of the i-th element (for getelementptr lookups).
    """
    if not isinstance(i.type, ir.VectorType) or not isinstance(i.type.element, ir.IntType):
        raise TypeError(i.type)
    self.pointee = self.element.pointee
    return self.element.pointee

def pointer_gep(self, i):
    """
    Resolve the type of the i-th element (for getelementptr lookups).
    """
    if isinstance(i.type, ir.VectorType):
        if not isinstance(i.type.element, ir.IntType):
            raise TypeError("vecdex")
        # self.pointee = self.element

    elif isinstance(i.type, ir.ArrayType):
        self.pointee = self.element

    elif not isinstance(i.type, ir.IntType):
        raise TypeError(i.type)

    if isinstance(self, ir.ArrayType):
        self.pointee = self.element

    return self.pointee

def irbuilder_vector_gep(self, *args, **kwargs):
    # print(args[:])
    # print(args[1].gep(args[2]))
    args[0].gep = vector_gep
    # print(args[1:])
    ret = self.old_gep(*args, **kwargs)
    if isinstance(args[0].type, ir.VectorType) and not isinstance(ret.type, ir.VectorType):
        count = args[0].type.count
        ret.type = ir.VectorType(ret.type, count)
    return ret


def bstruct_vector_cgep(self, i):
    raise ValueError('!!!')

ir.PointerType.gep = pointer_gep
ir.VectorType.gep = vector_gep
ir.ArrayType.gep = pointer_gep
ir.IRBuilder.old_gep = ir.IRBuilder.gep
ir.IRBuilder.gep = irbuilder_vector_gep
# ir.BaseStructType.gep = bstruct_vector_cgep

class EvalCodeGen:
    def __init__(self, ss, **kwargs): 
        options = kwargs
        target_triple = "unknown-unknown-unknown"
        
        if "target_triple" in options:
            target_triple = options["target_triple"]

        self._verbose = True
        if "verbose" in options:
            self._verbose["verbose"] = options['verbose']

        # Break the target triple up
        tt = target_triple.split('-')
        if len(tt) < 3:
            raise ValueError("Invalid target_triple '%s'!")

        self._arch, self._vendor, self._sys = tt[:3]

        self._target = target_triple
        self._known_target = False
        self.coset_decomposition = True

        # By default, setup the parameters for CPU based code
        if 'arm' in self._arch or 'aarch' in self._arch or 'x86' in self._arch or 'avr' in self._arch:
            self._known_target = True

        # Will this spline use texture fetches, if they're available?
        # texture fetches of at most dimension 3  are only supported my most GPUs
        self.has_texture_fetch = False

        # Will this spline use linear fetches on cartesian cosets?
        # by default this is on if we want to use texture fetches
        self.use_linear_fetch = False

        # Which type of memory fetches to use
        self.lookup_type = "ABSTRACT" # ABSTRACT, COSET, TEX, TEX_LIN

        # Whether or not to generate code with branching behaviour
        # this is only beneficial on the GPU
        self.branchless = False

        # Whether or not to use branch predication (if doing branchless)
        self.branch_predication = False

        # Break Horner evaluation up over different SIMD lanes,
        # these don't exist on GPUs, so we disable by default
        self.vectorize = True

        self.lane_count = 4       # Most CPUs have a vector extension with at least 4 elements

        # Should we refetch the transform on each lattice site read?
        # On CPUs - no, we should just store the transform on the stack, 
        # but on GPU we should always refetch --- this allows register
        # usage to be economized
        self.refetch_transform = False

        # How many memfetchest to queue before computing the polynomial
        self.group_size = floor(ss._ref_subregions[0].number_of_pts_per_reconstruction() / 4)
        self.group_size = ss._ref_subregions[0].number_of_pts_per_reconstruction() if self.group_size == 0 else self.group_size
        self.pipeline_depth = 0

        if self._arch in ['nvptx', 'nvptx64', 'amdil', 'amdil64', 'hsail', 'hsail64', 'spir', 'spir64', 'r600','amdgcn']:
            self._known_target = True

            self.has_texture_fetch = self._arch in ['nvptx', 'nvptx64']
            self.lookup_type = "TEX_LIN" if self.coset_decomposition and self.has_texture_fetch else "ABSTRACT"

            self.branchless = True
            self.branch_predication = True
            self.vectorize = False
            self.lane_count = 1
            self.refetch_transform = True

        if not self._known_target:
            self.log("Unknown architecture '%s' continuing with common PC parameters" % self._arch)

        if self.lookup_type == "TEX_LIN":
            self.group_size = 1

        self.reconstruction_name = "__reconstruct__"
        self.strict_branchless = False

        """
        Patch in options from the input options 
        """
        self.cache_path = None
        for k in options:
            if k == 'lane_count':
                self.lane_count = options[k]
            elif k == 'use_texture_fetch':
                self.use_texture_fetch = options[k]
            elif k == 'pipeline_depth':
                try:
                    self.pipeline_depth = int(options[k])
                except:
                    raise ValueError('pipeline_depth option must be convertable to an int!')
            elif k == 'use_linear_fetch':
                self.use_linear_fetch = options[k]
            elif k == 'branchless':
                self.branchless = options[k]
            elif k == 'strict_branchless':
                self.strict_branchless = options[k]
            elif k == 'branch_predication':
                self.branch_predication = options[k]
            elif k == 'reconstruction_name':
                if options[k] is not None:
                    self.reconstruction_name = options[k]
            elif k == 'refetch_transform':
                self.refetch_transform = options[k]
            elif k == 'group_size':
                try:
                    self.group_size = int(options[k])
                except:
                    raise ValueError('group_size option must be convertable to an int!')
            elif k == 'lookup_type':
                if options[k] not in ['ABSTRACT', 'COSET', 'COSET_RAW', 'TEX', 'TEX_LIN', 'CPU_LIN']:
                    raise ValueError("lookup_type must be one of ABSTRACT, COSET, TEX, TEX_LIN")
                self.lookup_type = options[k]
            elif k == 'cache_path':
                self.cache_path = options[k]
            elif k =='target_triple':
                pass
            else:
                raise ValueError('Unknown option "%s"' % k)

        """
        Double check that the parameters we have are good
        """ 
        if not self.branchless and self.branch_predication:
            raise ValueError("Cannot generate branch predicated code when using branching")

        if self.lookup_type in ['COSET', 'COSET_RAW', 'TEX', 'TEX_LIN', 'CPU_LIN'] and not self.coset_decomposition:
            raise ValueError("Cannot use a COSET, TEX or TEX_LIN if we don't know the coset decomposition for the "
                             "spline space!")

        if self.lookup_type in ['TEX', 'TEX_LIN'] and not self.has_texture_fetch:
            raise ValueError("Cannot generate code with texture fetches on architectures without texture units!")

        if (self.lookup_type == "TEX_LIN" or self.lookup_type == "CPU_LIN") and self.group_size != 1:
            raise ValueError("Cannot pipeline linear texture fetches")

        if not self.branch_predication and self.lane_count > 1:
            self.log("Warning: We can't jump when the lane count != 1, turning on branchless and branch_predication...")
            self.branchless = True
            self.branch_predication = True

        # Default variable types
        self.internal_precision = 'single'
        self.input_precision = 'single'
        self.output_precision = 'single'
        self.memfetch_precision = 'single'
        self.module_name = 'fastspline_codegen'

        self.log(f"Generating code for '{self._arch}':")
        self.log(f"   branchless: {self.branchless}")
        self.log(f"   branch_predication: {self.branch_predication}")
        self.log(f"   group_size: {self.group_size}")
        self.log(f"   pipeline_depth: {self.pipeline_depth}")
        self.log(f"   lane_count: {self.lane_count}")
        self.log(f"   refetch_transform: {self.refetch_transform}")
        self.log(f"   lookup_type: {self.lookup_type}")

        self._ss = ss
        self._dimension = ss._s
        self._process_sub_regions()

        self._bake_raw_fetches = True
        self._bake_ls_lookups = True
        self.ptr_size = int(64)
        self.test_stop = "" # "first_lindex" # "boundary0" #"poly_index"


    def _process_sub_regions(self):
        self._subregion_consistency = True

        # For each sub region, bin the lattice sites into their respective
        # coset structure, and do a greedy horner factorization on each 
        # subregion
        self.log("Processing reference subregions...")
        G, coset_vectors = self._ss.coset_vectors()
        for region_idx, region in enumerate(self._ss._ref_subregions):
            region_path = None
            is_linear_fetch = self.lookup_type == "TEX_LIN" or self.lookup_type == "CPU_LIN"
            if self.cache_path is not None:
                rt = "_tex" if is_linear_fetch else ""
                region_path = os.path.join(self.cache_path, 'planned_regions', f"region_{region_idx}{rt}")

                if region.load(region_path):
                    self.log("😊", ' ')
                    continue

            self.tick("Partitioning reference subregion %d into cosets... " % region_idx)
            region.partition_into_cosets(G, coset_vectors)

            self.tick("factoring...")
            region.factor_for_evaluation()

            # Use dancing links to cover each lookup, if necessary 
            self.tick("covering...")
            region.cover_lookups(is_linear_fetch)

            print(region._coverings)

            # Plan out the lookups --- if we use a texture fetch, 
            # then we can take advantage of texture locality
            self.tick("optimizing lookups...")
            region.order_lookups()

            if self.cache_path is not None:
                region.save(region_path)
            self.log("x" if self.cache_path is None else "", ' ')
            self.log("done!")

        # Check for consistency between sub-regions, that is, check to see 
        # that each sub region has the same number of accesses per coset,
        # this means we can avoid determining which coset a memory fetch is 
        # on every time we access a lattice site
        subregion_lk_sizes = [len(_) for _ in self._ss._ref_subregions[0]._coset_partition]
        for refsregion in self._ss._ref_subregions[1:]:
            _sizes = [len(_) for _ in refsregion._coset_partition]
            if any(_ != __ for (_, __) in zip(subregion_lk_sizes, _sizes)):
                self.log("Subregions are inconsistent -- turning off consistency optimization")
                self._subregion_consistency = False

        if self._subregion_consistency:
            self.log("Consistent subregions!")

    def _llvm_generate_vec_module(self):
        """
        This starts the LLVM code generation phase.

        :return: None
        """
        # Determine the width of calulation within the module
        calculation_type = ir.FloatType() if self.internal_precision == 'single' else ir.DoubleType()
        self._internal_type = ir.VectorType(calculation_type, self.lane_count) if self.lane_count > 1 else calculation_type

        # Setup the input and output types
        calculation_type = ir.FloatType() if self.input_precision == 'single' else ir.DoubleType()
        self._input_type = ir.VectorType(calculation_type, self.lane_count) if self.lane_count > 1 else calculation_type

        calculation_type = ir.FloatType() if self.output_precision == 'single' else ir.DoubleType()
        self._output_type = ir.VectorType(calculation_type, self.lane_count) if self.lane_count > 1 else calculation_type

        # Setup the type for memory fetches, note that this isn't vectorized
        self._memfetch_base = ir.FloatType() if self.memfetch_precision == 'single' else ir.DoubleType()
        self._memfetch_type = ir.VectorType(self._memfetch_base, self.lane_count) if self.lane_count > 1 else self._memfetch_base

        self._index_base = ir.IntType(int(32))
        self._index_type = ir.VectorType(self._index_base, self.lane_count) if self.lane_count > 1 else self._index_base

        # Setup the module
        self._llvm_module = ir.Module(name=self.module_name)
        self._llvm_module.triple = self._target

        # Okay, there's a bug with llvm ite where VectorTypes don't have an 'intrinsic_name' property, so we patch that
        # property to properly reflect the name of the vector type
        def setup_intrinsic_name(llvm_type):
            if hasattr(llvm_type, "intrinsic_name"):
                return llvm_type
            llvm_type.intrinsic_name = f"v{llvm_type.count}{llvm_type.element.intrinsic_name}"
            return llvm_type

        self._internal_type = setup_intrinsic_name(self._internal_type)
        self._input_type = setup_intrinsic_name(self._input_type)
        self._output_type = setup_intrinsic_name(self._output_type)

        # Declare intrinsics for a few important functions
        rfn_type = ir.FunctionType(self._internal_type, (self._internal_type,))
        self._round = self._llvm_module.declare_intrinsic("llvm.round", [self._internal_type], fnty=rfn_type)
        ffn_type = ir.FunctionType(self._internal_type, (self._internal_type,))
        self._floor = self._llvm_module.declare_intrinsic("llvm.floor", [self._internal_type], fnty=ffn_type)
        fafn_type = ir.FunctionType(self._internal_type, (self._internal_type,))
        self._fabs = self._llvm_module.declare_intrinsic("llvm.fabs", [self._internal_type], fnty=fafn_type)

        # For GPU code, we use texture lookups, we declare those intrinsics here
        self.constant_addrspace = 0
        if self.lookup_type == "TEX" or self.lookup_type == "TEX_LIN":
            assert self.lane_count == 1
            float4 = ir.LiteralStructType([ir.FloatType(), ir.FloatType(), ir.FloatType(), ir.FloatType()])
            fnty = ir.FunctionType(float4, (ir.IntType(int(64)), ir.FloatType(), ir.FloatType(), ir.FloatType()))
            self._tex3d = self._llvm_module.declare_intrinsic("llvm.nvvm.tex.unified.3d.v4f32", [ir.FloatType()], fnty)

            float4 = ir.LiteralStructType([ir.FloatType(), ir.FloatType(), ir.FloatType(), ir.FloatType()])
            fnty = ir.FunctionType(float4, (ir.IntType(int(64)), ir.FloatType(), ir.FloatType()))
            self._tex2d = self._llvm_module.declare_intrinsic("llvm.nvvm.tex.unified.2d.v4f32", [ir.FloatType()], fnty)

            float4 = ir.LiteralStructType([ir.FloatType(), ir.FloatType(), ir.FloatType(), ir.FloatType()])
            fnty = ir.FunctionType(float4, (ir.IntType(int(64)), ir.FloatType()))
            self._tex1d = self._llvm_module.declare_intrinsic("llvm.nvvm.tex.unified.1d.v4f32", [ir.FloatType()], fnty)

            # On the GPU constant address space refers to the actual 'constant' memory on chip, we need to explicitly
            # indicate this, here
            self.constant_addrspace = 4

        # Abstract Memory Fetch Type
        input_types = ([ir.IntType(int(32))] * self._dimension)  # The input is dim_s integer vectors
        self._abstract_memfetch_lk_type = ir.FunctionType(self._memfetch_base, input_types)

        # Coset Lookup Fetch type
        coset_levels = 1 if len(self._ss.get_L_cosets()) == 1 else 2
        input_types = ([ir.IntType(32)]*(self._dimension + coset_levels))
        self._cluptp = ir.FunctionType(self._memfetch_type, input_types)

        # Raw Coset Lookup type
        self._coset_info_struct = ir.LiteralStructType([
            ir.IntType(32).as_pointer(),
            self._memfetch_base.as_pointer()
        ])

        self._raw_lookup = ir.LiteralStructType([
            ir.IntType(32),
            ir.IntType(32),
            self._coset_info_struct.as_pointer()
        ])

        # linear interpolated lookups (test jig)
        input_types = ([ir.IntType(32)]+[self._internal_type] * self._dimension)
        self._clerptp = ir.FunctionType(self._memfetch_type, input_types)

        return self._llvm_module

    def llvm_generate_bake_function(self, module):
        print("Generating bake function")

        subregions = self._ss.get_subregions()
        num_cosets = self._ss.coset_count()

        for region in subregions:
            transform = region.get_transform()
            region = subregions[region.get_ref_region()]

            for coset in range(num_cosets):
                for ls in region.get_lattice_sites():
                    pass
                    # ls = ls + coset_shift(coset)
                    # ls = cartesianize(ls)
                    # offset = linearize(ls, boundaries)
                    # print(ls)

        #
        # # Build types
        # input_types = ([self._input_type] * self._dimension)
        # ntextures = len(self._ss.coset_vectors()[1]) * len(self._ss.get_L_cosets())
        #
        # # Create a new type for our function
        # fn_type = ir.FunctionType(self._output_type, tuple(input_types + mf))
        # self._function = ir.Function(self._llvm_module, fn_type, name=self.reconstruction_name)
        #
        # memfetch = self._function.args[self._dimension:]
        # xvars = self._function.args[0:self._dimension]
        # self._textures = self._function.args[self._dimension:]
        #
        # self._entry_blk = self._function.append_basic_block(name="entry")
        # builder = ir.IRBuilder(self._entry_blk)
        #
        # if self.lookup_type == 'COSET_RAW':
        #     # Generate reads for all the coset offset stuff
        #     print()
        #
        # return builder, xvars, memfetch
    def llvm_build_stack(self, builder):
        num_cosets = self._ss.coset_count()

        self._stack_xvars = [
            builder.alloca(self._internal_type, name='xvar_%d' % _) for _ in range(self._dimension)
        ]
        self._stack_txvars = [
            builder.alloca(self._internal_type, name='txvar_%d' % _) for _ in range(self._dimension)
        ]

        self._stack_rho = [
            builder.alloca(self._index_type, name='rho_%d' % _) for _ in range(self._dimension)
        ]

        self._stack_bounds = [[
            builder.alloca(self._index_type, name=f"bdry{coset_idx}_{_}") for _ in range(self._dimension)
        ] for coset_idx in range(num_cosets)]

        buffer_primitive = self._memfetch_base.as_pointer()

        self._stack_buffers = [
            builder.alloca(
                ir.VectorType(buffer_primitive, self.lane_count) if self.lane_count > 1 else buffer_primitive,
                name='coset_buffer_%d' % coset)
            for coset in range(num_cosets)
        ]

        if len(self._ss.get_subregions()) > 1:
            self._stack_xform = [
                [builder.alloca(self._internal_type, name=f"xmat_{i}_{j}") for i in range(self._dimension)]
                for j in range(self._dimension)
            ]

            self._stack_tshift = [
                builder.alloca(self._internal_type, name='shift_%d' % _) for _ in range(self._dimension)
            ]

        self._stack_pipeline = [
            builder.alloca(self._internal_type, name='p_%d' % _) for _ in range(self.pipeline_depth)
        ]
        self._stack_result = builder.alloca(self._internal_type, name='result')
        self._stack_coset_result = builder.alloca(self._internal_type, name='coset_result')
        self._stack_bspindex = builder.alloca(self._index_type, name='bspidx')
        if len(self._ss._ref_subregions) > 1:
            self._stack_polyindex = builder.alloca(self._index_type, name='polyidx')
        self._stack_cosetid = builder.alloca(self._index_type, name='coset_id')

        if self._ss._reflective:
            self._stack_xflip = [
                builder.alloca(self._internal_type, name='xflip_%d' % _) for _ in range(self._dimension)
            ]

        self._stack_gval = builder.alloca(self._internal_type, name='g_value')
        self._stack_g_pt = [
            builder.alloca(self._internal_type, name='g_pt_%d' % _) for _ in range(self._dimension)
        ]
        self._stack_rp_pt = [
            builder.alloca(self._internal_type, name='g_rp_%d' % _) for _ in range(self._dimension)
        ]

    def _llvm_generate_bsp_lookup_table(self, module):
        """
        Generates the lookup table from plane_index -> region_index.

        :param module:
        :return: False if no bsp lookup table is required, True
        """
        # If we only have one sub region, we don't need a bsp lookup
        if len(self._ss.get_subregions()) == 1:
            return False

        # We have the plane index with 32 bit wide
        plane_index_t = ir.IntType(int(32))
        plane_index_ta = ir.ArrayType(plane_index_t, len(self._ss.get_index()))

        # Convert the packed index to the correct type
        bsp_index_c = plane_index_ta([
            plane_index_t(_) for _ in self._ss.get_index()
        ])

        # Setup the global variables
        bsp_index = ir.GlobalVariable(module, plane_index_ta, "bsp_index", self.constant_addrspace)
        bsp_index.global_constant = True
        bsp_index.initializer = bsp_index_c

        self.bsp_index = bsp_index
        return True


    def _llvm_generate_xfrom_lookup_table(self, module):
        """
        Generates the lookup table for the sub-region transforms.

        :param module:
        :return:
        """
        if len(self._ss.get_subregions()) <= 1:
            return False

        # Take a look to ensure that all the sub regions are non-redundant
        for ref_idx in self._ss.get_index():
            if ref_idx >= self._ss.get_unique_subregion_count():
                raise ValueError("Index pointed to a redundant region!")

        # Take a look at all the transforms, and choose an appropriate type for them
        transform_type = ["integer", 32]
        shift_type = ["integer", 32]
        no_shift = True

        # We can just blast thru the non redundant sub regions
        for ridx in range(self._ss.get_unique_subregion_count()):
            x, t = self._ss.get_subregions()[ridx].get_transform()

            for value in x.list():
                if round(value) != value:
                    transform_type[0] = "float"

            for value in t.list():
                if round(value) == value:
                    if int(value) != 0:
                        no_shift = False
                else:
                    no_shift = False
                    shift_type[0] = "float"

        self._transform_type = ir.IntType(transform_type[1]) if transform_type[0] == "integer" else ir.FloatType()
        self._shift_type = ir.IntType(shift_type[1]) if shift_type[0] == "integer" else ir.FloatType()

        # Create the struct type
        msize = self._dimension
        msize2 = msize * msize
        msize = 0 if no_shift else msize

        self._xform_struct = ir.LiteralStructType(
            ([self._transform_type] * msize2) + ([self._shift_type] * msize),
            packed=True
        )
        self._xform_struct_array = ir.ArrayType(self._xform_struct, self._ss.get_unique_subregion_count())

        # Pack the LUT
        self.xform_lut = ir.GlobalVariable(module, self._xform_struct_array, "xform_lookup", self.constant_addrspace)

        sl = []
        for ridx in range(self._ss.get_unique_subregion_count()):
            x, t = self._ss.get_subregions()[ridx].get_transform()

            a = [self._transform_type(_) for _ in x.list()]
            if not no_shift:
                b = [self._shift_type(_) for _ in t.list()]
            else:
                b = []

            sl += [
                self._xform_struct(a + b)
            ]

        self.no_shift = no_shift
        self.xform_lut.global_constant = True
        self.xform_lut.initializer = self._xform_struct_array(sl)

    def _llvm_generate_ref_region_index_lookup_table(self, module):
        """
        Generates the polynomial index table

        :param module:
        :return: True if the table exists, false if it is unneccesary
        """
        if len(self._ss.get_ref_subregions()) <= 1:
            return False

        unique_sr_count = self._ss.get_unique_subregion_count()
        subregions = self._ss.get_subregions()

        poly_index_t = ir.IntType(int(32))
        poly_index_ta = ir.ArrayType(poly_index_t, unique_sr_count)

        # Convert the packed index to the correct type
        index = []
        for ridx in range(unique_sr_count):
            polyindex = subregions[ridx].get_ref_region()
            index += [poly_index_t(polyindex)]

        polyindex_inst = poly_index_ta(index)

        # Setup the global variables
        poly_index = ir.GlobalVariable(module, poly_index_ta, "poly_index", self.constant_addrspace)
        poly_index.global_constant = True
        poly_index.initializer = polyindex_inst

        self.poly_index_lut = poly_index

        # For GPU architectures, declare this in constant memory
        # bsp_index.addrspace = 0 # TODO  change to constant memory addr space

        return True

    def _llvm_generate_lattice_site_lookup_table(self, module):
        # TODO: Remove old/new lattice site lookup code
        subregions = self._ss.get_subregions()
        ref_regions = self._ss.get_ref_subregions()
        num_subregions = len(subregions)

        if num_subregions <= 1:
            return None

        pts_per_sr = max([_.number_of_pts_per_reconstruction() for _ in self._ss.get_ref_subregions()])

        # Setup the LLVM types for the index
        lsc_index_t = ir.IntType(32)
        ls_index_t = ir.ArrayType(lsc_index_t, self._dimension)
        lslst_index_t = ir.ArrayType(ls_index_t, pts_per_sr)
        lsidx_t = ir.ArrayType(lslst_index_t, num_subregions)
        #
        # self._transform_type = ir.IntType(transform_type[1]) if transform_type[0] == "integer" else ir.FloatType()
        # self._shift_type = ir.IntType(shift_type[1]) if shift_type[0] == "integer" else ir.FloatType()
        #
        # self._xform_struct = ir.LiteralStructType(
        #     ([self._transform_type] * msize2) + ([self._shift_type] * msize),
        #     packed=True
        # )
        # self._xform_struct_array = ir.ArrayType(self._xform_struct, self._ss.get_unique_subregion_count())

        # Pack the LUT
        # self.xform_lut = ir.GlobalVariable(module, self._xform_struct_array, "xform_lookup", self.constant_addrspace)


        lattice_site_int_type = ir.IntType(32)
        region_lsindex = ir.LiteralStructType([lattice_site_int_type for _ in range(pts_per_sr*self._dimension)],
                                               packed=True)
        ls_index_type = ir.ArrayType(region_lsindex, num_subregions)

        lattice_site_lookup = []
        unroll_ls_lookup = []
        for sridx, sub_region in enumerate(subregions):
            ref_region = sub_region.get_ref_region()
            ref_region = ref_regions[ref_region]
            region_list = []
            unroll_list = []

            transform, shift = sub_region.get_transform()
            for ls in ref_region.get_optimal_nn_lookups():
                ls = transform * vector(ls) + shift
                unroll_list += [lattice_site_int_type(_) for _ in ls]

                ls = ls_index_t([lsc_index_t(_) for _ in ls])
                region_list += [ls]


            region_list = lslst_index_t(region_list)
            lattice_site_lookup += [region_list]
            # print(len(unroll_list))
            unroll_ls_lookup += [region_lsindex(unroll_list)]

        lattice_site_lookup = lsidx_t(lattice_site_lookup)
        lattice_site_lookup_unroll = ls_index_type(unroll_ls_lookup)

        ls_index = ir.GlobalVariable(module, lsidx_t, "ls_index", self.constant_addrspace)
        ls_index.global_constant = True
        ls_index.initializer = lattice_site_lookup
        #
        ls_index_u = ir.GlobalVariable(module, ls_index_type, "ls_index_u", self.constant_addrspace)
        ls_index_u.global_constant = True
        ls_index_u.initializer = lattice_site_lookup_unroll

        self.ls_index = ls_index
        self.ls_index_u = ls_index_u
        self.ls_struct_type = region_lsindex

    def llvm_generate_lookup_tables(self, module):
        """
        Generates all the required tables for this interpolation space

        :param module:
        :return:
        """
        self._llvm_generate_bsp_lookup_table(module)
        self._llvm_generate_xfrom_lookup_table(module)
        self._llvm_generate_ref_region_index_lookup_table(module)
        self._llvm_generate_lattice_site_lookup_table(module)

    def llvm_generate_reconstruct_preamble(self, module):
        """

        :param module:
        :return:
        """
        # Build types
        input_types = ([self._input_type] * self._dimension)
        ntextures = len(self._ss.coset_vectors()[1])*len(self._ss.get_L_cosets())

        # mf (memory fetch) is an input type that depends on the type of code being generated (and architecture)
        if self.lookup_type == 'ABSTRACT':
            # Abstract lookups take a pointer to a function
            mf = [self._abstract_memfetch_lk_type.as_pointer()]
        elif self.lookup_type == 'COSET':
            # Same with coset (but the function call is different)
            mf = [self._cluptp.as_pointer()]
        elif self.lookup_type == 'COSET_RAW':
            # Coset raw unwraps a simple struct that contains all the info about a given lattice
            # instance
            mf = [self._raw_lookup.as_pointer()]
        elif self.lookup_type == 'CPU_LIN':
            mf = [self._clerptp.as_pointer()]  # Debug jig, a linear fetch emulated in CPU
        elif self.lookup_type == 'TEX' and 'nvptx64' in self._arch and self._vendor == 'nvidia':
            # Int type is the type textures are aliased to, so those are the input
            # to the mem lookup code
            mf = [ir.IntType(int(64)) for _ in range(ntextures)]
        elif self.lookup_type == 'TEX_LIN' and 'nvptx64' in self._arch and self._vendor == 'nvidia':
            mf = [ir.IntType(int(64)) for _ in range(ntextures)]
        else:
            raise ValueError("Invalid lookup_type?")

        # Create a new type for our function
        fn_type = ir.FunctionType(self._output_type, tuple(input_types + mf))
        self.function = ir.Function(module, fn_type, name=self.reconstruction_name)
        memfetch = self.function.args[self._dimension:]
        xvars = self.function.args[0:self._dimension]
        self._textures = self.function.args[self._dimension:]
        self.memfetch = memfetch
        #
        self._entry_blk = self.function.append_basic_block(name="entry")
        builder = ir.IRBuilder(self._entry_blk)
        # Emit the allocate instructions for the stack variables
        self.llvm_build_stack(builder)

        # Store the inputs in the stack variables
        [builder.store(x, sx) for x, sx in zip(xvars, self._stack_xvars)]

        # We need to store some info for the coset lookups
        if self.lookup_type == "COSET_RAW":

            # First we have to read all the lattice boundary info
            # addr = builder.gep(memfetch[0], [ir.IntType(32)(int(0)), ir.IntType(32)(int(0))])
            # input_coset_count = builder.load(addr)
            # addr = builder.gep(memfetch[0], [ir.IntType(32)(int(0)), ir.IntType(32)(int(1))])
            # input_dimension = builder.load(addr)
            coset = builder.gep(memfetch[0], [ir.IntType(32)(int(0)), ir.IntType(32)(int(2))])
            coset = builder.load(coset, name='lattice_info')

            coset_bounds = []
            coset_base = []
            num_cosets = self._ss.coset_count()
            for coset_idx in range(num_cosets):
                buffer = builder.gep(coset, [ir.IntType(32)(int(coset_idx)), ir.IntType(32)(int(1))])
                coset_base += [builder.load(buffer, name='buffer_base')]
                bounds = []
                for dim_idx in range(self._dimension):
                    bnd = builder.gep(coset, [ir.IntType(32)(int(coset_idx)), ir.IntType(32)(int(0))])
                    bnd = builder.load(bnd)
                    bnd = builder.gep(bnd, [ir.IntType(32)(int(dim_idx))])
                    bnd = builder.load(bnd, name='bound')

                    if self.lane_count > 1:
                        bcst = builder.insert_element(self._index_type(ir.Undefined),
                                                      bnd,
                                                      ir.IntType(32)(0), name='bnd_bcst')
                        bnd = builder.shuffle_vector(bcst, self._index_type(ir.Undefined), self._index_type(None))

                    bounds += [bnd]
                coset_bounds += [bounds]

            if self.lane_count > 1:
                bcst_bases = []
                for base in coset_base:
                    # ptr_type = ir.VectorType(self._internal_type.as_pointer(), self.lane_count)
                    ptr_size = self.ptr_size
                    index_type = ir.VectorType(ir.IntType(32), self.lane_count)
                    ptr_type = ir.VectorType(ir.IntType(ptr_size), self.lane_count)

                    bcst = builder.insert_element(ptr_type(ir.Undefined),
                                                  builder.ptrtoint(base, ir.IntType(ptr_size)),
                                                  ir.IntType(32)(0), name='bsp_bcast')
                    bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), index_type(None))
                    bcst = builder.inttoptr(bcst,
                                            ir.VectorType(self._memfetch_base.as_pointer(), self.lane_count)
                                            )
                    bcst_bases += [bcst]
                coset_base = bcst_bases

            # Write the base addresses
            [builder.store(cb, sb) for sb, cb in zip(self._stack_buffers, coset_base)]

            for coset_idx in range(num_cosets):
                [builder.store(cb, sb) for sb, cb in zip(self._stack_bounds[coset_idx], coset_bounds[coset_idx])]

        return builder

    def _llvm_generate_vec_poly_index(self, builder):
        """
        Returns a the bsp index

        Inputs: builder ir builder.
        x_variables
        """

        if len(self._ss._ref_subregions) == 1:
            return None

        ptr_size = self.ptr_size
        index_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(32), self.lane_count)
        ptr_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(ptr_size), self.lane_count)

        index = builder.load(self._stack_bspindex)

        if self.lane_count == 1:
            index_ptr = builder.gep(self.poly_index_lut, [ir.IntType(32)(0), index])
            index = builder.load(index_ptr)
            builder.store(index, self._stack_polyindex)
        else:
            # Broadcast the poly index table pointer ...
            bcst = builder.insert_element(ptr_type(ir.Undefined),
                                          self.poly_index_lut.ptrtoint(ir.IntType(ptr_size)),
                                          ir.IntType(32)(0), name='bsp_bcast')
            bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), index_type(None))
            bcst = builder.inttoptr(bcst, ir.VectorType(ir.IntType(32).as_pointer(), self.lane_count))

            index_ptr = builder.gep(bcst, [index])
            index = self.read_mem(index_ptr, builder)
            builder.store(index, self._stack_polyindex)
        return None

    def _llvm_generate_vec_bsp_index(self, builder):
        """
        Returns a the bsp index

        Inputs: builder ir builder.
        x_variables
        """
        ptr_size = self.ptr_size
        index_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(32), self.lane_count)
        ptr_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(ptr_size), self.lane_count)
        index = None

        # We only need to generate this test if we have more than one
        # sub-region
        if len(self._ss.get_subregions()) > 1:

            # Generate the plane test
            for i, (n, d) in enumerate(reversed(self._ss.get_plane_list())):
                running_sum = None

                for x_i, n_i in zip(self._stack_xvars, n):

                    x_i = builder.load(x_i)
                    partial = builder.fmul(x_i, self._internal_type(n_i), flags=['fast'])

                    if running_sum is not None:
                        running_sum = builder.fadd(partial, running_sum, flags=['fast'])
                    else:
                        running_sum = partial

                p = builder.fcmp_ordered(">=", running_sum, self._internal_type(float(-d)))

                idx = builder.select(p, index_type(int(1 << i)), index_type(0))

                if index is not None:
                    index = builder.load(self._stack_bspindex)
                    index = builder.or_(index, idx)
                    builder.store(index, self._stack_bspindex)
                else:
                    builder.store(idx, self._stack_bspindex)
                    index = idx

            # Do the modulus optimization
            index = builder.load(self._stack_bspindex)
            index = builder.urem(index, index_type(int(self._ss.get_modulus())))
            builder.store(index, self._stack_bspindex)

            if self.lane_count == 1:
                index = builder.load(self._stack_bspindex)
                index_ptr = builder.gep(self.bsp_index, [ir.IntType(32)(0), index])
                index = builder.load(index_ptr)
                builder.store(index, self._stack_bspindex)
            else:
                # Broadcast the index table pointer ...
                bcst = builder.insert_element(ptr_type(ir.Undefined),
                                              self.bsp_index.ptrtoint(ir.IntType(ptr_size)),
                                              ir.IntType(32)(0), name='bsp_bcast')
                bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), index_type(None))
                bcst = builder.inttoptr(bcst, ir.VectorType(ir.IntType(32).as_pointer(), self.lane_count))
                index_ptr = builder.gep(bcst, [index])
                index = self.read_mem(index_ptr, builder)
                builder.store(index, self._stack_bspindex)

    def _llvm_apply_transform(self, var_input, builder, target=None, transpose=False):
        product = [None] * self._dimension
        vec = [builder.load(_) for _ in var_input]
        for i in range(self._dimension):
            for j in range(self._dimension):
                if not transpose:
                    _ = builder.load(self._stack_xform[i][j])
                else:
                    _ = builder.load(self._stack_xform[j][i])
                _ = builder.fmul(_, vec[j], flags=['fast'])
                if product[i] is not None:
                    product[i] = builder.fadd(product[i], _, flags=['fast'])
                else:
                    product[i] = _

        if target is not None and len(target) == len(var_input):
            var_input = target

        for result, loc in zip(product, var_input):
            builder.store(result, loc)

    @classmethod
    def _llvm_apply_shift(cls, builder, vec_a, vec_b, subtract=False, target=None):
        if len(vec_a) != len(vec_b):
            raise ValueError("Lengths of vec_a and vec_b differ!")

        fn = builder.fsub if subtract else builder.fadd
        result = []
        for (a, b) in zip(vec_a, vec_b):
            a = builder.load(a)
            b = builder.load(b)
            result += [fn(a, b, flags=['fast'])]

        if target is not None and len(target) == len(vec_a):
            vec_a = target

        for res, loc in zip(result, vec_a):
            builder.store(res, loc)

    def _llvm_generate_transform_load(self, builder):
        """
        Fills the transform variables

        :param builder:
        :return:
        """

        if len(self._ss.get_subregions()) <= 1:
            return

        ptr_size = self.ptr_size
        index_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(32), self.lane_count)
        ptr_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(ptr_size), self.lane_count)

        sregion_index = builder.load(self._stack_bspindex)

        # Load the offset of the transform
        if self.lane_count == 1:
            subregion_transform = builder.gep(self.xform_lut, [ir.IntType(32)(0), sregion_index])
        else:
            bcst = builder.insert_element(ptr_type(ir.Undefined),
                                          self.xform_lut.ptrtoint(ir.IntType(ptr_size)),
                                          ir.IntType(32)(0), name='xform_map')

            bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), index_type(None))
            subregion_transform = builder.inttoptr(bcst,
                                                   ir.VectorType(self._xform_struct.as_pointer(), self.lane_count))

        for i in range(self._dimension):
            for j in range(self._dimension):
                offset = i * self._dimension + j
                if self.lane_count == 1:
                    mep = builder.gep(subregion_transform, [ir.IntType(32)(0), ir.IntType(32)(int(offset))])
                else:
                    mep = builder.gep(subregion_transform, [sregion_index, ir.IntType(32)(int(offset))])
                mel = self.read_mem(mep, builder)
                mel = self._int_to_float(mel, mel.type, self._internal_type, builder)
                builder.store(mel, self._stack_xform[i][j])

        if self.no_shift:
            return

        offset = self._dimension * self._dimension
        for index in range(self._dimension):
            if self.lane_count == 1:
                mp = builder.gep(subregion_transform, [ir.IntType(32)(0), ir.IntType(32)(int(offset + index))])
                sel = builder.load(mp)
                sel = self._llvm_convert_to_internal_type(builder, sel)
            else:
                mep = builder.gep(subregion_transform, [sregion_index, ir.IntType(32)(int(offset + index))])
                mel = self.read_mem(mep, builder)
                sel = self._int_to_float(mel, mel.type, self._internal_type, builder)

            builder.store(sel, self._stack_tshift[index])

        return None

    def _llvm_generate_rho_vec(self, builder):
        # Generate the pre-rho shift
        variables = {}
        expressions = []
        for i, sh in enumerate(self._ss._rho_shift):
            variables[f"ss{i}"] = self._internal_type(sh)
            variables[f"x{i}"] = self._stack_xvars[i]
            expressions += [f"x{i} - ss{i}"]

        for i, expression in enumerate(expressions):
            ast = self.parse_expression(expression, variables)
            self.generate_expression(ast, builder, variables, target=self._stack_txvars[i])

        # Some lattices have easy to express \rho functions, so we have these cases broken out here
        if self._ss._rho_type == "PP":
            # x_variables, i_variables = self._llvm_generate_rho_pp_vec(builder, x_variables)
            raise ValueError("Not Implemented")
        else:
            if self._dimension == 3:
                if self._ss.get_lattice_hash() == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                    self._llvm_generate_rho_bccind(builder, self._stack_txvars)

                elif self._ss.get_lattice_hash() == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                    self._llvm_generate_rho_fccind(builder, self._stack_txvars)

                elif self._ss.get_lattice_hash() == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                    self._llvm_generate_rho_cartesian(builder, self._stack_txvars)
                else:
                    # x_variables, i_variables = self._llvm_generate_rho_generic_vec(builder, x_variables)
                    raise ValueError("Not Implemented")

            elif self._ss.is_lattice_cartesian():
                self._llvm_generate_rho_cartesian(builder, self._stack_txvars)

            else:
                raise ValueError("Not Implemented")

        # Shift to the reference region
        for i in range(self._dimension):
            x = builder.load(self._stack_xvars[i])
            k = builder.load(self._stack_rho[i])
            k = builder.sitofp(k, self._internal_type)
            sx = builder.fsub(x, k, flags=['fast'])

            if self._ss._reflective:
                print('we should be using reflective symmetry?')
                cmp = builder.fcmp_ordered(">=", sx, self._internal_type(0))
                wt = builder.select(cmp, self._internal_type(1.0), self._internal_type(-1.0), name="refl_%d" % i)
                builder.store(wt, self._stack_xflip[i])
                sx = builder.fmul(wt, sx)

            builder.store(sx, self._stack_xvars[i])
        return None

    def llvm_generate_vec(self):
        # First we generate the module that hold all of our code
        self._llvm_generate_vec_module()


        # Then we generate all the lookup tables
        self.llvm_generate_lookup_tables(self._llvm_module)

        if self._bake_raw_fetches:
            self.llvm_generate_bake_function(self._llvm_module)

        builder = self.llvm_generate_reconstruct_preamble(self._llvm_module)
        builder.store(self._internal_type(0.0), self._stack_result)
        builder.store(self._index_type(0), self._stack_bspindex)

        if self.test_stop.startswith("input_x"):
            ret = builder.load(self._stack_xvars[0])
            builder.ret(ret)
            return str(self._llvm_module)

        if self.test_stop.startswith("boundary"):
            rho = int(self.test_stop[8:])
            print(f"stopping at {rho}")
            ret = builder.load(self._stack_bounds[0][rho])
            ret = builder.sitofp(ret, self._internal_type)
            builder.ret(ret)
            return str(self._llvm_module)

        if len(self._ss.get_L_cosets()) == 1:


            # Generate rho(x)
            self._llvm_generate_rho_vec(builder)
            if self.test_stop.startswith("rho"):
                rho = int(self.test_stop[3:])
                ret = builder.load(self._stack_rho[rho])
                ret = builder.sitofp(ret, self._internal_type)
                builder.ret(ret)
                return str(self._llvm_module)

            # Generate the BSP index
            self._llvm_generate_vec_bsp_index(builder)
            if self.test_stop == "bsp_index":
                ret = builder.load(self._stack_bspindex)
                ret = builder.sitofp(ret, self._internal_type)
                builder.ret(ret)
                return str(self._llvm_module)

            self._llvm_generate_transform_load(builder)

            # If we have multiple sub regions, then we need to generate the shift
            if len(self._ss.get_subregions()) > 1:
                if not self.no_shift:
                    xvars = self._llvm_apply_shift(builder,
                                                   self._stack_xvars,
                                                   self._stack_tshift,
                                                   target=self._stack_xvars,
                                                   subtract=True)

                self._llvm_apply_transform(self._stack_xvars, builder, target=self._stack_xvars, transpose=True)

                if len(self._ss._ref_subregions) > 1:
                    self._llvm_generate_vec_poly_index(builder)
                    if self.test_stop == "poly_index":
                        ret = builder.load(self._stack_polyindex)
                        ret = builder.sitofp(ret, self._internal_type)
                        builder.ret(ret)
                        return str(self._llvm_module)


            if self.lookup_type == "CPU_LIN" or self.lookup_type == "TEX_LIN":
                assert self.lane_count == 1
                if self.branchless:
                    self._llvm_generate_lerp_conv_sum_bpred(builder)
                    presult = builder.load(self._stack_coset_result)
                    result = builder.load(self._stack_result)
                    result = builder.fadd(presult, result, flags=['fast'])
                    builder.store(result, self._stack_result)
                else:
                    raise ValueError("Branchy code not supported for lin yet")
            else:
                print("GEN")
                self._llvm_generate_no_lerp_conv_sum_bpred(builder)
                presult = builder.load(self._stack_coset_result)
                result = builder.load(self._stack_result)
                result = builder.fadd(presult, result, flags=['fast'])
                builder.store(result, self._stack_result)
        else:
            # Create a basic block for the loop
            loop_blk = self.function.append_basic_block("coset_loop")
            exit_blk = self.function.append_basic_block("exit")

            # Make sure the first thing we do is jump into the loop
            builder.position_at_end(self._entry_blk)
            builder.branch(loop_blk)

            # Now let's build the loop preamble
            builder.position_at_end(loop_blk)

            # we have two
            loop_idx = builder.phi(ir.IntType(int(32)))
            loop_idx.add_incoming(ir.IntType(int(32))(0), self._entry_blk)

            sum_prev = builder.phi(self._internal_type)
            sum_prev.add_incoming(self._internal_type(0.0), self._entry_blk)

            raise ValueError('Verified up to this point...')

            # Now we can emit the interpolation code
            shifts = self._llvm_generate_outer_coset_index(builder, loop_idx, outer_coset_lut)
            xvars_shifted = [builder.fsub(x, s) for x, s in zip(xvars,shifts)]

            xvars_shifted, ivars = self._llvm_generate_rho(builder, xvars_shifted)
            srindex      = self._llvm_generate_bsp_index(builder, xvars_shifted)
            polyindex    = self._llvm_generate_poly_index(builder, srindex, refpp_lut)

            mat          = self._llvm_load_transform_matrix(builder, srindex, False)
            shift        = self._llvm_load_shift(builder, srindex)

            if not self.no_shift:
                xvars_shifted    = self._llvm_apply_shift(builder, xvars_shifted, shift, True)   # Subtract
            xvars_shifted        = self._llvm_apply_transform(builder, mat, xvars_shifted, True) # Transform with transposed matrix

            #########
            #
            # print("buidling thing")
            # builder.ret(polyindex)

            if self.lookup_type == "CPU_LIN" or self.lookup_type == "TEX_LIN":
                if self.branchless:
                    presult = self._llvm_generate_lerp_conv_sum_bpred_inconsistent(builder, xvars_shifted, ivars, mat, shift, srindex, polyindex, ls_lut, memfetch)
                else:
                    raise ValueError("Branchy code not supported for lin yet")
            else:
                presult = self._llvm_generate_no_lerp_conv_sum_bpred(builder, xvars_shifted, ivars, mat, shift, srindex, polyindex, ls_lut, memfetch, loop_idx, shifts)
            presult = self._llvm_convert_to_output_type(builder, presult)

            #######

            result = builder.fadd(sum_prev, presult)



            # Finally, do the branch back
            loop_inc = builder.add(loop_idx, ir.IntType(32)(1))
            cond = builder.icmp_unsigned("<", loop_inc, ir.IntType(int(32))(len(self._ss.get_L_cosets())))

            loop_idx.add_incoming(loop_inc, loop_blk)
            sum_prev.add_incoming(result, loop_blk)

            builder.cbranch(cond, loop_blk, exit_blk)
            builder.position_at_end(exit_blk)

        result = builder.load(self._stack_result)
        builder.ret(result)
        sout = str(self._llvm_module)
        # print(sout)
        return sout

    """
    ####################################
    # LLVM Code generation functions
    ####################################

    These are the workhorses of this class -- generally they don't
    depend on instance variables when they don't have to. This modular
    design allows me to test each one in isolation.
    """

    @staticmethod
    def _int_to_float(variable, type_from, type_to, builder, signed=True):
        if isinstance(type_from, ir.VectorType) != isinstance(type_to, ir.VectorType):
            raise ValueError('Tried to cast scalar<->vector or vector<->scalar!')

        element_from = type_from
        element_to = type_to
        if isinstance(type_from, ir.VectorType):
            element_from = type_from.element
            element_to = type_to.element

        # Casting to the same type is the identity mapping
        if element_from == element_to:
            return variable

        if signed:
            return builder.sitofp(variable, type_to)
        return builder.uitofp(variable, type_to)

    @staticmethod
    def _float_to_int(variable, type_from, type_to, builder, signed=True):
        if isinstance(type_from, ir.VectorType) != isinstance(type_to, ir.VectorType):
            raise ValueError('Tried to cast scalar<->vector or vector<->scalar!')

        element_from = type_from
        element_to = type_to
        if isinstance(type_from, ir.VectorType):
            element_from = type_from.element
            element_to = type_to.element

        # Casting to the same type is the identity mapping
        if element_from == element_to:
            return variable

        if signed:
            return builder.fptosi(variable, type_to)
        return builder.fptoui(variable, type_to)

    @staticmethod
    def _float_cast(variable, type_from, type_to, builder):
        if isinstance(type_from, ir.VectorType) != isinstance(type_to, ir.VectorType):
            raise ValueError('Tried to cast scalar<->vector or vector<->scalar!')

        element_from = type_from
        element_to = type_to
        if isinstance(type_from, ir.VectorType):
            element_from = type_from.element
            element_to = type_to.element

        # Casting to the same type is the identity mapping
        if element_from == element_to:
            return variable

        # If we're casting from something big to something small
        if element_from.get_abi_size() > element_to.get_abi_size():
            # Truncate
            if isinstance(variable, list):
                return [builder.fptrunc(_, type_to) for _ in variable]
            return builder.fptrunc(variable, type_to)

        # Otherwise, extend
        if isinstance(variable, list):
            return [builder.fpext(_, type_to) for _ in variable]
        return builder.fpext(variable, type_to)


    def _llvm_generate_rho_generic_vec(self, builder, x_variables):
        """
        Returns a tuple of x_variables, and i_variables (which is
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        D, offsets = self._ss._lattice.coset_structure()
        internal_type = self._internal_type
        best = None

        for i, offset in enumerate(offsets):

            # Shift by offset
            xt = [builder.fsub(x_i, internal_type(component), flags=['fast'])
                  for (x_i, component) in zip(x_variables, offset)]

            # Scale by 1/D
            xt = [builder.fdiv(x_i, internal_type(component), flags=['fast'])
                  for (x_i, component) in zip(xt, D.diagonal())]

            # Round -- note that on some architectures, a call to round() is generated -- I'd prefer
            # an inline version, so I inlined a 'round' implementation here.
            xt = [builder.fadd(x_i, internal_type(0.5), flags=['fast'])
                  for (x_i, component) in zip(xt, offset)]
            xt = [builder.fptosi(x_i, ir.IntType(32) if self.lane_count == 1
                  else ir.VectorType(ir.IntType(32), self.lane_count)) for x_i in xt]
            xt = [builder.sitofp(x_i, internal_type) for x_i in xt]

            # xt = [builder.call(self._round, [_]) for _ in xt]

            # Scale by D
            xt = [builder.fmul(x_i, internal_type(component), flags=['fast'])
                  for (x_i, component) in zip(xt, D.diagonal())]
            xt = [builder.fadd(x_i, internal_type(component), flags=['fast'])
                  for (x_i, component) in zip(xt, offset)]

            d = internal_type(0)
            for (a, b) in zip(x_variables, xt):
                t = builder.fsub(a, b, flags=['fast'])
                t = builder.fmul(t, t, flags=['fast'])
                d = builder.fadd(t, d, flags=['fast'])

            if i == 0:
                best = xt
                bestd = d
            else:
                p = builder.fcmp_ordered("<", d, bestd)
                bestd = builder.select(p, d, bestd)
                best = [builder.select(p, new, old) for (new, old) in zip(xt, best)]
        return x_variables, best

    def _llvm_generate_rho_pp_vec(self, builder, x_variables):
        """
        Returns a tuple of x_variables, and i_variables (which is
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """

        Ti = self._ss._rho_ppiped.inverse()
        T = self._ss._rho_ppiped

        # Inverse transform
        iv = []
        for i, xv in enumerate(x_variables):
            tmp = []
            for j, xw in enumerate(x_variables):
                tmp.append(builder.fmul(xw,  self._internal_type(Ti[i, j]), flags=['fast']))

            result = tmp[0]
            for k in tmp[1:]:
                result = builder.fadd(result, k)

            iv += [builder.call(self._floor, [result])]

        # Forward
        i_variables = []
        for i, xv in enumerate(x_variables):
            tmp = []
            for j, xw in enumerate(iv):
                tmp.append(builder.fmul(xw,  self._internal_type(T[i, j]), flags=['fast']))

            result = tmp[0]
            for l in tmp[1:]:
                result = builder.fadd(result, l)
            i_variables += [result]

        return x_variables, i_variables

    @staticmethod
    def read_mem(bcst, builder):
        # print(bcst.type)
        # if isinstance(bcst.type, ir.PointerType) and isinstance(bcst.type.pointee, ir.VectorType):
        #     print("!!", bcst.type.pointee)

        if isinstance(bcst.type, ir.VectorType):
            type = bcst.type.element.pointee

            ptrs = [builder.extract_element(bcst, ir.IntType(32)(_)) for _ in range(bcst.type.count)]
            loads = [builder.load(_) for _ in ptrs]
            buildvec = None
            for _ in range(bcst.type.count):
                input = ir.VectorType(type, bcst.type.count)(ir.Undefined) if buildvec is None else buildvec
                buildvec = builder.insert_element(input, loads[_], ir.IntType(32)(_))
            return buildvec
        return builder.load(bcst, 'read_mem')



    def _llvm_generate_outer_coset_table(self, module):
        coset_shifts = self._ss.get_L_cosets()
        if len(coset_shifts) == 1:
            return None

        # Determine the max bit-width we need for these
        bw = max([_bw(__) for __ in sum([list(_) for _ in coset_shifts], [])])


        ocs_index_t = ir.IntType(bw)
        ocsl_index_t = ir.ArrayType(ocs_index_t, self._ss._s)
        ocl_t = ir.ArrayType(ocsl_index_t, len(coset_shifts))

        outer_coset_lookup = []
        for coset_shift in coset_shifts:
            coset_offsets = []

            for x in coset_shift:
                coset_offsets += [ocs_index_t(x)]
            outer_coset_lookup += [ocsl_index_t(coset_offsets)]
        coset_lk = ocl_t(outer_coset_lookup)

        coset_index = ir.GlobalVariable(module, ocl_t, "ocs_index", self.constant_addrspace)
        coset_index.global_constant = True
        coset_index.initializer     = coset_lk

        return coset_index

    def _llvm_generate_rho_bccind(self, builder, stack_x):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        i_variables = [None]*3

        # Stage 1
        _ = builder.load(stack_x[0])
        _ = builder.fdiv(_, self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_, self._internal_type(0.5), flags=['fast'])
        _ = builder.fptosi(_, self._index_type)
        _ = builder.sitofp(_, self._internal_type)
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        i_variables[0] = _
        
        _ = builder.load(stack_x[1])
        _ = builder.fdiv(_, self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_, self._internal_type(0.5), flags=['fast'])
        _ = builder.fptosi(_, self._index_type)
        _ = builder.sitofp(_, self._internal_type)
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        i_variables[1] = _

        _ = builder.load(stack_x[2])
        _ = builder.fdiv(_, self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_, self._internal_type(0.5), flags=['fast'])
        _ = builder.fptosi(_, self._index_type)
        _ = builder.sitofp(_, self._internal_type)
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        i_variables[2] = _
        
        # Sd
        _ = builder.fsub(i_variables[0], builder.load(stack_x[0]), flags=['fast'])
        a = builder.fmul(_, _, flags=['fast'])
        
        _ = builder.fsub(i_variables[1], builder.load(stack_x[1]), flags=['fast'])
        b = builder.fmul(_, _, flags=['fast'])
        
        a = builder.fadd(a, b, flags=['fast'])
        _ = builder.fsub(i_variables[2], builder.load(stack_x[2]))
        _ = builder.fmul(_, _, flags=['fast'])
        sd = builder.fadd(a, _, flags=['fast'], name="sd")
        
        
        tx = [
            builder.fsub(builder.load(stack_x[0]), self._internal_type(1)),
            builder.fsub(builder.load(stack_x[1]), self._internal_type(1)),
            builder.fsub(builder.load(stack_x[2]), self._internal_type(1)),
        ]
        
        # Stage 3
        _ = tx[0]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.fadd(_, self._internal_type(0.5), flags=['fast'])
        _ = builder.fptosi(_, self._index_type)
        _ = builder.sitofp(_, self._internal_type)
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_,  self._internal_type(1.0), flags=['fast'])
        tx[0] = _
        
        _ = tx[1]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.fadd(_, self._internal_type(0.5), flags=['fast'])
        _ = builder.fptosi(_, self._index_type)
        _ = builder.sitofp(_, self._internal_type)
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_,  self._internal_type(1.0), flags=['fast'])
        tx[1] = _

        _ = tx[2]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.fadd(_, self._internal_type(0.5), flags=['fast'])
        _ = builder.fptosi(_, self._index_type)
        _ = builder.sitofp(_, self._internal_type)
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_,  self._internal_type(1.0), flags=['fast'])
        tx[2] = _

        # Sd
        _ = builder.fsub(tx[0], builder.load(stack_x[0]))
        a = builder.fmul(_, _, flags=['fast'])
        
        _ = builder.fsub(tx[1], builder.load(stack_x[1]))
        b = builder.fmul(_, _, flags=['fast'])
        
        a = builder.fadd(a, b, flags=['fast'])
        _ = builder.fsub(tx[2], builder.load(stack_x[2]))
        _ = builder.fmul(_, _, flags=['fast'])
        sd2 = builder.fadd(a,_, name="sd2")
        cmp_dist = builder.fcmp_ordered("<", sd, sd2)
        
        i_variables[0] = builder.select(cmp_dist, i_variables[0], tx[0])
        i_variables[1] = builder.select(cmp_dist, i_variables[1], tx[1])
        i_variables[2] = builder.select(cmp_dist, i_variables[2], tx[2])

        builder.store(builder.fptosi(i_variables[0], self._index_type), self._stack_rho[0])
        builder.store(builder.fptosi(i_variables[1], self._index_type), self._stack_rho[1])
        builder.store(builder.fptosi(i_variables[2], self._index_type), self._stack_rho[2])

    def _llvm_generate_rho_cartesian(self, builder, stack):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """

        variables = {}
        expressions = []

        # Shift
        for i in range(self._dimension):
            variables[f"ss{i}"] = self._internal_type(0.5)
            variables[f"x{i}"] = stack[i]
            expressions += [f"x{i} + ss{i}"]

        xvars = []
        for i, expression in enumerate(expressions):
            ast = self.parse_expression(expression, variables)
            xvars += [self.generate_expression(ast, builder, variables)]

        xt = [builder.fptosi(x_i, self._index_type, name='ccro') for x_i in xvars]
        [builder.store(ivar, rho) for ivar, rho in zip(xt, self._stack_rho)]


    def _llvm_generate_rho_fccind(self, builder, stack_x):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        i_variables = [None]*3
        x_variables = [None]*3

        x_variables[0] = builder.load(stack_x[0])
        x_variables[1] = builder.load(stack_x[1])
        x_variables[2] = builder.load(stack_x[2])

        i_i = builder.fadd(x_variables[0], self._internal_type(0.5))
        i_j = builder.fadd(x_variables[1], self._internal_type(0.5))
        i_k = builder.fadd(x_variables[2], self._internal_type(0.5))

        ai = builder.fptoui(i_i, self._index_type)
        aj = builder.fptoui(i_j, self._index_type)
        ak = builder.fptoui(i_k, self._index_type)

        i_i = builder.uitofp(ai, self._internal_type)
        i_j = builder.uitofp(aj, self._internal_type)
        i_k = builder.uitofp(ak, self._internal_type)

        trunc_type = ir.VectorType(ir.IntType(int(1)), self.lane_count) if self.lane_count > 1 else ir.IntType(int(1))
        onfcc = builder.and_(builder.add(builder.add(ai, aj), ak), self._index_type(1))
        onfcc = builder.trunc(onfcc, trunc_type)

        xx = builder.fsub(x_variables[0], i_i)
        yy = builder.fsub(x_variables[1], i_j)
        zz = builder.fsub(x_variables[2], i_k)

        xc = builder.fcmp_ordered(">=", xx, self._internal_type(0.0))
        yc = builder.fcmp_ordered(">=", yy, self._internal_type(0.0))
        zc = builder.fcmp_ordered(">=", zz, self._internal_type(0.0))

        xc = builder.select(xc, self._internal_type(1), self._internal_type(-1))
        yc = builder.select(yc, self._internal_type(1), self._internal_type(-1))
        zc = builder.select(zc, self._internal_type(1), self._internal_type(-1))

        xxa = builder.fmul(xc, xx)
        yya = builder.fmul(yc, yy)
        zza = builder.fmul(zc, zz)
        # xxa = builder.call(self._fabs, [xx])
        # yya = builder.call(self._fabs, [yy])
        # zza = builder.call(self._fabs, [zz])

        yyxx = builder.fcmp_ordered(">=", yya, xxa)
        yyzz = builder.fcmp_ordered(">=", yya, zza)
        zzxx = builder.fcmp_ordered(">=", zza, xxa)
        zzyy = builder.fcmp_ordered(">=", zza, yya)

        idx = builder.select(builder.and_(yyxx, yyzz), self._index_type(1), self._index_type(0))
        idx = builder.select(builder.and_(zzxx, zzyy), self._index_type(2), idx)

        xx = builder.select(builder.fcmp_ordered(">=", xx, self._internal_type(0.0)), self._internal_type(1), self._internal_type(-1))
        yy = builder.select(builder.fcmp_ordered(">=", yy, self._internal_type(0.0)), self._internal_type(1), self._internal_type(-1))
        zz = builder.select(builder.fcmp_ordered(">=", zz, self._internal_type(0.0)), self._internal_type(1), self._internal_type(-1))

        nonfcc = onfcc # builder.not_(onfcc, 'negate_onfcc')
        xx = builder.select(builder.and_(builder.icmp_unsigned("==", idx, self._index_type(0)), nonfcc), xx, self._internal_type(0))
        yy = builder.select(builder.and_(builder.icmp_unsigned("==", idx, self._index_type(1)), nonfcc), yy, self._internal_type(0))
        zz = builder.select(builder.and_(builder.icmp_unsigned("==", idx, self._index_type(2)), nonfcc), zz, self._internal_type(0))

        i_i = builder.fadd(xx,  i_i, flags=['fast'])
        i_j = builder.fadd(yy,  i_j, flags=['fast'])
        i_k = builder.fadd(zz,  i_k, flags=['fast'])

        builder.store(builder.fptosi(i_i, self._index_type), self._stack_rho[0])
        builder.store(builder.fptosi(i_j, self._index_type), self._stack_rho[1])
        builder.store(builder.fptosi(i_k, self._index_type), self._stack_rho[2])

        # return x_variables, i_variables

    def _llvm_convert_to_internal_type(self, builder, value):
        tp = value.type

        if str(tp)[0] == 'i':
            return builder.sitofp(value, self._internal_type)

        if tp == ir.DoubleType() and self._internal_type == ir.FloatType():
            return builder.fptrunc(value, self._internal_type)

        if tp == ir.FloatType() and self._internal_type == ir.DoubleType():
            return builder.fpext(value, self._internal_type)

        return value

    def _llvm_convert_to_output_type(self, builder, value):
        tp = value.type

        if str(tp)[0] == 'i':
            return builder.sitofp(value, self._output_type)

        if tp == ir.DoubleType() and self._output_type == ir.FloatType():
            return builder.fptrunc(value, self._output_type)

        if tp == ir.FloatType() and self._output_type == ir.DoubleType():
            return builder.fpext(value, self._output_type)

        return value

    def _llvm_generate_outer_coset_index(self, builder, index, outer_coset_lut):
        """
        REturns a count of self._s variables which correspond to the 
        outer coset shift.
        """
        shift = []
        if self.lane_count == 1:
            coset = builder.gep(outer_coset_lut, [ir.IntType(32)(0), index])
            for index in range(self._ss._s):
                mp = builder.gep(coset, [ir.IntType(32)(0), ir.IntType(32)(index)])
                sel = builder.load(mp)
                shift += [self._llvm_convert_to_internal_type(builder, sel)]
        else:
            ptr_size = int(64)
            ptr_type = ir.IntType(ptr_size)
            index_type = ir.VectorType(ir.IntType(32), self.lane_count)
            raise ValueError('Not done')
            # Broadcast the index table pointer ...
            bcst = builder.insert_element(ptr_type(ir.Undefined),
                                          self.bsp_index.ptrtoint(ir.IntType(ptr_size)),
                                          ir.IntType(32)(0), name='bsp_bcast')
            bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), index_type(None))
            self.log("Warning: Fix 2: -- intType should inherit from type of bsp")
            bcst = builder.inttoptr(bcst, ir.VectorType(ir.IntType(8).as_pointer(), self.lane_count))

            index_ptr = builder.gep(bcst, [index])
            sel = self.read_mem(index_ptr, builder)
            shift.append(self._int_to_float(sel, sel.type, self._internal_type))

        return shift

    def _llvm_generate_horner(self, builder, horner, variables, pipeline_map, callback=None, name=None):
        """
        Generate code for horner factorized polynomials
        builder: LLVMlite builder instance -- this is where the code gets emitted to
        horner:  Factorization of the polynomial
        variables: The input variables
        callback: a function that maps a constant to a mem fetch. I don't think this 
        is used anymore.
        """

        for _ in reversed(range(self._dimension)):
            v = var(f"x_{_}")
            pipeline_map[v] = self._stack_xvars[_]
            variables = [v] + variables

        def _build_power(power, constant):
            # TODO: be more clever with how we compute this
            # we could do some fast exponentiation thing, but
            # for low degree monomials, it might not be worth it
            if all([_ == 0 for _ in power]):
                try:
                    return self._internal_type(float(constant))
                except:
                    if callback is None:
                        raise ValueError("encountered memfetch in a polynomial with no callback!")
                    else:
                        return callback(constant)
            
            r = None
            for v_idx, ex in enumerate(power):
                for _ in range(int(ex)):
                    if r is None:
                        r =  builder.load(pipeline_map[variables[v_idx]])
                    else:
                        r = builder.fmul(r,  builder.load(pipeline_map[variables[v_idx]]), flags=['fast'])
            
            # Multiply in the constant if necessary
            if constant != 1.0 and constant != 1:
                try:
                    _ = self._internal_type(float(constant))
                except:
                    if callback is None:
                        raise ValueError("encountered memfetch in a polynomial with no callback!")
                    else:
                        _ = callback(constant)

                r = builder.fmul(r,  _, flags=['fast'])
            return r
        
        def _build_plist(plist):
            if len(plist) == 0:
                return self._internal_type(0)

            r = None

            for coeff, pwr in plist:
                if coeff == 0 or coeff == 0.0:
                    continue
                    
                mono = _build_power(pwr, coeff)
                if r is not None:
                    r = builder.fadd(r,  mono, flags=['fast'])
                else:
                    r = mono

            return r
        
        def _gen_horner_r(factorization):
            power = factorization['power']
            left = factorization['left']
            right = factorization['right']

            if left is None and power is None:
                return _build_plist(right)

            lp = _build_power(power, 1.0)
            l = _gen_horner_r(left)
            left = builder.fmul(lp,  l, flags=['fast'])
            
            if len(right) == 0:
                return builder.fmul(lp,  l, flags=['fast'])

            right = _build_plist(right)
            return builder.fadd(left,  right, flags=['fast'])
        
        # Build the horner eval and then scale the result
        scale = self._internal_type(horner['scale'])
        hf = _gen_horner_r(horner['horner'])
        return builder.fmul(hf,  scale, flags=['fast'])

    def _llvm_generate_coset_index(self, builder, ls):
        """
        Takes a lattice site, and determines which coset it belongs to. 
        This is the "coset index" (coset_index). It's an integer that
        indexes 
        """
        if self._dimension == 3:
            """
            These are some special cases which have easier to compute
            methods to determine which coset a point belongs to 
            """
            if self._ss.get_lattice_hash() == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                return self._llvm_generate_coset_index_bcc(builder, ls)

            if self._ss.get_lattice_hash() == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                return self._llvm_generate_coset_index_fcc(builder, ls)

            if self._ss.get_lattice_hash() == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                return self._index_type(0)

        if self._ss.is_lattice_cartesian():
            return self._index_type(0)

        raise ValueError("Unsupported lattice -- no coset decomposition")

    def _llvm_generate_coset_index_bcc(self, builder, ls):
        """
        Takes a lattice site, and determines which coset it belongs to. 
        This is the "coset index" (coset_index). It's an integer that
        indexes 

        This is a special case for the BCC lattice, where the parity
        of the Z component tells you which lattice a point belongs to
        IIRC, this also exists on the D_n lattice (or is it the D*_n)
        so that's another simplification that could be made down the line
        """
        # cs = builder.load(ls[2])
        cs = ls[2]
        cs = builder.and_(cs, self._index_type(1))
        return cs

    def _llvm_generate_coset_index_fcc(self, builder, ls):
        """
        Takes a lattice site, and determines which coset it belongs to. 
        This is the "coset index" (coset_index). It's an integer that
        indexes 

        This is a special case for the FCC lattice
        """
        # z = builder.load(ls[2])
        z = ls[2]
        z = builder.and_(z, self._index_type(1))

        # y = builder.load(ls[1])
        y = ls[1]
        y = builder.and_(y, self._index_type(1))

        z = builder.shl(z, self._index_type(1))
        return builder.or_(z, y)

    def _llvm_generate_coset_scale_bcc(self, builder, coset_idx, ls):
        is_floating_point = ls[0].type == self._internal_type

        if is_floating_point:
            ls = [builder.fsub(_, builder.sitofp(coset_idx, self._internal_type)) for _ in ls]
            ls = [builder.fmul(_, self._internal_type(0.5)) for _ in ls]
        else:
            ls = [builder.sub(_, coset_idx) for _ in ls]
            ls = [builder.lshr(_, self._index_type(1)) for _ in ls]
        return ls

    def _llvm_generate_coset_scale_fcc(self, builder, coset_idx, ls):

        is_floating_point = ls[0].type == self._internal_type
        y = builder.and_(coset_idx, self._index_type(1))
        z = builder.and_(coset_idx, self._index_type(2))
        z = builder.lshr(z, self._index_type(1))
        x = builder.xor(z, y)

        if is_floating_point:
            x = builder.fsub(ls[0], builder.sitofp(x, self._internal_type))
            y = builder.fsub(ls[1], builder.sitofp(y, self._internal_type))
            z = builder.fsub(ls[2], builder.sitofp(z, self._internal_type))

            x = builder.fmul(x, self._internal_type(0.5))
            y = builder.fmul(y, self._internal_type(0.5))
            z = builder.fmul(z, self._internal_type(0.5))

        else:
            x = builder.sub(ls[0], x)
            y = builder.sub(ls[1], y)
            z = builder.sub(ls[2], z)

            x = builder.lshr(x, self._index_type(1))
            y = builder.lshr(y, self._index_type(1))
            z = builder.lshr(z, self._index_type(1))

        return [x,y,z]

    def _llvm_generate_coset_scale_generic(self, builder, coset_idx, ls):
        """
        This takes the coset index, that is, which coset a point belongs to
        and shifts the input lattice site by that coset offset, and scales it
        to be unit cartesian. 

        REVIEW
        """
        d_matrix, v = self._ss.coset_vectors()

        if len(v) == 1:  # Special case, we're on a scaled cartesian grid
            return [builder.fdiv(k, self._internal_type(s)) for k, s in zip(ls, d_matrix.diagonal())]

        # Otherwise, we need to determine the shift, first        
        shift = [self._internal_type(_) for _ in v]

        for idx, coset_shift in v[1:]:
            comp = builder.icmp_unsigned("==", coset_idx, self._index_type(idx+1))
            shift = [builder.select(comp, self._internal_type(new), old) for new, old in zip(coset_shift, shift)]
        
        ls = [builder.fsub(k, s, flags=['fast']) for k, s in zip(ls, shift)]
        return [builder.fdiv(k, self._internal_type(s)) for k, s in zip(ls, d_matrix.diagonal())]

    def _llvm_generate_coset_scale(self, builder, coset_idx, ls):
        if self._dimension == 3:
            # BCC lattice indicator
            if self._ss.get_lattice_hash() == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                return self._llvm_generate_coset_scale_bcc(builder, coset_idx, ls)

            if self._ss.get_lattice_hash() == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                return self._llvm_generate_coset_scale_fcc(builder, coset_idx, ls)

            if self._ss.get_lattice_hash() == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                return ls

        if self._ss.is_lattice_cartesian():
            return ls
        
        raise ValueError("Unsupported lattice -- no coset decomposition")

    def _ls_coset_index_bcc(self, ls):
        return int(ls[2]) % 2

    def _ls_coset_index_fcc(self, ls):
        x = (ls[0] % 2) << 1
        y = ls[1] % 2
        return x | y

    def _ls_coset_index(self, ls):
        if self._dimension == 3:
            # BCC lattice indicator
            if self._ss.get_lattice_hash() == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                return self._ls_coset_index_bcc(ls)

            if self._ss.get_lattice_hash() == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                return self._ls_coset_index_fcc(ls)

            if self._ss.get_lattice_hash() == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                return int(0)

        if self._ss.is_lattice_cartesian():
            return int(0)
        
        raise ValueError("Unsupported lattice -- no coset decomposition")

    def lookup_lattice_shift(self, index, builder):
        subregions = self._ss.get_subregions()

        if len(subregions) <= 1:
            ls = subregions[0].get_optimal_nn_lookups()[index]
            return [self._index_type(_) for _ in ls]

        region_index = self.read_mem(self._stack_bspindex, builder)
        site_index = self._index_type(int(index))

        if self.lane_count == 1:
            ls = [builder.gep(self.ls_index, [self._index_type(0), region_index, site_index, self._index_type(_)])
                  for _ in range(self._dimension)]
            ls = [builder.load(_, name='load_site') for _ in ls]
        else:
            ptr_type = ir.VectorType(ir.IntType(self.ptr_size), self.lane_count)
            bcst = builder.insert_element(ptr_type(ir.Undefined),
                                          self.ls_index_u.ptrtoint(ir.IntType(self.ptr_size)),
                                          ir.IntType(32)(0), name='ls_index_u')
            bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), self._index_type(None))
            region_lsindex = builder.inttoptr(bcst, ir.VectorType(self.ls_struct_type.as_pointer(), self.lane_count))

            ls = [builder.gep(region_lsindex, [region_index, ir.IntType(32)(int(index * self._dimension + _))] )
                  for _ in range(self._dimension)]

            ls = [self.read_mem(_, builder) for _ in ls]
        return ls


    def get_lattice_site(self, index, coset):
        """
        This gets the lattice site at index 'index'. There's no need to worry about
        the coset, nor the

        :param index:
        :return:
        """

        # The most basic case is if we have one reference sub_region
        if len(self._ss._ref_subregions) == 1:
            # If we only have one sub region, just encode the
            # lattice site fetches in code
            return None

        # Otherwise we have more than one
        if self._subregion_consistency:
            # load polyindex
            # use polyindex to pick out where the ls read should come from
            pass

    def optimize_polynomial(self, site_polynomials):
        c_lookups = []
        c_to_pipe = {}
        c_polynomial = 0

        for cidx, (polynomial, pipeline_var) in enumerate(site_polynomials):
            # Generate a memfetch at ls, and store it in the c_lookups
            mf = var("c_%d" % cidx)
            c_to_pipe[mf] = pipeline_var
            c_lookups += [mf]
            c_polynomial += mf * polynomial

        # Factorize
        return horner_factor(c_polynomial, self._ss._s, c_lookups), c_lookups, c_to_pipe


    def emit_horner(self, builder, xvars, polyindex, hregs, c_fetches, result_current):
        """
        This little chunk of code gets reused in a bunch
        of sub-cases, so pull it out here to make it easier
        to change
        """
        xvars = list(xvars)
        c_fetches = list(c_fetches)

        # If we want to avoid branchy code (or we have only 1 case)
        # then we'll use branch predication
        if self.branchless or self.branch_predication or len(self._ss._ref_subregions) == 1:
            result = result_current
            for refregion_idx, __ in enumerate(self._ss._ref_subregions):

                # Generate the result
                chunk_result = self._llvm_generate_horner(builder, hregs[refregion_idx], xvars + c_fetches, None)
                _ = chunk_result

                if len(self._ss._ref_subregions) > 1:
                    # select result if polyindex == refregion_idx else zero
                    # add result to previous result
                    _ = builder.icmp_unsigned("==", polyindex, ir.IntType(32)(refregion_idx))
                    _ = builder.select(_, chunk_result, self._internal_type(0.0))
                result = builder.fadd(_, result, flags=['fast'])
            return result

        else:
            assert self.lane_count == 1
            # Get a reference to this bb
            bb = builder.basic_block

            # First, create a bunch of basic blocks for each case
            bblist = [builder.append_basic_block(name=_label_suffix(bb.name, '.case%d' % idx))
                      for idx, _ in enumerate(self._ss._ref_subregions)]

            # This last bb is where all the above join back together
            bbend = builder.append_basic_block(name=_label_suffix(bb.name, '.endcase'))

            # Ugh, polyindex might be a low bit type, cast to
            # int32 type --- this should be sufficient
            # I can't imagine splines with 2^32 sub-regions to be
            # practical
            pidx = builder.zext(polyindex, ir.IntType(int(32)))

            # Create the switch ... the default is the 0'th block
            c = builder.switch(pidx, bblist[0])

            # Add the different cases (effectively: add jumps to different basic blocks)
            for refregion_idx, __ in enumerate(self._ss._ref_subregions[1:]):
                c.add_case(ir.IntType(int(32))(refregion_idx + 1), bblist[refregion_idx + 1])

            # Generate code for the different blocks
            outputs = []
            for idx, block in enumerate(bblist):
                with builder.goto_block(block):
                    chunk_result = self._llvm_generate_horner(builder, hregs[idx], xvars + c_fetches, None)
                    outputs += [chunk_result]

                    # Branch back to a common point
                    builder.branch(bbend)

            # Resume serial branchless code
            builder.position_at_start(bbend)

            # Add a phi node to choose the correct value from the
            # above blocks
            phi = builder.phi(self._internal_type)
            for bb, val in zip(bblist, outputs):
                phi.add_incoming(val, bb)

            # Add the final result to the convolution sum
            result = builder.fadd(phi, result_current, flags=['fast'])
            return result

    def _llvm_read_lattice_memory(self, builder, ref_lattice_site, ls_idx, ll_coset_index, coset_index):
        rt_coset_index = self._ls_coset_index(ref_lattice_site)
        ls = None
        res = None
        if self.lookup_type == "ABSTRACT":
            ls = self.lookup_lattice_shift(ls_idx, builder)

            if self._ss._reflective:
                m = [builder.load(self._stack_xflip[d]) for d in range(self._dimension)]
                m = [builder.fptosi(_) for _ in m]
                ls = [builder.mul(val, m[d]) for d, val in enumerate(ls)]

            if self.lane_count == 1:
                res = builder.call(self.memfetch, ls)
            else:
                res = self._internal_type(0.)
                for lane in range(self.lane_count):
                    lsl = []
                    for d in range(self._dimension):
                        lsl += [builder.extract_element(ls[d], ir.IntType(int(32))(lane))]
                    _ = builder.call(self.memfetch[0], lsl)
                    res = builder.insert_element(res, _, ir.IntType(int(32))(lane))

        elif self.lookup_type in ["COSET_RAW", "COSET", "TEX"]:
            if (ll_coset_index is None or rt_coset_index != coset_index) or not self._subregion_consistency:
                ls = self.lookup_lattice_shift(ls_idx, builder)

                for i, s in enumerate(ls):
                    value = builder.load(self._stack_rho[i])
                    if self._ss._reflective:
                        m = builder.load(self._stack_xflip[i])
                        m = builder.fptosi(m, self._index_type)
                        s = builder.mul(m, s)
                    ls[i] = builder.add(s, value)
                ll_coset_index = self._llvm_generate_coset_index(builder, ls)
                coset_index = rt_coset_index

            if self.lookup_type in ["COSET", "TEX"]:

                # Load the lattice site
                if ls is None:
                    ls = self.lookup_lattice_shift(ls_idx, builder)
                    if self._ss._reflective:
                        m = [builder.load(self._stack_xflip[d]) for d in range(self._dimension)]
                        m = [builder.fptosi(_, self._index_type) for _ in m]
                        ls = [builder.mul(val, m[d]) for d, val in enumerate(ls)]

                    ls = [builder.add(val, builder.load(self._stack_rho[d])) for d, val in enumerate(ls)]

                # Then turn it into a cartesian site on ll_coset_index
                cart_addr = self._llvm_generate_coset_scale(builder, ll_coset_index, ls)

                if self.lookup_type == "TEX":
                    assert self.lane_count == 1

                    # Choose the texture to lookup from
                    texture = self._textures[0]
                    for idx_tex, tex in enumerate(self._textures[1:]):
                        p = builder.icmp_unsigned("==", ll_coset_index, ir.IntType(32)(idx_tex + 1))
                        texture = builder.select(p, tex, texture)

                    # Shift
                    pt = [builder.sitofp(_, self._internal_type) for _ in cart_addr]
                    pt = [builder.fadd(_,  self._internal_type(0.5), flags=['fast']) for _ in pt]
                    pt = [builder.fptrunc(_, ir.FloatType()) for _ in pt]

                    # if self.strict_branchless:
                        # Fetch
                    if self._ss._s == 1:
                        res = builder.call(self._tex1d, [texture]+pt)
                    elif self._ss._s == 2:
                        res = builder.call(self._tex2d, [texture]+pt)
                    elif self._ss._s == 3:
                        res = builder.call(self._tex3d, [texture]+pt)
                    else:
                        raise ValueError("Dimension not supported")

                    res = builder.extract_value(res, 0)

                elif self.lookup_type == "COSET":
                    raise ValueError("no time to debug rn")

            elif self.lookup_type == "COSET_RAW":
                buffer = builder.load(self._stack_buffers[0])

                # Choose the buffer that corresponds to the coset we're currently on
                for _ in range(1, self._ss.coset_count()):
                    tmp = builder.load(self._stack_buffers[_])
                    sel = builder.icmp_unsigned("==", ll_coset_index, self._index_type(_))
                    buffer = builder.select(sel, tmp, buffer)

                if self._precalc_fetch_offset:
                    raise ValueError("Optimization not supported yet")
                else:
                    if ls is None:
                        ls = self.lookup_lattice_shift(ls_idx, builder)
                        if self._ss._reflective:
                            m = [builder.load(self._stack_xflip[d]) for d in range(self._dimension)]
                            m = [builder.fptosi(_, self._index_type) for _ in m]
                            ls = [builder.mul(val, m[d]) for d, val in enumerate(ls)]

                        ls = [builder.add(val, builder.load(self._stack_rho[d])) for d, val in enumerate(ls)]

                    cart_addr = self._llvm_generate_coset_scale(builder, ll_coset_index, ls)

                    # Compute the linear index for the lattice site
                    # within the coset structure
                    lindex = cart_addr[-1]
                    for d in reversed(range(1, self._dimension)):
                        # Choose the correct boundary for the coset
                        bdry = builder.load(self._stack_bounds[0][d])
                        for _ in range(1, self._ss.coset_count()):
                            tmp = builder.load(self._stack_bounds[_][d])
                            sel = builder.icmp_unsigned("==", ll_coset_index, self._index_type(_))
                            bdry = builder.select(sel, tmp, bdry)

                        lindex = builder.mul(bdry, lindex)
                        lindex = builder.add(cart_addr[d - 1], lindex)

                    if self.test_stop == "first_lindex":
                        lindex = builder.sitofp(lindex, self._internal_type)
                        return lindex, ll_coset_index, coset_index

                    faddr = builder.gep(buffer, [lindex])
                    res = self.read_mem(faddr, builder)
            else:
                raise ValueError("Not re-implemented yet")
        return res, ll_coset_index, coset_index

    def _llvm_generate_no_lerp_conv_sum_bpred(self, builder, outer_coset_index=None, outer_coset_shifts=None):

        global polynomial
        result = self._internal_type(0.0)
        builder.store(result, self._stack_coset_result)

        ref_regions = self._ss.get_ref_subregions()
        rref_lookups = [_.get_optimal_nn_lookups() for _ in ref_regions]
        r0_lookups = rref_lookups[0]

        # Fill the pipeline up with empty slots
        pipeline = [None] * self.pipeline_depth
        pipeline_slot = 0

        # This is the max number of terms in the unwrapped conv sum
        weights_len = max([len(_.get_optimal_nn_lookups()) for _ in self._ss.get_ref_subregions()])
        polyplan = {}
        polyqueue = []
        idx = 0

        self._precalc_fetch_offset = False
        ll_current_coset, current_coset = (None, None)
        while idx < weights_len or len(polyqueue) > 0:

            # Try to fill up the pipeline
            while any([_ is None for _ in pipeline]) and idx < weights_len:
                ref_ls = r0_lookups[idx]

                # Generate a memory fetch
                res, ll_current_coset, current_coset = self._llvm_read_lattice_memory(
                    builder,
                    ref_ls,
                    idx,
                    ll_current_coset,
                    current_coset)

                # Store the mem fetch in the pipeline
                if self.test_stop == "first_fetch_mem" or self.test_stop == "first_lindex":
                    builder.store(res, self._stack_coset_result)
                    return None
                builder.store(res, self._stack_pipeline[pipeline_slot])

                # Map that pipeline slot to the correct polynomial
                pipeline[pipeline_slot] = (self._stack_pipeline[pipeline_slot], idx, res)
                polyplan[idx] = pipeline_slot
                polyqueue += [idx]

                # Increment the pipeline
                pipeline_slot = (pipeline_slot + 1) % self.pipeline_depth
                idx += 1

            # Evaluate the polynomial with data from the pipeline
            polynomial = [[] for _ in ref_regions]
            eval_result = self._internal_type(0.)

            for _ in range(self.group_size):
                if len(polyqueue) == 0:
                    break
                poly_idx = polyqueue.pop(0)
                pipeline_slot = polyplan[poly_idx]

                for ref_region_idx,_ in enumerate(ref_regions):
                    ls = rref_lookups[ref_region_idx][poly_idx]
                    poly = ref_regions[ref_region_idx].get_polynomial(ls)
                    polynomial[ref_region_idx] += [(poly, self._stack_pipeline[pipeline_slot])]

                if self.test_stop == "coeff_sum":
                    eval_result = builder.fadd(eval_result, builder.load(self._stack_pipeline[pipeline_slot]))

                # At the end of all the evaluation, clear out the pipeline element
                pipeline[pipeline_slot] = None

            if self.test_stop != "coeff_sum":
                if len(ref_regions) == 1 or self.branch_predication:
                    results = []
                    for poly_i in polynomial:
                        poly_eval, hregs, pl = self.optimize_polynomial(poly_i)
                        eval_result = self._llvm_generate_horner(builder, poly_eval, hregs, pl)
                        results += [eval_result]

                    eval_result = results[0]
                    if len(ref_regions) > 1:
                        poly_index = builder.load(self._stack_polyindex)
                        for ridx, result in enumerate(results[1:]):
                            ridx = ridx + 1
                            cmp = builder.icmp_unsigned("==", poly_index, self._index_type(ridx))
                            eval_result = builder.select(cmp, result, eval_result)
                else:
                    assert self.lane_count == 1
                    # Get a reference to this bb
                    bb = builder.basic_block

                    # First, create a bunch of basic blocks for each case
                    bb_list = [builder.append_basic_block(name=_label_suffix(bb.name, '.case%d' % idx))
                              for idx, _ in enumerate(self._ss._ref_subregions)]

                    # This last bb is where all the above join back together
                    bb_end = builder.append_basic_block(name=_label_suffix(bb.name, '.endcase'))

                    poly_index = builder.load(self._stack_polyindex)
                    # Create the switch ... the default is the 0'th block
                    c = builder.switch(poly_index, bb_list[0])

                    # Add the different cases (effectively: add jumps to different basic blocks)
                    for refregion_idx, __ in enumerate(self._ss._ref_subregions[1:]):
                        c.add_case(ir.IntType(int(32))(refregion_idx + 1), bb_list[refregion_idx + 1])

                    # Generate code for the different blocks
                    outputs = []
                    for idx_inner, block in enumerate(bb_list):
                        with builder.goto_block(block):
                            poly_eval, hregs, pl = self.optimize_polynomial(polynomial[idx_inner])
                            eval_result = self._llvm_generate_horner(builder, poly_eval, hregs, pl)
                            outputs += [eval_result]
                            # Branch back to a common point
                            builder.branch(bb_end)

                    # Resume serial branchless code
                    builder.position_at_start(bb_end)

                    # Add a phi node to choose the correct value from the
                    # above blocks
                    eval_result = builder.phi(self._internal_type)
                    for bb, val in zip(bb_list, outputs):
                        eval_result.add_incoming(val, bb)

            r = builder.load(self._stack_coset_result)
            r = builder.fadd(r, eval_result, flags=['fast'])
            builder.store(r, self._stack_coset_result)
        return None

    def _llvm_generate_lerp_conv_sum_bpred_consistent(self, builder, xvars, ivars, xform, shift, srindex, polyindex, ls_lut, memfetch):
        # This should be handled below, I think...
        pass

    def _llvm_generate_lerp_conv_sum_bpred(self, builder, outer_coset_index=None):

        # Clear out the result
        result = self._internal_type(0.0)
        builder.store(result, self._stack_coset_result)

        in_flight_q = []
        G, coset_vectors = self._ss.coset_vectors()
        current_coset = (None, None, None)
        debug_cnt = {}

        # For every fetch
        for fetch_index, rrchunks in enumerate(zip_longest(*[_.get_optimal_lin_lookups(True) for _ in self._ss._ref_subregions])):
            coset_var = None

            # Iterate over all of the reference regions
            for ref_idx, grp in enumerate(rrchunks):

                if ref_idx not in debug_cnt:
                    debug_cnt[ref_idx] = 0

                if grp is not None:
                    debug_cnt[ref_idx] += len(grp)

                if self._ss._s == 1:
                    raise ValueError("Dimension=1 not supported yet")
                elif self._ss._s == 2:
                    raise ValueError("Dimension=2 not supported yet")
                elif self._ss._s == 3:
                    if grp is None:
                        pt = [self._internal_type(0) for _ in range(self._ss._s)]
                        rp = pt[:]
                        g = self._internal_type(0.0)

                    # Generate the parameters for the lookup
                    elif len(grp) == 1:

                        # Simulate a nearest neighbor fetch
                        A = grp[0]
                        pt = [self._internal_type(_) for _ in A]
                        rp = pt[:]

                        h = self._ss._ref_subregions[ref_idx]._factored_terms[A]
                        g = self._llvm_generate_horner(builder, h, [], {})

                    elif len(grp) == 2:
                        A, B = grp
                        rp = [self._internal_type(float(_)) for _ in A]

                        hw0 = self._ss._ref_subregions[ref_idx]._factored_terms[A]
                        hw1 = self._ss._ref_subregions[ref_idx]._factored_terms[B]

                        w0 = self._llvm_generate_horner(builder, hw0, [], {})
                        w1 = self._llvm_generate_horner(builder, hw1, [], {})

                        g = builder.fadd(w0, w1, flags=['fast'])
                        t = builder.fdiv(w1, g)

                        terr = builder.call(self._fabs, [g])
                        c = builder.fcmp_unordered("<=", g, self._internal_type(0.0000000001))
                        t = builder.select(c, self._internal_type(0), t)

                        pt = [
                            builder.fadd(
                                builder.fmul(t, self._internal_type(b - a)),
                                self._internal_type(a)) for (a, b) in zip(A, B)
                        ]
                        # g = self._internal_type(2.0)

                    elif len(grp) == 4:

                        A, B, C, D = grp
                        rp = [self._internal_type(float(_)) for _ in A]

                        hw0 = self._ss._ref_subregions[ref_idx]._factored_terms[A]
                        hw1 = self._ss._ref_subregions[ref_idx]._factored_terms[B]
                        hw2 = self._ss._ref_subregions[ref_idx]._factored_terms[C]
                        hw3 = self._ss._ref_subregions[ref_idx]._factored_terms[D]

                        w0 = self._llvm_generate_horner(builder, hw0, [], {})
                        w1 = self._llvm_generate_horner(builder, hw1, [], {})
                        w2 = self._llvm_generate_horner(builder, hw2, [], {})
                        w3 = self._llvm_generate_horner(builder, hw3, [], {})

                        a = builder.fadd(w0, w3, flags=['fast'])
                        b = builder.fadd(w1, w2, flags=['fast'])

                        g = builder.fadd(a, b, flags=['fast'])
                        tx = builder.fdiv(b, g)

                        a = builder.fadd(w1, w3, flags=['fast'])
                        ty = builder.fdiv(a, g)

                        terr = builder.call(self._fabs, [g])
                        c = builder.fcmp_unordered("<=", g, self._internal_type(0.0000000001))

                        # Select the
                        tx = builder.select(c, self._internal_type(0), tx)
                        ty = builder.select(c, self._internal_type(0), ty)

                        pt = [self._internal_type(0) for _ in range(self._ss._s)]

                        d0 = vector(C) - vector(A)
                        d1 = vector(D) - vector(A)
                        d0 = [i for i, x in enumerate(d0) if x != 0][0]
                        d1 = [i for i, x in enumerate(d1) if x != 0][0]

                        pt[d0] = tx
                        pt[d1] = ty

                        pt = [builder.fmul(_, self._internal_type(G[i, i])) for i, _ in enumerate(pt)]
                        pt = [builder.fadd(a, b, flags=['fast']) for a, b in zip(pt, rp)]


                    elif len(grp) == 8:
                        rp = [self._internal_type(float(_)) for _ in grp[0]]

                        hw0 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[0]]
                        hw1 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[4]]
                        hw2 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[2]]
                        hw3 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[6]]
                        hw4 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[1]]
                        hw5 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[5]]
                        hw6 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[3]]
                        hw7 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[7]]

                        w0 = self._llvm_generate_horner(builder, hw0, [], {}, name='w0')
                        w1 = self._llvm_generate_horner(builder, hw1, [], {}, name='w1')
                        w2 = self._llvm_generate_horner(builder, hw2, [], {}, name='w2')
                        w3 = self._llvm_generate_horner(builder, hw3, [], {}, name='w3')
                        w4 = self._llvm_generate_horner(builder, hw4, [], {}, name='w4')
                        w5 = self._llvm_generate_horner(builder, hw5, [], {}, name='w5')
                        w6 = self._llvm_generate_horner(builder, hw6, [], {}, name='w6')
                        w7 = self._llvm_generate_horner(builder, hw7, [], {}, name='w7')

                        # g = w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7
                        # hz = (w4 + w5 + w6 + w7)
                        # hy = (w2 + w3 + w6 + w7)
                        # hx = (w1 + w3 + w5 + w7)

                        a = builder.fadd(w0, w1, flags=['fast'])
                        b = builder.fadd(w2, w3, flags=['fast'])
                        c = builder.fadd(w4, w5, flags=['fast'])
                        d = builder.fadd(w6, w7, flags=['fast'])

                        _a = builder.fadd(a, b, flags=['fast'])
                        _b = builder.fadd(c, d, flags=['fast'])

                        _c = builder.fadd(w1, w3, flags=['fast'])
                        _d = builder.fadd(w5, w7, flags=['fast'])

                        g = builder.fadd(_a, _b, flags=['fast'], name='g_out')

                        hz = builder.fadd(c, d, flags=['fast'])
                        hy = builder.fadd(b, d, flags=['fast'])
                        hx = builder.fadd(_c, _d, flags=['fast'])

                        tz = builder.fdiv(hz, g)
                        ty = builder.fdiv(hy, g)
                        tx = builder.fdiv(hx, g)

                        terr = builder.call(self._fabs, [g])
                        c = builder.fcmp_unordered("<=", g, self._internal_type(0.0000000001))

                        tz = builder.select(c, self._internal_type(0), tz)
                        ty = builder.select(c, self._internal_type(0), ty)
                        tx = builder.select(c, self._internal_type(0), tx)

                        # These are our scales
                        pt = [tx, ty, tz]

                        # Scale them for the coset lattice

                        pt = [builder.fmul(_, self._internal_type(G[i, i])) for i, _ in enumerate(pt)]
                        pt = [builder.fadd(a, b, flags=['fast']) for a, b in zip(pt, rp)]
                        # pt[0] = w0
                        # pt[1] = w1
                        # pt[2] = w2

                    else:
                        raise ValueError("Cannot group into higher than 8!")

                    if ref_idx == 0:
                        builder.store(g, self._stack_gval)

                        for _ in range(self._dimension):
                            builder.store(pt[_], self._stack_g_pt[_])
                            builder.store(rp[_], self._stack_rp_pt[_])

                    else:
                        # Load the current polyindex
                        polyindex = builder.load(self._stack_polyindex)
                        c = builder.icmp_unsigned("==", polyindex, ir.IntType(32)(ref_idx))

                        # predicate G(x)
                        g_old = builder.load(self._stack_gval)
                        g_pred = builder.select(c, g, g_old)
                        builder.store(g_pred, self._stack_gval)

                        for _ in range(self._dimension):
                            old_pt = builder.load(self._stack_g_pt[_])
                            pt_pred = builder.select(c, pt[_], old_pt)
                            builder.store(pt_pred, self._stack_g_pt[_])

                            old_rp = builder.load(self._stack_rp_pt[_])
                            rp_pred = builder.select(c, rp[_], old_rp)
                            builder.store(rp_pred, self._stack_rp_pt[_])

            # Transform the tri-lin lookup point
            if len(self._ss.get_subregions()) > 1:
                self._llvm_apply_transform(self._stack_g_pt, builder, target=self._stack_g_pt)
                self._llvm_apply_shift(builder, self._stack_g_pt, self._stack_tshift, target=self._stack_g_pt)

                if coset_var is None:
                    self._llvm_apply_transform(self._stack_rp_pt, builder, target=self._stack_rp_pt)
                    self._llvm_apply_shift(builder, self._stack_rp_pt, self._stack_tshift, target=self._stack_rp_pt)

            # If we have reflective symmetry, we need to apply that now
            if self._ss._reflective:
                for d in range(self._dimension):
                    m = builder.load(self._stack_xflip[d])
                    pt = builder.load(self._stack_g_pt[d])
                    builder.store(builder.fmul(pt, m), self._stack_g_pt[d])
                    if coset_var is None:
                        rp = builder.load(self._stack_rp_pt[d])
                        builder.store(builder.fmul(rp, m), self._stack_rp_pt[d])

            # self._llvm_apply_shift(builder, self._stack_g_pt, self._stack_rho, target=self._stack_g_pt)
            for d in range(self._dimension):
                f = builder.load(self._stack_g_pt[d])
                ro = builder.load(self._stack_rho[d])
                fro = builder.sitofp(ro, self._internal_type)
                ro = builder.fadd(fro, f)

                builder.store(ro, self._stack_g_pt[d])

            if coset_var is None:
                ls = []
                for d in range(self._dimension):
                    f = builder.load(self._stack_rp_pt[d])
                    i = builder.fptosi(f, ir.IntType(32))
                    ro = builder.load(self._stack_rho[d])

                    ls += [builder.add(i, ro)]

                # self._llvm_apply_shift(builder, self._stack_rp_pt, self._stack_rho, target=self._stack_rp_pt)
                # rp = [builder.fptosi(builder.load(_), ) for _ in self._stack_rp_pt]
                coset_var = self._llvm_generate_coset_index(builder, ls)

            if self.lookup_type == "TEX_LIN":
                _ = self._textures[0]

                if outer_coset_index is not None:
                    n_ocosets = len(self._ss.get_L_cosets())
                    offset = builder.mul(outer_coset_index, ir.IntType(32)(n_ocosets))
                    coset = builder.add(offset, coset_var)
                else:
                    coset = coset_var

                # Given the coset, select the correct texture.
                texture = self._textures[0]
                for idx_tex, tex in enumerate(self._textures[1:]):
                    p = builder.icmp_unsigned("==", coset, ir.IntType(32)(idx_tex + 1))
                    texture = builder.select(p, tex, texture)

                # Shift
                pt = [builder.load(_) for _ in self._stack_g_pt]
                int_pt = [builder.fptosi(_, ir.IntType(32)) for _ in pt]
                pt = self._llvm_generate_coset_scale(builder, coset_var, pt)
                pt = [builder.fadd(_,  self._internal_type(0.5), flags=['fast']) for _ in pt]
                pt = [builder.fptrunc(_, ir.FloatType()) for _ in pt]
                for d in range(self._dimension):
                    builder.store(pt[d], self._stack_g_pt[d])

                if self.strict_branchless:
                    # Fetch
                    if self._ss._s == 1:
                        mf = builder.call(self._tex1d, [texture]+pt)
                    elif self._ss._s == 2:
                        mf = builder.call(self._tex2d, [texture]+pt)
                    elif self._ss._s == 3:
                        mf = builder.call(self._tex3d, [texture]+pt)
                    else:
                        raise ValueError("Dimension not supported")

                    if self.pipeline_depth >= 1:
                        if len(in_flight_q) >= self.pipeline_depth:
                            # We can't queue another read, so we consume a read before
                            # we queue
                            fetch, gf = in_flight_q.pop(0)
                            fetch = builder.extract_value(fetch, 0)
                            prd = builder.fmul(fetch,  gf, flags=['fast'])
                            result = builder.fadd(prd,  result, flags=['fast'])

                        # Queue up the last
                        g = builder.load(self._stack_gval)
                        in_flight_q += [(mf, g)]
                    else:
                        mf = builder.extract_value(mf, 0)
                        g = builder.load(self._stack_gval)
                        _ = builder.fmul(mf,  g, flags=['fast'])
                        result = builder.fadd(_,  result, flags=['fast'])
                else:
                    if self.pipeline_depth > 1:
                        raise ValueError("Invalid combo")

                    g = builder.load(self._stack_gval)
                    pred = builder.fcmp_unordered("!=", g, self._internal_type(0.0))
                    c = self._internal_type(0.0)

                    bb = builder.basic_block
                    bbif = builder.append_basic_block(name=_label_suffix(bb.name, '.if'))
                    bbend = builder.append_basic_block(name=_label_suffix(bb.name, '.endif'))
                    br = builder.cbranch(pred, bbif, bbend)

                    with builder.goto_block(bbif):
                        if self._ss._s == 1:
                            output = builder.call(self._tex1d, [texture]+pt)
                        elif self._ss._s == 2:
                            output = builder.call(self._tex2d, [texture]+pt)
                        elif self._ss._s == 3:
                            output = builder.call(self._tex3d, [texture] + pt)
                        else:
                            raise ValueError("Dimension not supported")
                        output = builder.extract_value(output, 0)
                        builder.branch(bbend)

                    builder.position_at_start(bbend)

                    mf = builder.phi(self._internal_type)
                    mf.add_incoming(c, bb)
                    mf.add_incoming(output, bbif)

                    _ = builder.fmul(mf,  g, flags=['fast'])
                    result = builder.fadd(_,  result, flags=['fast'])

        # We need to take care of any left over fetches
        for fetch, g in in_flight_q:
            if self.lookup_type == "TEX_LIN":
                fetch = builder.extract_value(fetch, 0)
            prd = builder.fmul(fetch,  g, flags=['fast'])
            result = builder.fadd(prd,  result, flags=['fast'])

        r = builder.load(self._stack_coset_result)
        r = builder.fadd(r, result, flags=['fast'])
        builder.store(r, self._stack_coset_result)
        return None

    def _llvm_generate_lerp_conv_sum_bpred_inconsistent(self,
                                                        builder,
                                                        xvars,
                                                        ivars,
                                                        xform,
                                                        shift,
                                                        srindex,
                                                        polyindex,
                                                        ls_lut,
                                                        memfetch,
                                                        outer_coset_index=None,
                                                        outer_coset_shifts=None):
        """
        This is a fairly lengthy method, let me break down the method call name
        _llvm_generate_(no_lerp)_conv_sum_(bpred)
        no_lerp => No linear interpolation
        bpred => uses branch predication
        basically breaks down the convulatuon sup into chunks (big outer loop)
        """
        # Need to think about how the coset lookup works here....
        # We're assured that, if the _subregion_consistency is set, 
        # then we know that we only need to generate a coset identification
        # when any lattice site breaks across the coset boundary
        ls_start = 0
        result = self._internal_type(0.0)
        reg0 = self._ss._ref_subregions[0]

        in_flight_q = []

        G, coset_vectors = self._ss.coset_vectors()

        current_coset = (None, None, None)

        
        debug_cnt = {}

        for fetch_index, rrchunks in enumerate(zip_longest(*[_.get_optimal_lin_lookups(True) for _ in self._ss._ref_subregions])):

            # list_chunk = zip(*list_chunk)
            coset_var = None

            pt_pred = None 
            rp_pred = None
            g_pred  = None

            # If we need to, generate a load for the transform registers,
            # then do it now
            if self.refetch_transform:
                xform = self._llvm_load_transform_matrix(builder, srindex, False)
                if not self.no_shift:
                    shift = self._llvm_load_shift(builder, srindex)

            for ref_idx, grp in enumerate(rrchunks):

                if ref_idx not in debug_cnt:
                    debug_cnt[ref_idx] = 0
    
                if grp is not None:
                    debug_cnt[ref_idx] += len(grp)

                if self._ss._s == 1:
                    raise ValueError("Dimension=1 not supported yet")
                elif self._ss._s == 2:
                    raise ValueError("Dimension=2 not supported yet")
                elif self._ss._s == 3:

                    if grp is None:
                        pt = [self._internal_type(0) for _ in range(self._ss._s)]
                        rp = pt[:]
                        g = self._internal_type(0.0)

                    # Generate the parameters for the lookup
                    elif len(grp) == 1:

                        # Simulate a nearest neighbor fetch
                        A = grp[0]
                        pt = [self._internal_type(_) for _ in A]
                        rp = pt[:]

                        h = self._ss._ref_subregions[ref_idx]._factored_terms[A]
                        g = self._llvm_generate_horner(builder, h, xvars,  None)

                    elif len(grp) == 2:
                        A,B = grp
                        rp = [self._internal_type(float(_)) for _ in A]

                        hw0 = self._ss._ref_subregions[ref_idx]._factored_terms[A]
                        hw1 = self._ss._ref_subregions[ref_idx]._factored_terms[B]

                        w0 = self._llvm_generate_horner(builder, hw0, xvars, None)
                        w1 = self._llvm_generate_horner(builder, hw1, xvars, None)

                        g = builder.fadd(w0,  w1, flags=['fast'])
                        t = builder.fdiv(w1, g)

                        terr = builder.call(self._fabs, [g])
                        c    = builder.fcmp_unordered("<=", g, self._internal_type(0.0000000001))
                        t    = builder.select(c, self._internal_type(0), t)

                        pt = [
                            builder.fadd(
                                builder.fmul(t, self._internal_type(b-a)), 
                                self._internal_type(a)) for (a,b) in zip(A,B)
                            ] 
                        # g = self._internal_type(2.0)

                    elif len(grp) == 4:
                                                
                        A,B,C,D = grp
                        
                        hw0 = self._ss._ref_subregions[ref_idx]._factored_terms[A]
                        hw1 = self._ss._ref_subregions[ref_idx]._factored_terms[B]
                        hw2 = self._ss._ref_subregions[ref_idx]._factored_terms[C]
                        hw3 = self._ss._ref_subregions[ref_idx]._factored_terms[D]

                        w0 = self._llvm_generate_horner(builder, hw0, xvars, None)
                        w1 = self._llvm_generate_horner(builder, hw1, xvars, None)
                        w2 = self._llvm_generate_horner(builder, hw2, xvars, None)
                        w3 = self._llvm_generate_horner(builder, hw3, xvars, None)
                        
                        a = builder.fadd(w0,  w3, flags=['fast'])
                        b = builder.fadd(w1,  w2, flags=['fast'])

                        g  = builder.fadd(a, b, flags=['fast'])
                        tx = builder.fdiv(b,g)

                        a = builder.fadd(w1,  w3, flags=['fast'])
                        ty = builder.fdiv(a, g)

                        terr = builder.call(self._fabs, [g])
                        c    = builder.fcmp_unordered("<=", g, self._internal_type(0.0000000001))

                        # Select the 
                        tx  = builder.select(c, self._internal_type(0), tx)
                        ty  = builder.select(c, self._internal_type(0), ty)


                        pt = [self._internal_type(0) for _ in range(self._ss._s)]
                        
                        d0 = vector(C) - vector(A)
                        d1 = vector(D) - vector(A)
                        d0 = [i for i,x  in enumerate(d0) if x != 0][0]        
                        d1 = [i for i,x  in enumerate(d1) if x != 0][0]
                        
                        pt[d0] = tx
                        pt[d1] = ty

                        pt = [builder.fmul(_, self._internal_type(G[i,i]) ) for i,_ in enumerate(pt)]
                        pt = [builder.fadd(a,  b, flags=['fast']) for a,b in zip(pt, rp)]

    
                    elif len(grp) == 8:
                        rp = [self._internal_type(float(_)) for _ in grp[0]]        

                        hw0 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[0]]
                        hw1 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[4]]
                        hw2 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[2]]
                        hw3 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[6]]
                        hw4 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[1]]
                        hw5 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[5]]
                        hw6 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[3]]
                        hw7 = self._ss._ref_subregions[ref_idx]._factored_terms[grp[7]]

                        w0 = self._llvm_generate_horner(builder, hw0, xvars, None)
                        w1 = self._llvm_generate_horner(builder, hw1, xvars, None)
                        w2 = self._llvm_generate_horner(builder, hw2, xvars, None)
                        w3 = self._llvm_generate_horner(builder, hw3, xvars, None)
                        w4 = self._llvm_generate_horner(builder, hw4, xvars, None)
                        w5 = self._llvm_generate_horner(builder, hw5, xvars, None)
                        w6 = self._llvm_generate_horner(builder, hw6, xvars, None)
                        w7 = self._llvm_generate_horner(builder, hw7, xvars, None)

                        # g = w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7
                        # hz = (w4 + w5 + w6 + w7)
                        # hy = (w2 + w3 + w6 + w7)
                        # hx = (w1 + w3 + w5 + w7)


                        a = builder.fadd(w0,  w1, flags=['fast'])
                        b = builder.fadd(w2,  w3, flags=['fast'])
                        c = builder.fadd(w4,  w5, flags=['fast'])
                        d = builder.fadd(w6,  w7, flags=['fast'])

                        _a = builder.fadd(a, b, flags=['fast'])
                        _b = builder.fadd(c, d, flags=['fast'])

                        _c = builder.fadd(w1,  w3, flags=['fast'])
                        _d = builder.fadd(w5,  w7, flags=['fast'])

                        g = builder.fadd(_a, _b, flags=['fast'])
                        
                        hz = builder.fadd(c, d, flags=['fast'])
                        hy = builder.fadd(b, d, flags=['fast'])
                        hx = builder.fadd(_c, _d, flags=['fast'])

                        tz = builder.fdiv(hz, g)
                        ty = builder.fdiv(hy, g)
                        tx = builder.fdiv(hx, g)

                        terr = builder.call(self._fabs, [g])
                        c    = builder.fcmp_unordered("<=", g, self._internal_type(0.0000000001))
                        
                        tz  = builder.select(c, self._internal_type(0), tz)
                        ty  = builder.select(c, self._internal_type(0), ty)
                        tx  = builder.select(c, self._internal_type(0), tx)

                        # These are our scales
                        pt = [tx, ty, tz]

                        # Scale them for the coset lattice

                        pt = [builder.fmul(_, self._internal_type(G[i,i]) ) for i,_ in enumerate(pt)]
                        pt = [builder.fadd(a,  b, flags=['fast']) for a,b in zip(pt, rp)]
                        # pt[0] = w0
                        # pt[1] = w1
                        # pt[2] = w2

                    else:
                        raise ValueError("Cannot group into higher than 8!")

                    if ref_idx == 0:
                        pt_pred = pt
                        rp_pred = rp
                        g_pred  = g
                    else:
                        c = builder.icmp_unsigned("==", polyindex, ir.IntType(32)(ref_idx))

                        # Do the predication thing...
                        pt_pred = [
                            builder.select(c, newv, current) for
                            newv, current in zip(pt, pt_pred)
                        ]

                        rp_pred = [
                            builder.select(c, newv, current) for
                            newv, current in zip(rp, rp_pred)
                        ]

                        g_pred = builder.select(c, g, g_pred)
                        # raise ValueError('need to predicate')

            pt = pt_pred
            rp = rp_pred
            g  = g_pred

            # Apply the transform to the lookup 
            pt = self._llvm_apply_transform(builder, xform, pt)
            if coset_var is None:
                rp = self._llvm_apply_transform(builder, xform, rp)

            pt = self._llvm_apply_shift(builder, pt, shift)
            if coset_var is None:
                rp = self._llvm_apply_shift(builder, rp, shift)

            # If we have reflective symmetry, we need to apply that now
            if self._ss._reflective:
                pt = [builder.fmul(_,  self._reflective[i], flags=['fast']) for i, _ in enumerate(pt)]
                if coset_var is None:
                    rp = [builder.fmul(_,  self._reflective[i], flags=['fast']) for i, _ in enumerate(rp)]


            # Finally, add in the reference lattice sites
            pt = [builder.fadd(_,  ivars[i], flags=['fast']) for i, _ in enumerate(pt)]

            if coset_var is None:
                rp = [builder.fadd(_,  ivars[i], flags=['fast']) for i, _ in enumerate(rp)]
                coset_var = self._llvm_generate_coset_index(builder, rp)

            if self.lookup_type == "CPU_LIN":
                pt     = self._llvm_generate_coset_scale(builder, coset_var, pt)

                if self.strict_branchless:
                    _ = builder.call(memfetch, [coset_var] + pt)
                else:
                    pred = builder.fcmp_unordered("!=", g, self._internal_type(0.0))
                    c    = self._internal_type(0.0) 

                    bb = builder.basic_block
                    bbif = builder.append_basic_block(name=_label_suffix(bb.name, '.if'))
                    bbend = builder.append_basic_block(name=_label_suffix(bb.name, '.endif'))
                    br = builder.cbranch(pred, bbif, bbend)

                    with builder.goto_block(bbif):
                        output = builder.call(memfetch, [coset_var] + pt)
                        builder.branch(bbend)
                        
                    builder.position_at_start(bbend)

                    _ = builder.phi(self._internal_type)
                    _.add_incoming(c, bb)
                    _.add_incoming(output, bbif)

                if self.pipeline_depth >= 1:
                    if len(in_flight_q) >= self.pipeline_depth:
                        # We can't queue another read, so we consume a read before 
                        # we queue
                        fetch, gf = in_flight_q.pop(0)
                        prd = builder.fmul(fetch,  gf, flags=['fast'])
                        result = builder.fadd(prd,  result, flags=['fast'])
                    # Queue up the last 
                    in_flight_q += [(_, g)]
                else:
                    _ = builder.fmul(_,  g, flags=['fast'])
                    result = builder.fadd(_,  result, flags=['fast'])

            elif self.lookup_type == "TEX_LIN":
                _ = self._textures[0]

                if outer_coset_index is not None:
                    n_ocosets = len(self._ss.get_L_cosets())
                    offset = builder.mul(outer_coset_index, ir.IntType(32)(n_ocosets))
                    coset = builder.add(offset, coset_var)
                else:
                    coset = coset_var

                # predicate texture 
                for idx_tex, tex in enumerate(self._textures[1:]):
                    idx_tex += 1
                    p = builder.icmp_unsigned("==", coset, ir.IntType(32)(idx_tex))
                    _ = builder.select(p, tex, _)

                texture = _ 

                # Shift 
                pt = self._llvm_generate_coset_scale(builder, coset_var, pt)
                pt = [builder.fadd(_,  self._internal_type(0.5), flags=['fast']) for _ in pt]
                pt = [builder.fptrunc(_, ir.FloatType()) for _ in pt]

                if self.strict_branchless:

                    # Fetch 
                    if self._ss._s == 1:
                        mf = builder.call(self._tex1d, [texture]+pt)
                    elif self._ss._s == 2:
                        mf = builder.call(self._tex2d, [texture]+pt)
                    elif self._ss._s == 3:
                        mf = builder.call(self._tex3d, [texture]+pt)
                    else:
                        raise ValueError("Dimension not supported")

                    if self.pipeline_depth >= 1:
                        if len(in_flight_q) >= self.pipeline_depth:
                            # We can't queue another read, so we consume a read before 
                            # we queue
                            fetch, gf = in_flight_q.pop(0)
                            fetch = builder.extract_value(fetch, 0)
                            prd = builder.fmul(fetch,  gf, flags=['fast'])
                            result = builder.fadd(prd,  result, flags=['fast'])
                        # Queue up the last 
                        in_flight_q += [(mf, g)]
                    else:
                        mf = builder.extract_value(mf, 0)
                        _ = builder.fmul(mf,  g, flags=['fast'])
                        result = builder.fadd(_,  result, flags=['fast'])
                else:
                    if self.pipeline_depth >= 1:
                        raise ValueError("Invalid combo")

                    pred = builder.fcmp_unordered("!=", g, self._internal_type(0.0))
                    c    = self._internal_type(0.0) 

                    bb = builder.basic_block
                    bbif = builder.append_basic_block(name=_label_suffix(bb.name, '.if'))
                    bbend = builder.append_basic_block(name=_label_suffix(bb.name, '.endif'))
                    br = builder.cbranch(pred, bbif, bbend)

                    with builder.goto_block(bbif):
                        if self._ss._s == 1:
                            output = builder.call(self._tex1d, [texture]+pt)
                        elif self._ss._s == 2:
                            output = builder.call(self._tex2d, [texture]+pt)
                        elif self._ss._s == 3:
                            output = builder.call(self._tex3d, [texture] + pt)
                        else:
                            raise ValueError("Dimension not supported")
                        output = builder.extract_value(output, 0)
                        builder.branch(bbend)
                        
                    builder.position_at_start(bbend)

                    mf = builder.phi(self._internal_type)
                    mf.add_incoming(c, bb)
                    mf.add_incoming(output, bbif)

                    _ = builder.fmul(mf,  g, flags=['fast'])
                    result = builder.fadd(_,  result, flags=['fast'])

        # We need to take care of any left over fetches   
        for fetch, g in in_flight_q:
            if self.lookup_type == "TEX_LIN":
                fetch = builder.extract_value(fetch, 0)
            prd = builder.fmul(fetch,  g, flags=['fast'])
            result = builder.fadd(prd,  result, flags=['fast'])

        return result

    def _llvm_generate_lerp_conv_sum_branchy(self, xvars, xform, shift, srindex, polyindex, lsindex):
        pass

    def tick(self, indicator='.', end=' '):
        """
        A simple function to show a some sort of 
        tick when progress is made. Broken off into
        """
        if self._verbose:
            print(indicator, end=end)
        sys.stdout.flush()

    def log(self, strng, end='\n'):
        """
        Logs all the text to 
        """
        if self._verbose:
            print(strng,end=end)
        sys.stdout.flush()


    def setup_parser(self):
        import sys
        import pyparsing

        ppc = pyparsing.pyparsing_common

        pyparsing.ParserElement.enablePackrat()
        sys.setrecursionlimit(3000)

        integer = ppc.integer
        variable = pyparsing.Word(pyparsing.alphas,  pyparsing.alphanums)
        operand = integer | variable

        expop = pyparsing.Literal("^")
        signop = pyparsing.oneOf("+ -")
        multop = pyparsing.oneOf("* /")
        plusop = pyparsing.oneOf("+ -")
        factop = pyparsing.Literal("!")
        self.infix_parser = pyparsing.infixNotation(
            operand,
            [
                ("!", 1, pyparsing.opAssoc.LEFT),
                ("^", 2, pyparsing.opAssoc.RIGHT),
                (signop, 1, pyparsing.opAssoc.RIGHT),
                (multop, 2, pyparsing.opAssoc.LEFT),
                (plusop, 2, pyparsing.opAssoc.LEFT),
            ],
        )

    def parse_expression(self, str_expression, variables):
        """

        :param str_expression:
        :param stack_vars:
        :return:
        """
        self.setup_parser()
        L = self.infix_parser.parseString(str_expression).asList()
        valid_ops = "*","/","-","+","^","|","&"

        def recursive_check(expression_in):
            for _ in expression_in:
                if isinstance(_, list):
                    recursive_check(_)
                else:
                    if _ not in valid_ops and _ not in variables:
                        raise ValueError(f"Expression '{str_expression}' contains invalid variable '{_}'" )
            return True
        if recursive_check(L):
            return L
        return None

    def generate_expression(self, expression, builder, variables, target=None):
        """
        Generates a "sensible" sequence of LLVM for a high level expression.
        :param expression:
        :return:
        """
        global output_type
        valid_ops = "*","/","-","+","^","|","&"
        output_type = None
        if target is not None:
            output_type = target.type

        def load_var(variable):
            if isinstance(variable, ir.instructions.AllocaInstr):
                return builder.load(variable)

            if isinstance(variable, ir.values.Constant):
                return variable
            return variable
            # elif isinstance(variable, ir.instructions)

        def convert_to_type(a, type, builder):
            if a.type == type:
                return a
            return a

        def generate_instruction(token, a, b, type, builder):
            a = convert_to_type(a, type, builder)
            b = convert_to_type(b, type, builder)

            if token == '+':
                res = builder.fadd(a, b, flags=['fast'])
            elif token == "-":
                res = builder.fsub(a, b, flags=['fast'])
            elif token == "*":
                res = builder.fmul(a, b, flags=['fast'])
            else:
                raise ValueError(f"")
            return res

        def recursive_generate(expr):
            global output_type
            # First, recurse on all lower expression
            expr = [recursive_generate(token) if isinstance(token, list) else token for token in expr]
            result = None
            next_op = None
            for token in expr:
                if token not in valid_ops and (token in variables or isinstance(token, ir.Instruction)):
                    if token in variables:
                        value = load_var(variables[token])
                    else:
                        value = token

                    if result is None:
                        result = value

                    if output_type is None:
                        output_type = value.type

                    if next_op is not None:
                        result = generate_instruction(next_op, result, value, output_type, builder)
                        next_op = None

                elif token in valid_ops:
                    next_op = token
                    continue

                elif result is None:
                    result = token

            return result

        retval = recursive_generate(expression)
        if target is not None:
            builder.store(retval, target)
        return retval

    def _llvm_load_transform_matrix_vec(self, builder, srindex, transpose=False):
        """
        Returns a 2D array of matrix elements [i][j]
        """
        if srindex is None:
            return None, []
        mat = [[None for _ in range(self._dimension)] for _ in range(self._dimension)]

        ptr_size = int(64)
        index_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(32), self.lane_count)
        ptr_type = ir.IntType(32) if self.lane_count == 1 else ir.VectorType(ir.IntType(ptr_size), self.lane_count)

        # Load the offset of the transform
        if self.lane_count == 1:
            subregion_transform = builder.gep(self.xform_lut, [ir.IntType(32)(0), srindex])
        else:
            bcst = builder.insert_element(ptr_type(ir.Undefined),
                                          self.xform_lut.ptrtoint(ir.IntType(ptr_size)),
                                          ir.IntType(32)(0), name='xform_map')

            bcst = builder.shuffle_vector(bcst, ptr_type(ir.Undefined), index_type(None))
            subregion_transform = builder.inttoptr(bcst,
                                                   ir.VectorType(self._xform_struct.as_pointer(), self.lane_count))

        for i in range(self._dimension):
            for j in range(self._dimension):
                offset = i * self._dimension + j

                if self.lane_count == 1:
                    mep = builder.gep(subregion_transform, [ir.IntType(32)(0), ir.IntType(32)(int(offset))])
                else:
                    mep = builder.gep(subregion_transform, [srindex, ir.IntType(32)(int(offset))])
                mel = self.read_mem(mep, builder)
                mel = self._int_to_float(mel, mel.type, self._internal_type, builder)

                if transpose:
                    mat[j][i] = mel
                else:
                    mat[i][j] = mel

        if self.no_shift:
            return mat, []

        offset = self._dimension * self._dimension
        shift = []
        for index in range(self._dimension):
            if self.lane_count == 1:
                mp = builder.gep(subregion_transform, [ir.IntType(32)(0), ir.IntType(32)(int(offset + index))])
                sel = builder.load(mp)
                sel = self._llvm_convert_to_internal_type(builder, sel)
            else:
                mep = builder.gep(subregion_transform, [srindex, ir.IntType(32)(int(offset + index))])
                mel = self.read_mem(mep, builder)
                sel = self._int_to_float(mel, mel.type, self._internal_type, builder)

            shift += [sel]

        return mat, shift
