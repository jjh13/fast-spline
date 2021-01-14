#!/usr/bin/env python

"""
Main class that implements part I of the paper. Basically all of the analysis 
happens here, anything that needs a one-time computation is intended to be in 
part I.

Copyright © 2016, 2020 Joshua Horacsek

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


"""
Overall todo:
Check for address spaces
"""

class EvalCodeGen:

    def __init__(self, ss, **kwargs): 
        options = kwargs

        target_triple="unknown-unknown-unknown"
        
        if "target_triple" in options:
            target_triple = options["target_triple"]

        self._verbose = True
        if "verbose" in options:
            self._verbose["verbose"]


        # Break the target triple up
        tt = target_triple.split('-')
        if len(tt) < 3:
            raise ValueError("Invalid target_triple '%s'!")

        self._arch, self._vendor, self._sys = tt[:3]

        self._target       = target_triple
        self._known_target = False
        self.coset_decomposition   = True


        # By default, setup the parameters for CPU based code
        if 'arm' in self._arch or 'aarch' in self._arch or 'x86' in self._arch or 'avr' in self._arch:
            self._known_target = True

        #@ Will this spline use texture fetches, if they're available?
        # texture fetches of at most dimenion 3  are only supported my most GPUs
        self.has_texture_fetch  = False

        #@ Will this spline use linear fetches on cartesian cosets? 
        # by default this is on if we want to use texture fetches
        self.use_linear_fetch   = False

        #@ Which type of memory fetches to use
        self.lookup_type        = "COSET" if self.coset_decomposition else "ABSTRACT" # ABSTRACT, COSET, TEX, TEX_LIN

        #@ Whether or not to generate code with branching behaviour
        # this is only beneficial on the GPU
        self.branchless         = False

        #@ Whether or not to use branch predication (if doing branchless)
        self.branch_predication = False

        #@ Break horner evaluation up over different SIMD lanes, 
        # these don't exist on GPUs, so we disable by default
        self.vectorize          = True 

        self.lane_count         = 4       # Most CPUs have a vector extension with at least 4 elements

        # Should we refetch the transform on each lattice site read?
        # On CPUs - no, we should just store the transform on the stack, 
        # but on GPU we should always refetch --- this allows registers 
        # to be saved
        self.refetch_transform  = False

        #@ How many memfetchest to queue before computing the polynomial
        self.group_size     = floor(ss._ref_subregions[0].number_of_pts_per_reconstruction() / 4)
        self.group_size     = ss._ref_subregions[0].number_of_pts_per_reconstruction() if self.group_size == 0 else self.group_size
        self.pipeline_depth  = 0


        if self._arch in ['nvptx', 'nvptx64', 'amdil', 'amdil64', 'hsail', 'hsail64', 'spir', 'spir64', 'r600','amdgcn']:
            self._known_target      = True

            self.has_texture_fetch  = self._arch in ['nvptx', 'nvptx64'] 
            self.lookup_type        = "TEX_LIN" if self.coset_decomposition and self.has_texture_fetch else "ABSTRACT"

            self.branchless         = True
            self.branch_predication = True
            self.vectorize          = False
            self.lane_count         = 1
            self.refetch_transform  = True

        if not self._known_target:
            self.log("Unknown architecture '%s' continuing with common PC parameters" % self._arch)

        if self.lookup_type == "TEX_LIN":
             self.group_size = 1

        self.reconstruction_name = "__reconstruct__"
        self.strict_branchless   = False

        """
        Patch in options from the input options
        """
        for k in options:
            if k == 'use_texture_fetch':
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
                self.reconstruction_name = options[k]
            elif k == 'refetch_transform':
                self.refetch_transform = options[k]
            elif k == 'group_size':
                try:
                    self.group_size = int(options[k])
                except:
                    raise ValueError('group_size option must be convertable to an int!')
            elif k == 'lookup_type':
                if options[k] not in ['ABSTRACT','COSET','TEX','TEX_LIN', 'CPU_LIN']:
                    raise ValueError("lookup_type must be one of ABSTRACT, COSET, TEX, TEX_LIN")
                self.lookup_type = options[k]

            else:
                raise ValueError('Unknown option "%s"' % k)


        """
        Double check that the parameters we have are good
        """ 
        if not self.branchless and self.branch_predication:
            raise ValueError("Cannot generate branch predicated code when using branching")

        if self.lookup_type in ['COSET','TEX','TEX_LIN', 'CPU_LIN'] and not self.coset_decomposition:
            raise ValueError("Cannot use a COSET, TEX or TEX_LIN if we don't know the coset decomposition for the spline space!")

        if self.lookup_type in ['TEX','TEX_LIN'] and not self.has_texture_fetch:
            raise ValueError("Cannot generate code with texture fetches on architectures without texture units!")

        if (self.lookup_type == "TEX_LIN" or self.lookup_type == "CPU_LIN") and self.group_size != 1:
            raise ValueError("Cannot pipeline linear texture fetches")


        self.log("Generating code for '%s':" % self._arch)
        self.log("   branchless: " + str(self.branchless))
        self.log("   branch_predication: " + str(self.branch_predication))
        self.log("   group_size: " + str(self.group_size))
        self.log("   lane_count: " + str(self.lane_count))
        self.log("   refetch_transform: " + str(self.refetch_transform))
        self.log("   lookup_type: " + str(self.lookup_type))


        self._ss = ss
        self._process_sub_regions()

    def _process_sub_regions(self):
        self._subregion_consistency = True

        # For each sub region, bin the lattice sites into their respective
        # coset structure, and do a greedy horner factorization on each 
        # subregion
        self.log("Processing reference subregions...")
        G, coset_vectors = self._ss.coset_vectors()
        for region_idx, region in enumerate(self._ss._ref_subregions):

            self.tick("Partitioning reference subregion %d into cosets... " % region_idx)
            region.partition_into_cosets(G, coset_vectors)

            self.tick("factoring...")
            region.factor_for_evaluation()

            # Use dancing links to cover each lookup, if necessary 
            self.tick("covering...")
            region.cover_lookups(self.lookup_type == "TEX_LIN" or self.lookup_type == "CPU_LIN")

            # Plan out the lookups --- if we use a texture fetch, 
            # then we can take advantage of texture locality
            self.tick("optimizing lookups...")
            region.order_lookups()

            self.log("done!")

        # Check for consistency between sub-regions, that is, check to see 
        # that each sub region has the same number of accesses per coset,
        # this means we can avoid determining which coset a memory fetch is 
        # on every time we access a lattice site
        subregion_lk_sizes = [len(_) for _ in self._ss._ref_subregions[0]._coset_partition]
        for refsregion in self._ss._ref_subregions[1:]:
            _sizes = [len(_) for _ in refsregion._coset_partition]
            if any(_ != __ for (_,__) in zip(subregion_lk_sizes, _sizes)):
                self.log("Subregions are inconsistent -- turning off consistency optimization")
                self._subregion_consistency = False

        if self._subregion_consistency:
            self.log("Consistent subregions!")


    """
    ####################################
    # LLVM Code generation functions
    ####################################

    These are the workhorses of this class -- generally they don't
    depend on instance variables when they don't have to. This modular
    design allows me to test each one in isolation.
    """

    def _llvm_generate_module(self):
        """
        Sets up the module, takes values from the init 
        to create all the type definitions we need. Also
        declares and sets up intrinsics for this platform
        """
        self._intp = self._internal_type = ir.FloatType()
        self._iptp = self._input_type    = ir.FloatType()
        self._optp = self._output_type   = ir.FloatType()
        self._mftp = self._memfetch_type = ir.FloatType()

        self._llvm_module = ir.Module(name="fastspline_codegen")
        self._llvm_module.triple = self._target

        # Declare intrinsics for the round and floor functions
        rfn_type = ir.FunctionType(self._internal_type, (self._internal_type,))
        self._round = self._llvm_module.declare_intrinsic("llvm.round", [self._internal_type], fnty=rfn_type)
        
        ffn_type = ir.FunctionType(self._internal_type, (self._internal_type,))
        self._floor = self._llvm_module.declare_intrinsic("llvm.floor", [self._internal_type], fnty=ffn_type)

        ffn_type = ir.FunctionType(self._internal_type, (self._internal_type,))
        self._fabs = self._llvm_module.declare_intrinsic("llvm.fabs", [self._internal_type], fnty=ffn_type)
        
        #####
        # Declare the memory lookup call types
        #####
        
        ####################
        # Abstract lookups
        # When using abstract lookups, the reconstruction function will expect 
        # to be paseed a pointer to a function that returns the 'memfetch' type
        # and takes 'dimension_s' arguments of integers --- i.e. lattice sites
        input_types = ([ir.IntType(int(32))]*self._ss._s)
        self._amftp = ir.FunctionType(self._mftp, input_types)
        
        ####################
        # Coset lookups
        # When using coset lookups, the reconstruction function will expect 
        # to be paseed a pointer to a function that returns the 'memfetch' type
        # and takes one integer which tells the function which coset to lookup
        # and  'dimension_s' arguments of integers --- i.e. lattice sites on the coset
        #
        # In reality, this complexity could be pused off into the abstract lookup function, 
        # but this is a good middle-ground for testing as it somewhat simulates the case with
        # nn texture fetches
        coset_levels = 1 if len(self._ss.get_L_cosets()) == 1 else 2
        input_types = ([ir.IntType(32)]*(self._ss._s+coset_levels))
        self._cluptp = ir.FunctionType(self._mftp, input_types)
        
        # linear interpolated lookups (test jig)
        input_types = ([ir.IntType(32)]+[self._internal_type]*(self._ss._s))
        self._clerptp = ir.FunctionType(self._mftp, input_types)

        self.constant_addrspace = 0
        # TODO: Cuda intrinsics for texture fetches
        if self.lookup_type == "TEX" or self.lookup_type == "TEX_LIN":
            float4 = ir.LiteralStructType([ir.FloatType(),ir.FloatType(),ir.FloatType(),ir.FloatType()])
            fnty = ir.FunctionType(float4, (ir.IntType(int(64)), ir.FloatType(), ir.FloatType(), ir.FloatType()))
            self._tex3d =  self._llvm_module.declare_intrinsic("llvm.nvvm.tex.unified.3d.v4f32", [ir.FloatType()], fnty)

            float4 = ir.LiteralStructType([ir.FloatType(),ir.FloatType(),ir.FloatType(),ir.FloatType()])
            fnty = ir.FunctionType(float4, (ir.IntType(int(64)), ir.FloatType(), ir.FloatType()))
            self._tex2d =  self._llvm_module.declare_intrinsic("llvm.nvvm.tex.unified.2d.v4f32", [ir.FloatType()], fnty)

            float4 = ir.LiteralStructType([ir.FloatType(),ir.FloatType(),ir.FloatType(),ir.FloatType()])
            fnty = ir.FunctionType(float4, (ir.IntType(int(64)), ir.FloatType()))
            self._tex1d =  self._llvm_module.declare_intrinsic("llvm.nvvm.tex.unified.1d.v4f32", [ir.FloatType()], fnty)
            self.constant_addrspace = 4

        return self._llvm_module

    def _llvm_generate_function(self, module):
        """
        Generates the function and its preamble (i.e, allocates space on the stack for any variables)
        that "should" presist -- i.e. if any passes assume stack variables are somehow imporant will
        get the proper cue.

        Arguments: 
        module - the llvm module

        What's left:
        """

        # Build types
        input_types = ([self._input_type]*self._ss._s)
        ntextures = len(self._ss.coset_vectors()[1])*len(self._ss.get_L_cosets())

        # mf  (memory fetch) is an input type that depends on th
        if self.lookup_type == 'ABSTRACT':
            mf = [self._amftp.as_pointer()] # Abstract lookups take a pointer to a function
        elif self.lookup_type == 'COSET':
            mf = [self._cluptp.as_pointer()] # Same with coset (but the function call is different)
        elif self.lookup_type == 'CPU_LIN':
            mf = [self._clerptp.as_pointer()] # Debug jig, a linear fetch emulated in CPU
        elif self.lookup_type == 'TEX' and 'nvptx64' in self._arch and self._vendor == 'nvidia':
            # Int type is the type textures are aliased as, so those are the input
            # To the mem lookup code
            mf = [ir.IntType(int(64)) for _ in range(ntextures)]
        elif self.lookup_type == 'TEX_LIN' and 'nvptx64' in self._arch and self._vendor == 'nvidia':
            mf = [ir.IntType(int(64)) for _ in range(ntextures)]

        # Create a new type for our function
        fn_type = ir.FunctionType(self._output_type, tuple(input_types + mf))

        self._function = ir.Function(self._llvm_module, fn_type, name=self.reconstruction_name)

        # TODO: Other memfetch types
        # DOING: I don't think this applies anymore
        # there aren't any other memfetch types? Maybe AMD?
        memfetch = self._function.args[-1] # <- NVM this, this needs to be fixed
        xvars    = self._function.args[0:self._ss._s]

        try:
            self._textures = self._function.args[self._ss._s:] 
        except:
            pass # TODO: Uh

        self._entry_blk = self._function.append_basic_block(name="entry")

        return  ir.IRBuilder(self._entry_blk), xvars, memfetch


    def _llvm_generate_bsp_lookup_table(self, module):
        """
        Emits the BSP lookup table, (i.e. the sub-region table)
        to the module

        Returns None if no table is needed, otherwise a 
        reference to the table

        """  
        # If we only have one sub region, we don't need a bsp lookup
        if len(self._ss._subregions) <= 0:
            return
        
        # Check how many sub regions we have so we can restict to an
        # appropriate bit width
        bw = int(_bw(len(self._ss._subregions)))
        
        plane_index_t = ir.IntType(bw)
        plane_index_ta = ir.ArrayType(plane_index_t, len(self._ss._sregion_index))
        
        # Convert the packed index to the correct type
        bsp_index_c = plane_index_ta([
            plane_index_t(_) for _ in self._ss._sregion_index 
        ])
        
        # Setup the global variables
        bsp_index = ir.GlobalVariable(module, plane_index_ta, "bsp_index", self.constant_addrspace)
        bsp_index.global_constant = True
        bsp_index.initializer = bsp_index_c

        self.bsp_index = bsp_index
        return bsp_index


    def _llvm_generate_xfrom_lookup_table(self, module):
        """
        Emits the transform lookup table to the module

        Returns None if no table is needed, otherwise a 
        reference to the table

        TODO:
          - Nothing
        """  
        if len(self._ss._subregions) <= 0:
            return None
        
        # Take a look to ensure that all the sub regions are non-redundant
        # this shol
        for ref_idx in self._ss._sregion_index:
            if ref_idx >= self._ss._non_redundant_srgn_cnt:
                raise ValueError("Index pointed to a redundant region!")

        # Take a look at all the transforms, and choose an 
        # appropriate type for them
        T_type = ["integer", 8]
        t_type = ["integer", 8]
        no_shift = True
        
        # We can just blast thru the non redundant sub regions 
        for ridx in range(self._ss._non_redundant_srgn_cnt):
            T,t = self._ss._subregions[ridx]._transform
            
            for value in T.list():
                if round(value) == value:
                    if _bw(abs(int(value))) > T_type[1]:
                        T_type[1] = _bw(abs(int(value)))
                else:
                    T_type[0] = "internal_float"
                    
            for value in t.list():
                if round(value) == value:
                    if int(value) != 0:
                        no_shift = False
                    if _bw(abs(int(value))) > t_type[1]:
                        t_type[1] = _bw(abs(int(value)))
                    
                else:
                    no_shift = False
                    t_type[0] = "internal_float"
        
        no_shift = False
        # Construct the types
        if T_type[0] == "integer":
            self._transform_type = ir.IntType(T_type[1])
        else:
            self._transform_type = self._internal_type
                
        if t_type[0] == "integer":
            self._shift_type = ir.IntType(t_type[1])
        else:
            self._shift_type = self._internal_type
            
        # Create the strct type
        msize   = self._ss._s
        msize2  = msize*msize
        
        if no_shift:
            msize = 0
        
        self._xform_struct = ir.LiteralStructType(
            ([self._transform_type]*msize2) + ([self._shift_type]*msize), 
            packed=True
        )
        
        self._xform_struct_array = ir.ArrayType(
            self._xform_struct, 
            self._ss._non_redundant_srgn_cnt
        )

        # Pack the LUT
        self.xform_lut = ir.GlobalVariable(module, self._xform_struct_array, "xform_lookup", self.constant_addrspace)

        sl = []
        for ridx in range(self._ss._non_redundant_srgn_cnt):
            T,t = self._ss._subregions[ridx]._transform
            
            a = [self._transform_type(_) for _ in T.list()]
            if not no_shift:
                b = [self._shift_type(_) for _ in t.list()]
            else:
                b = []
                
            sl += [
                self._xform_struct(a+b)
            ]
            
        self.no_shift = no_shift
        self.xform_lut.global_constant = True
        self.xform_lut.initializer = self._xform_struct_array(sl)

    def _llvm_generate_ref_region_index_lookup_table(self, module):
        if len(self._ss._ref_subregions) <= 1:
            return None

        bw = int(_bw(len(self._ss._ref_subregions)))
        poly_index_t  = ir.IntType(bw)
        poly_index_ta = ir.ArrayType(poly_index_t, self._ss._non_redundant_srgn_cnt)
        
        # Convert the packed index to the correct type
        index = []
        for ridx in range(self._ss._non_redundant_srgn_cnt):
            polyindex = self._ss._subregions[ridx]._ref_region
            index += [poly_index_t(polyindex)]

        polyindex_inst = poly_index_ta(index)
        
        # Setup the global variables
        poly_index                 = ir.GlobalVariable(module, poly_index_ta, "poly_index", self.constant_addrspace)
        poly_index.global_constant = True
        poly_index.initializer     = polyindex_inst
        
        # For GPU architectures, declare this in constant memory
        # bsp_index.addrspace = 0 # TODO  change to constant memory addr space
    
        return poly_index

    def _llvm_generate_lattice_site_lookup_table(self, module):
        if len(self._ss._ref_subregions) <= 1:
            return None

        # 
        rsr_count  = len(self._ss._ref_subregions)
        pts_per_sr = self._ss._ref_subregions[0].number_of_pts_per_reconstruction() 

        # Go thru each sub region, and determine the widest type we need
        bw = max([_bw(abs(_)) for _ in sum(sum([[list(k) for k in rsr._weighted_sum_terms.keys()] for rsr in self._ss._ref_subregions],[]), [])])

        lsc_index_t    = ir.IntType(bw)
        ls_index_t     = ir.ArrayType(lsc_index_t, self._ss._s)
        lslst_index_t  = ir.ArrayType(ls_index_t, pts_per_sr)
        lsidx_t        = ir.ArrayType(lslst_index_t, rsr_count)

        lattice_site_lookup = []
        for rsr in self._ss._ref_subregions:
            srls_index = []
            for ls in rsr.get_optimal_nn_lookups():
                srls_index += [ls_index_t([lsc_index_t(_) for _ in ls])]

            lattice_site_lookup += [lslst_index_t(srls_index)]
        lslk = lsidx_t(lattice_site_lookup)

        ls_index                 = ir.GlobalVariable(module, lsidx_t, "ls_index", self.constant_addrspace)
        ls_index.global_constant = True
        ls_index.initializer     = lslk
        return ls_index


    def _llvm_generate_outer_coset_table(self, module):
        coset_shifts = self._ss.get_L_cosets()
        if len(coset_shifts) == 1:
            return None

        # Determine the max bit-width we need for these
        bw = max([_bw(__) for __ in sum([list(_) for _ in coset_shifts], [])])


        ocs_index_t    = ir.IntType(bw)
        ocsl_index_t   = ir.ArrayType(ocs_index_t, self._ss._s)
        ocl_t          = ir.ArrayType(ocsl_index_t, len(coset_shifts))

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
        #_llvm_generate_outer_coset_index

        # raise ValueError("Not implemented")

    def _llvm_generate_rho_bccind(self, builder, x_variables):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        i_variables = [None]*3

        # Stage 1
        _ = x_variables[0]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.call(self._round, [_])
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        i_variables[0] = _
        
        _ = x_variables[1]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.call(self._round, [_])
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        i_variables[1] = _

        _ = x_variables[2]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.call(self._round, [_])
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        i_variables[2] = _
        
        # Sd
        _ = builder.fsub(i_variables[0], x_variables[0])
        a = builder.fmul(_, _, flags=['fast'])
        
        _ = builder.fsub(i_variables[1], x_variables[1])
        b = builder.fmul(_, _, flags=['fast'])
        
        a = builder.fadd(a, b, flags=['fast'])
        _ = builder.fsub(i_variables[2], x_variables[2])
        _ = builder.fmul(_, _, flags=['fast'])
        sd = builder.fadd(a,_, name="sd")
        
        
        tx = [
            builder.fsub(x_variables[0], self._internal_type(1)),
            builder.fsub(x_variables[1], self._internal_type(1)),
            builder.fsub(x_variables[2], self._internal_type(1)),
        ]
        
        # Stage 3
        _ = tx[0]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.call(self._round, [_])
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_,  self._internal_type(1.0), flags=['fast'])
        tx[0] = _
        
        _ = tx[1]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.call(self._round, [_])
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_,  self._internal_type(1.0), flags=['fast'])
        tx[1] = _

        _ = tx[2]
        _ = builder.fdiv(_, self._internal_type(2.0))
        _ = builder.call(self._round, [_])
        _ = builder.fmul(_,  self._internal_type(2.0), flags=['fast'])
        _ = builder.fadd(_,  self._internal_type(1.0), flags=['fast'])
        tx[2] = _

        # Sd
        _ = builder.fsub(tx[0], x_variables[0])
        a = builder.fmul(_, _, flags=['fast'])
        
        _ = builder.fsub(tx[1], x_variables[1])
        b = builder.fmul(_, _, flags=['fast'])
        
        a = builder.fadd(a, b, flags=['fast'])
        _ = builder.fsub(tx[2], x_variables[2])
        _ = builder.fmul(_, _, flags=['fast'])
        sd2 = builder.fadd(a,_, name="sd2")
        cmp_dist = builder.fcmp_ordered("<", sd, sd2)
        
        i_variables[0] = builder.select(cmp_dist, i_variables[0], tx[0])
        i_variables[1] = builder.select(cmp_dist, i_variables[1], tx[1])
        i_variables[2] = builder.select(cmp_dist, i_variables[2], tx[2])

        return x_variables, i_variables
        

    def _llvm_generate_rho_cartesian(self, builder, x_variables):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        i_variables = [
            builder.call(self._round, [xv]) for xv in x_variables
        ]
        return x_variables, i_variables

    def _llvm_generate_rho_fccind(self, builder, x_variables):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        i_variables = [None]*3

        i_variables[0] = builder.call(self._round, [x_variables[0]])
        i_variables[1] = builder.call(self._round, [x_variables[1]])
        i_variables[2] = builder.call(self._round, [x_variables[2]])

        ai = builder.fptoui(i_variables[0], ir.IntType(int(32)))
        aj = builder.fptoui(i_variables[1], ir.IntType(int(32)))
        ak = builder.fptoui(i_variables[2], ir.IntType(int(32)))

        onfcc = builder.and_(builder.add(builder.add(ai, aj), ak), ir.IntType(int(32))(1))
        onfcc = builder.trunc(onfcc, ir.IntType(int(1)))

        xx = builder.fsub(x_variables[0], i_variables[0])
        yy = builder.fsub(x_variables[1], i_variables[1])
        zz = builder.fsub(x_variables[2], i_variables[2])

        xxa = builder.call(self._fabs, [xx])
        yya = builder.call(self._fabs, [yy])
        zza = builder.call(self._fabs, [zz])

        yyxx = builder.fcmp_ordered(">=", yya, xxa)
        yyzz = builder.fcmp_ordered(">=", yya, zza)
        zzxx = builder.fcmp_ordered(">=", zza, xxa)
        zzyy = builder.fcmp_ordered(">=", zza, yya)

        idx = builder.select(builder.and_(yyxx, yyzz), ir.IntType(8)(1), ir.IntType(8)(0))
        idx = builder.select(builder.and_(zzxx, zzyy), ir.IntType(8)(2), idx)

        xx = builder.select(builder.fcmp_ordered(">=", xx, self._internal_type(0.0)), self._internal_type(1), self._internal_type(-1))
        yy = builder.select(builder.fcmp_ordered(">=", yy, self._internal_type(0.0)), self._internal_type(1), self._internal_type(-1))
        zz = builder.select(builder.fcmp_ordered(">=", zz, self._internal_type(0.0)), self._internal_type(1), self._internal_type(-1))

        nonfcc = onfcc # builder.not_(onfcc, 'negate_onfcc')
        xx = builder.select(builder.and_(builder.icmp_unsigned("==", idx, ir.IntType(8)(0)), nonfcc), xx, self._internal_type(0))
        yy = builder.select(builder.and_(builder.icmp_unsigned("==", idx, ir.IntType(8)(1)), nonfcc), yy, self._internal_type(0))
        zz = builder.select(builder.and_(builder.icmp_unsigned("==", idx, ir.IntType(8)(2)), nonfcc), zz, self._internal_type(0))

        i_variables[0] = builder.fadd(xx,  i_variables[0], flags=['fast'])
        i_variables[1] = builder.fadd(yy,  i_variables[1], flags=['fast'])
        i_variables[2] = builder.fadd(zz,  i_variables[2], flags=['fast'])


        return x_variables, i_variables

    def _llvm_generate_rho_pp(self, builder, x_variables):
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
                tmp += [
                    builder.fmul(xw,  self._internal_type(Ti[i,j]), flags=['fast'])
                ]

            result = tmp[0]
            for l in tmp[1:]:
                result = builder.fadd(result, l)

            iv += [builder.call(self._floor, [result])]

        # Forward 
        i_variables = []
        for i, xv in enumerate(x_variables):
            tmp = []
            for j, xw in enumerate(iv):
                tmp += [
                    builder.fmul(xw,  self._internal_type(T[i,j]), flags=['fast'])
                ]

            result = tmp[0]
            for l in tmp[1:]:
                result = builder.fadd(result, l)

            i_variables += [result]

        return x_variables, i_variables



    def _llvm_generate_rho_generic(self, builder, x_variables, vectorize=False):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        D, offsets = self._ss._lattice.coset_structure()
        _it = self._internal_type

        for i, offset in enumerate(offsets):

            # Shift by offset
            xt = [builder.fsub(x_i, _it(component)) for (x_i, component) in zip(x_variables,offset)]

            # Scale by 1/D
            xt = [builder.fdiv(x_i, _it(component)) for (x_i, component) in zip(xt,D.diagonal())]

            # Round
            xt = [builder.call(self._round, [_]) for _ in xt]

            # Scale by D
            xt = [builder.fmul(x_i, _it(component)) for (x_i, component) in zip(xt,D.diagonal())]
            xt = [builder.fadd(x_i, _it(component)) for (x_i, component) in zip(xt,offset)]


            d = _it(0)
            for (a,b) in zip(x_variables, xt):
                t = builder.fsub(a,b)
                t = builder.fmul(t,t)
                d = builder.fadd(t,d)

            if i == 0:
                best = xt
                bestd = d
            else:
                p = builder.fcmp_ordered("<", d, bestd)
                bestd = builder.select(p, d, bestd)
                best = [builder.select(p, new, old) for (new,old) in zip(xt, best)]
        return x_variables, best

    def _llvm_generate_rho(self, builder, x_variables):
        """
        Returns a tuple of x_variables, and i_variables (which is 
        the lattice shift \\rho).

        Inputs: builder ir builder.
        x_variables
        """
        i_variables = [None]*self._ss._s
        i_refs      = None
        
        # Extend bit if needed
        if self._internal_type == ir.DoubleType() and self._input_type == ir.FloatType():
            x_variables = [builder.fpext(_, self._internal_type) for _ in x_variables]

        # or truncate
        if self._internal_type == ir.FloatType() and self._input_type == ir.DoubleType():
            x_variables = [builder.fptrunc(_, self._internal_type) for _ in x_variables]

        # store the original xvariables 
        cxvar = x_variables[:]
        
        # shift by the rho shift
        x_variables = [
                builder.fsub(x_variables[i], self._internal_type(sh)) if sh != 0 else x_variables[i]
            for i, sh in enumerate(self._ss._rho_shift)]

        # TODO: This doesn't account for PP???
        if self._ss._rho_type == "PP":
            x_variables, i_variables = self._llvm_generate_rho_pp(builder, x_variables)
        else:
            if self._ss._s == 3:

                if self._ss._lattice.hash(True) == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                    x_variables, i_variables = self._llvm_generate_rho_bccind(builder, x_variables)

                elif self._ss._lattice.hash(True) == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':                
                    x_variables, i_variables = self._llvm_generate_rho_fccind(builder, x_variables)

                elif self._ss._lattice.hash(True) == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                    x_variables, i_variables = self._llvm_generate_rho_cartesian(builder, x_variables)
                else:
                    x_variables, i_variables = self._llvm_generate_rho_generic(builder, x_variables)

            elif self._ss._lattice_matrix == matrix.identity(self._ss._s):
                x_variables, i_variables = self._llvm_generate_rho_cartesian(builder, x_variables)

            else:
                x_variables, i_variables = self._llvm_generate_rho_generic(builder, x_variables)

        # Shift to the reference region
        cxvar = [builder.fsub(x,  i, flags=['fast']) for (x,i) in zip(cxvar, i_variables)]
        
        # Use reflective symmetry hack
        if self._ss._reflective:
            self._reflective = []
            for i,x_i in enumerate(cxvar):
                self._reflective += [builder.select(
                    builder.fcmp_ordered(">=", x_i, self._internal_type(0)),
                    self._internal_type(1.0), self._internal_type(-1.0), name="refl_%d" % i)]
                
            cxvar = [builder.fmul(_,  self._reflective[i], flags=['fast']) for i, _ in enumerate(cxvar)]
            
        # Generate the shift to the reference region, 
        # Return x_0... and i_0...
        return cxvar, i_variables

    def _llvm_generate_bsp_index(self, builder, x_variables):
        """
        Returns a the bsp index

        Inputs: builder ir builder.
        x_variables
        """
        index = ir.IntType(32)(int(0))

        # We only need to generate this test if we have more than one 
        # sub-region
        if len(self._ss._subregions) > 1:

            # Generate the plane test
            for i, (n,d) in enumerate(reversed(self._ss._plist)):
                running_sum = None

                for x_i, n_i in zip(x_variables, n):
                    # Do the dot product
                    _ = builder.fmul(x_i,  self._internal_type(n_i), flags=['fast'])
                    if running_sum is not None:
                        running_sum = builder.fadd(_,  running_sum, flags=['fast'])
                    else:
                        running_sum = _
                      
                p = builder.fcmp_ordered(">=", running_sum, self._internal_type(float(-d)))
                idx = builder.select(p, ir.IntType(int(32))(int(1 << i)), ir.IntType(int(32))(0))

                if index is not None:
                    index = builder.or_(index, idx)
                else:
                    index = idx

            index = builder.urem(index, ir.IntType(32)(int(self._ss._sregion_modulus)))
            index_ptr =  builder.gep(self.bsp_index, [ir.IntType(32)(0), index])
            index = builder.load(index_ptr)

        return index


    def _llvm_generate_poly_index(self, builder, srindex, poly_tbl):
        """
        Returns a variable that holds the identity (an integer) of the polynomial
        based on the sub-region index (srindex)

        Inputs: builder ir builder.
        x_variables
        """
        if len(self._ss._ref_subregions) <= 1:
            return ir.IntType(int(32))(0)

        ptptr = builder.gep(poly_tbl, [ir.IntType(32)(0), srindex])
        return builder.load(ptptr)

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

    def _llvm_load_transform_matrix(self, builder, srindex, transpose = False):
        """
        Returns a 2D array of matrix elements [i][j] 
        """
        mat = [[None for _ in range(self._ss._s)] for _ in range(self._ss._s)]

        # Load the offset of the transform 
        subregion_Tt = builder.gep(self.xform_lut, [ir.IntType(32)(0), srindex])

        for i in range(self._ss._s):
            for j in range(self._ss._s):
                offset = i * self._ss._s + j
                
                mep = builder.gep(subregion_Tt, [ir.IntType(32)(0), ir.IntType(32)(int(offset))])
                mel = builder.load(mep)
                mel = self._llvm_convert_to_internal_type(builder, mel)

                if transpose:
                    mat[j][i] = mel
                else:
                    mat[i][j] = mel

        return mat

    def _llvm_apply_transform(self, builder, mat, vec, transpose = False):
        product = [None] * self._ss._s

        for i in range(self._ss._s):
            for j in range(self._ss._s):
                if not transpose:
                    _ = builder.fmul(mat[i][j], vec[j])
                else:
                    _ = builder.fmul(mat[j][i], vec[j])

                if product[i] is not None:
                    product[i] = builder.fadd(product[i], _)
                else:
                    product[i] = _
        return product


    def _llvm_load_shift(self, builder, srindex):
        if self.no_shift:
            return None

        subregion_Tt = builder.gep(self.xform_lut, [ir.IntType(32)(0), srindex])
        offset = self._ss._s * self._ss._s
        
        shift = []
        for index in range(self._ss._s):
            mp  = builder.gep(subregion_Tt, [ir.IntType(32)(0), ir.IntType(32)(int(offset + index))])
            sel = builder.load(mp)
            shift += [self._llvm_convert_to_internal_type(builder, sel)]

        return shift


    def _llvm_generate_outer_coset_index(self, builder, index, outer_coset_lut):
        """
        REturns a count of self._s variables which correspond to the 
        outer coset shift.
        """
        coset = builder.gep(outer_coset_lut, [ir.IntType(32)(0), index])
        shift = []
        for index in range(self._ss._s):
            mp = builder.gep(coset, [ir.IntType(32)(0), ir.IntType(32)(index)])
            sel = builder.load(mp)
            shift += [self._llvm_convert_to_internal_type(builder, sel)]
        return shift

    def _llvm_apply_shift(self, builder, vec_a, vec_b, subtract=False):
        if len(vec_a) != len(vec_b):
            raise ValueError("Lengths of vec_a and vec_b differ!")

        if subtract:
            return [builder.fsub(a, b, flags=['fast']) for (a,b) in zip(vec_a, vec_b)]
        return [builder.fadd(a, b, flags=['fast']) for (a,b) in zip(vec_a, vec_b)]

    def _llvm_generate_horner(self, builder, horner, variables, callback=None):
        """
        Generate code for horner factorized polynomials
        builder: LLVMlite builder instance -- this is where the code gets emitted to
        horner:  Factorization of the polynomial
        variables: The input variables
        callback: a function that maps a constant to a mem fetch. I don't think this 
        is used anymore.
        """
        def _build_power(power, constant):
            # TODO: be more clever with how we compute this
            # we could do some fast exponentiation thing, but
            # for low degree monomials, it might not be worth it
            if all([_==0 for _ in power]):
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
                        r = variables[v_idx]
                    else:
                        r = builder.fmul(r,  variables[v_idx], flags=['fast'])
            
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
            left  = factorization['left']
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
        if self._ss._s == 3:
            """
            These are some special cases which have easier to compute
            methods to determine which coset a point belongs to 
            """
            if self._ss._lattice.hash(True) == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                return self._llvm_generate_coset_index_bcc(builder, ls)

            if self._ss._lattice.hash(True) == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                return self._llvm_generate_coset_index_fcc(builder, ls)

            if self._ss._lattice.hash(True) == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                return ir.IntType(32)(0)

        if self._ss._lattice_matrix == matrix.identity(self._ss._s):
            return ir.IntType(32)(0)

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
        cs = builder.fptosi(ls[2], ir.IntType(32))
        cs = builder.and_(cs, ir.IntType(32)(1))
        return cs

    def _llvm_generate_coset_index_fcc(self, builder, ls):
        """
        Takes a lattice site, and determines which coset it belongs to. 
        This is the "coset index" (coset_index). It's an integer that
        indexes 

        This is a special case for the FCC lattice
        """
        z = builder.fptosi(ls[2], ir.IntType(32))
        z = builder.and_(z, ir.IntType(32)(1))

        y = builder.fptosi(ls[1], ir.IntType(32))
        y = builder.and_(y, ir.IntType(32)(1))

        z = builder.shl(z, ir.IntType(32)(1))
        return builder.or_(z,y)


    def _llvm_generate_coset_scale_bcc(self, builder, coset_idx, ls):
        return [builder.fdiv(builder.fsub(_, builder.uitofp(coset_idx, self._internal_type)), self._internal_type(2.)) for _ in ls]

    def _llvm_generate_coset_scale_fcc(self, builder, coset_idx, ls):
        y = builder.and_(coset_idx, ir.IntType(32)(1))

        z = builder.and_(coset_idx, ir.IntType(32)(2))
        z = builder.lshr(z, ir.IntType(32)(1))

        x = builder.xor(z, y)

        x = builder.uitofp(x, self._internal_type)
        y = builder.uitofp(y, self._internal_type)
        z = builder.uitofp(z, self._internal_type)

        return [builder.fdiv(builder.fsub(a, b, flags=['fast']), self._internal_type(2.0)) for (a,b) in zip(ls, [x,y,z])]

    def _llvm_generate_coset_scale_generic(self, builder, coset_idx, ls):
        """
        This takes the coset index, that is, which coset a point belongs to
        and shifts the input lattice site by that coset offset, and scales it
        to be unit cartesian. 

        REVIEW
        """

        D, v = self._ss.coset_vectors()

        if len(v) == 1: # Special case, we're on a scaled cartesian grid
            return [builder.fdiv(l,self._internal_type(s)) for l,s in zip(ls, D.diagonal())]

        # Otherwise, we need to determine the shift, first        
        shift = [self._internal_type(_) for _ in v]

        for idx, coset_shift in v[1:]:
            comp = builder.icmp_unsigned("==", coset_idx, ir.IntType(32)(idx+1))
            shift = [builder.select(comp, self._internal_type(new), old) for new,old in zip(coset_shift, shift)]
        
        ls = [builder.fsub(l,s, flags=['fast']) for l,s in zip(ls, shift)]
        return [builder.fdiv(l,self._internal_type(s)) for l,s in zip(ls, D.diagonal())]


    def _llvm_generate_coset_scale(self, builder, coset_idx, ls):
        if self._ss._s == 3:
            # BCC lattice indicator
            if self._ss._lattice.hash(True) == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                return self._llvm_generate_coset_scale_bcc(builder, coset_idx, ls)

            if self._ss._lattice.hash(True) == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                return self._llvm_generate_coset_scale_fcc(builder, coset_idx, ls)

            if self._ss._lattice.hash(True) == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                return ls

        if self._ss._lattice_matrix == matrix.identity(self._ss._s):
            return ls
        
        raise ValueError("Unsupported lattice -- no coset decomposition")

    def _ls_coset_index_bcc(self, ls):
        return int(ls[2]) % 2

    def _ls_coset_index_fcc(self, ls):
        x = (ls[0] % 2) << 1
        y = ls[1] % 2
        return x | y

    def _ls_coset_index(self, ls):
        if self._ss._s == 3:
            # BCC lattice indicator
            if self._ss._lattice.hash(True) == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
                return self._ls_coset_index_bcc(ls)

            if self._ss._lattice.hash(True) == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
                return self._ls_coset_index_fcc(ls)

            if self._ss._lattice.hash(True) == '569e5cf8c93c16fe49c1702381418d0d86455c36694e53eb1c5cc4954195010b':
                return int(0)

        if self._ss._lattice_matrix == matrix.identity(self._ss._s):
            return int(0)
        
        raise ValueError("Unsupported lattice -- no coset decomposition")

    def _llvm_generate_no_lerp_conv_sum_bpred(
        self, builder, xvars, ivars, xform, shift, srindex, polyindex, ls_lut, memfetch, 
        outer_coset_index=None,
        outer_coset_shifts=None):
        """
        This is a fairly lengthy method, let me break down the method call name
        _llvm_generate_(no_lerp)_conv_sum_(bpred)
        no_lerp => No linear interpolation
        bpred => uses branch predication

        basically breaks down the convulatuon sup into chunks (big outer loop)
        """

        # First, a not so trivial helper function
        def emit_horner(builder, xvars, polyindex, hregs, c_fetches, result_current):
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
                for refregion_idx,__ in enumerate(self._ss._ref_subregions):

                    # Generate the result 
                    chunk_result = self._llvm_generate_horner(builder, hregs[refregion_idx], xvars + c_fetches,  None)
                    _ = chunk_result

                    if len(self._ss._ref_subregions) > 1:
                        # select result if polyindex == refregion_idx else zero
                        # add result to previous result
                        _ = builder.icmp_unsigned("==", polyindex, ir.IntType(32)(refregion_idx))
                        _ = builder.select(_, chunk_result, self._internal_type(0.0))
                    result = builder.fadd(_,  result, flags=['fast'])
                return result

            else: 

                # Get a reference to this bb
                bb = builder.basic_block

                # First, create a bunch of basic blocks for each case
                bblist = [builder.append_basic_block(name=_label_suffix(bb.name, '.case%d' % idx)) 
                            for idx,_ in enumerate(self._ss._ref_subregions)]

                # This last bb is where all the above join back together
                bbend = builder.append_basic_block(name=_label_suffix(bb.name, '.endcase'))

                # Ugh, polyindex might be a low bit type, cast to
                # int32 type --- this should be sufficient 
                # I can't imagine splines with 2^32 sub-regions to be
                # practical
                pidx = builder.zext(polyindex, ir.IntType(int(32)))

                # Create the switch ... the default is the 0'th block
                c = builder.switch(pidx,  bblist[0])
                
                # Add the different cases (effectively: add jumps to different basic blocks)
                for refregion_idx,__ in enumerate(self._ss._ref_subregions[1:]):
                    c.add_case(ir.IntType(int(32))(refregion_idx + 1),  bblist[refregion_idx+1])

                # Generate code for the different blocks 
                outputs  = []
                for idx, block in enumerate(bblist):
                    with builder.goto_block(block):
                        chunk_result = self._llvm_generate_horner(builder, hregs[idx], xvars + c_fetches,  None)
                        outputs += [chunk_result]

                        # Branch back to a common point
                        builder.branch(bbend)

                # Resume serial branchless code
                builder.position_at_start(bbend)

                # Add a phi node to choose the correct value from the 
                # above blocks
                phi = builder.phi(self._internal_type)
                for bb,val in zip(bblist, outputs):
                    phi.add_incoming(val, bb)

                # Add the final result to the convolution sum
                result = builder.fadd(phi,  result_current, flags=['fast'])
                return result

        ls_start = 0
        result = self._internal_type(0.0)
        current_coset = (None, None, None)
        r0_lookups = self._ss._ref_subregions[0].get_optimal_nn_lookups()

        in_flight_q = []
        """
        Group all the memory reads in chunks of group_size, returns a list
        of elements, each of size group_size
        """
        for list_chunk in grouper(self.group_size, zip(*[_.get_optimal_nn_lookups() for _ in self._ss._ref_subregions])):

            # If we need to, generate a load for the transform registers,
            # then do it now
            if self.refetch_transform:
                xform = self._llvm_load_transform_matrix(builder, srindex, False)
                if not self.no_shift:
                    shift = self._llvm_load_shift(builder, srindex)

            #####################################
            # Generate the memory fetches
            #####################################
            c_fetches = []
            for idx, chunk in enumerate(list_chunk):
                # group_list pads, so ignore padding
                if chunk is None:
                    continue

                # Do memfetch from chunk 
                if len(self._ss._ref_subregions) == 1:
                    # If we only have one sub region, just encode the 
                    # lattice site fetches in code
                    ls = [ir.IntType(32)(_) for _ in chunk[0]]
                else:
                    # Otherwise we have to fetch them based on the polynomial
                    _ = builder.gep(ls_lut, [ir.IntType(32)(0), polyindex, ir.IntType(32)(idx + ls_start)])
                    ptrs = [builder.gep(_, [ir.IntType(32)(0), ir.IntType(32)(i)]) for i,__ in enumerate(chunk[0])]
                    ls = [builder.sext(builder.load(_), ir.IntType(32)) for _ in ptrs]

                ls = [self._llvm_convert_to_internal_type(builder, _) for _ in ls]
                # Apply the transform
                ls = self._llvm_apply_transform(builder, xform, ls)
                ls = self._llvm_apply_shift(builder, ls, shift)

                # If we have reflective symmetry, we need to apply that now
                if self._ss._reflective:
                    ls = [builder.fmul(_,  self._reflective[i], flags=['fast']) for i, _ in enumerate(ls)]

                # Finally, add in the reference lattice sites
                ls = [builder.fadd(_,  ivars[i], flags=['fast']) for i, _ in enumerate(ls)]

                if self.lookup_type == "ABSTRACT":
                    if outer_coset_index is None:
                        c_fetches += [builder.call(memfetch, [builder.fptosi(_, ir.IntType(32)) for _ in ls])]
                    else:
                        ls_shift = [builder.fadd(l,s, flags=['fast']) for l,s in zip(ls, outer_coset_shifts)]
                        c_fetches += [builder.call(memfetch, [builder.fptosi(_, ir.IntType(32)) for _ in ls_shift])]

                if self.lookup_type == "COSET":
                    if not self._subregion_consistency:
                        coset = self._llvm_generate_coset_index(builder, ls)
                        ls = self._llvm_generate_coset_scale(builder, coset, ls)
                        tmp = ([outer_coset_index] if outer_coset_index is not None else []) + [coset]
                        c_fetches += [builder.call(memfetch, tmp+[builder.fptosi(_, ir.IntType(32)) for _ in ls])]

                    else:
                        cidx = self._ls_coset_index(r0_lookups[ls_start + idx])
                        if current_coset[0] is None or cidx != current_coset[1]:
                            coset = self._llvm_generate_coset_index(builder, ls)
                            current_coset = (coset, cidx)
                        ls = self._llvm_generate_coset_scale(builder, current_coset[0], ls)
                        tmp = ([outer_coset_index] if outer_coset_index is not None else []) + [current_coset[0]]
                        c_fetches += [builder.call(memfetch, tmp+[builder.fptosi(_, ir.IntType(int(32))) for _ in ls])]
                        
                if self.lookup_type == "TEX":
                    cidx = self._ls_coset_index(r0_lookups[ls_start + idx])
                    if current_coset[0] is None or cidx != current_coset[1]:
                        coset = self._llvm_generate_coset_index(builder, ls)

                        if outer_coset_index is not None:
                            n_ocosets = len(self._ss.get_L_cosets())
                            offset = builder.mul(outer_coset_index, ir.IntType(32)(n_ocosets))
                            coset = builder.add(offset, coset)

                        _ = self._textures[0]
                        for idx_tex, tex in enumerate(self._textures[1:]):
                            idx_tex += 1
                            p = builder.icmp_unsigned("==", coset, ir.IntType(32)(idx_tex))
                            _ = builder.select(p, tex, _)

                        current_coset = (coset, cidx, _)

                    ls = self._llvm_generate_coset_scale(builder, current_coset[0], ls)

                    ls = [builder.fadd(_,  self._internal_type(0.5), flags=['fast']) for _ in ls]
                    ls = [builder.fptrunc(_, ir.FloatType()) for _ in ls]

                    if self._ss._s == 1:
                        mf = builder.call(self._tex1d, [current_coset[2]]+ls)
                    elif self._ss._s == 2:
                        mf = builder.call(self._tex2d, [current_coset[2]]+ls)
                    elif self._ss._s == 3:
                        mf = builder.call(self._tex3d, [current_coset[2]]+ls)
                    else:
                        raise ValueError("Dimension not supported")

                    c_fetches += [ mf]


            ls_start += len(list_chunk)
            list_chunk = zip(*[_ for _ in list_chunk if _ is not None])

            ########################################################
            # End of memfetch code
            ##################################################
            """
            Generate the horner factorizations, but don't emit code
            """
            hregs = []
            for refregion_idx, chunk in enumerate(list_chunk):
                c_lookups    = []
                c_polynomial = 0

                for cidx, ls in enumerate(chunk):
                    # Generate a memfetch at ls, and store it in the c_lookups
                    mf = var("c_%d" % cidx)
                    c_lookups    += [mf]
                    c_polynomial += mf * self._ss._ref_subregions[refregion_idx]._weighted_sum_terms[ls]

                # Factorize
                horner = horner_factor(c_polynomial, self._ss._s, c_lookups)
                hregs += [horner]


            if self.pipeline_depth == 0:
                if self.lookup_type == "TEX":
                    c_fetches = [builder.extract_value(_,0 ) for _ in c_fetches]
                c_fetches = [self._llvm_convert_to_internal_type(builder, _) for _ in c_fetches]

                result = emit_horner(builder, xvars, polyindex, hregs, c_fetches, result)
            else:
                if len(in_flight_q) >= self.pipeline_depth:
                    fetches, hr = in_flight_q.pop(0)

                    # Grab all the fetches
                    if self.lookup_type == "TEX":
                        fetches = [builder.extract_value(_,0 ) for _ in fetches]
                    fetches = [self._llvm_convert_to_internal_type(builder, _) for _ in fetches]

                    result = emit_horner(builder, xvars, polyindex, hr, fetches, result)
                # Now add the current fetch to the queue
                in_flight_q += [(c_fetches, hregs)]

                 
        # Empty out anything left in the queue
        for fetches, hr  in in_flight_q:
            if self.lookup_type == "TEX":
                fetches = [builder.extract_value(_,0 ) for _ in fetches]
            fetches = [self._llvm_convert_to_internal_type(builder, _) for _ in fetches]
            result = emit_horner(builder, xvars, polyindex, hr, fetches, result)

        return result    

    def _llvm_generate_lerp_conv_sum_bpred_consistent(self, builder, xvars, ivars, xform, shift, srindex, polyindex, ls_lut, memfetch):
        # This should be handled below, I think...
        pass

    def _llvm_generate_lerp_conv_sum_bpred_inconsistent(
        self, 
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

    def llvm_generate(self):

        self._llvm_generate_module()

        # Generate all lookup tables
        bsp_lut         = self._llvm_generate_bsp_lookup_table(self._llvm_module)
        xform_lut       = self._llvm_generate_xfrom_lookup_table(self._llvm_module)
        refpp_lut       = self._llvm_generate_ref_region_index_lookup_table(self._llvm_module)
        ls_lut          = self._llvm_generate_lattice_site_lookup_table(self._llvm_module)
        outer_coset_lut = self._llvm_generate_outer_coset_table(self._llvm_module)
        
        # Generate the function
        builder, xvars, memfetch = self._llvm_generate_function(self._llvm_module) 

        # Generate the standard preamble
        if len(self._ss.get_L_cosets()) == 1:
            xvars, ivars = self._llvm_generate_rho(builder, xvars)
            srindex      = self._llvm_generate_bsp_index(builder, xvars)
            polyindex    = self._llvm_generate_poly_index(builder, srindex, refpp_lut)

            mat          = self._llvm_load_transform_matrix(builder, srindex, False)
            shift        = self._llvm_load_shift(builder, srindex)

            if not self.no_shift:
                xvars    = self._llvm_apply_shift(builder, xvars, shift, True)   # Subtract
            xvars        = self._llvm_apply_transform(builder, mat, xvars, True) # Transform with transposed matrix

                 
            if self.lookup_type == "CPU_LIN" or self.lookup_type == "TEX_LIN":
                if self.branchless:
                    result = self._llvm_generate_lerp_conv_sum_bpred_inconsistent(builder, xvars, ivars, mat, shift, srindex, polyindex, ls_lut, memfetch)
                else:
                    raise ValueError("Branchy code not supported for lin yet")
            else:
                result = self._llvm_generate_no_lerp_conv_sum_bpred(builder, xvars, ivars, mat, shift, srindex, polyindex, ls_lut, memfetch)
            result = self._llvm_convert_to_output_type(builder, result)

        else:
            # Create a basic block for the loop
            loop_blk = self._function.append_basic_block("coset_loop")
            exit_blk = self._function.append_basic_block("exit")

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


            # Now we can emit the interpolation code
            shifts = self._llvm_generate_outer_coset_index(builder, loop_idx, outer_coset_lut)
            xvars_shifted = [builder.fsub(x,s) for x,s in zip(xvars,shifts)]

            xvars_shifted, ivars = self._llvm_generate_rho(builder, xvars_shifted)
            srindex      = self._llvm_generate_bsp_index(builder, xvars_shifted)
            polyindex    = self._llvm_generate_poly_index(builder, srindex, refpp_lut)

            mat          = self._llvm_load_transform_matrix(builder, srindex, False)
            shift        = self._llvm_load_shift(builder, srindex)

            if not self.no_shift:
                xvars_shifted    = self._llvm_apply_shift(builder, xvars_shifted, shift, True)   # Subtract
            xvars_shifted        = self._llvm_apply_transform(builder, mat, xvars_shifted, True) # Transform with transposed matrix

            #########3

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




        builder.ret(result)

        return str(self._llvm_module)

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