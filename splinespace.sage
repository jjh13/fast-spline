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

load("subregion.sage")
load("lattice.sage")
load("helpers.sage")

from os import path
import time
import sys
from sage.structure.element import is_Matrix

class SplineSpace:

    """ 
    Initialize spline space

    spline: A BoxSpline object, or a list of  (polyhedron, polynomial pairs)
    lattice: A lattice spec, either the abbreviation of the lattice, or the generating matrix
    kwargs:
    coset_matrix: the matrix G in the paper
    rho: Either a matrix that describes 
    """
    def __init__(self, spline, lattice, **kwargs):

        # Unpack kwargs
        coset_matrix= None
        if "coset_matrix" in kwargs:
            coset_matrix = kwargs["coset_matrix"]

        rho="indicator" 
        if "rho" in kwargs:
            rho = kwargs["rho"]

        rho_shift=None
        if "rho_shift" in kwargs:
            rho_shift = kwargs["rho"]

        self._allow_reflective = True
        if "enable_reflective" in kwargs:
            self._allow_reflective = kwargs['enable_reflective']

        self._id = None
        if "cache_id" in kwargs:
            self._id = kwargs["cache_id"]
        
        self._verbose=True
        if "verbose" in kwargs:
            self._verbose = kwargs["verbose"]

        # Setup the initial lattice and do some sanity checks
        if isinstance(lattice, str):
            lattice_matrix = lattice_from_str(lattice)
            if lattice_matrix is None:
                raise ValueError("Invalid lattice spec '%s'" % lattice)
            self._s = len(lattice_matrix.rows()[0])

        elif isinstance(lattice, IntegerLattice):
            lattice_matrix = matrix(ZZ,lattice.get_basis()).transpose()
            self._s        = lattice.dimension()

        elif isinstance(lattice, sage.matrix.matrix_integer_dense.Matrix_integer_dense):
            lattice_matrix = lattice
            self._s = len(lattice.rows()[0])
        else:
            raise ValueError("Invalid input lattice! Should be an integer matrix, IntegerLattice object, or Lattice Spec")

        # Ensure that the lattice matrix is square
        if lattice_matrix.dimensions()[0] != lattice_matrix.dimensions()[1]:
            raise ValueError("Input lattice matrix must be square")

        # Check to make sure that the coset matrix is valid
        if coset_matrix is None:
            coset_matrix = matrix.identity(ZZ, self._s)
        elif isinstance(coset_matrix, sage.matrix.matrix_integer_dense.Matrix_integer_dense):
            if coset_matrix.dimensions()[0] != coset_matrix.dimensions()[1]:
                raise ValueError("Coset matrix should be square")
            if coset_matrix.dimensions()[0] != lattice_matrix.dimensions()[0]:
                raise ValueError("Coset matrix should have the same dimension as the lattice")
        else:
            raise ValueError("Coset matrix should be a square integer matrix")

        # Construct the lattice objects
        orig_lattice = IntegerLattice(lattice_matrix)
        self._lattice_matrix = lattice_matrix*coset_matrix
        self._lattice        = IntegerLattice(self._lattice_matrix) 

        # Do the coset decomposition
        p = Polyhedron(vertices=self._lattice.get_parallelpiped(), base_ring=AA)

        pp_inv = matrix(self._lattice.get_basis()).transpose().inverse()

        cosets = []
        for pt in p.integral_points():
            if orig_lattice.is_lattice_site(pt) and not self._lattice.is_lattice_site(pt):
                if all([0 <= _ < 1 for _ in pp_inv*vector(pt)]):
                    cosets += [vector(pt)]

        self._cosets = [vector([0]*self._s)] + cosets

        """ 
        Setup the region of evaluation
        """
        self._rho_ppiped = matrix.identity(self._s)
        self._rho_shift  = vector([0]*self._s) if rho_shift is None else rho_shift
        self._rho_type   = "PP"

        if isinstance(rho, str):
            if rho.lower() == "indicator":
                self._rho_type = "IND"
            else:
                raise ValueError("Invalid $\\rho$ specified")
        elif is_Matrix(rho): # isinstance(rho, sage.matrix.matrix.Matrix):
            self._rho_ppiped = rho
        else:
            raise ValueError("Invalid $\\rho$ specified --- should be a matrix or 'indicator'")
        self._build_roe()


        self._xvars = [var('x_%d' % _) for _ in range(self._s)]

        """
        Initialize the rest of the spline space, box splines have a slightly less involved setup
        but collections of peicewise polynomials need extra work
        """
        if isinstance(spline, BoxSpline):
            self.__init_box_spline__(spline=spline, lattice=lattice, rho=rho, rho_shift=rho_shift)
        elif isinstance(spline, (list, tuple)):
            self.__init_pp_spline__(spline=spline, lattice=lattice, rho=rho, rho_shift=rho_shift)
        else:
            raise ValueError("Invalid spline specified --- should be a BoxSpline object or a list of (polyhedron, polynomial) pairs")

    def get_L_cosets(self):
        """
        returns the cosets of the lattice G*L
        """
        return self._cosets

    def tick(self, indicator='.',end=' '):
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

    def _check_cache(self, name):
        """ 
        Checks to see if an object exists in the cache, and 
        loads it. If the user hasn't specified an identity for 
        the space, then this function always returns None. 
        """
        if self._id is None:
            return 

        try:
            _path = "cache/ss/%s_%s" % (self._id, name)

            if path.exists(_path+".sobj"):
                self.log("Loading %s from cache..." % (name), end = '')
            obj = load("cache/ss/%s_%s" % (self._id, name))
            self.log("Ok!")
            return obj
        except:
            self.log("Cache miss or load fail at %s!" % name)
            
        return None

    def _cache_obj(self, obj, name):
        if self._id is None:
            return 

        if not os.path.exists('cache'):
            os.makedirs('cache')

        if not os.path.exists('cache/ss'):
            os.makedirs('cache/ss')

        self.log("Caching %s..." % name, end='')
        t0 = time.time()
        _path = "cache/ss/%s_%s" % (self._id, name)
        save(obj, _path)
        t1 = time.time()
        self.log("done in %s seconds!" % (t1-t0))

    def __init_box_spline__(self, spline, lattice, rho, rho_shift):
        """
        This sub-routine implements all of the analysis 
        in part-I of the paper --- it's specialized to 
        box-splines, since the box-spline class has some
        helper routines that expidite the process.
        """

        # We need the support of the spline for various things, grab it
        self._support = spline.get_polyhedron()

        # then do the symmetry analysis
        self._analyze_symmetry_box_spline(spline)

        # Collect all the planes of the finer mesh
        # This is a tuple of (normal, [list of d values])
        planes = spline._get_grouped_planes_for_eval() 

        # Determine which lattice shifts contribute to the ROE
        mink = self._roe.minkowski_sum(self._support)

        # Filter out points not on the lattice (and not in the mink supp)
        pts = [pt for pt in mink.integral_points() if
            self._lattice.is_lattice_site(vector(pt)) 
        ]

        # Decompose the region of evaluation into its sub regions
        self._decompose_sub_regions(pts, planes)

        # Distribute the polynomial to the subregions 
        self._pp_regions = [(_['polyhedron'], _['polynomial']) for _ in spline.get_pp_regions()]

        self._distribute_polynomials(pts, self._pp_regions)

        self._subregion_symmetry_analysis()
        self._collect_sr_planes()
        self._pack_regions()

    def __init_pp_spline__(self, spline, lattice, rho, rho_shift):

        # We need the support of the spline for various things, construct it        
        self._support = Polyhedron(vertices=sum([list(polyhedron.vertices()) for polyhedron,_ in spline],[]))


        # then do the symmetry analysis
        self._analyze_symmetry_pp(spline)

        # Concatenate all the planes together
        planes = []
        for polyhedron, _ in spline:
            for plane in polyhedron.Hrepresentation():
                d = plane[0]
                normal = vector(plane[1:])

                planes += [(d,normal)]

        # Group and consistently orient them
        plane_set = {}
        for d, n in planes:
            tn = tuple(n)
            if tn not in plane_set:
                tn = tuple(-vector(n))

                if tn not in plane_set:
                    plane_set[tuple(n)] = set([d])
                else:
                    plane_set[tn] = plane_set[tn].union(set([-d]))
            else:
                plane_set[tn] = plane_set[tn].union(set([d]))

        for p in plane_set:
            plane_set[p] = sorted(list(plane_set[p]))

        grouped_planes =  list([(n,plane_set[n]) for n in plane_set])
        self._gp = grouped_planes

        # Determine which lattice shifts contribute to the ROE
        mink = self._roe.minkowski_sum(self._support)

        pts = [pt for pt in mink.integral_points() if
            self._lattice.is_lattice_site(vector(pt))
        ]

        self._decompose_sub_regions(pts, grouped_planes)
        self._distribute_polynomials(pts, spline)
        self._subregion_symmetry_analysis()
        self._collect_sr_planes()
        self._pack_regions()

    def _rho_check(self):
        """
            Checks that the region of evaluation is valid.
            Validity means that it tesselates space if 
            det|spline lattice| = 
        """
        pass
        
    def _build_roe(self):
        """
            Builds the explicit region of evaluation for a 
            spline space. This is a polyhedron shifted by roe_shift
        """
        self.tick("Builing region of evaluation...")
        self._roe = None
        if self._rho_type == "PP":
            self._roe = Polyhedron(vertices=build_ppiped(self._rho_ppiped, False)[0])
        elif self._rho_type == "IND":
            # pass
            self._roe = self._lattice.voronoi_region()

        self._roe = shift_polyhedron(self._roe, self._rho_shift)
        self.log("done!")
        return self._roe

    def _decompose_sub_regions(self, lsites, gplanes):
        """
        This takes the region of evaluation, and breaks it into 
        the sub regions of evaluation
        """
        regions = [self._roe]
        new_regions = []
        aug_planes = []

        aug_planes = self._check_cache("aug_planes")
        if aug_planes is None:
            aug_planes = []
            # Augment the grouped planes with the shifted lattice sites
            self.log("Decomposing sub-regions - augmenting planes to create finer mesh...")

            for normal, dlist in gplanes:
                new_d_set = []

                for site in lsites:
                    offset = vector(normal)*vector(site)

                    for d in dlist:
                        new_d_set += [offset + d]

                aug_planes += [(normal, sorted(list(set(new_d_set))))]
            self._cache_obj(aug_planes, "aug_planes")

        self.aug_planes = aug_planes
        self._subregions = self._check_cache("sub_region1")
        if self._subregions is None: 
            # Dice the region of evaluation by all the planes above
            self.log("Decomposing sub-regions - slicing ROE by fine mesh...")

            for idx,(normal, dlist) in enumerate(aug_planes):
                for d in dlist:
                    new_regions = []
                    for r in regions:
                        A,B = split_polyhedron(r, normal, d)

                        if A is None or B is None:
                            new_regions += [r]
                        else:
                            new_regions += [A,B]
                    regions = new_regions
                    self.tick()
                self.tick('%d/%d' % (idx+1, len(aug_planes)))

                # if self._verbose:
                #     print ".",
                # for region in regions:
                #     print ",",
                #     new_regions += split_polyhedron_by_group(region, normal, dlist)
                # regions     = new_regions
                # new_regions = []
            self.log("done!")
            self.log("Constructing sub-regions...")

            self._subregions = [SubRegion(_, self) for _ in regions]
            self._cache_obj(self._subregions, "sub_region1")

        self.log("Decomposing sub regions - checking for reflective symmetry in ROE...")

        # First we check to see if we have reflective 
        # Symmetry in our symmetry group, this lets us 
        # reduce the number of reference regions
        self._reflective = self._allow_reflective
        for i in range(self._s):
            X = matrix.identity(self._s)
            X[i,i] = -1
            roe_prime = Polyhedron(vertices=[X*vector(v) for v in self._roe.vertices()], base_ring=AA)

            if X not in self._symmetry or self._roe != roe_prime:
                self._reflective = False
                break

        if  self._reflective:
            self.log("SplineSpace has reflective symmetry!")
        else:
            self.log("SplineSpace doesn't have reflective symmetry :<")

        # If we only have one sub region, then we can turn off
        # the reflective symmetry hack
        if len(self._subregions) == 1:
            self._reflective = False

        # If we have reflective symmetry, we can prune away regions
        # that aren't strictly positive
        if self._reflective:

            srp = self._check_cache("sub_region_pruned")
            if srp is None:

                self.log("Subregion symmetry analysis - Pruning redundant regions")

                positive = Polyhedron(ieqs=[ 
                                tuple([1 if i - 1 == j else 0 for i in range(self._s+1)])
                                for j in range(self._s)
                            ], base_ring=AA)

                new_sregions = []
                for sr in self._subregions:
                    intersect = positive.intersection(sr._region)
                    if intersect.volume() == 0 or intersect.dim() < self._s:
                        continue
                    new_sregions += [sr]
                
                self._subregions = new_sregions
                self._cache_obj(new_sregions, "sub_region_pruned")
            else:
                self._subregions = srp


        return self._subregions

    def get_subregions(self):
        return self._subregions

    def _analyze_symmetry_box_spline(self, bs):
        self._symmetry = self._check_cache("symmetry")
        if self._symmetry is not None:
            return
        self._symmetry = []

        # Center the support, just in case
        center = self._support.center()
        supp = Polyhedron(vertices=[vector(_) - center for _ in self._support.vertices()])

        self.log("Symmetry analysis - signed permutation matrices...", end='')

        #  check for  rotational symmetry 
        for X in signed_permutation_matrices(self._s):
            supp_t = Polyhedron(vertices=[X*vector(_)  for _ in supp.vertices()], base_ring=AA)

            if supp == supp_t:
                self._symmetry += [X]
                self.tick('↻')
            else:
                self.tick('x')

        self.log('done!')
        self.log("Symmetry analysis - reflections...", end='') 

        # Next, we check if we have reflective symmetry along the knot planes
        # and Xi direction vectors
        planes = [vector(_).normalized() for _ in bs.Xi_]
 
        H1, H2 = bs.calc_knotplanes()
        for d, normal in H1 + H2:
            planes += [normal]

        planes = set([tuple(_) for _ in planes])

        for n in planes:
            P = matrix.identity(self._s) - 2*matrix(n).transpose()*matrix(n)

            if P in self._symmetry:
                self.tick("✔")
                continue

            supp_t = Polyhedron(vertices=[P*vector(AA,_)  for _ in supp.vertices()], base_ring=AA)

            if supp == supp_t:
                self._symmetry += [P]
                self.tick("😊")
            else:
                self.tick("☠")

        self.log("done!")
        self._cache_obj(self._symmetry, "symmetry")
        return 

        
    def _analyze_symmetry_pp(self, spline):
        self._symmetry = self._check_cache("symmetry")
        if self._symmetry is not None:
            return

        self._symmetry = []

        # Center the support, just in case
        center = self._support.center()
        supp = Polyhedron(vertices=[vector(_) - center for _ in self._support.vertices()])

        self.log("Symmetry analysis - signed permutation matrices...", end='') 


        #  check for  rotational symmetry 
        for X in signed_permutation_matrices(self._s):
            supp_t = Polyhedron(vertices=[X*vector(_)  for _ in supp.vertices()])

            if supp == supp_t:
                self._symmetry += [X]
                self.tick("↻")
            else:
                self.tick('x')

        self.log("done!")
        self.log("Symmetry analysis - reflections...", end='')

        planes = {}

        for polyhedron, _ in spline:
            for plane in polyhedron.Hrepresentation():
                d = plane[0]
                normal = vector(plane[1:].normalized())

                if tuple(normal) not in planes:
                    planes[tuple(normal)] = 0
                planes[tuple(normal)] += 1

        planes = planes.keys()

        for n in planes:
            P = matrix.identity(self._s) - 2*matrix(n).transpose()*matrix(n)

            if P in self._symmetry:
                self.tick('✔')
                continue

            supp_t = Polyhedron(vertices=[P*vector(_)  for _ in supp.vertices()])

            if supp == supp_t:
                self._symmetry += [P]
                self.tick('😊')
            else:
                self.tick('☠')

        self.log('done!')
        self._cache_obj(self._symmetry, "symmetry")
        return 


    def _distribute_polynomials(self, ls, pp_regions):
        # Name lattice sites for sub regions
        index = 0

        chk = self._check_cache("site_mapping")
        if chk is None:
            cvar_to_site = {}
            site_to_cvar = {}

            self.log("Unrolling convolution sum - naming lattice sites... ")

            for site in ls:
                p = shift_polyhedron(self._support, site)
                p = p.intersection(self._roe)
                if p.dim() < self._s or p.volume() == 0:
                    continue
                    
                site_to_cvar[tuple(site)] = var('c_%d' % index)
                cvar_to_site[var('c_%d' % index)] = tuple(site)
                index += 1
                
            self._cache_obj([cvar_to_site, site_to_cvar], "site_mapping")
        else:
            cvar_to_site, site_to_cvar = chk
        self._site_to_cvar = site_to_cvar
        self._cvar_to_site = cvar_to_site



        chk = self._check_cache("distributed_pp")
        if chk is None:

            # Distribute weighted polynomials
            regions = pp_regions
            hits = 0

            # We use polyhedra over the RDF to speed things up here ---
            # we only need intersection tests, everything else doesn't really matter
            approx_roe = approx_polyhedron(self._roe)
            self.log("Unrolling convolution sum - distributing polynomials... ", end='')
            self._dbg = ls

            for s in site_to_cvar:
                self.tick()
                for polyhedron, pp in pp_regions:
                    polyhedron = shift_polyhedron(polyhedron,s)
                    # print "-",

                    # If the ROE touches the shifted polynomial from the spline 
                    # if polyhedron.intersection(approx_roe).volume() > 0:
                    if polyhedron.center() in self._roe or any([_ in approx_roe for _ in polyhedron.vertices()]):

                        # The we should distribute it into the sub regions
                        for idx, subregion in enumerate(self._subregions):
                            # print "+",

                            # If we have any overlap between the subregion
                            if subregion._region.center() in polyhedron:
                            # if subregion.intersection(polyhedron).volume() > 0:
                                polynomial = pp.subs({
                                        var('x_%d'%i): var('x_%d'%i) - s[i] for i in range(self._s)
                                    })
                                subregion.distribute(site_to_cvar[tuple(s)], tuple(s), polynomial)
            self._cache_obj(self._subregions, "distributed_pp")
            self.log("done!")
        else:
            self._subregions = chk




    def _subregion_symmetry_analysis(self):
        import itertools

        check = self._check_cache("ref_subregions")
        if check is not None:
            self._ref_subregions = check[0]
            self._subregions = check[1]

        else:
            # Analyze the symmetry of each region with respect to the set of 
            # reference regions
            self._ref_subregions = [self._subregions[0]]

            
            pss = self._check_cache("partial_sym_subregions")
            ignore = -1

            if pss is not None:
                ignore = pss[0]
                self._subregions = pss[1]
                self._ref_subregions = pss[2]

            self.log("Subregion symmetry analysis - cataloging reference regions...", end='')

            for sidx, sregion in enumerate(self._subregions):
                if sidx < ignore + 1:
                    self.tick("C")
                    continue

                best_transform = None
                bridx = 0 
                ignore_shift_search = True

                for (ridx, ref_sr), X  in itertools.product(enumerate(self._ref_subregions), self._symmetry):
                    transform = self._compare_regions_sym(ref_sr, sregion, X)

                    if transform is not None:
                        X, shift = transform
                        best_transform = transform
                        bridx = ridx

                        if all([_==0 for _ in shift]) or ignore_shift_search:
                            break

                # if we found a transform, catalogue it
                if best_transform is not None:
                    sregion._ref_region = bridx
                    sregion._transform = best_transform
                    self.tick('↷')
                    if sidx % 5 == 0:
                        self.tick("(%d/%d)" % (sidx, len(self._subregions)))
                else:
                    # Otherwise add new reference region with 
                    # an identity transform
                    sregion._ref_region = len(self._ref_subregions)
                    sregion._transform = matrix.identity(self._s), vector([0]*self._s)
                    self._ref_subregions += [sregion]
                    if self._verbose:
                        self.tick('⌧')
                        self.tick("(%d/%d)" % (sidx, len(self._subregions)))

                if sidx % 3 == 0 and sidx > 0:
                    self._cache_obj([sidx, self._subregions, self._ref_subregions], "partial_sym_subregions")


            self.log("done!")
            self._cache_obj([self._ref_subregions, self._subregions], "ref_subregions")

    def _compare_regions_sym(self, A, B, X):
        import itertools

        if not isinstance(A, SubRegion) or not isinstance(B, SubRegion):
            raise ValueError("A and B need to be SubRegion objects")

        if X.det() == 0:
            raise ValueError("X must be a transformation matrix")

        if A.number_of_pts_per_reconstruction() != B.number_of_pts_per_reconstruction():
            self.log("Warning: Two regions had a different amount of points per reconstruction, this shouldn't happen, your approximation space is invalid?")
            return None

        A_center = vector(A.get_center())
        B_center = vector(B.get_center())


        # Construct the transforms
        xt_a = vector(self._xvars)
        xt_b = X*(vector(self._xvars) - A_center ) + B_center

        G = BipartiteGraph()
        for _ in range(A.number_of_pts_per_reconstruction()):
            G.add_vertex("A_%d" % _, left=True)
            G.add_vertex("B_%d" % _, right=True)

        A_list = A.get_ws_list()
        B_list = B.get_ws_list()

        Xi = X.inverse()
        A_list = [(ls_a, pp_a.subs({a:b for a, b in zip(self._xvars, xt_a)}).expand()) for (ls_a, pp_a) in A_list]
        B_list = [(ls_b, pp_b.subs({a:b for a, b in zip(self._xvars, xt_b)}).expand()) for (ls_b, pp_b) in B_list]

        # Check all pairs of terms in each weighted sum form
        for (i, (ls_a, a_polynomial)), (j, (ls_b, b_polynomial)) in itertools.product(enumerate(A_list), enumerate(B_list)):
            if tuple(Xi*(vector(ls_b) - B_center) + A_center) == ls_a:
                if (a_polynomial.expand()-b_polynomial.expand()).expand() == 0:
                    G.add_edge("A_%d" % i, "B_%d" % j)

        # Check for a perfect matching
        if int(G.matching(True)) != int(A.number_of_pts_per_reconstruction()):
            return None

        # Check if the transformation is rigid
        sd = {}
        t = B_center - X*A_center
        for ls_a,pp_a in A.get_ws_list():
            sd[tuple(ls_a)] = True

        for ls_b,pp_a in B.get_ws_list():
            ls_a = X*vector(ls_b) + t
            if tuple(ls_a) not in sd:
                return None

        # return the transformation
        return X, t

    def _collect_sr_planes(self):
        self.log("Trimming ROE to only relevant subregions...")
        roe_trimmed = Polyhedron(vertices=sum([list(_._region.vertices()) for _ in self._subregions],[]), base_ring=AA)

        # Collect all the planes 
        pdict = {}

        # Check the cache first
        self._plist = self._check_cache("plist")
        if self._plist is not None:
            return self._plist
        
        self.log("Collecting planes for evaluation...")
        
        candidates = []
        # Go thru all sub regions, and collect the
        # planes that bound the region --- we'll remove duplicates
        # at a later step
        for subregion in self._subregions:
            for plane in subregion._region.Hrepresentation():
                d = plane[0]
                normal = plane[1:]

                A,B = split_polyhedron(roe_trimmed, normal, -d)
                
                # Discard planes on the exterior of the ROE
                if A is None or B is None:
                    continue

                candidates += [(tuple(normal),d)]

        # remove obvious duplicates
        candidates = list(set(candidates))

        pdict = {}
        # Iterate thru the planes we collected, 
        # and check to see if any have the same direction
        # then check for duplicates
        for idx,(pa,d) in enumerate(candidates):
            found = False
            for n in pdict:
                nv  = vector(n).normalized()
                pav = vector(pa).normalized()
                dp = nv*pav
                alpha = vector(pa)*vector(n)/(vector(n)*vector(n))

                if abs(dp) == 1:
                    pdict[n] += [d/alpha]
                    # We have the same direction
                    found = True
                    break

            if not found:
                pdict[tuple(pa)] = [d]

        for n in pdict:
            pdict[n] = list(set(pdict[n]))

        # Try to orient planes so that their distance values are all positive
        pprime = {}
        self.log("Reorienting planes...")
        for k in pdict:
            if all([_ <= 0 for _ in pdict[k]]):
                pprime[tuple([-_ for _ in k])] = sorted([- _ for  _ in pdict[k]])
            else:
                pprime[k] = sorted(pdict[k])
        
        # Collect all the planes in a list
        plist = []
        for k in pprime:
            for d in pprime[k]:
                l = vector(k).length()
                plist+= [(vector(k)/l,d/l)]

        # Cache n' return         
        self._plist = plist
        self._cache_obj(plist, "plist")

        return self._plist


    def _pack_regions(self):
        """
        Pack all regions into a BSP index
        """

        # First, look for any regions that share the EXACT same 
        # polynomial --- this means they could have been joined,
        # however, it's easier to run them thru the pipeline if they're
        # split.
        # When we build the index, if we have identical regions, then we just 
        # reference the original index. 
        aliased_sregions = {}
        sregions = []
        redundant_sregions = []
        rmap = {}

        self.log("Reording regions for redundancy removal...", end='')

        # A redundant region is a region whose function and lookups
        # are equivalent to some (already seen) subregion. Really,
        # These two regions just need to share the same polynomial reference
        for (i, sr_a) in enumerate(self._subregions):
            if sr_a._polynomial.expand() not in aliased_sregions:
                self.tick('+')
                sregions += [sr_a]
                aliased_sregions[sr_a._polynomial.expand()] = True
            else:
                self.tick('-')
                redundant_sregions += [sr_a]

        # Save the re-ordered sub-regions
        self._subregions = sregions + redundant_sregions
        aliased_sregions = {}

        # Create a map from redundant regions to sregions
        for idx, sr in enumerate(sregions):
            rmap[idx] = idx
            aliased_sregions[sr._polynomial.expand()] = idx

        for idx, sr in enumerate(redundant_sregions):
            idx += len(sregions)
            rmap[idx] = aliased_sregions[sr._polynomial.expand()]

        self.log('done!')
        self.log('Creating BSP index...', end='')

        packed_regions = []
        for ridx, r in enumerate(self._subregions):
            index = []

            for n,d in self._plist:

                A,B = split_polyhedron(r._region, n, -d)

                if A is None or B is None:
                    index += ['1' if vector(n)*(r._region.center()) >= -d else '0']
                else:
                    index += ['X']

            packed_regions += [(index, rmap[ridx])]

            self.log(''.join(index), end=' ')

        region_index = [-1 for _ in range(2**len(self._plist))]
        nz_indices = 0

        for bits, rindex in packed_regions:
            indices = [0]
            for bit_index, bt in enumerate(reversed(bits)):
                if bt == '1':
                    new_indices = [_ | (1 << bit_index) for _ in indices]
                    indices = new_indices
                elif bt == 'X':
                    new_indices = [_ | (1 << bit_index) for _ in indices]
                    indices = new_indices + indices
            for idx in indices:
                # if region_index[idx] != -1:
                    # print "WHAT?", region_index[idx], rindex
                region_index[idx] = rindex
                nz_indices += 1
                
        self.log('done!')
        self.log('Refining index...', end='')

        # Things can get weird around knot planes,
        # so check to make sure they're being mapped to
        # a region
        ticker = 0
        check = []
        for ridx, region in enumerate(self._subregions):            
            for d in reversed(range(1, self._s)):
                for a in region._region.faces(d):
                    v = a.vertices()
                    c = sum([vector(_) for _ in v])/len(v)
                    check += [c]

        verts = []

        for ridx, region in enumerate(self._subregions):
            verts += list(region._region.vertices())

        check = list(set([tuple(_) for _ in check]))
        verts = list(set([tuple(_) for _ in verts]))
        # print region_index
        for x in check + verts:
            # print x
            index = 0
            for i, (n,d) in enumerate(reversed(self._plist)):
                index |= (1 if vector(n)*vector(x) >= -d else 0) << i
                
            if region_index[index] == -1:
                self.tick('+')
                region_index[index] = rmap[ridx]
            else:
                if ticker % 4 == 0:
                    self.tick('-')
                ticker += 1
        # The initial BSP lookup table can be quite big, so we try to pack 
        # it mod some number

        self.log('done!')
        self.log('Packing BSP index...', end='')

        # brute force search 
        for limit in range(max(nz_indices-1, 1), 2^len(self._plist) + 1):
            packed_rindex = [None]*limit
            failed = False
            for index,rindex in enumerate(region_index):
                pindex = index % limit
                if packed_rindex[pindex] is None and rindex != -1:
                    packed_rindex[pindex] = rindex
                elif rindex != -1:
                    if packed_rindex[pindex] != rindex:
                        failed = True
                        break
                else:
                    packed_rindex[pindex] = 0
            if failed:
                self.tick("(%d, ✖)" % limit)
            else:
                self.log("✓ success with modulus %d, finished packing BSP!" % limit)
                break

        self._non_redundant_srgn_cnt = len(sregions)
        self._sregion_modulus = limit
        self._sregion_index = packed_rindex

    def eval_test(self, x, mmap=None, verbose=False, slow=True):
        x_orig = vector(x)
        if mmap == None:
            mmap = lambda __: 1 if all([_ == 0 for _ in __]) else 0
        
        cumsum = 0
        for shift in self._cosets:
            x = x_orig - vector(shift)

            # Map to the region of evaluation
            ls_ref = self.rho(x)
            # if verbose:
            #     print ls_ref
            x = vector(x) - vector(ls_ref)
            # if verbose:
            #     print x

            # Exploit the reflective symmetry, if we have it
            if self._reflective:
                Xr = matrix.diagonal([1 if _ >= 0 else -1 for _ in x])
            else:
                Xr = matrix.identity(self._s)
            x = Xr * x
            
            if not slow:
                # Do the BSP test to figure out our primary index
                index = 0
                for i, (n,d) in enumerate(reversed(self._plist)):
                    index |= (1 if vector(n)*vector(x) >= -d else 0) << i


                index %= self._sregion_modulus

                # Use that index to reference the correct sub_region
                subregion = self._sregion_index[index]

                # grab the proper region
                subregion = self._subregions[subregion]
            else:
                for idx, subregion in enumerate(self._subregions):
                    if x in subregion._region:
                        break

                # if verbose:
                #     print idx
                result = 0
                for ls, pp in subregion.get_ws_list():
                    pw = pp.subs({a:b for a,b in zip(self._xvars, x)})
                    p = Xr*vector(ls) + vector(ls_ref)
                    c = mmap(p + vector(shift))
                    # if verbose:
                    #     print c, ls, ls_ref,p
                    result +=  c * pw
                cumsum += result
                continue


            refregion = self._ref_subregions[subregion._ref_region]
            T, t = subregion._transform
            # if verbose:
            #     print T,t
            x_prime = T.transpose()*(vector(x) - t) 
            
            result = 0
            for ls, pp in refregion.get_ws_list():
                # Transform the points to the reference region
                ls = T*vector(ls) + t
                
                # Evaluate within the reference region
                pw = pp.subs({a:b for a,b in zip(self._xvars, x_prime)})
                p = Xr*ls + vector(ls_ref)
                c = mmap(p + vector(shift))
                # if verbose:
                #     print c, ls, ls_ref,p
                result +=  c * pw
            cumsum += result
        return cumsum

    def rho(self, x):
        x = vector(x)
        x = x - self._rho_shift

        if self._rho_type == "PP":
            return self._rho_ppiped *vector([floor(_) for _ in self._rho_ppiped.inverse()*x])

        # Handle any cartesian lattice
        if self._lattice_matrix == matrix.identity(self._s):
            return tuple([int(round(_)) for _ in x])

        # BCC lattice
        elif self._lattice.hash(True) == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
            return rho_bcc(x[0], x[1], x[2])

        # FCC lattice
        elif self._lattice.hash(True) == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
            return rho_fcc(x[0], x[1], x[2])

        # If we don't have a fast method for rho, we simply
        # find all the nearest lattice sites on the different cartesian cosets

        D, offsets = self._lattice.coset_structure()
        pts = []
        # D = D.transpose()
        for i, offset in enumerate(offsets):
            xt = x - vector(offset)
            xt = D.solve_right(xt)
            xt = (D*vector([round(_) for _ in xt])) + vector(offset)
            pts += [xt]

        # Then take the nearest of those
        dst = lambda a,b: (vector(a) - vector(b))*(vector(a) - vector(b))
        best  = pts[0]
        bestd = dst(best, x)
        for pt in pts[1:]:
            d = dst(pt, x)
            if d < bestd:
                best = pt
                bestd = d

        return best


    def coset_vectors(self):
        """
         Returns the cartesian coset structure of this spline space's lattice
         in the form (D, [shifts]).
        """


        # We have some special cases here, since we'd like to ensure that 
        # the coset vectors have a particular order. Additionally, I wrote
        # this prior to the general lattice code, so I don't want to break
        # something that depends on this (and I think the codegen
        # depends on the order here)
        if self._lattice_matrix == matrix.identity(self._s):
            return (matrix.identity(self._s), [vector([0,0,0])])

        # BCC lattice
        elif self._lattice.hash(True) == '3791d26e4028bb9279729f6280b07b615a80c12fe59a6f247e0f48cf28e2aee9':
            return (2*matrix.identity(self._s), [vector([0,0,0]), vector([1,1,1])])

        # FCC lattice
        elif self._lattice.hash(True) == 'c0050f04a3529c5bb8591d5e3fdc93fb6d5bd2a0c9597cf9c5caf194d1520b7e':
            return (2*matrix.identity(self._s), [vector((0,0,0)), vector((1,1,0)), vector([1,0,1]), vector((0,1,1))])

        return self._lattice.coset_structure()
    """
    Handy getters
    """
    def getRegionOfEvaluation(self):
        return self._roe

    def getLatticeMatrix(self):
        return self._lattice_matrix