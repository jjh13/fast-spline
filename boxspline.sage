#!/usr/bin/env python
"""
boxspline.sage

This is the bread & butter of the paper "Fast and exact evaluation of box
splines via the PP-form". It's where the actual set decomposition code is,
conveniently supplied in a ``nice'' BoxSpline class. See the README for
a few quick examples on how to use this class.

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

from sage.all import *
from itertools import product, combinations, chain
from operator import mul
import logging, sys

load("./helpers.sage")

class BoxSpline:
    def __init__(self, Xi, centered=True, shift = None, weight=1,verbose=False):
        """
        Sets the default parameters for a box spline.

        Arguments:
        Xi -- A list of direction tuples that form the direction matrix, for
        example:
            Xi = [(1,0),(0,1),(1,1)] # Courant Element
            Xi = [(1,0),(0,1)]*4 # Cubic 2d-tensor product b-spline
            Xi = [(1,)]*4 #  1d-tensor cubic b-spline
            Xi = [(-1,1,1),(1,-1,1),(1,1,-1),(-1,-1,-1)] # BCC Linear spline
        centered -- True or False, automatically shift the box spline so that
            its center is at the origin.

        shift -- A vector by which to shift the spline, example: vector([-1,-1])

        weight -- This scales the internal polynomials so that the output values
            are scaled by this factor.
            
        verbose -- Print status messages as the decomposition proceeds

        """
        # Get simple parameters for the box spline
        self.s_ = len(Xi[0])
        self.n_ = len(Xi)

        # Define the fourier and spatial variables for this BS
        self.w_ = [var('w_%d' % i) for i in range(self.n_)]
        self.x_ = [var('x_%d' % i) for i in range(self.s_)]

        factor = 0 if not centered else 1
        self.c_xi = factor * sum([vector(x) for x in Xi])/2
        self.centered = centered

        if shift:
            self.c_xi -= vector(shift)

        self.Xi_ = Xi[:]
        self.kerXi_ = list(matrix(SR,Xi).kernel().basis())
        self.weight = weight

        # Setup caches for each of these objects
        self.greens_cache = None
        self.differ_cache = None
        self.polytope_cache = None
        self.cached_regions = None
        self.polyhedron_cache = None
        self.polyterm_cache = None
        self.polytermx_cache = None
        self.gt_cache = None
        self._verbose = verbose
        
    def log(self, message):
        if self._verbose:
            print(message)

    def warmup(self):
        """
        Builds the internal caches for many of the operations used internally.
        """
        _ = self.get_difference_operator()
        _ = self.get_polyhedron()
        _ = self.decompose_greens()
        _ = self._get_grouped_planes_for_eval()
        _ = self.get_polyterms_w_xform()

    def cleanup(self):
        """
        Forcibly delete the caches. Note that this doesn't reclaim memory, since
        most of our computations are via sage varaibles, which consistently
        leak memory : /. This IS useful if you're monkey patching or poking
        about in the internals of the class, and want to invalidate the cached
        objects.
        """
        del self.greens_cache
        del self.differ_cache
        del self.polytope_cache
        del self.cached_regions

    def get_fourier_form(self):
        """
        Returns the Fourier form of the box spline. This is a function of the
        variables w_0, w_1, ... w_n. Where n is the dimension of the input
        space.
        """
        w = vector([var('w_%d' % i) for i in range(self.s_)])
        phi = prod([
                    (1 - exp(-I* w.dot_product(vector(xi))))/(I*w.dot_product(vector(xi)))
                    for xi in self.Xi_
        ])
        ret =  phi * exp(I * self.c_xi.dot_product(w)) * self.weight
        return ret

    def get_centered_fourier(self, use_sinc=False):
        """
        Returns the Fourier form of the box spline. This is a function of the
        variables w_0, w_1, ... w_n. Where n is the dimension of the input
        space.

        This form is automatically centered.

        Arguments:
            use_sinc -- forces the spline to use the symbolic sinc function
                which is not part of sage. This is nice if you want to
                numerically evaluate these functions.
        """
        w = vector([var('w_%d' % i) for i in range(self.s_)])
        if use_sinc:
            phi = prod([
                w.dot_product(vector(xi))
                for xi in self.Xi_
            ])
        else:
            phi = prod([
                        sin(w.dot_product(vector(xi))/2)/(w.dot_product(vector(xi))/2)
                        for xi in self.Xi_
            ])* self.weight
        return phi

    def get_polyhedron(self):
        """
        Returns a Polyhedron object that represents the support for this
        box spline.
        """
        if self.polyhedron_cache:
            return self.polyhedron_cache

        diff = self.get_difference_operator()
        self.polyhedron_cache = Polyhedron(vertices=list(set([tuple(v) for v,_ in diff])), base_ring=AA)
        return self.polyhedron_cache

    def get_difference_operator(self):
        """
        Returns a liat of tuples of the form (lattice_site, weight) that
        specify the difference operator.

        This is mainly to be used internally by the class.
        """
        if self.differ_cache:
            return self.differ_cache

        D  = {tuple([0]*self.s_):1}

        for xi in self.Xi_:
            Dp = {}
            for v in D:
                p = tuple(vector(v) + vector(xi))
                Dp[p] = -D[v]

            for v in D:
                if v in Dp:
                    Dp[v] += D[v]
                else:
                    Dp[v] = D[v]
            D = Dp

        self.differ_cache = [(vector(k)-self.c_xi, D[k]) for k in D if D[k] != 0]
        return self.differ_cache

    def simplify_and_split_from_ker(self, v, isolate, gfunc):
        """
        This is basically equivalent to the r vector function in the paper.

        This is mainly to be used internally by the class.
        """
        # A term in the greens function (fourier) can be represented as
        # (coeffecient, [deg(w_1), ..., deg(w_n)])

        # If we simplify with a vector from the kernel
        # v = [v_1, v_2, ..., v_n]
        # choosing the i'th term as the focal point
        # v_i =/= 0
        #
        # Then the simplification proceeds as
        # for each v_k in v such that v_k != 0
        # (-coeffecient*v_k/v_j, [deg(w_1), deg(w_2), ..., deg(w_k)-1, ..., deg(w_i)+1, ..., deg(w_n)])

        (coeffecient, w) = gfunc
        i = isolate
        olv_v = len([x for x in w if x != 0])
        new_terms = []

        for k, v_k in [(kk, v_kk) for kk, v_kk in enumerate(v) if v_kk != 0 and kk != i]:
            w_prime = w[:]
            w_prime[k] -= 1
            w_prime[i] += 1
            new_terms += [(-coeffecient*v[k]/v[i], w_prime)]

        return new_terms

    def search_nullspace(self, zeros):
        """
        Constructs the \\nu vector as in the paper.

        This is mainly to be used internally by the class.
        """
        K = matrix(SR,self.kerXi_).transpose()
        N = []

        for i,_ in [(i,z) for i,z in enumerate(zeros) if z != 0 ]:
            N.append(K[i])

        C = matrix(SR, N).transpose().kernel().basis()[:]
        lc = vector(C[0])
        return matrix(SR, self.kerXi_).transpose() * lc

    def decompose_term(self, gfd):
        """
        This is the actual recursive form for S, in the paper. It calls the
        two helper methods ``search_nullspace'', and
        ``simplify_and_split_from_ker'', that constuct \\nu and simplify the
        terms in S.

        This is mainly to be used internally by the class.
        """
        (coeff, w) = gfd
        if len([x for x in w if x != 0]) == self.s_:
            return [(gfd)]

        # These are the terms that MUST be zero
        constraint = map(lambda x: 1 if x == 0 else 0, w)
        simp = self.search_nullspace(constraint)

        # find any appropriate term to simplify
        sterm = [i for i,v in enumerate(simp) if v != 0][0]
        decomposition = []
        for term in  self.simplify_and_split_from_ker(simp, sterm, gfd):
            decomposition += self.decompose_term(term)

        # Collect like terms
        compacted = {}
        for v,k in decomposition:
            k = tuple(k)
            compacted[k] = v if k not in compacted else (v + compacted[k])

        return [(compacted[k], list(k)) for k in compacted]

    def decompose_greens(self):
        """
        Decompose the actual greens function into sepeable terms, this Returns
        a list of (w_vars, coeffecient) such that w_vars is a vector with
        s non-zero terms.
        """
        if self.greens_cache:
            return self.greens_cache

        # The starting point for the simplification is the fourier domain greens function
        nu_alpha = [1] * self.n_


        # If we find repetitions of the vector, change the initial representation
        xi_index = {}
        for i, xi in enumerate(self.Xi_):
            xi = tuple(xi)
            if xi not in xi_index:
                xi_index[xi] = i

        for j, xi in enumerate(self.Xi_):
            xi = tuple(xi)
            if xi in xi_index:
                i = xi_index[xi]
                nu_alpha[i] += 1
                nu_alpha[j] -= 1

        gf = (self.weight, nu_alpha)

        if sum([1 if _ != 0 else 0 for _ in nu_alpha]) == self.s_:
            return [gf]

        # pick any vector from the kernel to simplify about
        if len([i for i in nu_alpha if i == 0]) == 0:
            # This is the base case
            simp = self.kerXi_[0]
        else:
            constraint = map(lambda x: 1 if x == 0 else 0, nu_alpha)
            simp = self.search_nullspace(constraint)


        # pick a non zero term in that vector
        term = [i for i,v in enumerate(simp) if v != 0][0]

        # recursively decompose terms
        final = []
        for term in  self.simplify_and_split_from_ker(simp, term, gf):
            final += self.decompose_term(term)

        # Do one more pass collecting like terms
        compacted = {}
        for v,k in final:
            k = tuple(k)
            compacted[k] = v if k not in compacted else (v + compacted[k])

        self.greens_cache = [(compacted[k], list(k)) for k in compacted]
        return self.greens_cache


    def calc_knotplanes(self):
        """
        Calculates the knot planes that touch the support of the box spline.

        Returns two lists, the first is a list of planes that intersect the open
           set of the support of the spline (i.e, only the internal planes), and
           the second list is the planes that touch only the exterior of the
           spline.

           These lists have elements of the form (d,n) where n.dot(x,y,z) -d = 0
        """
        # Calculate the set of knot planes at the origin
        if self.s_ == 1:
            H = list(set(self.Xi_))
        else:
            tmp = set([ncross_product(nt) for nt in combinations(set(self.Xi_), self.s_ - 1)])
            H = [x for x in tmp if len([y for y in x if y != 0]) > 0]
        H = [vector(v).normalized() for v in H]
        #
        Hprime = []
        Hshell = []

        for plane in H:
            d_list = [0]

            for v in self.Xi_:
                dlp = set(d_list[:])
                d = vector(v).dot_product(vector(plane))
                for dp in dlp:
                    d_list.append(d + dp)

            d_list = list(set(d_list))
            d_list.sort()

            min_d, max_d = d_list.pop(0), d_list.pop()
            Hshell += [(min_d, plane), (max_d, plane)]

            for d in d_list:
                Hprime.append((d, plane))

        Hprime = [(d- vector(self.c_xi)*vector(n), n)  for (d, n) in Hprime]
        Hshell = [(d- vector(self.c_xi)*vector(n), n)  for (d, n) in Hshell]

        return Hprime, Hshell

    def _group_planes(self, planes):
        """
        Groups planes by normal. Takes in a list of (d,n), and returns a list
        of (n, [d_1, d_2, ...]).

        Private method.
        """
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

        return  list([(n,plane_set[n]) for n in plane_set])

    def _recursive_split(self, L, P, depth = 0):
        """
        Recursively split P by the planes in L.

        Private method.
        """

        if depth == 0:
            self.pp_num = 0

        if len(L) == 0:
            self.pp_num += 1
            return [P]

        p, D = L[0]
        mid = len(D)//2

        D_A = D[0:mid]
        D_B = D[mid+1:]
        d = D[mid]

        result = []

        B,A = split_polyhedron(P, p, d)

        # Left
        if A is not None:
            Lnew = L[1:]
            if len(D_A) > 0:
                Lnew += [(p, D_A)]
            result += self._recursive_split(Lnew, A, depth+1)
        # Right
        if B is not None:
            Lnew = L[1:]
            if len(D_B) > 0:
                Lnew += [(p, D_B)]
            result += self._recursive_split(Lnew, B, depth+1)

        return result

    def _get_grouped_planes_for_eval(self):
        """
        Groups all the planes that touch the support of the spline.

        Private method.
        """
        if self.gt_cache:
            return self.gt_cache

        Hprime, _ = self.calc_knotplanes()
        A = set([(d, tuple(a)) for (d,a) in Hprime])
        B = set([(d, tuple(a)) for (d,a) in _])
        self.gt_cache = self._group_planes(A.union(B))
        return self.gt_cache


    def get_regions(self):
        """
        Returns a list of all Polyhedron objects that correspond to the distinct
        regions of evaulation for this box spline.
        """

        differential = self.get_difference_operator()
        Hprime, _ = self.calc_knotplanes()
        L = self._group_planes(set([(d, tuple(a)) for (d,a) in Hprime]))
        poly = Polyhedron(vertices = map(lambda x: list(x[0]), differential), base_ring=AA)
        return self._recursive_split(L, poly)


    def poly_term(self, term):
        """
        Returns a ``polyterm'' for a given element of the set S. That is, when
        we have a separable term from the greens frunction, this gives the
        truncated form of that polynomial.

        i.e. it returns a tuple of (polynomial, heaviside expression)
        """

        def _hside(self, x, parent=None, algorithm=None):
            return 0 if x < 0 else (1 if x > 0 else 1/2)
        H =  function_factory('H', 1, '\\text{H}', evalf_func=_hside)
        
        (coeffecient, w) = term
        v = [var('v_%d' % i) for i in range(self.s_)]
        def make_term(i, k):
            return v[i]**(k-1)/factorial(k-1)
        def make_heavy(i):
            return H(v[i])

        itr = [(s,m) for s,m in enumerate(w) if m != 0]
        xi_sigma1 = matrix(SR, [self.Xi_[i] for i,_ in itr]).transpose().inverse()

        transform = coeffecient * prod([make_term(i,m) for i,(s,m) in enumerate(itr)], 1)*abs(xi_sigma1.det())
        heavy = prod([make_heavy(i) for i,(s,m) in enumerate(itr)], 1)

        subs = xi_sigma1 * vector(self.x_)
        for idx, (sigma, mu) in enumerate(itr):
            transform = transform.substitute(v[idx] == subs[idx])
            heavy = heavy.substitute(v[idx] == subs[idx])
        return (transform,heavy)

    def poly_term_w_xi(self, term):
        """
        Returns a ``polyterm'' for a given element of the set S. That is, when
        we have a separable term from the greens frunction, this gives the
        truncated form of that polynomial.

        i.e. it returns a tuple of (polynomial, heaviside expression, Xi)
        where Xi is the transform for the heaviside expression.
        """
        (coeffecient, w) = term
        v = [var('v_%d' % i) for i in range(self.s_)]
        def make_term(i, k):
            return v[i]**(k-1)/factorial(k-1)
        def make_heavy(i):
            return H(v[i])

        itr = [(s,m) for s,m in enumerate(w) if m != 0]
        xi_sigma1 = matrix(SR, [self.Xi_[i] for i,_ in itr]).transpose().inverse()

        transform = coeffecient * prod([make_term(i,m) for i,(s,m) in enumerate(itr)], 1)*abs(xi_sigma1.det())
        heavy = prod([make_heavy(i) for i,(s,m) in enumerate(itr)], 1)

        subs = xi_sigma1 * vector(self.x_)
        for idx, (sigma, mu) in enumerate(itr):
            transform = transform.substitute(v[idx] == subs[idx])
            heavy = heavy.substitute(v[idx] == subs[idx])
        return (transform, heavy, matrix(AA, xi_sigma1))

    def get_polyterms(self):
        """
        Returns the polyterms for each element in the decomposed greens set.
        See ``poly_term_w_xi''. This is basically a list of polynomials,
        heaviside function products, and transforms to a local space for the
        truncated polynomial.
        """
        if self.polyterm_cache:
            return self.polyterm_cache
        greens = self.decompose_greens()
        self.polyterm_cache = []
        for (pp,hs) in [self.poly_term(t) for t in greens]:
            self.polyterm_cache += [(pp.full_simplify(), hs)]

        return self.polyterm_cache

    def get_polyterms_w_xform(self):
        """
        Returns the polyterms for each element in the decomposed greens set.
        See ``poly_term''. This is basically a list of polynomials and heaviside
        function products.
        """
        if self.polytermx_cache:
            return self.polytermx_cache
        greens = self.decompose_greens()
        self.polytermx_cache = []
        for (pp,hs,xi) in [self.poly_term_w_xi(t) for t in greens]:
            self.polytermx_cache += [(pp.full_simplify(), hs, xi)]

        return self.polytermx_cache

    def get_pp_regions(self):
        """
        Returns a list of dictionaries that correspond to the polynomial regions
        of the spline. Each dictionary has the form
            {
               'polynomial': Symbolic polynomial for this region,
               'center': a vector that specifies the center of this region,
               'polyhedron': A polyhedron object that specifies the actual
                   hyper-volume of the region.
            }
        """
        if self.cached_regions:
            return self.cached_regions

        # Get auxilary info
        regions = self.get_regions()
        greens = self.decompose_greens()
        differential = self.get_difference_operator()
        polyterms = self.get_polyterms_w_xform()

        # Compute all the polynomials in each region
        summary = []

        # Stack all the transforms for each region
        stack_xi = polyterms[0][2]
        for (pp, hs, xi) in polyterms[1:]:
            stack_xi = stack_xi.stack(xi)

        # We can easily use that to quickly evaluate xi now

        for idx, polyhedron in enumerate(regions):
            c = polyhedron.center().n()

            # Compute the pp form of the polynomial
            PP = 0
            for jdx, (pos, val) in enumerate(differential):
                posr = vector(pos).n()
                bpp = 0

                pts = stack_xi * (c-vector(pos)).n()
                for idx, p in enumerate(grouper(self.s_, pts)):
                    if any([x <= 0. for x in p]):
                        continue
                    pp, _, _ = polyterms[idx]
                    bpp += pp.subs({self.x_[i]:  self.x_[i] - pos[i] for i,_ in enumerate(pos)})
                PP += val*bpp

            summary.append({
                    'polynomial': PP,
                    'center': c,
                    'polyhedron': polyhedron
                })
        self.cached_regions = summary[:]
        return summary
    
    def evaluate(self,pt):
        """
        A simple evaluation function that takes the input point,
        checks which mesh-region it belongs to, then evaluates
        the polynomial within the region.
        
        pt -- input vector 
        """
        pt = vector(pt)
        preg = self.get_pp_regions()
        poly = self.get_polyhedron()
        pt = vector(pt)

        if pt not in poly:
            return 0

        for region in preg:
            if vector(pt) in region['polyhedron']:
                return region['polynomial'].subs({x:v for (x,v) in zip(bs.x_, pt)})
        return 0
    
    def stable_eval(self, pt):
        """
        This method evaluates a box spline for a given point. This uses a binary
        search over the planes to deduce what region a point is in (on a knot
        plane it arbitrarily makes a decision), then the convolution sum is
        performed over that region.

        Ex:
        bs = BoxSpline([(1,0), (0,1), (1,1), (1,-1)], False)
        bs.stable_eval((0,1))


        This is a method not found in the paper, and it isn't super fast, but
        is still useful at times.
        """

        if not self.get_polyhedron().interior_contains(pt):
            return 0

        pt = vector(pt)
        planes = self._get_grouped_planes_for_eval()

        planes_for_polyhedron = []
        planes_to_split = []

        for n, dlist in planes:
            n = vector(n)
            psign = None
            dot_prod = n.dot_product(pt)

            region_idx = [0, len(dlist)]
            left_sign = sign(dot_prod - dlist[0])
            on_plane = False

            # Do a binary search over the dlist
            while region_idx[1] - region_idx[0] > 1:
                new_idx = sum(region_idx)//2
                d = dlist[new_idx]

                if dot_prod - d == 0:  # Add both sides
                    planes_for_polyhedron += [[-dlist[new_idx - 1]] + list(n)]
                    planes_for_polyhedron += [[dlist[new_idx + 1]] + list(-n)]
                    planes_to_split += [(n, d)]
                    on_plane = True
                    break

                if sign(dot_prod - d) == left_sign:
                    region_idx[0] = new_idx
                else:
                    region_idx[1] = new_idx

            if not on_plane:
                planes_for_polyhedron += [[-dlist[region_idx[0]]] + list(n)]
                planes_for_polyhedron += [[ dlist[region_idx[1]]] + list(-n)]

        poly = Polyhedron(ieqs=planes_for_polyhedron, base_ring=AA)

        for n,d in planes_to_split:
            poly, _ = split_polyhedron(poly, n, d)
            if poly is None:
                poly = _

        differential = self.get_difference_operator()

        c = poly.center()
        PP = 0
        for pos, val in differential:
            bpp = 0
            for (pp, hs) in self.get_polyterms():
                hs = hs.subs({self.x_[i]: c[i] - pos[i] for i,_ in enumerate(pos)})
                if hs.n() > 0.5:
                    bpp += pp.subs({self.x_[i]:  pt[i] - pos[i] for i,_ in enumerate(pos)})
            PP += val*bpp
        return PP
