#!/usr/bin/env python
"""
lattice.sage

A simple helper class to break down an integer lattice into components we care 
about. Specifically, module will compute the Voronoi cell and the Cartesian coset 
structure of the lattice.

Copyright © 2020 Joshua Horacsek

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

class IntegerLattice:
    def __init__(self, generating_matrix):
        """
        Inits the lattice based on the input generating
        matrix. This corresponds to L in part I of the 
        paper.
        """
        if isinstance(generating_matrix, str):
            generating_matrix = lattice_from_str(generating_matrix)
            
        # Check the base field to ensure that we're coming 
        # from Z.
        if generating_matrix.base_ring() != ZZ:
            try: 
                generating_matrix = matrix(ZZ, generating_matrix)
            except:
                raise ValueError("Input generating_matrix is not in ZZ")

        self._d, _ = generating_matrix.dimensions()

        if self._d != _:
            raise ValueError("Input generating_matrix is not square")

        if self._d != generating_matrix.rank():
            raise ValueError("Input generating_matrix is not full rank")

        self._matrix = generating_matrix

        # Setup variables to be memoized
        self._coset       = None
        self._coset_scale = None
        self._voronoi     = None
        
    def hash(self, asstr=False):
        """
        Returns a hash of this space, basically takes an
        ordering of the vectors of the voronoi region, then
        hashes that. This should be unique provided the given 
        lattice isn't a rotated/scaled version of another.
        """
        import hashlib
        u = "%s,%s" % (str(self._d), ','.join([','.join([str(__) for __ in _]) for _ in self.get_basis()]))
        if asstr:
            return hashlib.sha256(u.encode()).hexdigest()
        return hashlib.sha256(u.encode()).digest()

    def dimension(self):
        """
        Returns the dimension of the matrix
        """
        return self._matrix.dimensions()[0]

    def voronoi_region(self):
        """
        Uses the diamond cutting algorithm to obtain
        the voronoi region for this lattice.

        Cached on first call

        returns a Polyhedron object
        """
        if self._voronoi is None:
            from sage.modules.diamond_cutting import calculate_voronoi_cell
            v = calculate_voronoi_cell(self._matrix.transpose())

            # # First, let's try re-casting the polyhedron as one in ZZ
            # try:
            #     v = Polyhedron(vertices=[vector(_)  for _ in v.vertices()], base_ring=ZZ)
            # except:
            #     # If that didn't work, take it as QQ
            #     try:
            #         v = Polyhedron(vertices=[vector(_)  for _ in v.vertices()], base_ring=QQ)
            #     except:
            #         # Just pass and let it be in AA
            #         pass
            self._voronoi = v #Polyhedron(vertices=[vector(_)  for _ in v.vertices()], base_ring=AA)
        return self._voronoi

    def is_lattice_site(self, ls):
        """ 
        @param: ls is the input lattice site, either a vector or list
                of integers

        Determines whether or not the the given lattice site 'ls' 
        lies on the lattice.
        """
        P, cs = self.coset_structure()
        ls = vector(ls)
        for l in cs:
            pt = ls - vector(l)
            if all([x % m == 0 for m,x in zip(P.diagonal(), pt)]):
                return True
        return False 

    def coset_structure(self):
        """
        Returns the Cartesian coset structure of the lattice. Returns
        two values D and cs. 'D' is a diagonal integer matrix that encodes
        axial scales of the cosets. 'cs' is a list of coset shifts, one
        for each coset. For example, the D_n lattices should produce

        D = 2 * matrix.identity(n)
        cs = [ (0,0,0,....,0), (1,1,1,....,1)]

         """
        from itertools import product as p
        frac = lambda x: ceil(x)-x

        if self._coset is None:
            D = matrix.identity(ZZ, self._d)
            maximum = self.voronoi_region().volume()

            # This is sort of naive, but it works
            # Basically step along each axis until 
            # we find a lattice site
            for i in range(self._d):
                for j in range(1, maximum*2):
                    l = [0]*self._d
                    l[i] = j
                    z = self._matrix.solve_right(vector(l))
                    if all([frac(_) == 0 for _ in z]):
                        D[i,i] = j
                        break
            self._coset_scale = D


            Q = Polyhedron(vertices = [self._matrix.solve_right(D*vector(k)) for k in p(*([[0,1]]*self._d))])
            pts = [self._matrix*vector(_) for _ in Q.integral_points()]
            self._coset = []
            for pt in pts:
                touched = False
                for j in range(self._d):
                    if pt[j] == D[j,j]:
                        touched = True
                        break
                if not touched:
                    self._coset += [pt]

            self._coset = sorted(self._coset, reverse=False)

        return self._coset_scale, self._coset

    def count_01(self, scale):
        """
        Counts the number of integer lattice sites
        in the hyper-cube [0,1]^s when the lattice is 
        scaled by 'scale'
        """
        D, c = self.coset_structure()
        try:
            extent = floor(1/scale)
        except:
            return None
        e = vector(ZZ, [extent]*len(c[0]))
        cnt = 0
        for offset in c:
            c_extent = e - vector(offset)
            counts = [floor(ce/d) + 1 for (ce,d) in zip(c_extent, D.diagonal())]
            cnt += prod(counts)

        return cnt

    def get_parallelpiped(self):
        """
        Returns the fundamental parallelpiped 
        for the shortest basis for this lattice. 
        """
        return build_ppiped(self.get_basis())[0]
    def get_basis(self):
        """
        Returns a basis such that all the vectors are minimal,
        then ordered.
        """

        # The basic idea for this method is to enumerate 
        # all the k-facets for the voronoi region
        verts = [_.as_polyhedron().center()*2 for _ in self.voronoi_region().faces(self._d - 1)] 
        

        verts = sorted(verts)
        set_a = []
        set_b = []
        for v in verts:
            if all([_>=0 for _ in v]):
                set_a += [v]
            else:
                set_b += [v]

        verts = set_a + set_b

        basis = [verts[0]]

        for p in verts[1:]:
            basis_prime = basis[:] + [p]

            if matrix(AA, basis).transpose().rank() < matrix(AA, basis_prime).transpose().rank():
                basis = basis_prime
            if matrix(AA, basis).transpose().rank() == self._d:
                break

        return basis
