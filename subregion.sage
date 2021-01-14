#!/usr/bin/env python


"""

 
"""


class SubRegionTransform:
    def __init__(self):
        self._translate = None
        self._rotate    = None


    def invert_transform(self):
        pass



class SubRegion:
    def __init__(self, region, splinespace):

        self._xv = splinespace._xvars

        # geometric region of the sub-region
        self._region     = region

        # Polynomial in terms of (x_0, x_1, ..., x_{n-1}, c_0, ..., c_k)
        self._polynomial = 0

        self._parent_spline = splinespace

        # Contains the list of (lattice site offset, polynomial)
        self._weighted_sum_terms = {}

        self._factored_terms = {}

        # Index into the re
        self._ref_region = None
        self._transform = None
        self._coset_partition = None

    def partition_into_cosets(self, G, coset_vectors):
        cosets = [[] for _ in range(len(coset_vectors))]

        for coset_index, coset in enumerate(coset_vectors):
            for ls, __ in self.get_ws_list():
                if all([frac(_) == 0 for _ in G.inverse()*(vector(ls)-vector(coset))]):
                    cosets[coset_index] += [ls]
        self._coset_partition = cosets
        return self._coset_partition

    def factor_for_evaluation(self):
        for ls, pp in self.get_ws_list():
            f = horner_factor(pp, self._parent_spline._s)
            self._factored_terms[ls] = f

    def intersection(self, polyhedron, approx=True):
        if approx:
            return approx_polyhedron(self._region).intersection(polyhedron)
        return self._region.intersection(polyhedron)

    def distribute(self, coeff, ls, polynomial):
        self._polynomial += (coeff * polynomial).expand()
        self._weighted_sum_terms[ls] = polynomial.expand()

    def number_of_pts_per_reconstruction(self):
        return len(self.get_ws_list())

    def get_center(self):
        center = None

        for ls, pp in self.get_ws_list():
            if center is None:
                center = vector(ls)
            else:
                center += vector(ls)

        return (1/len(self.get_ws_list())) * center

    def get_ws_list(self):
        return [(k,self._weighted_sum_terms[k]) for k in self._weighted_sum_terms]

    def get_optimal_nn_lookups(self):
        """
        We make the assumption that the optimal 
        nn lookup pattern is the same as the
        optimal linear lookup pattern, not strictly true,
        but can be changed in the future
        """
        v = []
        for coset_idx, covering in enumerate(self._coverings):
            for chunk in covering:
                for pt in chunk:
                    v += [pt]

        return v
    
    def get_optimal_lin_lookups(self, linearize=False):
        if linearize:
            return sum(self._coverings, [])
        return self._coverings

    def order_lookups(self):
        global debug, gg

        debug = None
        gg = None
        T, coset_vectors = self._parent_spline.coset_vectors()

        for coset_idx, covering in enumerate(self._coverings):
            if len(covering) == 1:
                continue

            # Construct the matrix for the TSP
            G = matrix(QQ, len(covering), len(covering))
            for i, cx in enumerate(covering):
                for j, cy in enumerate(covering):
                    if i == j:
                        continue
                    G[i,j] = self.cache_cost_common_cube(cx,cy, T)
            graph = DiGraph(G, format='weighted_adjacency_matrix')

            TSP = graph.traveling_salesman_problem(True)
            edges = {}

            for vin, vout, w in TSP.edges():
                if vin not in edges:
                    edges[vin] = None
                    
                if vout not in edges:
                    edges[vout] = None
                edges[vin] = vout

            num_lookups = len(edges) - 1

            # the heuristic we use is to start at the mem read closest to the origin
            origin = vector([0] *self._parent_spline._s)
            closest_idx = None
            distance = None

            for vertex_idx, vertex_grp in enumerate(covering):
                
                for vertex in vertex_grp:
                    d = origin * vector(vertex)
                    if distance is None or d < distance:
                        closest_idx = vertex_idx
                        distance = d
                        
            vertex = closest_idx
            nc = []
            for _ in range(num_lookups):
                nc += [covering[vertex]]
                vertex = edges[vertex]

            # Bake the covering back into
            self._coverings[coset_idx]

    def cache_cost_intersect(self, g0, g1, G):
        def conv(A, B):
            d = {}
            for a,b in itertools.product(A,B):
                k = tuple(vector(a) + vector(b))
                if k not in d:
                    d[k] = 0
                d[k] += 1
            return d.keys()
        dimension = self._parent_spline._s

        # Construct the stencil of nearest neighbors 
        nn = [
           tuple([1 if i==j else 0 for j in range(dimension)]) for i in range(dimension)
        ] + [
           tuple([-1 if i==j else 0 for j in range(dimension)]) for i in range(dimension)
        ] + [([0]*dimension)]

        nn = [tuple(G*vector(_)) for _ in nn]

        A = conv(nn, g0)
        B = conv(nn, g1)
        overlap = {}

        for pt in A+B:
            if pt not in overlap:
                overlap[pt] = 0
            overlap[pt] += 1

        # print len(B) - sum([1 for k in overlap if overlap[k] > 1])


    def cache_cost_common_cube(self, g0, g1, G):
        min_vertex  = min([vector(_) for _ in g0])
        common_cube = product(*([[0,1]]*self._parent_spline._s))
        common_cube = set([tuple(G*vector(_) + min_vertex) for _ in common_cube])
        return len(common_cube) - sum([1 if pt in common_cube else 0 for pt in g1])

    def cover_lookups(self, linear_fetch):
        G, coset_vectors = self._parent_spline.coset_vectors()


        self._coverings = [[]  for _ in self._coset_partition]
        
        for coset_index, coset in enumerate(self._coset_partition):
            
            # Group the lattice sites based on what the spline allows
            groupings = [[_] for _ in  coset[:]]
            
            site_functions= {ls: self._weighted_sum_terms[ls] for ls in coset}
            if linear_fetch == False:
                pass

            elif self._parent_spline._s == 1:
                self._parent_spline.log("\nWarning: Dimension 1 not supported yet, defaulting to identity covering")
            elif self._parent_spline._s == 2:
                self._parent_spline.log("\nWarning: Dimension 2 not supported yet, defaulting to identity covering")
            elif self._parent_spline._s == 3:
                groupings = enum_8_groups_3d(site_functions, G)
                groupings += enum_4_groups_3d(site_functions, G)
                groupings += enum_2_groups_3d(site_functions, G)
                groupings += enum_1_groups_3d(site_functions, G)

            # Construct the dancing links matrix
            dlm = [[i+1,[j+1 for j, c in enumerate(site_functions) if tuple(c) in g]] for i, g in enumerate(groupings)]
            for i, g in enumerate(groupings):
                dlm += [[i+1,[j+1 for j, c in enumerate(site_functions) if tuple(c) in g]]]
            DLXM = DLXMatrix(dlm)
            
            # Dance, and take the best solution
            min_sol  = len(groupings)+1
            solution = None
            for _ in range(1000000):
                try:
                    x = DLXM.next()
                    if len(x) < min_sol:
                        min_sol = len(x)
                        solution = x
                        self._parent_spline.tick("+")  
                except:
                    break
            
            # There's always a solution, so we don't need to check
            # for a none value
            self._coverings[coset_index] = [groupings[idx-1] for idx in solution]

        # print self._coverings
        return self._coverings