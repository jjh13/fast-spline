#!/usr/bin/env python
"""
This file contains numerous helpers for the box-spline decomposition code.
Mostly, the helper functions in this file are to do with splitting polyhedral
regions, but there are some hacky combinatorial functions, and some function
definitions in here as well.

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

from sage.all_cmdline import *
from itertools import product
from sage.symbolic.function_factory import function_factory
import sympy

import six

if six.PY2:
    from itertools import izip_longest
else:
    from itertools import zip_longest

"""
Define some useful symbolic functions, namely a symbolic heaviside and sinc
function.
"""
def _hside(self, x, parent=None, algorithm=None):
    return 0 if x < 0 else (1 if x > 0 else 1/2)
H =  function_factory('H', 1, '\\text{H}', evalf_func=_hside)
sympy.sinc = sympy.Function("sinc")
def _sinc_fn__(self, x, parent=None, algorithm=None):
    return 2*sin(x/2)/x if x != 0 else 1
sinc = function_factory('sinc', 1, '\\text{sinc}', evalf_func=_sinc_fn__,
        conversions= {'mathematica':'sinc', 'sympy':'sinc'})
"""
 Rho for BCC lattice
"""
def rho_bcc(x,y,z):
    i = round(x/2) * 2;
    j = round(y/2) * 2;
    k = round(z/2) * 2;
    
    sd = (i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z)
    
    _tx = (round((x-1)/2) * 2)+1;
    _ty = (round((y-1)/2) * 2)+1;
    _tz = (round((z-1)/2) * 2)+1;
    pred =  sd < (_tx-x)*(_tx-x) + (_ty-y)*(_ty-y) + (_tz-z)*(_tz-z)
    i = i if pred else _tx
    j = j if pred else _ty
    k = k if pred else _tz
    return (i,j,k)


def rho_fcc(x,y,z):
    i = int(round(x))
    j = int(round(y))
    k = int(round(z))
    
    on_fcc = False
    on_fcc |= (abs(i) % 2 == 0 and abs(j) % 2 == 0 and abs(k) % 2 == 0)
    on_fcc |= (abs(i-1) % 2 == 0 and abs(j-1) % 2 == 0 and abs(k) % 2 == 0)
    on_fcc |= (abs(i-1) % 2 == 0 and abs(j) % 2 == 0 and abs(k-1) % 2 == 0)
    on_fcc |= (abs(i) % 2 == 0 and abs(j-1) % 2 == 0 and abs(k-1) % 2 == 0)
    
    xx = x - i
    yy = y - j
    zz = z - k

    idx = 0 
    idx = 1 if abs(yy) >= abs(xx) and abs(yy) >= abs(zz) else idx
    idx = 2 if abs(zz) >= abs(xx) and abs(zz) >= abs(yy) else idx

    xx = 1 if xx >= 0 else -1
    yy = 1 if yy >= 0 else -1
    zz = 1 if zz >= 0 else -1
    i = i + xx if idx == 0 and not on_fcc else i
    j = j + yy if idx == 1 and not on_fcc else j
    k = k + zz if idx == 2 and not on_fcc else k
        
    return (int(i),int(j),int(k))



def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    if six.PY2:
        return izip_longest(*args, fillvalue=fillvalue)
    else:
        return zip_longest(*args, fillvalue=fillvalue)

def plane_hits_polyhedron(poly, n, d, s = 2):
    """
    Checks whether a plane intersects a polyhedron, it also discards
    intersections at a vertex.
    """
    n = vector(n)
    plane = Polyhedron(eqns=[[-d]+list(n)], base_ring=AA)
    return plane.intersection(poly).dimension() >= s - 1

def n_choose_rho(n, rho):
    """
    Multivatiate combinations:

    return: n!/rho!
    """
    return factorial(n)/prod([factorial(r) for r in rho],1)

def int_solutions_to(n, terms):
    """
    Returns all the integer solutions to the equation
    x_1 + x_2 + ... + x_terms = n
    such that x_i >= 0 for all 1<=i<=terms.

    note that this is easy to count, but less
    trivial to construct

    """
    def recurse(n, max_add, remaining, level = 1):
         if level == n+1:
             return [remaining*[0]]
         results = []
         add_s = min(remaining, max_add)
         for i in range(add_s + 1):
             chunk = [[0]*i]
             results += list(product(chunk, recurse(n, max_add, remaining - i, level + 1)))
         return results
    def unpack(v):
        val = v[0]
        rest = v[1]

        if isinstance(rest, tuple):
            return [val] + unpack(rest)
        return [val, rest]
    def zeros_for_soln(n, m):
        z = recurse(n, m-n, m-n)
        return [unpack(s) for s in z]
    def partition(number):
        answer = set()
        answer.add((number, ))
        for x in range(1, number):
            for y in partition(number - x):
                answer.add((x, ) + y)
        return answer
    psols = [s for s in partition(n) if len(s) <= terms]
    results = []
    for s in psols:
        s = [[i] for i in s] + [[]]
        z = zeros_for_soln(len(s)-1, terms)
        for zp in z:
            results += [sum([y+x for (x,y) in zip(s,zp)], [])]
    return results

def signed_permutation_matrices(n = 3):
    if n <= 1:
        yield matrix([[ 1]])
        yield matrix([[-1]])
    else:
        # Construct matrix out of sub matrices
        for i in range(n):
            for A in signed_permutation_matrices(n-1):
                P = A[0:i]
                Q = A[i:n]
                P = matrix(P.rows() + [vector([0]*(n-1))] + Q.rows())
                Q = matrix([1  if _ == i else 0 for _ in range(n) ]).transpose()
                yield matrix(block_matrix([Q,P], ncols=2).rows())
                Q[i,0] = -1
                yield matrix(block_matrix([Q,P], ncols=2).rows())

def has_negative_element(m):
    for r in m.rows():
        for e in r:
            if e < 0: return True
    return False


def approx_polyhedron(p):
    return Polyhedron(vertices=[vector(RDF, _) for _ in p.vertices()])


def is_degenerate(P):
    """
    Returns wether or not a polyhedron is degenerate.
    """
    return not P.is_compact() or len(P.equations()) > 0 or P.is_empty()

def split_polyhedron(P, plane, d):
    """
    Splits a polyhedron along the plane defined by plane.dot(x,y,z) = d,
    and returns two polyhedra A, B -- one of which may be None.
    """
    plane = vector(plane)
    
    L = Polyhedron(ieqs=[[-d] + list(plane)], base_ring=AA)
    R = Polyhedron(ieqs=[[d] + list(-plane)], base_ring=AA)

    try:
        P_l = P.intersection(L)
        if P_l.dim() < P.dim():
            P_l = None
    except: # This is a lazy way to catch degenerate polytopes
        P_l = None

    try:
        P_r = P.intersection(R)
        if P_r.dim() < P.dim():
            P_r = None
    except:
        P_r = None
    return (P_l, P_r)

def ncross_product(vectors):
    """
    Generialized cross product. Takes in a list of n-1 vectors in n-dimensions,
    and returns an orthogonal vector.
    """
    vectors = list(vectors)
    dim = len(vectors[0])
    if len(vectors) != dim - 1:
        return None
    rv = [0]*dim
    for i in range(dim):
        v = [0]*dim
        v[i] = 1
        rv[i] = matrix(vectors + [v]).det()
    return tuple(rv)

def is_same_plane(p1, p2, flip=True):
    """
    Checks whether two planes are the same (normalizes and compares the
    direction).
    """
    p1 = list(p1)
    p2 = list(p2)

    p1[0] = p1[0]*(-1 if flip else 1)

    n1 = ([i for i in p1 if i != 0][0])
    n2 = [i for i in p2 if i != 0][0]

    return all([ i/n1 == j/n2 for i,j in zip(p1, p2)])

def lattice_sites_in(polyhedron, lattice_f = None):
    """
    Returns the lattice sites within a polyhedron
    for a given lattice
    """

    # If there's no lattice
    if lattice_f is None:
        lattice_f = lambda _: true

    nint = lambda v: sign(v)*ceil(abs(v))
    s = polyhedron.dim()
    v = matrix([vector(v) for v in polyhedron.vertices()]).transpose()

    # Get lattice sites that touch the support
    lattice = [x for x in
        product(*[range(nint(min(list(v)[i])), nint(max(list(v)[i]))+1 )
        for i in range(s)])
        if vector(x) in polyhedron and lattice_f(x)
    ]
    return lattice

def lattice_from_str(lattice):
    if lattice.upper() == "CP": # Cartesian planar
        return matrix.identity(ZZ, 2)

    elif lattice.upper() == "QC":
        return matrix(ZZ, [[1,1],[-1,1]])

    elif lattice.upper() == "CC":
        return matrix.identity(ZZ, 3)
    elif lattice.upper() == "BCC":
        return matrix(ZZ, [[-1,1,1],[1,-1,1],[1,1,-1]])
    elif lattice.upper() == "FCC":
        return matrix(ZZ, [[1,1,0],[1,0,1],[0,1,1]])

    if lattice.upper()[0] == "C":
        dim = int(lattice[1:])
        return matrix.identity(ZZ, dim)

    if lattice.upper()[0] == "D" and lattice.upper()[-1] == 'S':
        dim = int(lattice[1:-1])
        return None

    if lattice.upper()[0] == "D":
        dim = int(lattice[1:])

        return matrix(ZZ, [
            [2 if i == j and i < dim-1 else (1 if i == dim -1 else 0) for i in range(dim)] for j in range(dim)
            ])
    return None


def build_ppiped(P, c=False):
    """
    Builds the parallelpiped defined 
    from the direction vectors in P.

    Returns a set of 2^s points and 
    the center of the parallel piped
    """
    d = len(P[0])
    P = [vector(_) for _ in P]
    S = [vector([0]*d)]
    Sp = []
    for p in P:
        for s in S:
            Sp += [s + p, s]
        S = Sp
        Sp = []
    c_xi = vector([0]*d)
    if c:
        c_xi = sum(S)/len(S)
        S = [_-c_xi for _ in S]
    return S, c_xi


def _label_suffix(label, suffix):
    """Returns (label + suffix) or a truncated version if it's too long.
    Parameters
    ----------
    label : str
        Label name
    suffix : str
        Label suffix
    """
    if len(label) > 50:
        nhead = 25
        return ''.join([label[:nhead], '..', suffix])
    else:
        return label + suffix

def shift_polyhedron(p, x):
    """
    Returns a the polyhedron defined by p, but shifted
    by x
    """
    x = vector(x)
    return Polyhedron(vertices=[vector(_)+x for _ in p.vertices()])

def _bw(integer):
    if integer == 0 or integer == 0.0:
        return 8
    
    r = int(ceil(log(float(integer))))
    if r <= 8:
        return 8
    if r <= 16:
        return 16
    if r <= 32:
        return 32
    return 64


def enum_8_groups_3d(sf_s, G):
    cubes = []
    cond1 = [
        [[(1, 1, 0), (0, 0, 1)], 1],
        [[(1, 1, 0), (1, 0, 1)], 1],
        [[(1, 1, 0), (0, 1, 1)], 1],
        [[(0, 0, 0), (1, 1, 1)], -1],
        [[(1, 0, 0), (1, 1, 1)], -1],
        [[(0, 1, 0), (1, 1, 1)], -1],
    ]
    cond2 = [
        [[(0, 1, 0), (1, 0, 1)], 1],
        [[(1, 1, 0), (1, 0, 1)], 1],
        [[(1, 0, 1), (0, 1, 1)], 1],
        [[(0, 0, 0), (1, 1, 1)], -1],
        [[(1, 0, 0), (1, 1, 1)], -1],
        [[(0, 0, 1), (1, 1, 1)], -1],
    ]
    cond3 = [
        [[(1, 0, 0), (0, 1, 1)], 1],
        [[(1, 1, 0), (0, 1, 1)], 1],
        [[(1, 0, 1), (0, 1, 1)], 1],
        [[(0, 0, 0), (1, 1, 1)], -1],
        [[(0, 1, 0), (1, 1, 1)], -1],
        [[(0, 0, 1), (1, 1, 1)], -1],
    ]
    cond4 = [
        [[(1, 0, 0), (0, 1, 0)], -1],
        [[(0, 0, 0), (1, 1, 0)], 1],
        [[(1, 1, 0), (0, 0, 1)], 1],
        [[(0, 1, 0), (1, 0, 1)], -1],
        [[(1, 0, 0), (0, 1, 1)], -1],
        [[(1, 0, 1), (0, 1, 1)], -1],
        [[(0, 0, 0), (1, 1, 1)], 1],
        [[(0, 0, 1), (1, 1, 1)], 1],
    ]
    cond5 = [
        [[(1, 0, 0), (0, 0, 1)], -1],
        [[(1, 1, 0), (0, 0, 1)], -1],
        [[(0, 0, 0), (1, 0, 1)], 1],
        [[(0, 1, 0), (1, 0, 1)], 1],
        [[(1, 0, 0), (0, 1, 1)], -1],
        [[(1, 1, 0), (0, 1, 1)], -1],
        [[(0, 0, 0), (1, 1, 1)], 1],
        [[(0, 1, 0), (1, 1, 1)], 1],
    ]
    cond6 = [
        [[(0, 1, 0), (0, 0, 1)], -1],
        [[(1, 1, 0), (0, 0, 1)], -1],
        [[(0, 1, 0), (1, 0, 1)], -1],
        [[(1, 1, 0), (1, 0, 1)], -1],
        [[(0, 0, 0), (0, 1, 1)], 1],
        [[(1, 0, 0), (0, 1, 1)], 1],
        [[(0, 0, 0), (1, 1, 1)], 1],
        [[(1, 0, 0), (1, 1, 1)], 1],
    ]
    cond7 = [
        [[(1, 0, 0), (0, 1, 0), (0, 0, 1)], 1],
        [[(1, 0, 0), (1, 1, 0), (0, 0, 1)], 1],
        [[(0, 1, 0), (1, 1, 0), (0, 0, 1)], 1],
        [[(1, 1, 0), (1, 1, 0), (0, 0, 1)], 1],
        [[(1, 0, 0), (0, 1, 0), (1, 0, 1)], 1],
        [[(1, 0, 0), (1, 1, 0), (1, 0, 1)], 1],
        [[(0, 1, 0), (1, 1, 0), (1, 0, 1)], 1],
        [[(1, 1, 0), (1, 1, 0), (1, 0, 1)], 1],
        [[(0, 1, 0), (0, 0, 1), (1, 0, 1)], 1],
        [[(1, 1, 0), (0, 0, 1), (1, 0, 1)], 1],
        [[(0, 1, 0), (1, 0, 1), (1, 0, 1)], 1],
        [[(1, 1, 0), (1, 0, 1), (1, 0, 1)], 1],
        [[(1, 0, 0), (0, 1, 0), (0, 1, 1)], 1],
        [[(1, 0, 0), (1, 1, 0), (0, 1, 1)], 1],
        [[(0, 1, 0), (1, 1, 0), (0, 1, 1)], 1],
        [[(1, 1, 0), (1, 1, 0), (0, 1, 1)], 1],
        [[(1, 0, 0), (0, 0, 1), (0, 1, 1)], 1],
        [[(1, 1, 0), (0, 0, 1), (0, 1, 1)], 1],
        [[(1, 0, 0), (1, 0, 1), (0, 1, 1)], 1],
        [[(0, 1, 0), (1, 0, 1), (0, 1, 1)], 1],
        [[(1, 1, 0), (1, 0, 1), (0, 1, 1)], 2],
        [[(0, 0, 1), (1, 0, 1), (0, 1, 1)], 1],
        [[(1, 0, 1), (1, 0, 1), (0, 1, 1)], 1],
        [[(1, 0, 0), (0, 1, 1), (0, 1, 1)], 1],
        [[(1, 1, 0), (0, 1, 1), (0, 1, 1)], 1],
        [[(1, 0, 1), (0, 1, 1), (0, 1, 1)], 1],
        [[(0, 0, 0), (0, 0, 0), (1, 1, 1)], -1],
        [[(0, 0, 0), (1, 0, 0), (1, 1, 1)], -2],
        [[(1, 0, 0), (1, 0, 0), (1, 1, 1)], -1],
        [[(0, 0, 0), (0, 1, 0), (1, 1, 1)], -2],
        [[(1, 0, 0), (0, 1, 0), (1, 1, 1)], -1],
        [[(0, 1, 0), (0, 1, 0), (1, 1, 1)], -1],
        [[(0, 0, 0), (1, 1, 0), (1, 1, 1)], -2],
        [[(1, 0, 0), (1, 1, 0), (1, 1, 1)], -1],
        [[(0, 1, 0), (1, 1, 0), (1, 1, 1)], -1],
        [[(0, 0, 0), (0, 0, 1), (1, 1, 1)], -2],
        [[(1, 0, 0), (0, 0, 1), (1, 1, 1)], -1],
        [[(0, 1, 0), (0, 0, 1), (1, 1, 1)], -1],
        [[(0, 0, 1), (0, 0, 1), (1, 1, 1)], -1],
        [[(0, 0, 0), (1, 0, 1), (1, 1, 1)], -2],
        [[(1, 0, 0), (1, 0, 1), (1, 1, 1)], -1],
        [[(1, 1, 0), (1, 0, 1), (1, 1, 1)], 1],
        [[(0, 0, 1), (1, 0, 1), (1, 1, 1)], -1],
        [[(0, 0, 0), (0, 1, 1), (1, 1, 1)], -2],
        [[(0, 1, 0), (0, 1, 1), (1, 1, 1)], -1],
        [[(1, 1, 0), (0, 1, 1), (1, 1, 1)], 1],
        [[(0, 0, 1), (0, 1, 1), (1, 1, 1)], -1],
        [[(1, 0, 1), (0, 1, 1), (1, 1, 1)], 1],
        [[(0, 0, 0), (1, 1, 1), (1, 1, 1)], -2],
        [[(1, 0, 0), (1, 1, 1), (1, 1, 1)], -1],
        [[(0, 1, 0), (1, 1, 1), (1, 1, 1)], -1],
        [[(0, 0, 1), (1, 1, 1), (1, 1, 1)], -1]
    ]

    for site in sf_s:
        # First we check to see if each
        cube = [tuple(G * vector(_)) for _ in product(*([[0, 1]] * 3))]
        if any([tuple(vector(site) + vector(v)) not in sf_s for v in cube]):
            continue
        satisfies = True

        for c in [cond1, cond2, cond3, cond4, cond5, cond6, cond7]:
            t = 0
            for clist, w in c:
                clist = [tuple(vector(site) + G * vector(_)) for _ in clist]
                t += w * prod([sf_s[_] for _ in clist])
            if t.expand() != 0:
                satisfies = False
                break

        if satisfies:
            cubes += [[tuple(vector(site) + vector(_)) for _ in cube]]

    return cubes


def enum_4_groups_3d(sf_s, G):
    squares = []
    for site in sf_s:

        for e in range(3):
            # Create a square
            square = [[0, 0], [1, 1], [1, 0], [0, 1]]
            for s in square:
                s.insert(e, 0)
            square = [tuple(G * vector(_) + vector(site)) for _ in square]

            if any([v not in sf_s for v in square]):
                continue

            if (sf_s[square[0]] * sf_s[square[1]] - sf_s[square[2]] * sf_s[square[3]]).expand() == 0:
                squares += [[tuple(vector(_)) for _ in square]]

            # todo: check conditions
    return squares


def enum_2_groups_3d(sf_s, G):
    pairs = []

    for site in sf_s:
        site = vector(site)
        for e in range(3):
            s2 = vector(site[:]) + G * vector([1 if e == i else 0 for i in range(3)])
            if tuple(s2) in sf_s:
                pairs += [[tuple(site), tuple(s2)]]
    return pairs


def enum_1_groups_3d(sf_s, scale):
    return [[_] for _ in sf_s.keys()]


def poly_coeffs(p, cvar_to_site):
    p = p.expand()
    term_map = {}
    for term in p.operands():
        # Grab the coeff, if we donk out, that's because term was only a coeff
        try:
            coeff = term.subs({v: 1 for v in term.variables()})
        except:
            coeff = term

        # Grab c
        c = 0
        for v in term.variables():
            if v in cvar_to_site:
                c = v

        term = term / (c * coeff)

        # Figure out the degree of the term
        # Constant term
        if term == 1:
            index = (0,) * dim_s
        else:
            index = tuple([term.degree(var('x_%d' % v)) for v in xrange(dim_s)])

        if index not in term_map:
            term_map[index] = [
                0,
                []
            ]
        term_map[index][0] += coeff * c
        term_map[index][1] += [coeff * c]

    new_poly = 0
    for d in term_map:
        x = prod([var('x_%d' % i) ^ deg for i, deg in enumerate(d)])
        new_poly += term_map[d][0] * x

    return new_poly, term_map


def consistent_vars(p, cvar_to_site):
    pvar_map = {}
    for k in p:
        for v in p[k][0].variables():
            if v not in pvar_map:
                pvar_map[v] = True

    pvars = []
    for c in cvar_to_site:
        if c in pvar_map:
            pvars += [c]

    return pvars


def construct_matrix(p, deg_order, dim_s, c_var_to_site, R=None):
    """
       p is a dict with key polynomial degree (in (0,2,3) form), maps to (coeff, c_list)
       deg_order is a list of degrees (so we have a consistent ordering between matrices)

       returns a matrix whose rows correspond to a degree, and cols correspond to the the
       weight of a coeff in the linear sum of coeffs
    """
    if R is None:
        R = matrix.identity(dim_s)

    pvars = consistent_vars(p, c_var_to_site)
    pp = {}
    try:
        for i, deg in enumerate(p):
            sgn = prod([sign(_) ** abs(_) for _ in vector(deg) * R])
            dnew = tuple([abs(_) for _ in vector(deg) * R])
            pp[dnew] = sgn * p[deg][0]
    except:
        return None

    A = matrix(QQ, len(p), len(pvars))
    for i, k in enumerate(deg_order):
        for j, v in enumerate(pvars):
            A[i, j] = pp[k].coefficient(v)
    return A


def match_coeffs(B, A, q, p, c_var_to_site):
    pvars = []
    qvars = []

    # Check for row consistency
    for (a, b) in zip(A.columns(), B.columns()):
        if any([c != d for (c, d) in zip(sorted(a), sorted(b))]):
            Acdf = matrix(CDF, A)
            Bcdf = matrix(CDF, B)

            _, As, _ = Acdf.SVD()
            _, Bs, _ = Bcdf.SVD()
            closeness = (As - Bs).norm()
            print( "inconsistent cols!", closeness)
            return None

    Acdf = matrix(CDF, A)
    Bcdf = matrix(CDF, B)

    _, As, _ = Acdf.SVD()
    _, Bs, _ = Bcdf.SVD()
    closeness = (As - Bs).norm()

    if closeness > 1e-10:
        print("SVD OPT")
        return None
    else:
        print("close?", closeness)

    n = A.nrows()
    P = matrix(QQ, n, n)
    K = matrix(QQ, n, n)

    P = recursive_find_fasty(A, B, P, K, n, -1, -1)

    if P is None or A != P * B:
        return None

    pvars = consistent_vars(p, c_var_to_site)
    qvars = consistent_vars(q, c_var_to_site)

    return {t: d for (t, d) in zip(vector(pvars) * P, vector(qvars))}


def find_P_and_t(packed_info_p, packed_info_q, dim_s, c_var_to_site):
    # unpack
    p_shape, pp, pc, c_p = packed_info_p
    q_shape, pq, qc, c_q = packed_info_q

    __, pvarc = poly_coeffs(pc, c_var_to_site)
    __, pvar = poly_coeffs(pp, c_var_to_site)

    __, qvarc = poly_coeffs(qc, c_var_to_site)
    __, qvar = poly_coeffs(pq, c_var_to_site)
    __, qvar = poly_coeffs(pq, c_var_to_site)

    # Choose a consistent ordering for the coeffs
    deg_orderc = pvarc.keys()
    deg_order = pvar.keys()

    # Precompute the B matrix
    Bc = construct_matrix(pvarc, deg_orderc, dim_s, c_var_to_site).transpose()
    B = construct_matrix(pvar, deg_order, dim_s, c_var_to_site).transpose()

    t = vector(p_shape.center()) - vector(q_shape.center())

    best_solution = None

    # find a rotation/reflection that works
    for P in signed_permutation_matrices(dim_s):
        # Do the change of variables
        X = vector([var('x_%d' % i) for i in range(dim_s)])
        subs_x = {a: b for (a, b) in zip(X, P * X)}

        Ac = construct_matrix(qvarc, deg_orderc, P, dim_s, c_var_to_site).transpose()
        A = construct_matrix(qvar, deg_order, P, dim_s, c_var_to_site).transpose()

        if Ac is None or A is None:
            print("Something inconsistent")

        if Ac is None and A is None:
            print("Inconsistent Renaming...")
            continue

        c_subs = None
        c_subsc = None

        if A is not None:
            c_subs = match_coeffs(B, A, pvar, qvar, c_var_to_site)

        if Ac is not None:
            c_subsc = match_coeffs(Bc, Ac, pvarc, qvarc, c_var_to_site)

        if c_subs is not None:
            print("Best solution")
            if pq.subs(subs_x).subs(c_subs).expand() == pp.expand():
                print("and it's valid...")
                best_solution = (P, (0,) * dim_s, c_subs)

        elif c_subsc is not None:
            print("Okay solution")
            if qc.subs(subs_x).subs(c_subsc).expand() == pc.expand():
                best_solution = (P, t, c_subsc)

        if not has_negative_element(P) and best_solution is not None:
            break

    if best_solution is not None:
        P, t, c_subs = best_solution
        subs_x = {a: b for (a, b) in zip(X, P * X)}
        # print
        # pp.expand()
        # print
        # pq.subs(subs_x).subs(c_subs).expand()

    return best_solution

import time
import copy

def recursive_find_fasty(A, B, P, K, n, i, j, Bmap=None, depth=0):
    P = copy(P)
    K = copy(K)

    t0 = time.time()
    if i >= 0 and j >= 0:
        if K[i, j] == 0 and P[i, j] != 1:
            P[i, j] = 1
            for ii in range(n):
                K[ii, j] = -1
                K[i, ii] = -1

        else:
            return None

    changed = True
    while changed:
        changed = False
        for a, b in zip(A.columns(), B.columns()):
            for i, p in enumerate(a):

                count = 0
                jdx = -1

                for j, q in enumerate(b):
                    if p == q and K[i, j] == 0:
                        count += 1
                        if count == 1:
                            jdx = j

                if count == 1 and jdx >= 0:
                    if K[i, jdx] == -1 and P[i, jdx] != 1:
                        print("(%d,%d) violated constraint" % (i, jdx))
                        return None

                    if P[i, jdx] == 1:
                        continue

                    changed = True
                    P[i, jdx] = 1
                    for ii in range(n):
                        K[ii, jdx] = -1
                        K[i, ii] = -1

    if P.rank() < n:
        for a, b in zip(A.columns(), B.columns()):
            for i, p in enumerate(a):
                for j, q in enumerate(b):
                    if p == q and K[i, j] == 0:
                        # make sure that this choice is consistent...
                        if A[i, :] != B[j, :]:
                            continue
                        Pq = recursive_find_fasty(A, B, P, K, n, i, j, Bmap, depth + 1)
                        if Pq is not None:
                            return Pq
        return None

    if A != P * B:
        return None

    return P
