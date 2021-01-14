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
        for i in xrange(s)])
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
