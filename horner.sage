#!/usr/bin/env python
from sage.symbolic.operators import add_vararg

def pp_to_list(variables, pp):

    # Try to get the leading coeff, if this fails, we have a constant
    try:
        pp.leading_coefficient(variables[0])
    except:
        return [(pp, tuple([0]*len(variables)))]

    # Expand the polynomial
    p = pp.expand()
    plist = []
    
    # Edge case: single term
    operands = p.operands()
    if p.operator() != add_vararg:
        operands = [p]
        
    
    for term in operands:
        # Grab the coeff, if we donk out, that's because term was only a coeff
        try:
            coeff = term.subs({v: 1 for v in variables})
        except:
            coeff = term

        # print coeff

        t = term/coeff
        pwrs = [0]*len(variables)
        for _ in range(len(variables)):
            s = {v: (1 if idx != _ else exp(1)) for idx, v in enumerate(variables)}
            try:
                pwrs[_] = ln(t.subs(s)).simplify()
            except: 
                pass
        plist += [(coeff, tuple(pwrs))]
    return plist

def _greedy_horner(plist, dimension_s):
    factor = greedy_factorization(plist, dimension_s)
    
    # Pulling a factor out would be useless
    if factor is None:
        return {'power':None, 'left': None, 'right': plist}
    
    right = []
    left = []
    
    # Iterate thru the list of terms, and partition it 
    # into left and right sets based on what we can factor out
    common_denom = 1
    
    for coeff, power in plist:
        diff = vector(power) - vector(factor)
        
        if(all([_>= 0 for _ in diff])):
            left += [(coeff, tuple(diff))]
        else:
            right += [(coeff, power)]
    
    return {'power':factor, 'left':_greedy_horner(left, dimension_s), 'right':right}

def horner_factor(polynomial, dimension_s, aux_variables = []):

    plist = pp_to_list(list(var(','.join(['x_%d' % _  for _ in range(dimension_s)]))) + aux_variables, polynomial)
    
    scale = 1
    for coeff, power in plist:
        try:
            if len(coeff.variables()) > 0:
                coeff = coeff.subs({_:1 for _ in coeff.variables()})
        except:
            pass

        try:
            scale = gcd(coeff, scale)
        except:
            scale = scale
    
    plist = [(coeff/scale,power) for (coeff,power) in plist]
    
    # print plist
    
    return {'scale': scale, 'horner': _greedy_horner(plist, dimension_s + len(aux_variables))}
    
def greedy_factorization(plist, dimension_s):
    factor = [0] * dimension_s
    max_index = 0

    # Count the amount of terms pulling out x_index
    # would affect
    for index in range(dimension_s):
        for coeff, power in plist:
            if power[index] > 0:
                factor[index]+=1    
    
    # If we only can only factor a monomial out of one term
    # it's not worth continuing
    if all([_ <= 1 for _ in factor]):
        return None
    
    # Take the maximum
    max_index = factor.index(max(factor))
    factor = [0] * dimension_s
    factor[max_index] = 1
    factor = tuple(factor)
    return factor
    
def power_str(power, variables = None):
    if all([ _ == 0 for _ in power]):
        return '1'
    
    if variables is None:
        variables = ['x_%d' % i for i,_ in enumerate(power)]
        
    r = []
    for term_idx, ex in enumerate(power):
        if ex == 0:
            continue
        r += ['*'.join([variables[term_idx]]*int(ex))]
        
    return '*'.join(r)

def coeff_str(coeff, approx=False):
    if approx:
        return str(RDF(coeff))
    
    try:
        if QQ(coeff).denom() == 1 or QQ(coeff).denom() == -1:
            return "%s" % str(coeff)
        return "(%s)" % str(coeff)
    except:
        pass
    return "%s" % str(coeff) 

def plist_str(plist, variables=None, approx=False):
    if len(plist) == 0:
        return '0'
    
    r = []
    
    for coeff, pwr in plist:
        if coeff == 0:
            continue
        ps = power_str(pwr, variables)
        
        if ps == '1':
            s = coeff_str(coeff, approx)
        elif coeff == 1:
            s = ps
        else:
            s = '%s*%s' % ( coeff_str(coeff,approx), power_str(pwr, variables))
            
        r += [s]
        
    return '+'.join(r)
            
    
def horner_str(h, variables = None, approx=False):
    return "%s*(%s)" % (coeff_str(h['scale'], approx), _horner_str(h['horner'], variables, approx) )

def _horner_str(factorization, variables = None, approx=False):
    power = factorization['power']
    left  = factorization['left']
    right = factorization['right']
    
    if left is None and power is None:
        return plist_str(right,variables, approx)
    
    if len(right) == 0:
        return '%s*(%s)' % (power_str(power,variables), _horner_str(left,variables, approx))
    return '%s*(%s)+%s' % (power_str(power,variables), _horner_str(left,variables,approx), plist_str(right,variables,approx))
    
    