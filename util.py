"""
This module aims at providing utility functions for other modules.
"""

from sage.all import *

SMALL_PRIMES = (3, 5, 7, 11, 13, 17)
# Starts at 2 as y^2 = x^3 + 1 cannot have prime order
COEFFICIENT_RANGE = range(2, 20)
RHO_SECURITY = 125
EXTENSION_SECURITY = 120
TWIST_SECURITY = 100
MIN_EMBEDDING_DEGREE = 200
# For Pollard-Rho security analysis
PI_12 = (pi/12).numerical_approx()


def gcd_small_primes(p):
    return (r for r in SMALL_PRIMES if gcd(p-1, r) == 1)


####################################################
# BLS UTILITY FUNCTIONS
####################################################

def bls12_scalar(x):
    return Integer(cyclotomic_value(12, x))


def bls12_base(x, r=0):
    if r == 0:
        tmp = (x-1)**2 * cyclotomic_value(12, x) / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    else:
        tmp = (x-1)**2 * r / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    return Integer(0)


def bls24_scalar(x):
    return Integer(cyclotomic_value(24, x))


def bls24_base(x, r=0):
    if r == 0:
        tmp = (x-1)**2 * cyclotomic_value(24, x) / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    else:
        tmp = (x-1)**2 * r / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    return Integer(0)


def bls48_scalar(x):
    return Integer(cyclotomic_value(48, x))


def bls48_base(x, r=0):
    if r == 0:
        tmp = (x-1)**2 * cyclotomic_value(48, x) / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    else:
        tmp = (x-1)**2 * r / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    return Integer(0)


def twoadicity(x, base=0, limit=256):
    return max(i for i in range(base, limit) if ((x-1) % (1 << i) == 0))


def threeadicity(x, base=0, limit=128):
    return max(i for i in range(base, limit) if ((x-1) % (3**i) == 0))


def h_weight(x):
    return x.binary().count('1')


####################################################
# CURVE SECURITY FUNCTIONS
####################################################

def curve_security(p, q, main_factor=0):
    sys.stdout.write('!')
    sys.stdout.flush()
    r = main_factor if main_factor != 0 else factor(q)[-1][0]
    return (log(PI_12 * r, 4), embedding_degree(p, r))


def embedding_degree(p, r):
    assert gcd(p, r) == 1
    Z_q = Integers(r)
    u = Z_q(p)
    d = r-1
    V = factor(d)
    for (v, k) in V:
        while d % v == 0:
            if u**(d/v) != 1:
                break
            d /= v

    return Integer(d)


def twist_security(p, q):
    return curve_security(p, 2*(p+1) - q)


def twist_security_ignore_embedding_degree(p, q):
    r = factor(2*(p+1) - q)[-1][0]
    return log(PI_12 * r, 4)


####################################################
# FIELD UTILITY FUNCTIONS
####################################################

def repr_field_element(n, nb_hex=64, output_hex=True):
    assert nb_hex % 16 == 0
    n = str(hex(Integer(n)))[2:]
    while len(n) < nb_hex:
        n = "0" + n
    num_list = reversed(range(nb_hex//16))
    res = []
    for i in num_list:
        if output_hex:
            res.append("0x" + n[i*16:i*16+16])
        else:
            res.append(Integer("0x" + n[i*16:i*16+16]))
    return res


def pretty_element_repr(n, nb_hex=64, output_hex=True):
    num_list = repr_field_element(n)
    res = "\n[\n"
    for i in num_list:
        res += "    %s,\n" % i
    res += "]"
    return res


def repr_field_element_bytes(n, nb_bytes=32, output_hex=False):
    assert nb_bytes % 16 == 0
    n = str(hex(Integer(n)))[2:]
    while len(n) < nb_bytes * 2:
        n = "0" + n
    res = []
    for i in range(nb_bytes-1, -1, -1):
        if output_hex:
            res.append("0x" + n[i*2:i*2+2])
        else:
            res.append(Integer("0x" + n[i*2:i*2+2]))
    return res
