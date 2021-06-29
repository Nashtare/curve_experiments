from sage.all import *

SMALL_PRIMES = (3, 5, 7, 11, 13, 17)
# Starts at 2 as y^2 = x^3 + 1 cannot have prime order
COEFFICIENT_RANGE = range(2, 500)

def gcd_small_primes(p):
    return (r for r in SMALL_PRIMES if gcd(p-1, r) == 1)

# Outputs a BLS scalar field prime order, given its generator
def bls_scalar(x):
    return Integer(cyclotomic_value(12, x))

# Outputs a BLS base field prime order, given its generator
def bls_base(x, r = 0):
    if r == 0:
        tmp = (x-1)**2 * cyclotomic_value(12, x) / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    else:
        tmp = (x-1)**2 * r / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    return Integer(0)

# Outputs 2-adicity of a number
def twoadicity(x, base=0, limit=64):
    return max(i for i in range(base, limit) if ((x-1) % (1<<i) == 0))

# Outputs Hamming weight of a number
def h_weight(x):
    return x.binary().count('1')