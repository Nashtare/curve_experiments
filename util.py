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

def curve_security(p, q, no_endo=False, main_factor=0):
    sys.stdout.write('!')
    sys.stdout.flush()
    r = main_factor if main_factor != 0 else ecm.factor(q)[-1]
    if no_endo:
        return (log(PI_12 * 3 * r, 4), embedding_degree(p, r))
    else:
        return (log(PI_12 * r, 4), embedding_degree(p, r))


def embedding_degree(p, r):
    assert gcd(p, r) == 1
    Z_q = Integers(r)
    u = Z_q(p)
    d = r-1
    V = ecm.factor(d)
    for v in V:
        while d % v == 0:
            if u**(d/v) != 1:
                break
            d /= v

    return Integer(d)


def twist_security(p, q, no_endo=False):
    return curve_security(p, 2*(p+1) - q, no_endo)


def twist_security_ignore_embedding_degree(p, q, no_endo=False):
    sys.stdout.write('^')
    sys.stdout.flush()
    r = ecm.factor(2*(p+1) - q)[-1]
    if no_endo:
        return log(PI_12 * 3 * r, 4)
    else:
        return log(PI_12 * r, 4)


####################################################
# FIELD UTILITY FUNCTIONS
####################################################

def make_finite_field(k):
    r""" Return the finite field isomorphic to this field.

    INPUT:

    - ``k`` -- a finite field

    OUTPUT: a tuple `(k_1,\phi,\xi)` where `k_1` is a 'true' finite field,
    `\phi` is an isomorphism from `k` to `k_1` and `\xi` is an isomorphism
    from `k_1` to `k`.

    This function is useful when `k` is constructed as a tower of extensions
    with a finite field as a base field.

    Adapted from https://github.com/MCLF/mclf/issues/103.

    """

    assert k.is_field()
    assert k.is_finite()
    # TODO: partially solved sage9.4 issue but still failing for higher extensions (wrong isomorphic field)
    if k.base_ring().is_prime_field():
        return k, k.hom(k.gen(), k), k.hom(k.gen(), k)
    else:
        k0 = k.base_field()
        G = k.modulus()
        assert G.parent().base_ring() is k0
        k0_new, phi0, _ = make_finite_field(k0)
        G_new = G.map_coefficients(phi0, k0_new)
        k_new = k0_new.extension(G_new.degree())

        alpha = G_new.roots(k_new)[0][0]
        Pk0 = k.cover_ring()
        Pk0_new = k0_new[Pk0.variable_name()]
        psi1 = Pk0.hom(phi0, Pk0_new)
        psi2 = Pk0_new.hom(alpha, k_new)
        psi = psi1.post_compose(psi2)
        # psi: Pk0 --> k_new
        phi = k.hom(Pk0.gen(), Pk0, check=False)
        phi = phi.post_compose(psi)

        k_inv = k0.base_ring()
        phi0_inv = k_inv.hom(k_inv.gen(), k_inv)
        G_new_inv = k_new.modulus().map_coefficients(phi0_inv, k0_new)
        alpha_inv = G_new_inv.roots(k)[0][0]
        phi_inv = k_new.hom(alpha_inv, k)

        return k_new, phi, phi_inv


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
    num_list = repr_field_element(n, nb_hex, output_hex)
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


def poly_weight(poly, p):
    return sum(dist(t, p) for t in poly.coefficients())


def dist(n, p):
    if n > p//2:
        return Integer(p-n)
    else:
        return Integer(n)


def find_sparse_irreducible_poly(ring, degree, use_root=False, max_coeff=10):
    r"""Return an irreducible polynomial of the form X^k - j with smallest j
    in absolute value below max_coeff if any, or 0.

    INPUT:

    - ``ring`` -- a polynomial ring
    - ``degree`` -- the degree of the irreducible polynomial
    - ``use_root`` -- boolean indicating whether using only the ring base field elements as coefficients
                      or using also an element not belonging to the base field (default False)
    - ``max_coeff`` -- maximum absolute value for polynomial coefficients

    OUTPUT: an irreducible polynomial of the form X^k - j with smallest j
    in absolute value below max_coeff if any, or 0.

    """

    x = ring.gen()

    for j in range(1, max_coeff + 1):
        poly = x ** degree - j
        if poly.is_irreducible():
            return poly

    if use_root:
        root = ring.base().gen()
        for j in range(1, max_coeff + 1):
            poly = x ** degree - root*j
            if poly.is_irreducible():
                return poly

    return 0


def find_irreducible_poly(ring, degree, use_root=False, max_coeff=3, output_all=False):
    r"""Return a list of irreducible polynomials with small and few coefficients.

    INPUT:

    - ``ring`` -- a polynomial ring
    - ``degree`` -- the degree of the irreducible polynomial
    - ``use_root`` -- boolean indicating whether using only the ring base field elements as coefficients
                      or using also an element not belonging to the base field (default False)
    - ``max_coeff`` -- maximum absolute value for polynomial coefficients
    - ``output_all`` -- boolean indicating whether outputting only one polynomial or all (default False)

    OUTPUT: a list of irreducible polynomials.

    The default behaviour, to return a single polynomial, still outputs a list of length 1 to keep the
    function output consistent when `output_all == True`.

    """

    x = ring.gen()

    set_coeffs_1 = set(combinations_with_replacement(
        range(-max_coeff, max_coeff), degree))
    set_coeffs_2 = set(combinations_with_replacement(
        reversed(range(-max_coeff, max_coeff)), degree))
    set_coeffs = set_coeffs_1.union(set_coeffs_2)

    list_poly = []
    for coeffs in set_coeffs:
        p = x ** degree
        for n in range(len(coeffs)):
            p += coeffs[n]*x ** n
        if p.is_irreducible():
            list_poly.append(p)

    if use_root:
        root = ring.base().gen()
        for regular_coeffs in set_coeffs:
            p = x ** degree
            for n in range(len(regular_coeffs)):
                p += regular_coeffs[n]*x ** n
            for special_coeffs in set_coeffs:
                q = p
                for n in range(len(special_coeffs)):
                    q += root * special_coeffs[n]*x ** n
                if q.is_irreducible():
                    list_poly.append(q)
                    # Exhaustive search usually becomes too heavy with this,
                    # hence stop as soon as one solution is found
                    if not output_all:
                        return [min(list_poly, key=lambda t: len(t.coefficients()))]

    if output_all or list_poly == []:
        return list_poly
    else:
        return [min(list_poly, key=lambda t: len(t.coefficients()))]
