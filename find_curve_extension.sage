import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations_with_replacement

from util import *

if sys.version_info[0] == 2: range = xrange

# Adapted from https://github.com/MCLF/mclf/issues/103
def make_finite_field(k):
    r""" Return the finite field isomorphic to this field.

    INPUT:

    - ``k`` -- a finite field

    OUTPUT: a tuple `(k_1,\phi,\xi)` where `k_1` is a 'true' finite field,
    `\phi` is an isomorphism from `k` to `k_1` and `\xi` is an isomorphism
    from `k_1` to `k`.

    This function is useful when `k` is constructed as a tower of extensions
    with a finite field as a base field.

    """

    assert k.is_field()
    assert k.is_finite()
    k0 = k.base_field()
    G = k.modulus()
    assert G.parent().base_ring() is k0
    # Won't work with higher tower extensions but previous version
    # failing with Sage 9.4, investigate how to solve it
    k0_new = k.base_field()
    phi0 = k0.hom(k0.gen(), k0)
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


def find_curve(extension, max_cofactor, wid = 0, processes = 1):
    a = extension.primitive_element()
    for i in range(wid + 1, 1000000000, processes):
        sys.stdout.write(".")
        sys.stdout.flush()

        coeff_a = 1
        coeff_b = a**i

        E = EllipticCurve(extension, [coeff_a, coeff_b])

        n = E.count_points()
        prime_order = n.factor()[-1][0]
        cofactor = n // prime_order
        if cofactor > max_cofactor:
            continue

        sys.stdout.write("o")
        sys.stdout.flush()
        g = E.random_element()
        while g.order() != prime_order:
            g = E.random_element()
            
        (rho_sec, k) = curve_security(extension.base().cardinality(), n, prime_order)

        if k.nbits() < MIN_EMBEDDING_DEGREE:
            continue

        sys.stdout.write("+")
        sys.stdout.flush()

        if rho_sec < RHO_SECURITY:
            continue

        yield (extension, E, g, prime_order, cofactor, i, coeff_a, coeff_b, rho_sec, k)


# Outputs parameters of valid curves over an extension of a small field
def print_curve(prime, extension_degree, max_cofactor, wid = 0, processes = 1):
    Fp = GF(prime)
    K.<x> = Fp[]
    factors = list(extension_degree.factor())
    for n in range(len(factors)):
        degree = factors[n][0]
        for i in range(factors[n][1]): # multiplicity
            poly = find_low_poly(K, degree)
            if poly == []:
                poly = find_low_poly(K, degree, use_root=True)
                if poly == []:
                    raise ValueError('Could not find an irreducible polynomial with specified parameters.')
            Fp = Fp.extension(poly, "u_{0}".format(n))
            K.<x> = Fp[]

    extension, phi, phi_inv = make_finite_field(Fp)

    for (extension, E, g, order, cofactor, index, coeff_a, coeff_b, rho_security, embedding_degree) in find_curve(extension, max_cofactor, wid, processes):
        output = "\n\n\n"
        output += "E(GF((%s)^%s)) : y^2 = x^3 + x + %s (b == a^%s)\n" % (extension.base_ring().order().factor(), extension.degree(), coeff_b, index)
        output += "E(GF((%s)^%s)) : y^2 = x^3 + x + %s\n" % (Fp.base_ring().order().factor(), Fp.degree(), phi_inv(coeff_b))
        output += "E generator point: %s\n" % g
        output += "E prime order: %s (%s bits)\n" % (order, order.nbits())
        output += "E cofactor: %s\n" % cofactor
        output += "E security (Pollard-Rho): %s\n" % rho_security
        output += "E embedding degree: %s\n" % embedding_degree
        print(output)
    return

def find_low_poly(ring, degree, use_root = False, max_coeff = 5, output_all = False):
    x = ring.gen()
    
    set_coeffs_1 = set(combinations_with_replacement(range(-max_coeff, max_coeff), degree))
    set_coeffs_2 = set(combinations_with_replacement(reversed(range(-max_coeff, max_coeff)), degree))
    set_coeffs = set_coeffs_1.union(set_coeffs_2)

    list_poly = []
    for coeffs in set_coeffs:
        p = x^degree
        for n in range(len(coeffs)):
            p += coeffs[n]*x^n
        if p.is_irreducible():
            list_poly.append(p)
    
    if use_root:
        root = ring.base().gen()
        for regular_coeffs in set_coeffs:
            p = x^degree
            for n in range(len(regular_coeffs)):
                p += regular_coeffs[n]*x^n
            for special_coeffs in set_coeffs:
                q = p
                for n in range(len(special_coeffs)):
                    q += root * special_coeffs[n]*x^n
                if q.is_irreducible():
                    list_poly.append(q)
                    # Exhaustive search usually becomes too heavy with this,
                    # hence stop as soon as one solution is found
                    if not output_all:
                        return min(list_poly, key = lambda t: len(t.coefficients()))

    if output_all or list_poly == []:
        return list_poly
    else:
        return min(list_poly, key = lambda t: len(t.coefficients()))

########################################################################

def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    strategy = print_curve
    help = True if "--help" in args else False
    args = [arg for arg in args if not arg.startswith("--")]

    if help:
        print("""
Cmd: sage find_curve_extension.sage [--sequential] <prime> <extension_degree> <max_cofactor>

Args:
    --sequential        Uses only one process
    <prime>             A prime number, default 2^62 - 111 * 2^39 + 1
    <extension_degree>  The extension degree of the prime field, default 6
    <max_cofactor>      Maximum cofactor of the curve, default 64
""")
        return

    prime = int(args[0]) if len(args) > 0 else 2^62 - 111 * 2^39 + 1
    extension_degree = int(args[1]) if len(args) > 1 else 6
    max_cofactor = int(args[2]) if len(args) > 2 else 64

    if processes == 1:
        strategy(prime, extension_degree, max_cofactor)
    else:
        print("Using %d processes." % processes)
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(worker, (strategy, prime, extension_degree, max_cofactor, wid, processes))

            while True:
                sleep(1000)
        except (KeyboardInterrupt, SystemExit):
            pass
        finally:
            pool.terminate()

def worker(strategy, *args):
    res = []
    try:
        res = real_worker(strategy, *args)
    except (KeyboardInterrupt, SystemExit):
        pass
    except:
        print_exc()
    finally:
        return res

def real_worker(strategy, *args):
    return strategy(*args)

main()
