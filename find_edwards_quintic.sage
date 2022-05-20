"""
This module aims at finding curves defined over finite field extensions.
"""

import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations_with_replacement

from util import *
from util_hashtocurve import OptimizedSSWU

if sys.version_info[0] == 2:
    range = xrange


def find_curve(extension, wid=0, processes=1):
    r"""Yield curve constructed over a prime field extension.

    INPUT:

    - ``extension`` -- the field direct extension
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    OUTPUT:

    - ``extension`` -- the field extension
    - ``E`` -- the curve definition
    - ``g`` -- a generator of the large prime order subgroup
    - ``prime_order`` -- the prime order of the large subgroup generated by g
    - ``cofactor`` -- the cofactor of the curve
    - ``i`` -- the exponent of the extension primitive element defining the b coefficient
    - ``coeff_a`` -- the a coefficient of the curve in short Weierstrass form (always 1)
    - ``rho_sec`` -- the Pollard-Rho security of the curve
    - ``k`` -- the embedding degree of the curve
    - ``twist_rho_sec`` -- the Pollard-Rho security of the twist

    """

    ext_gen = extension.gen()
    a = -1
    for i in range(wid + 1, 1000000000, processes):
        coeff_a = ext_gen*4*i + 2
        if (coeff_a**2 - 4).is_square():
            continue

        d = - (coeff_a - 2) / (coeff_a + 2)
        if d.is_square():
            continue
        # Required for Ristretto encoding, (ad - 1)^2 == (a - d)^2 for a == -1
        if not (a*d - 1).is_square():
            continue

        sys.stdout.write(".")
        sys.stdout.flush()

        E = EllipticCurve(extension, [0, coeff_a, 0, 1, 0])

        n = E.count_points()
        prime_order = list(ecm.factor(n))[-1]
        cofactor = n // prime_order
        if cofactor != 4 or cofactor != 8:
            continue

        sys.stdout.write("o")
        sys.stdout.flush()

        # We generate a point on the curve with the SSWU hash-to-curve algorithm.
        # If the point is not in the prime-order subgroup, we multiply it by the cofactor.
        curve_sswu = OptimizedSSWU(extension, coeff_a, 1)
        bin = BinaryStrings()
        sswu_bin_encoding = bin.encoding("Falcon")
        sswu_int = extension(int(str(sswu_bin_encoding), 2))
        g = curve_sswu.map_to_curve(sswu_int)

        if prime_order * g != E(0, 1, 0):
            g = cofactor * g

        (rho_sec, k) = curve_security(
            extension.cardinality(), n, True, prime_order)

        if k.nbits() < MIN_EMBEDDING_DEGREE:
            continue

        sys.stdout.write("+")
        sys.stdout.flush()

        if rho_sec < RHO_SECURITY:
            continue

        twist_rho_sec = twist_security_ignore_embedding_degree(
            extension.cardinality(), n, True)

        if twist_rho_sec < TWIST_SECURITY:
            continue

        yield (extension, E, g, prime_order, cofactor, i, coeff_a, rho_sec, k, twist_rho_sec)


def print_curve(prime, extension_degree, wid=0, processes=1):
    r"""Print parameters of curves defined over a prime field extension

    INPUT:

    - ``prime`` -- the base prime defining Fp
    - ``extension_degree`` -- the targeted extension degree, defining Fp^n on which the curves will be constructed
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    """

    Fp = GF(prime)
    if wid == 0:
        info = f"\n{Fp}.\n"
    Fpx = Fp['x']
    poly = find_sparse_irreducible_poly(Fpx, extension_degree, use_root=True)
    if poly == 0:
        poly_list = find_irreducible_poly(
            Fpx, extension_degree, output_all=True)
        if poly_list == []:
            poly_list = find_irreducible_poly(
                Fpx, extension_degree, use_root=True, output_all=True)
        if poly_list == []:
            raise ValueError(
                'Could not find an irreducible polynomial with specified parameters.')
        poly_list.sort(key=lambda e: poly_weight(e, prime))
        poly = poly_list[0]  # extract the polynomial from the list
    Fp = Fp.extension(poly, "u")
    assert(Fp(-1).is_square())
    if wid == 0:
        info += f"Modulus: {poly}.\n"

    if wid == 0:
        info += f"Looking for Montgomery curves with cofactor 4 or 8.\n"
        print(info)

    for (extension, E, g, order, cofactor, index, coeff_a, rho_security, embedding_degree, twist_rho_security) in find_curve(Fp, wid, processes):
        output = "\n\n\n"
        output += f"E(GF(({extension.base_ring().order().factor()})^{extension.degree()})) : y^2 = x^3 + ({coeff_a})x^2 + x (b == a*4*{index} + 2)\n"
        output += f"\t\twith a = {extension.primitive_element()}\n"
        output += f"E generator point (from SSWU on 'Falcon'): {g}\n"
        output += f"Curve prime order: {order} ({order.nbits()} bits)\n"
        output += f"Curve cofactor: {cofactor}"
        if cofactor > 4:
            output += f" ( = {cofactor % 4} % 4 )"
        output += f"\nCurve security (Pollard-Rho): {'%.2f'%(rho_security)}\n"
        output += f"Curve embedding degree: {embedding_degree} (>2^{embedding_degree.nbits()-1}) \n"
        output += f"Twist security (Pollard-Rho): {'%.2f'%(twist_rho_security)}\n"
        if extension_degree == 6:
            output += f"Curve extension security: ≥ {'%.2f'%(extension_security)}\n\n"
        print(output)
    return


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
        p = x ^ degree
        for n in range(len(coeffs)):
            p += coeffs[n]*x ^ n
        if p.is_irreducible():
            list_poly.append(p)

    if use_root:
        root = ring.base().gen()
        for regular_coeffs in set_coeffs:
            p = x ^ degree
            for n in range(len(regular_coeffs)):
                p += regular_coeffs[n]*x ^ n
            for special_coeffs in set_coeffs:
                q = p
                for n in range(len(special_coeffs)):
                    q += root * special_coeffs[n]*x ^ n
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


########################################################################


def main():
    """Main function"""
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    strategy = print_curve
    help = "--help" in args
    args = [arg for arg in args if not arg.startswith("--")]

    if help:
        print("""
Cmd: sage find_direct_curve_extension.sage [--sequential] [--small-order] <prime> <extension_degree> <max_cofactor> <sswu_string>

Args:
    --sequential        Uses only one process
    <prime>             A prime number, default 2^64 - 2^32 + 1
""")
        return

    prime = int(args[0]) if len(
        args) > 0 else 18446744069414584321  # 2^64 - 2^32 + 1
    extension_degree = 5

    if processes == 1:
        strategy(prime, extension_degree)
    else:
        print(f"Using {processes} processes.")
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(
                    worker, (strategy, prime, extension_degree, wid, processes))

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