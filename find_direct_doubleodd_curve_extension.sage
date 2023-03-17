"""
This module aims at finding double-odd curves defined over finite field extensions.
"""

from sage.misc.verbose import set_verbose

import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations_with_replacement

from util import *

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
    - ``c`` -- the a4 coefficient of the curve in long Weierstrass form
    - ``rho_sec`` -- the Pollard-Rho security of the curve
    - ``k`` -- the embedding degree of the curve
    - ``twist_rho_sec`` -- the Pollard-Rho security of the twist

    """

    p = extension.base_ring().order()
    for i in range(wid + 1, 1000000000, processes):
        sys.stdout.write(".")
        sys.stdout.flush()

        c = i
        z = extension.gen()
        for j in range(1, extension.degree()):
            zi = z**j
            try:
                E = EllipticCurve(extension, [0, 2, 0, c*zi, 0])
            except:
                continue

            n = E.count_points()
            prime_order = list(ecm.factor(n))[-1]
            cofactor = n // prime_order
            if cofactor != 2:
                continue

            sys.stdout.write("o")
            sys.stdout.flush()

            g = E.random_point()

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

            yield (extension, E, g, prime_order, c*zi, rho_sec, k, twist_rho_sec)

        for j in range(1, extension.degree()):
            zi = z**j
            try:
                E = EllipticCurve(extension, [0, 2, 0, -c*zi, 0])
            except:
                continue

            n = E.count_points()
            prime_order = list(ecm.factor(n))[-1]
            cofactor = n // prime_order
            if cofactor != 2:
                continue

            sys.stdout.write("o")
            sys.stdout.flush()

            g = E.random_point()

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

            yield (extension, E, g, prime_order, -c*zi, rho_sec, k, twist_rho_sec)


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
    if wid == 0:
        info += f"Modulus: {poly}.\n"

    if wid == 0:
        print(info)

    for (extension, E, g, order, coeff_c, rho_security, embedding_degree, twist_rho_security) in find_curve(Fp, wid, processes):
        output = "\n\n\n"
        output += f"E(GF(({extension.base_ring().order().factor()})^{extension.degree()})) : y^2 = x.(x^2 + 2x + {coeff_c})\n"
        output += f"E generator point: {g}\n"
        output += f"Curve prime order: {order} ({order.nbits()} bits)\n"
        output += f"\nCurve security (Pollard-Rho): {'%.2f'%(rho_security)}\n"
        output += f"Curve embedding degree: {embedding_degree} (>2^{embedding_degree.nbits()-1}) \n"
        output += f"Twist security (Pollard-Rho): {'%.2f'%(twist_rho_security)}\n"
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
Cmd: sage find_direct_curve_extension.sage [--sequential] <prime> <extension_degree>

Args:
    --sequential        Uses only one process
    --small-order       Looks for curves with 252-bit to 255-bit prime order (overrides cofactor)
    <prime>             A prime number, default 2^64 - 2^32 + 1
    <extension_degree>  The extension degree of the prime field, default 5
""")
        return

    prime = int(args[0]) if len(
        args) > 0 else 2**64 - 2**32 + 1
    extension_degree = int(args[1]) if len(args) > 1 else 5

    # Silence warnings about groebner_basis toy implementation
    set_verbose(-1)

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
