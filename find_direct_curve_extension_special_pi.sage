"""
This module aims at finding curves defined over finite field extensions.
"""

import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations_with_replacement

from util import *

if sys.version_info[0] == 2:
    range = xrange

B5 = "3141592653589793238"
B4 = "4626433832795028841"
B3 = "9716939937510582097"
B2 = "4944592307816406286"
B1 = "2089986280348253421"
B0 = "1706798214808651328"
PI_STRING = "314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328"


def find_curve(extension, min_cofactor, max_cofactor, small_order, wid=0, processes=1):
    r"""Yield curve constructed over a prime field extension.

    INPUT:

    - ``extension`` -- the field direct extension
    - ``min_cofactor`` -- the minimum cofactor for the curve order
    - ``max_cofactor`` -- the maximum cofactor for the curve order
    - ``small_order`` -- boolean indicating whether to look for small orders (252/255 bits).
            Overrides `min_cofactor` and `max_cofactor` if set to `True`.
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
    - ``coeff_b`` -- the b coefficient of the curve in short Weierstrass form
    - ``rho_sec`` -- the Pollard-Rho security of the curve
    - ``k`` -- the embedding degree of the curve
    - ``extension_sec`` -- the extension-specific attacks security of the curve
    - ``twist_rho_sec`` -- the Pollard-Rho security of the twist

    """

    a = extension.primitive_element()
    p = extension.base_ring().order()
    for length in range(wid + 18, 19 * 6 + 1, processes):
        sys.stdout.write(".")
        sys.stdout.flush()

        for coeff_a in range(-3, 4):
            string = PI_STRING[0:length]
            while len(string) < 19 * 6:
                string = "0" + string
            coeff_b = extension([int(string[5*19:6*19]), int(string[4*19:5*19]), int(
                string[3*19:4*19]), int(string[2*19:3*19]), int(string[19:2*19]), int(string[0:19])])

            E = EllipticCurve(extension, [coeff_a, coeff_b])

            n = E.count_points()
            prime_order = list(ecm.factor(n))[-1]
            cofactor = n // prime_order
            if small_order:
                if prime_order.nbits() < 252 or prime_order.nbits() > 255:
                    continue
            elif cofactor < min_cofactor:
                continue
            elif cofactor > max_cofactor:
                continue

            sys.stdout.write("o")
            sys.stdout.flush()

            # TODO: use proper hash-to-curve algorithm
            bin = BinaryStrings()
            gen_x_bin = bin.encoding("Topos")
            gen_x = extension(int(str(gen_x_bin), 2))
            gen_y2 = (gen_x ^ 3 + coeff_a * gen_x + coeff_b)
            while True:
                if gen_y2.is_square():
                    g = E((gen_x, gen_y2.sqrt()))
                    if cofactor * g != E(0, 1, 0):  # ord(g) >= prime_order
                        sys.stdout.write("@")
                        sys.stdout.flush()
                        break
                gen_x += 1
                gen_y2 = (gen_x ^ 3 + coeff_a * gen_x + coeff_b)

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

            if extension.degree() == 6:
                sys.stdout.write("~")
                sys.stdout.flush()

                extension_sec = degree_six_security(
                    extension, p, E, n)
                if extension_sec < EXTENSION_SECURITY:
                    continue
            else:
                extension_sec = 0

            twist_rho_sec = twist_security_ignore_embedding_degree(
                extension.cardinality(), n, True)

            if twist_rho_sec < TWIST_SECURITY:
                continue

            yield (extension, E, g, prime_order, cofactor, i, coeff_a, coeff_b, rho_sec, k, extension_sec, twist_rho_sec)


def print_curve(prime, extension_degree, min_cofactor, max_cofactor, small_order, wid=0, processes=1):
    r"""Print parameters of curves defined over a prime field extension

    INPUT:

    - ``prime`` -- the base prime defining Fp
    - ``extension_degree`` -- the targeted extension degree, defining Fp^n on which the curves will be constructed
    - ``min_cofactor`` -- the minimum cofactor for the curve order
    - ``max_cofactor`` -- the maximum cofactor for the curve order
    - ``small_order`` -- boolean indicating whether to look for small orders (252/255 bits).
            Overrides `min_cofactor` and `max_cofactor` if set to `True`.
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
        if small_order:
            info += f"Looking for curves with 252-bit to 255-bit prime order.\n"
        else:
            info += f"Looking for curves with max cofactor: {max_cofactor}.\n"
        print(info)

    for (extension, E, g, order, cofactor, index, coeff_a, coeff_b, rho_security, embedding_degree, extension_security, twist_rho_security) in find_curve(Fp, min_cofactor, max_cofactor, small_order, wid, processes):
        output = "\n\n\n"
        output += f"E(GF(({extension.base_ring().order().factor()})^{extension.degree()})) : y^2 = x^3 + {coeff_a}x + {coeff_b} (b == a^{index})\n"
        output += f"\t\twith a = {extension.primitive_element()}\n"
        output += f"E generator point: {g}\n"
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


def degree_six_security(field, p, E, curve_order):
    assert field.order() == p ^ 6

    # Construct a tower extension isomorphic to field
    Fp = GF(p)
    Fpx = Fp["x"]
    poly = find_irreducible_poly(Fpx, 2, use_root=True)[0]
    Fp = Fp.extension(poly, "a1")
    Fpx = Fp["x"]
    poly = find_irreducible_poly(Fpx, 3)[0]
    Fp = Fp.extension(poly, "a2")

    field_bis, _, psi2 = make_finite_field(Fp)
    psi1 = field.Hom(field_bis)[0]
    assert(psi1.is_injective())
    assert(psi1.is_surjective())
    psi = psi1.post_compose(psi2)
    a = E.a4()
    b = E.a6()

    K = Fp["x"]
    x = K.gen()
    curve_polynomial = K(x ^ 3 + psi(a)*x + psi(b))

    # Sieving/Decomp. on Jac_H(\mathbb{F}_{p^2}), g = 3

    #   - hyperelliptic case
    if curve_order % 4 == 0:
        return p.nbits() * 5.0/3

    #   - non-hyperelliptic case
    # Heavier cost than Weil descent attack, so discarded here

    # Decomp. on Jac_H(\mathbb{F}_{p^3}), g = 2

    if curve_order % 2 == 1 and True:  # no j-invariant check as we can always convert to Sholten
        return p.nbits() * 12.0/7  # form with the current construction if curve_order % 2 == 1
    elif E.two_torsion_rank() == 2:
        return p.nbits() * 12.0/7

    # Ind. calc. on Jac_C(\mathbb{F}_p), d = 10
    # More than 2^128 operations: see https://eprint.iacr.org/2014/346.pdf

    # GHS method, with either Fp2 or Fp3 as basefield
    roots = curve_polynomial.roots(multiplicities=false)
    if roots != []:
        for root in roots:
            if (root ** (p**2) in roots) or (root ** (p**3) in roots):
                return p.nbits() * 8.0/3

    # Ind. calc. on Jac_H(\mathbb{F}_{p^2}), g = 3
    # Heavier cost than Weil descent attack, so discarded here

    # Weil descent attack with genus 2
    g = 2
    return log((g ^ 2 * log(p, 2) ^ 3) * factorial(g) * p + (g ^ 2 * log(p, 2)) * p ^ 2, 2)

########################################################################


def main():
    """Main function"""
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    small_order = "--small-order" in args
    strategy = print_curve
    help = "--help" in args
    args = [arg for arg in args if not arg.startswith("--")]

    if help:
        print("""
Cmd: sage find_curve_extension.sage [--sequential] [--small-order] <prime> <extension_degree> <max_cofactor>

Args:
    --sequential        Uses only one process
    --small-order       Looks for curves with 252-bit to 255-bit prime order (overrides cofactor)
    <prime>             A prime number, default 2^62 + 2^56 + 2^55 + 1
    <extension_degree>  The extension degree of the prime field, default 6
    <max_cofactor>      Maximum cofactor of the curve, default 64
""")
        return

    prime = int(args[0]) if len(
        args) > 0 else 4719772409484279809  # 2^62 + 2^56 + 2^55 + 1
    extension_degree = int(args[1]) if len(args) > 1 else 6
    min_cofactor = 0
    max_cofactor = int(args[2]) if len(args) > 2 else 64

    if processes == 1:
        strategy(prime, extension_degree,
                 min_cofactor, max_cofactor, small_order)
    else:
        print(f"Using {processes} processes.")
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(
                    worker, (strategy, prime, extension_degree, min_cofactor, max_cofactor, small_order, wid, processes))

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
