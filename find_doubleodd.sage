"""
This module aims at finding doubleodd curves.
Based from Daira Hopwood's curvesearch repository
https://github.com/daira/curvesearch
and adapted for doubleodd curves over primefields.
"""

import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations

from util import *


if sys.version_info[0] == 2:
    range = xrange


def find_curve(p, wid=0, processes=1):
    r"""Yield parameters defining a doubleodd curve over primefield Fp

    INPUT:

    - ``p`` -- a prime number, modulus of Fp

    OUTPUT:

    - ``p`` -- a prime number, modulus of Fp
    - ``q`` -- the prime order of the large subgroup (curve order = 2q)
    - ``a`` -- the a coefficient of the curve defined over Fp in fractional representation, always 2 here
    (so that alpha = 2 and beta = 0)
    - ``b`` -- the b coefficient of the curve defined over Fp in fractional representation
    - ``rho_sec`` -- the Pollard-Rho security of the curve
    - ``k`` -- the embedding degree of the curve
    - ``twist_sec`` -- the Twist security of the curve
    - ``twist_k`` -- the embedding degree of the twist of the curve

    All found curves have equations of the form E: y^2 = x(x^2 + 2x + b).

    """

    Fp = GF(p)
    a = 2
    for b in range(wid + 1, 1000000, processes):
        if Fp(b).is_square():
            continue
        if Fp(a**2 - 4*b).is_square():
            continue
        sys.stdout.write(".")
        sys.stdout.flush()
        # Weierstrass coefficients
        A = Fp((3*b - a^2) / 3)
        B = Fp(a*(2*a^2 - 9*b) / 27)
        try:
            Ep = EllipticCurve(GF(p), [A, B])
        except:
            continue
        q = Ep.count_points()
        if q % 2 != 0:
            continue
        q = q/2
        if not q.is_prime():
            continue
        
        sys.stdout.write("#")
        sys.stdout.flush()

        (rho_sec, k) = curve_security(p, q)

        if k.nbits() < MIN_EMBEDDING_DEGREE:
            continue

        if rho_sec < RHO_SECURITY:
            continue

        sys.stdout.write("o")
        sys.stdout.flush()

        (twist_sec, twist_k) = twist_security(p, q)

        if twist_sec < TWIST_SECURITY:
            continue

        yield (p, q, a, b, rho_sec, k, twist_sec, twist_k)
    return

def print_curve(p, wid=0, processes=1):
    """Iterate over BLS generators to print cycles of curves including a BLS scalar field

    INPUT:

    - ``p`` -- the prime characteristic of the basefield
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    """

    for (p, q, a, b, rho_sec, k, twist_sec, twist_k) in find_curve(p, wid, processes):
        output = "\n\n\n"
        output += f"p = {p} - ({p.n_bits()} bits)\n"
        output += f"p = 0b{p.binary()} - ({twoadicity(p)} 2-adicity)\n"
        output += f"q = {q} - ({q.nbits()} bits)\n"
        output += f"q = {q.binary()})\n"
        output += f"E/Fp : y^2 = x(x^2 + {a}x + {b})\n"
        output += f"E/Fp : y^2 = x^3 + {E.a4}x + {E.a6} (short Weierstrass form)\n"
        output += f"Embedding degree of E/Fp:  {k} ( > 2^{k.nbits() - 1})\n"
        output += f"E/Fp Pollard Rho security: {rho_sec}\n"
        output += f"E/Fp Twist security:       {twist_sec}\n"
        print(output)
    return

########################################################################

def main():
    """Main function"""
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    help = "--help" in args
    args = [arg for arg in args if not arg.startswith("--")]

    if (len(args) < 1) or help:
        print("""
Cmd: sage find_doubleodd.sage [--sequential] prime

Args:
    --sequential    Uses only one process
    <prime>         The prime used to define the prime basefield of the targeted curve
""")
        return

    p = args[0]

    if processes == 1:
        return print_curve(p)
    else:
        print(f"Using {processes} processes.")
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(worker, (p, wid, processes))

            while True:
                sleep(1000)
        except (KeyboardInterrupt, SystemExit):
            pass
        finally:
            pool.terminate()


def worker(*args):
    res = []
    try:
        res = real_worker(*args)
    except (KeyboardInterrupt, SystemExit):
        pass
    except:
        print_exc()
    finally:
        return res


def real_worker(*args):
    return print_curve(*args)


main()
