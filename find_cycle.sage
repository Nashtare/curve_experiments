"""
This module aims at finding cycle of curves.
Based from Daira Hopwood's curvesearch repository
https://github.com/daira/curvesearch
and adapted for BLS curves
"""

import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations

from util import *

if sys.version_info[0] == 2:
    range = xrange


def find_curves(p, q):
    """Yield parameters defining a cycle of curves between two prime fields Fp and Fq

    INPUT:

    - ``p`` -- a prime number, modulus of Fp
    - ``q`` -- a prime number, modulus of Fq

    OUTPUT:

    - ``p`` -- a prime number, modulus of Fp
    - ``b_p`` -- the b coefficient of the curve defined over Fp in short Weierstrass form
    - ``rho_sec_p`` -- the Pollard-Rho security of the curve defined over Fp
    - ``k_p`` -- the embedding degree of the curve defined over Fp
    - ``twist_sec_p`` -- the Twist security of the curve defined over Fp
    - ``twist_k_p`` -- the embedding degree of the twist of the curve defined over Fp
    - ``q`` -- a prime number, modulus of Fq
    - ``b_q`` -- the b coefficient of the curve defined over Fq in short Weierstrass form
    - ``rho_sec_q`` -- the Pollard-Rho security of the curve defined over Fq
    - ``k_q`` -- the embedding degree of the curve defined over Fq
    - ``twist_sec_q`` -- the Twist security of the curve defined over Fq
    - ``twist_k_q`` -- the embedding degree of the twist of the curve defined over Fq

    All found curves have equations in short Weierstrass form E: y^2 = x^3 + b (i.e. a coefficient is zero).

    """

    Fp = GF(p)
    for p_coeff_b in COEFFICIENT_RANGE:
        if Fp(p_coeff_b).is_square():
            # if b is square, the order cannot be prime
            continue
        Ep = EllipticCurve(GF(p), [0, p_coeff_b])
        if q == Ep.count_points():
            sys.stdout.write("#")
            sys.stdout.flush()

            Fq = GF(q)
            for q_coeff_b in COEFFICIENT_RANGE:
                if Fq(q_coeff_b).is_square():
                    continue
                Eq = EllipticCurve(GF(q), [0, q_coeff_b])
                if Eq.count_points() == p:
                    if Mod(p_coeff_b, p).multiplicative_order() != p-1:
                        continue
                    if Mod(q_coeff_b, q).multiplicative_order() != q-1:
                        continue

                    (rho_sec_p, k_p) = curve_security(p, q)
                    (rho_sec_q, k_q) = curve_security(q, p)

                    if k_p.nbits() < MIN_EMBEDDING_DEGREE or k_q.nbits() < MIN_EMBEDDING_DEGREE:
                        continue

                    if rho_sec_p < RHO_SECURITY or rho_sec_q < RHO_SECURITY:
                        continue

                    (twist_sec_p, twist_k_p) = twist_security(p, q)
                    (twist_sec_q, twist_k_q) = twist_security(q, p)

                    # TODO: check if this twist security flaw is by design or we are just unlucky
                    # if twist_sec_p < TWIST_SECURITY:
                    #     continue

                    if twist_sec_q < TWIST_SECURITY:
                        continue

                    yield (p, p_coeff_b, rho_sec_p, k_p, twist_sec_p, twist_k_p, q, q_coeff_b, rho_sec_q, k_q, twist_sec_q, twist_k_q)


def find_cycle(generator_list, is_bls12, wid=0, processes=1):
    """Iterate over BLS generators to print cycles of curves including a BLS scalar field

    INPUT:

    - ``generator_list`` -- list of BLS generators `x` defining a scalar and a base fields
    - ``is_bls12`` -- boolean indicating whether we are using generators for BLS12 or BLS24 curves
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    BLS48 curves are currently not supported.

    """

    for (x, x_form, p, q, V, T, q_form) in solve_CM(generator_list, is_bls12, wid, processes):
        sys.stdout.write("o")
        sys.stdout.flush()
        for (p, p_coeff_b, rho_sec_p, k_p, twist_sec_p, twist_k_p, q, q_coeff_b, rho_sec_q, k_q, twist_sec_q, twist_k_q) in find_curves(p, q):
            output = "\n\n\n"
            output += f"x = {x}\n"
            output += f"x = {x_form}\n"
            output += f"p = {p} - ({p.n_bits()} bits)\n"
            output += f"p = 0b{p.binary()} - ({twoadicity(p)} 2-adicity)\n"
            output += f"q = {q} - ({q.nbits()} bits)\n"
            output += f"q = {q.binary()} - ({twoadicity(q)} 2-adicity)\n"
            output += f"q = {q_form}\n"
            output += f"Ep/Fp : y^2 = x^3 + {p_coeff_b}\n"
            output += f"Embedding degree of Ep/Fp:  {k_p} ( > 2^{k_p.nbits() - 1})\n"
            output += f"Ep/Fp Pollard Rho security: {rho_sec_p}\n"
            output += f"Ep/Fp Twist security:       {twist_sec_p}\n"
            output += f"Eq/Fq : y^2 = x^3 + {q_coeff_b}\n"
            output += f"Embedding degree of Eq/Fq:  {k_q} ( > 2^{k_q.nbits() - 1})\n"
            output += f"Eq/Fq Pollard Rho security: {rho_sec_q}\n"
            output += f"Eq/Fq Twist security:       {twist_sec_q}\n\n"
            print(output)
    return


def solve_CM(generator_list, is_bls12=True, wid=0, processes=1):
    """Iterate over BLS generators and solve the complex multiplication equation

    INPUT:

    - ``generator_list`` -- list of BLS generators `x` defining a scalar and a base fields
    - ``is_bls12`` -- boolean indicating whether we are using generators for BLS12 or BLS24 curves
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    BLS48 curves are currently not supported.

    """

    scalar_func = bls12_scalar if is_bls12 else bls24_scalar
    if len(generator_list) > 1:
        V_var, T_var = var('V,T', domain=ZZ)
        for i in range(wid, len(generator_list), processes):
            sys.stdout.write('.')
            sys.stdout.flush()
            x = generator_list[i][0]
            p = scalar_func(x)
            solutions = solve([4*p == 3*V_var ^ 2 + T_var ^ 2], V_var, T_var)
            # To prevent duplicates with T1 == T2 and V1 == -V2 or reciprocally
            q_set = set()
            for (V, T) in solutions:
                V = Integer(V)
                T = Integer(T)
                q = Integer(p + 1 - T)
                if q in q_set:
                    continue
                q_set.add(q)
                if q != p and q.is_pseudoprime():  # q == p if T == 1
                    yield(x, generator_list[i][1], p, q, V, T, "p+1-T")
                q = Integer(p + 1 + (T - 3*V) // 2)
                if q in q_set:
                    continue
                q_set.add(q)
                if q.is_pseudoprime():
                    yield(x, generator_list[i][1], p, q, V, T, "p+1+(T-3V)/2")
    else:
        p = bls12_scalar(generator_list[0])
        p_len = p.nbits()
        assert p_len <= 255
        adicity = twoadicity(p)
        V_bit_size = (p_len-1)//2
        Vbase = 1 << V_bit_size
        trailing_zeros = adicity+1
        p4 = 4*p
        for w in range(wid, V_bit_size-trailing_zeros, processes):
            sys.stdout.write('.')
            sys.stdout.flush()
            for Vc in combinations(range(trailing_zeros, V_bit_size), w):
                V = Vbase + sum([1 << i for i in Vc]) + 1
                assert ((V-1)/2) % (1 << adicity) == 0
                T2 = p4 - 3*V**2
                if T2 > 0:
                    T = sqrt(T2)
                    if T in ZZ:
                        print(p, T, V)
    # TODO: Handle end signal in main()
    while True:
        sys.stdout.write('x')
        sys.stdout.flush()
        sleep(30)


########################################################################

def main():
    """Main function"""
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    jubjub = "--jubjub" in args
    strategy = find_cycle
    help = "--help" in args
    args = [arg for arg in args if not arg.startswith("--")]

    if (not jubjub and len(args) < 1) or help:
        print("""
Cmd: sage find_cycle.sage [--jubjub] [--sequential]
                          [<filename>]

Args:
    --jubjub        Tries finding a cycle with jubjub
    --sequential    Uses only one process
    <filename>      File listing BLS generators

""")
        return

    if jubjub and len(args) >= 1:
        print("Invalid arguments. Type `sage find_cycle.sage --help` for help")

    if jubjub:
        # BLS generator yielding Jubjub base field
        list_x = [-0xd201000000010000]
    else:
        file_name = str(args[0])
        is_bls12 = "bls12" in file_name
        with open(file_name, 'r') as f1:
            list_x = list(csv.reader(f1))

        for i in range(len(list_x)):
            tmp = list_x[i]
            list_x[i] = [Integer(int(tmp[0])), tmp[2]]

    if processes == 1:
        strategy(list_x)
    else:
        print(f"Using {processes} processes.")
        print(f"Generators list size: {len(list_x)}")
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(
                    worker, (strategy, list_x, is_bls12, wid, processes))

            while True:
                sleep(1000)
        except (KeyboardInterrupt, SystemExit):
            pass
        finally:
            pool.terminate()


def worker(*args):
    try:
        real_worker(*args)
    except (KeyboardInterrupt, SystemExit):
        pass
    except:
        print_exc()


def real_worker(strategy, list_x, is_bls12, wid, processes):
    return strategy(list_x, is_bls12, wid, processes)


main()
