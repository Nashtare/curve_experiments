# Based from Daira Hopwood's curvesearch repository
# https://github.com/daira/curvesearch
# and adapted for BLS curves

import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations

from util import *

if sys.version_info[0] == 2: range = xrange

# Find cycle of curves for given prime
# `a` coefficients can be set to 0 when using CM method
def find_curves(p, q):
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

                    #TODO: check if this twist security flaw is by design or we are just unlucky
                    # if twist_sec_p < TWIST_SECURITY:
                    #     continue

                    if twist_sec_q < TWIST_SECURITY:
                        continue

                    yield (p, p_coeff_b, rho_sec_p, k_p, twist_sec_p, twist_k_p, q, q_coeff_b, rho_sec_q, k_q, twist_sec_q, twist_k_q)

# Iterates over BLS generators to try finding cycle of curves including a BLS scalar
def find_cycle(generator_list, wid = 0, processes = 1):
    for (x, x_form, p, q, V, T, q_form) in solve_CM(generator_list, wid, processes):
        sys.stdout.write("o")
        sys.stdout.flush()
        for (p, p_coeff_b, rho_sec_p, k_p, twist_sec_p, twist_k_p, q, q_coeff_b, rho_sec_q, k_q, twist_sec_q, twist_k_q) in find_curves(p, q):
            output = "\n\n\n"
            output += "x = %s\n" % x
            output += "x = %s\n" % x_form
            output += "p = %s - (%d bits)\n" % (p, p.nbits())
            output += "p = %s - (%d 2-adicity)\n" % ("0b" + p.binary(), twoadicity(p))
            output += "q = %s - (%d bits)\n" % (q, q.nbits())
            output += "q = %s - (%d 2-adicity)\n" % ("0b" + q.binary(), twoadicity(q))
            output += "q = %s\n" % q_form
            output += "Ep/Fp : y^2 = x^3 + %d\n" % p_coeff_b
            output += "Embedding degree of Ep/Fp:  %s ( > 2^%d)\n" % (k_p, k_p.nbits() - 1)
            output += "Ep/Fp Pollard Rho security: %s\n" % rho_sec_p
            output += "Ep/Fp Twist security:       %s\n" % twist_sec_p
            output += "Eq/Fq : y^2 = x^3 + %d\n" % q_coeff_b
            output += "Embedding degree of Eq/Fq:  %s ( > 2^%d)\n" % (k_q, k_q.nbits() - 1)
            output += "Eq/Fq Pollard Rho security: %s\n" % rho_sec_q
            output += "Eq/Fq Twist security:       %s\n\n" % twist_sec_q
            print(output)
    return

# Tries solving CM method for BLS scalar fields of high 2-adicity
# Unlikely to work as it is, should be reworked
def solve_CM(generator_list, wid = 0, processes = 1):
    if len(generator_list) > 1:
        V_var,T_var = var('V,T', domain=ZZ)
        for i in range(wid, len(generator_list), processes):
            sys.stdout.write('.')
            sys.stdout.flush()
            x = generator_list[i][0]
            p = bls_scalar(x)
            solutions = solve([4*p == 3*V_var^2 + T_var^2], V_var,T_var)
            # To prevent duplicates with T1 == T2 and V1 == -V2 or reciprocally
            q_set = set()
            for (V,T) in solutions:
                V = Integer(V)
                T = Integer(T)
                q = Integer(p + 1 - T)
                if q in q_set:
                    continue
                q_set.add(q)
                if q != p and q.is_pseudoprime(): # q == p if T == 1
                    yield(x, generator_list[i][1], p, q, V, T, "p+1-T")
                q = Integer(p + 1 + (T - 3*V) // 2)
                if q in q_set:
                    continue
                q_set.add(q)
                if q.is_pseudoprime():
                    yield(x, generator_list[i][1], p, q, V, T, "p+1+(T-3V)/2")
        return
    else:
        p = bls_scalar(generator_list[0])
        L = p.nbits()
        assert(L <= 255)
        adicity = twoadicity(p)
        V_bit_size = (L-1)//2
        Vbase = 1 << V_bit_size
        trailing_zeros = adicity+1
        p4 = 4*p
        for w in range(wid, V_bit_size-trailing_zeros, processes):
            sys.stdout.write('.')
            sys.stdout.flush()
            for Vc in combinations(range(trailing_zeros, V_bit_size), w):
                V = Vbase + sum([1 << i for i in Vc]) + 1
                assert ((V-1)/2) % (1<<adicity) == 0
                T2 = p4 - 3*V**2
                if T2 > 0:
                    T = sqrt(T2)
                    if T in ZZ:
                        print(p,T,V)
        return


def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    jubjub = True if "--jubjub" in args else False
    strategy = find_cycle
    help = True if "--help" in args else False
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
        list_x = [-0xd201000000010000] # BLS generator yielding jubjub base field
    else:
        file_name = str(args[0])
        with open(file_name, 'r') as f1:
            list_x = list(csv.reader(f1))

        for i in range(len(list_x)):
            tmp = list_x[i]
            list_x[i] = [Integer(int(tmp[0])), tmp[2]]

    if processes == 1:
        strategy(list_x)
    else:
        print("Using %d processes." % processes)
        print("Generators list size: %d" % len(list_x))
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(worker, (strategy, list_x, wid, processes))

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

def real_worker(strategy, list_x, wid, processes):
    return strategy(list_x, wid, processes)

main()
