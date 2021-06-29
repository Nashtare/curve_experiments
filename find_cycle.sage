import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations

from util import *

if sys.version_info[0] == 2: range = xrange

# Find cycle of curves for given prime
# `a` coefficients can be set to 0 when using CM method
def find_curves(p):
    Fp = GF(p)
    for p_coeff_b in COEFFICIENT_RANGE:
        if Fp(p_coeff_b).is_square():
            # if b is square, the order cannot be prime
            continue
        Ep = EllipticCurve(GF(p), [0, p_coeff_b])
        q = Ep.count_points()
        if not is_pseudoprime(q):
            continue

        sys.stdout.write("#")
        sys.stdout.flush()

        Fq = GF(q)
        for q_coeff_b in COEFFICIENT_RANGE:
            if Fp(q_coeff_b).is_square():
                continue
            Eq = EllipticCurve(GF(p), [0, q_coeff_b])
            if Eq.count_points() % p == 0:
                yield (p, p_coeff_b, q, q_coeff_b)

# Iterates over BLS generators to try finding cycle of curves including a BLS scalar
def find_cycle(generator_list, wid = 0, processes = 1):
    end = len(generator_list)
    for index in range(0 + wid, end, processes):
        sys.stdout.write(".")
        sys.stdout.flush()
        x = generator_list[index]
        p = bls_scalar(x)
        for (p, p_coeff_b, q, q_coeff_b) in find_curves(p):
            output = "\n"
            output += "x = %s\n" % x
            output += "p = %s\n" % p
            output += "q = %s\n" % q
            output += "Ep/Fp : y^2 = x^3 + %d\n" % p_coeff_b
            output += "Eq/Fq : y^2 = x^3 + %d\n" % q_coeff_b
            print(output)
    return

# Tries solving CM method for BLS scalar fields of high 2-adicity
# Unlikely to work as it is, should be reworked
def solve_CM(generator_list, wid = 0, processes = 1):
    if len(generator_list) > 1:
        for i in range(wid, len(generator_list), processes):
            sys.stdout.write('.')
            sys.stdout.flush()
            x = generator_list[i]
            p = bls_scalar(x)
            L = p.nbits()
            adicity = twoadicity(p)
            V_bit_size = (L-1)//2
            Vbase = 1 << V_bit_size
            trailing_zeros = adicity+1
            p4 = 4*p
            for w in range(0, V_bit_size-trailing_zeros):
                for Vc in combinations(range(trailing_zeros, V_bit_size), w):
                    V = Vbase + sum([1 << i for i in Vc]) + 1
                    assert ((V-1)/2) % (1<<adicity) == 0
                    T2 = p4 - 3*V**2
                    if T2 > 0:
                        T = sqrt(T2)
                        if T in ZZ:
                            print(p,T,V)
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

def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    jubjub = True if "--jubjub" in args else False
    strategy = solve_CM if "--cm-method" in args else find_cycle
    help = True if "--help" in args else False
    args = [arg for arg in args if not arg.startswith("--")]

    if (not jubjub and len(args) < 1) or help:
        print("""
Cmd: sage find_cycle.sage [--jubjub] [--solve-cm] [--sequential]
                          [<filename>]

Args:
    --jubjub        Tries finding a cycle with jubjub
    --solve-cm      Tries solving the CM equation
    --sequential    Uses only one process
    <filename>      File listing BLS generators

--jubjub can be combined with --solve-cm to try solving the CM
equation with jubjub, or given as standalone for exhaustive search.
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
            list_x[i] = Integer(int(tmp[0]))

    if processes == 1:
        strategy(list_x)
    else:
        print("Using %d processes." % (processes,))
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
