import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations

if sys.version_info[0] == 2: range = xrange

LOW_PRIMES = (3, 5, 7, 11, 13, 17)
COEFFICIENT_RANGE = range(1, 20)
ACCEPTABLE_PRIMES = Primes()
ISOGENY_DEGREE_MAX = 3

# BLS scalar field
def bls_scalar(x):
    return cyclotomic_polynomial(12)(x)

# BLS base field
def bls_base(x, r = 0):
    if r == 0:
        tmp = (x-1)^2 * cyclotomic_polynomial(12)(x) / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    else:
        tmp = (x-1)^2 * r / 3 + x
        if tmp.is_integer():
            return Integer(tmp)
    return 0

def find_curves(p):
    for p_coeff_b in COEFFICIENT_RANGE:
        for p_coeff_a in COEFFICIENT_RANGE:
            Ep = EllipticCurve(GF(p), [p_coeff_a, p_coeff_b])
            q = Ep.count_points()
            if not is_pseudoprime(q):
                continue
            sys.stdout.write("#")
            sys.stdout.flush()

            for q_coeff_b in COEFFICIENT_RANGE:
                for q_coeff_a in COEFFICIENT_RANGE:
                    Eq = EllipticCurve(GF(p), [q_coeff_a, q_coeff_b])
                    if Eq.count_points() == p:
                        yield (p_coeff_a, p_coeff_b, q, q_coeff_a, q_coeff_b)

def find_cycle(generator_list, wid = 0, processes = 1):
    end = len(generator_list)
    for index in range(0 + wid, end, processes):
        sys.stdout.write(".")
        sys.stdout.flush()
        x = generator_list[index]
        p = bls_scalar(x)
        for (p_coeff_a, p_coeff_b, q, q_coeff_a, q_coeff_b) in find_curves(p):
            output = "\n"
            output += "x = %s\n" % x
            output += "p = %s\n" % p
            output += "q = %s\n" % q
            output += "Ep/Fp : y^2 = x^3 + %dx + %d\n" % (p_coeff_a, p_coeff_b)
            output += "Eq/Fq : y^2 = x^3 + %dx + %d\n" % (q_coeff_a, q_coeff_b)
            print(output)
    return

def twoadicity(x, base=0, limit=64):
    return max(i for i in range(base, limit) if ((x-1) % (1<<i) == 0))

def solve_CM(generator_list, wid = 0, processes = 1):
    for i in range(wid, len(generator_list), processes):
        sys.stdout.write('.')
        sys.stdout.flush()
        x = generator_list[i]
        p = bls_scalar(x)
        L = p.nbits()
        adicity = twoadicity(p-1)
        Vlen = (L-1)//2
        Vbase = 1 << Vlen
        trailing_zeros = adicity+1
        p4 = 4*p
        for w in range(0, Vlen-trailing_zeros):
            for Vc in combinations(range(trailing_zeros, Vlen), w):
                V = Vbase + sum([1 << i for i in Vc]) + 1
                assert ((V-1)/2) % (1<<adicity) == 0
                T2 = p4 - 3*V^2
                if T2 > 0:
                    T = sqrt(T2)
                    if T in ZZ:
                        print(p,T,V)

def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1 or help:
        print("""
Cmd: sage find_cycle_given_prime.sage [--sequential]
                                      <filename>

Args:
    --sequential    Uses only one process
    <filename>      File listing BLS generators
""")
        return

    file_name = str(args[0])

    with open(file_name, 'r') as f1:
        list_x = list(csv.reader(f1))

    for i in range(len(list_x)):
        tmp = list_x[i]
        list_x[i] = Integer(int(tmp[0]))

    if processes == 1:
        result_list = solve_CM(list_x)
    else:
        print("Using %d processes." % (processes,))
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(worker, (list_x, wid, processes))
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

def real_worker(*args):
    return solve_CM(*args)

main()
