import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from math import ceil
from itertools import combinations, combinations_with_replacement

from util import *

if sys.version_info[0] == 2: range = xrange

def find_curve(extension, wid = 0, processes = 1):
    a = extension.primitive_element()
    for i in range(wid + 1, 10000, processes):
        sys.stdout.write(".")
        sys.stdout.flush()

        coeff = a**i

        E = EllipticCurve(extension, [1, coeff])

        n = E.count_points()
        prime_order = n.factor()[-1][0]
        if prime_order.nbits() < 150:
            continue

        sys.stdout.write("o")
        sys.stdout.flush()
        g = E.random_element()
        while g.order() != prime_order:
            g = E.random_element()
            
        # (rho_sec, k) = curve_security(p, q)

        # if k.nbits() < MIN_EMBEDDING_DEGREE:
        #     continue
        # sys.stdout.write("X")
        # sys.stdout.flush()

        # if rho_sec < RHO_SECURITY:
        #     continue

        return (extension, E, g, n)


# Outputs parameters of valid curves over an extension of F62
def print_curve(prime = 2^62 - 111 * 2^39 + 1, extension_degree = 4, wid = 0, processes = 1):
    extension.<a> = GF(prime^extension_degree, modulus="primitive")
    for (extension, E, g, order) in find_curve(extension, wid, processes):
        output = "\n\n\n"
        output += "E(%s) : y^2 = x^3 + x + %d\n" % (extension, E.a6())
        output += "E generator point: %s\n" % g
        output += "E prime order: %s\n" % order
        output += "E prime order: %s bits\n" % order.nbits()
        print(output)
    return

########################################################################

def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    strategy = print_curve
    help = True if "--help" in args else False
    args = [arg for arg in args if not arg.startswith("--")]

    if help:
        print("""
Cmd: sage find_curve_extension.sage [--sequential] <prime> <extension_degree>

Args:
    --sequential        Uses only one process
""")
        return

    prime = int(args[0]) if len(args) > 0 else 2^62 - 111 * 2^39 + 1
    extension_degree = int(args[1]) if len(args) > 1 else 4
    if processes == 1:
        strategy(prime, extension_degree)
    else:
        print("Using %d processes." % processes)
        pool = Pool(processes=processes)

        try:
            for wid in range(processes):
                pool.apply_async(worker, (strategy, prime, extension_degree, wid, processes))

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
