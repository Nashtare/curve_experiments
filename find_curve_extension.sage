import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from math import ceil

from util import *

if sys.version_info[0] == 2: range = xrange

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

        if rho_sec < RHO_SECURITY - 5:
            continue

        yield (extension, E, g, prime_order, cofactor, i, coeff_a, coeff_b, rho_sec, k)


# Outputs parameters of valid curves over an extension of F62
def print_curve(prime = 2^62 - 111 * 2^39 + 1, extension_degree = 4, max_cofactor = 256, wid = 0, processes = 1):
    extension.<a> = GF(prime^extension_degree, modulus="primitive")
    for (extension, E, g, order, cofactor, index, coeff_a, coeff_b, rho_security, embedding_degree) in find_curve(extension, max_cofactor, wid, processes):
        output = "\n\n\n"
        output += "E(%s) : y^2 = x^3 + %s * x + %s (b == a^%s)\n" % (extension.order().factor(), coeff_a, coeff_b, index)
        output += "E generator point: %s\n" % g
        output += "E prime order: %s (%s bits)\n" % (order, order.nbits())
        output += "E cofactor: %s\n" % cofactor
        output += "E security (Pollard-Rho): %s\n" % rho_security
        output += "E embedding degree: %s\n" % embedding_degree
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
Cmd: sage find_curve_extension.sage [--sequential] <prime> <extension_degree> <max_cofactor>

Args:
    --sequential        Uses only one process
    <prime>             A prime number, default 2^62 - 111 * 2^39 + 1
    <extension_degree>  The extension degree of the prime field, default 4
    <max_cofactor>      Maximum cofactor of the curve, default 256
""")
        return

    prime = int(args[0]) if len(args) > 0 else 2^62 - 111 * 2^39 + 1
    extension_degree = int(args[1]) if len(args) > 1 else 4
    max_cofactor = int(args[2]) if len(args) > 2 else 256

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
