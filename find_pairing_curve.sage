"""
This module aims at finding pairing-friendly curves defined over finite field extensions.
"""

import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc

from util import *

if sys.version_info[0] == 2:
    range = xrange


def find_curve(p, m=6, embedding_degree=12, wid=0, processes=1):
    q = Integer(p ** m)
    Rq = Integers(q)

    for D in range(-3 - wid, -1000000000, -processes):
        if -D >= 4*q:
            break
        if (D % 4 == 2) or (D % 4 == 3):
            continue

        # D = 0 or 1 mod 4
        h = kronecker(D, p)
        if h == -1:
            continue

        v = Rq(D).sqrt()
        if v not in Rq:
            continue
        if (ZZ(v) % 2 != ZZ(D) % 2):
            v = q - v

        a = ZZ(2*q)
        b = ZZ(v)
        l = floor(2 * q.sqrt())

        while b > l:
            a, b = b, a % b

        c = (4*q - b ** 2) // abs(D)
        if ((4*q - b**2) % abs(D) != 0) or Rq(c).sqrt() not in Rq:
            continue

        t, y = b, c

        sys.stdout.write(".")
        sys.stdout.flush()
        phi = cyclotomic_polynomial(embedding_degree)(t-1)
        if phi.is_prime():
            r = phi
        else:
            r = ecm.factor(phi)[-1]
        if (q + 1 - t) % ZZ(r) == 0:
            yield(p, q, m, t, r, embedding_degree)


def print_curve(prime, m=6, embedding_degree=12, wid=0, processes=1):
    for (p, q, m, t, r, k) in find_curve(prime, m, embedding_degree, wid, processes):
        output = "\n\n\n"
        output += f"p = {p}\n"
        output += f"q = p^{m}\n"
        output += f"t = {t}\n"
        output += f"r = {r}\n"
        output += f"k = {k}\n"
        print(output)
    return

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
Cmd: sage find_pairing_curve.sage [--sequential] <prime> <extension_degree> <embedding_degree>

""")
        return

    prime = int(args[0]) if len(
        args) > 0 else 18446744069414584321
    extension_degree = int(args[1]) if len(args) > 1 else 6
    embedding_degree = int(args[1]) if len(args) > 1 else 12

    if processes == 1:
        strategy(prime, extension_degree)
    else:
        print(f"Using {processes} processes.")
        pool = Pool(processes=processes)

        # Simpler to manage for 2n parallel jobs
        assert(processes % 2 == 0)
        wid_list = [0, 1]
        for i in range((processes // 2) - 1):
            wid_list.append(wid_list[2*i] + 4)
            wid_list.append(wid_list[2*i + 1] + 4)

        try:
            for wid in wid_list:
                pool.apply_async(
                    worker, (strategy, prime, extension_degree, embedding_degree, wid, processes))

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
