"""
This module aims at finding curves defined over non highly 2-adic fields but with
a smooth subgroup, for ECFFT applications.
"""

import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from itertools import combinations_with_replacement

from util import *
from util_hashtocurve import OptimizedSSWU

if sys.version_info[0] == 2:
    range = xrange


def find_curve(field, adicity=15, wid=0, processes=1):
    r"""Yield curve with highly 2-adic group.

    INPUT:

    - ``field`` -- the basefield
    - ``adicity`` -- the required adicity of the curve's order
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    OUTPUT:

    - ``E`` -- the curve definition
    - ``L`` -- the factorisation list of the curve group order

    """

    coeff_a = 1

    for i in range(wid + 1, 1000000000, processes):
        sys.stdout.write(".")
        sys.stdout.flush()

        coeff_b = field(i)

        E = EllipticCurve(field, [coeff_a, coeff_b])

        n = E.count_points()
        if n % 2**adicity != 0:
            continue

        sys.stdout.write("o")
        sys.stdout.flush()

        yield (E, n.factor())


def print_curve(prime, adicity, curve_name, wid=0, processes=1):
    r"""Print parameters of curves defined over a prime field extension

    INPUT:

    - ``prime`` -- the base prime defining Fp
    - ``adicity`` -- the required adicity of the curve's order
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)

    """

    if wid == 0:
        if curve_name != "":
            print(
                f"Searching curves over {curve_name} basefield with curve's order 2-adicity at least {adicity}.")
        else:
            print(
                f"Searching curves over basefield of characteristic {prime} with curve's order 2-adicity at least {adicity}.")

    Fp = GF(prime)

    for (E, L) in find_curve(Fp, adicity, wid, processes):
        output = f"\n{E}\n"
        output += f"Curve's order factorisation: {L}\n"
        print(output)
    return


def main():
    """Main function"""
    args = sys.argv[1:]
    processes = cpu_count()

    prime = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787  # BLS12-381 basefield
    curve_name = "bls12_381"

    if "--secp256k1" in args:
        prime = 115792089237316195423570985008687907853269984665640564039457584007908834671663
        curve_name = "secp256k1"
    elif "--secp256r1" in args:
        prime = 115792089210356248762697446949407573530086143415290314195533631308867097853951
        curve_name = "secp256r1"
    elif "--curve25519" in args:
        prime = 57896044618658097711785492504343953926634992332820282019728792003956564819949
        curve_name = "curve25519"
    elif "--altbn_128" in args:
        prime = 21888242871839275222246405745257275088696311157297823662689037894645226208583
        curve_name = "altbn_128"

    strategy = print_curve
    help = "--help" in args
    args = [arg for arg in args if not arg.startswith("--")]

    if help:
        print("""
Cmd: sage find_ecfft.sage [--secp256k1] [--curve25519] [--altbn128] <2adicity> <prime>

Args:
    --secp256k1         Uses secp256k1 basefield. Overrides ``prime``.
    --secp256r1         Uses secp256r1 basefield. Overrides ``prime``.
    --curve25519        Uses curve25519 basefield. Overrides ``prime``.
    --altbn_128         Uses altbn_128 basefield. Overrides ``prime``.
    <2adicity>          The desired 2-adicity of the curve's order, default 15.
    <prime>             A prime number, default BLS12-381 basefield.
""")
        return

    adicity = int(args[0]) if len(args) > 0 else 15

    if len(args) > 1:
        prime = int(args[1])
        curve_name = ""

    assert prime.is_prime()

    print(f"Using {processes} processes.")
    pool = Pool(processes=processes)

    try:
        for wid in range(processes):
            pool.apply_async(
                worker, (strategy, prime, adicity, curve_name, wid, processes))

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
