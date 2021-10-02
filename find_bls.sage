import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from math import ceil
from itertools import combinations, combinations_with_replacement

from util import *

if sys.version_info[0] == 2:
    range = xrange


def find_BLS12_curve(adicity, weight_start=2, weight_end=8, conservative=False, wid=0, processes=1, extended=False, verbose=True):
    """Return parameters of valid BLS12 curves

    INPUT:

    - ``adicity`` -- minimum 2-adicity of the scalar field
    - ``weight_start`` -- minimum Hamming weight for generator `x` (default 2)
    - ``weight_end`` -- maximum Hamming weight for generator `x` (default 8)
    - ``conservative`` -- boolean indicating whether selecting conservative parameters or not (default False)
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)
    - ``extended`` -- boolean indicating whether to use negative powers of two in addition to the binary one (default False)
    - ``verbose`` -- boolean for more verbose outputs (default True)

    OUTPUT: a list of tuples `(x, w, repr, a)` where `x` is the generator of the BLS fields,
    `w` is the Hamming weight of `x`, repr is x representation in either regular or extended form
    and `a` is the adicity of the scalar field of the BLS.

    Non-conservative parameters will yield base field primes at around 384 bits,
    and conservative ones will yield base field primes at 448 bits.

    """

    assert(weight_start > 0)
    adicity = adicity // 2
    L = []
    if conservative:
        limit = 74  # Base field size < 448bits
    else:
        limit = 63  # Base field size < 384bits
    assert(adicity <= limit - weight_end)

    for weight in range(weight_start-1 + wid, weight_end, processes):
        count = 0
        List_wx = list(combinations(range(adicity, limit), weight))
        if extended:
            Set_sign_wx_1 = set(
                combinations_with_replacement(range(0, 2), weight))
            Set_sign_wx_2 = set(combinations_with_replacement(
                reversed(range(0, 2)), weight))
            # Actually still a Set, maybe there's a better way to remove duplicates
            List_sign_wx = Set_sign_wx_2 - Set_sign_wx_1
            List_sign_wx = list(Set_sign_wx_1) + list(List_sign_wx)
        else:
            List_sign_wx = [[0 for i in range(weight)]]
        output = "Weight %s\n" % (weight+1)
        total = len(List_wx) * len(List_sign_wx)
        output += "\tTotal cases: %s\n" % total
        for item in List_wx:
            for sign in List_sign_wx:
                x = 1 << (limit)  # to start already at desired size for r
                for i in range(len(sign)):
                    if sign[i] == 0:
                        x += 1 << item[i]
                    else:
                        x -= 1 << item[i]
                r = bls12_scalar(x)
                if (conservative and r.nbits() > 297) or (not conservative and r.nbits() > 255):
                    continue
                if gcd_small_primes(r) == None:
                    continue
                if r.is_pseudoprime():
                    w = len(item) + 1
                    ad = twoadicity(r)
                    p = bls12_base(x, r)
                    bin_x = "2^%s" % limit
                    for i in reversed(range(0, len(sign))):
                        bin_x += " %s 2^%s" % ("+" if sign[i]
                                               == 0 else "-", item[i])
                    if p.is_pseudoprime():
                        if extended:
                            L.append([x, len(item) + 1, bin_x, ad])
                        else:
                            L.append([x, w, bin_x, ad])
                        count += 1
                    p = bls12_base(-x, r)
                    if p.is_pseudoprime():
                        bin_x = "-(" + bin_x + ")"
                        if extended:
                            L.append([-x, len(item) + 1, bin_x, ad])
                        else:
                            L.append([-x, w, bin_x, ad])
                        count += 1
        output += "\tValid cases: %s (%s" % (count,
                                             round(1.0 * count / total, 4))
        output += " %)\n"
        if verbose:
            print(output)
    return L


def find_BLS24_curve(adicity, weight_start=2, weight_end=8, _conservative=False, wid=0, processes=1, extended=False, verbose=True):
    """Return parameters of valid BLS24 curves

    INPUT:

    - ``adicity`` -- minimum 2-adicity of the scalar field
    - ``weight_start`` -- minimum Hamming weight for generator `x` (default 2)
    - ``weight_end`` -- maximum Hamming weight for generator `x` (default 8)
    - ``conservative`` -- boolean indicating whether selecting conservative parameters or not (default False)
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)
    - ``extended`` -- boolean indicating whether to use negative powers of two in addition to the binary one (default False)
    - ``verbose`` -- boolean for more verbose outputs (default True)

    OUTPUT: a list of tuples `(x, w, repr, a)` where `x` is the generator of the BLS fields,
    `w` is the Hamming weight of `x`, repr is x representation in either regular or extended form
    and `a` is the adicity of the scalar field of the BLS.

    """

    assert(weight_start > 0)
    adicity = adicity // 4
    L = []
    limit = 31  # Scalar field size < 256bits
    assert(adicity <= limit - weight_end)

    for weight in range(weight_start-1 + wid, weight_end, processes):
        count = 0
        List_wx = list(combinations(range(adicity, limit), weight))
        if extended:
            Set_sign_wx_1 = set(
                combinations_with_replacement(range(0, 2), weight))
            Set_sign_wx_2 = set(combinations_with_replacement(
                reversed(range(0, 2)), weight))
            # Actually still a Set, maybe there's a better way to remove duplicates
            List_sign_wx = Set_sign_wx_2 - Set_sign_wx_1
            List_sign_wx = list(Set_sign_wx_1) + list(List_sign_wx)
        else:
            List_sign_wx = [[0 for i in range(weight)]]
        output = "Weight %s\n" % (weight+1)
        total = len(List_wx) * len(List_sign_wx)
        output += "\tTotal cases: %s\n" % total
        for item in List_wx:
            for sign in List_sign_wx:
                x = 1 << (limit)  # to start already at desired size for r
                for i in range(len(sign)):
                    if sign[i] == 0:
                        x += 1 << item[i]
                    else:
                        x -= 1 << item[i]
                r = bls24_scalar(x)
                if r.nbits() > 255:
                    continue
                if gcd_small_primes(r) == None:
                    continue
                if r.is_pseudoprime():
                    w = len(item) + 1
                    ad = twoadicity(r)
                    p = bls24_base(x, r)
                    bin_x = "2^%s" % limit
                    for i in reversed(range(0, len(sign))):
                        bin_x += " %s 2^%s" % ("+" if sign[i]
                                               == 0 else "-", item[i])
                    if p.is_pseudoprime():
                        if extended:
                            L.append([x, len(item) + 1, bin_x, ad])
                        else:
                            L.append([x, w, bin_x, ad])
                        count += 1
                    p = bls24_base(-x, r)
                    if p.is_pseudoprime():
                        bin_x = "-(" + bin_x + ")"
                        if extended:
                            L.append([-x, len(item) + 1, bin_x, ad])
                        else:
                            L.append([-x, w, bin_x, ad])
                        count += 1
        output += "\tValid cases: %s (%s" % (count,
                                             round(1.0 * count / total, 4))
        output += " %)\n"
        if verbose:
            print(output)
    return L


def find_BLS48_curve(adicity, weight_start=2, weight_end=8, _conservative=False, wid=0, processes=1, extended=False, verbose=True):
    """Return parameters of valid BLS48 curves

    INPUT:

    - ``adicity`` -- minimum 2-adicity of the scalar field
    - ``weight_start`` -- minimum Hamming weight for generator `x` (default 2)
    - ``weight_end`` -- maximum Hamming weight for generator `x` (default 8)
    - ``conservative`` -- boolean indicating whether selecting conservative parameters or not (default False)
    - ``wid`` -- current job id (default 0)
    - ``processes`` -- number of concurrent jobs (default 1)
    - ``extended`` -- boolean indicating whether to use negative powers of two in addition to the binary one (default False)
    - ``verbose`` -- boolean for more verbose outputs (default True)

    OUTPUT: a list of tuples `(x, w, repr, a)` where `x` is the generator of the BLS fields,
    `w` is the Hamming weight of `x`, repr is x representation in either regular or extended form
    and `a` is the adicity of the scalar field of the BLS.

    """

    assert(weight_start > 0)
    adicity = adicity // 8
    L = []
    limit = 31  # Scalar field size < 512bits
    assert(adicity <= limit - weight_end)

    for weight in range(weight_start-1 + wid, weight_end, processes):
        count = 0
        List_wx = list(combinations(range(adicity, limit), weight))
        if extended:
            Set_sign_wx_1 = set(
                combinations_with_replacement(range(0, 2), weight))
            Set_sign_wx_2 = set(combinations_with_replacement(
                reversed(range(0, 2)), weight))
            # Actually still a Set, maybe there's a better way to remove duplicates
            List_sign_wx = Set_sign_wx_2 - Set_sign_wx_1
            List_sign_wx = list(Set_sign_wx_1) + list(List_sign_wx)
        else:
            List_sign_wx = [[0 for i in range(weight)]]
        output = "Weight %s\n" % (weight+1)
        total = len(List_wx) * len(List_sign_wx)
        output += "\tTotal cases: %s\n" % total
        for item in List_wx:
            for sign in List_sign_wx:
                x = 1 << (limit)  # to start already at desired size for r
                for i in range(len(sign)):
                    if sign[i] == 0:
                        x += 1 << item[i]
                    else:
                        x -= 1 << item[i]
                r = bls48_scalar(x)
                if r.nbits() > 511:
                    continue
                if gcd_small_primes(r) == None:
                    continue
                if r.is_pseudoprime():
                    w = len(item) + 1
                    ad = twoadicity(r)
                    p = bls48_base(x, r)
                    bin_x = "2^%s" % limit
                    for i in reversed(range(0, len(sign))):
                        bin_x += " %s 2^%s" % ("+" if sign[i]
                                               == 0 else "-", item[i])
                    if p.is_pseudoprime():
                        if extended:
                            L.append([x, len(item) + 1, bin_x, ad])
                        else:
                            L.append([x, w, bin_x, ad])
                        count += 1
                    p = bls48_base(-x, r)
                    if p.is_pseudoprime():
                        bin_x = "-(" + bin_x + ")"
                        if extended:
                            L.append([-x, len(item) + 1, bin_x, ad])
                        else:
                            L.append([-x, w, bin_x, ad])
                        count += 1
        output += "\tValid cases: %s (%s" % (count,
                                             round(1.0 * count / total, 4))
        output += " %)\n"
        if verbose:
            print(output)
    return L


########################################################################

def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    strategy = find_BLS24_curve if "--bls24" in args else (
        find_BLS48_curve if "--bls48" in args else find_BLS12_curve)
    conservative = True if "--conservative" in args else False
    extended = True if "--extended" in args else False
    sortadicity = True if "--sort_adicity" in args else False
    silent = True if "--silent" in args else False
    save = True if "--save" in args else False
    help = True if "--help" in args else False
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1 or help:
        print("""
Cmd: sage find_bls.sage [--sequential] [--conservative] [--sort_adicity] [--silent] [--save] [--bls24] [--bls48]
                          <min-2adicity> [<min-weight> [<max-weight>]]

Args:
    --sequential        Uses only one process
    --bls24             Searches BLS24 curves instead of BLS12. Cannot be combined with conservative
    --bls48             Searches BLS48 curves instead of BLS12. Cannot be combined with conservative
    --conservative      Uses conservative sizes for BLS security
    --extended          Uses pos/neg combinations of powers of 2 for the BLS generator
    --sort_adicity      Sorts results by decreasing 2-adicity rather than increasing Hamming weight
    --silent            Ignores stats at each search step
    --save              Saves list into csv file
    <min-2adicity>      Minimum two-adicity of the scalar field of BLS
    <min-weight>        Minimum Hamming weight to start lookup
    <max-weight>        Maximum Hamming weight to end lookup
""")
        return

    adicity = int(args[0])
    min_weight = int(args[1]) if len(args) > 1 else 2
    max_weight = int(args[2]) if len(args) > 2 else 6

    result_list = []
    async_result_list = []

    def collect_result(result):
        async_result_list.append(result)

    if processes == 1 or min_weight == max_weight:
        result_list = strategy(
            adicity, min_weight, max_weight, conservative, 0, 1, extended, not silent)
    else:
        if not silent:
            print("Using %d processes." % (processes,))
        pool = Pool(processes=processes)

        for wid in range(processes):
            pool.apply_async(worker, (strategy, adicity, min_weight, max_weight, conservative,
                             wid, processes, extended, not silent), callback=collect_result)

        pool.close()
        pool.join()
        for i in async_result_list:
            result_list += i

    if len(result_list) > 1:
        if extended:
            if sortadicity:
                # sorting by 2-adicity first
                result_list.sort(key=lambda x: (x[3], -x[1]))
            else:
                # sorting by weight first
                result_list.sort(key=lambda x: (-x[1], x[3]))
        else:
            if sortadicity:
                # sorting by 2-adicity first
                result_list.sort(key=lambda x: (x[3], -x[1]))
            else:
                # sorting by weight first
                result_list.sort(key=lambda x: (-x[1], x[3]))
        result_list.reverse()

    if save:
        if strategy == find_BLS12_curve:
            filename = "bls12_generators/generator_list_%d_%d_%d" % (
                adicity, min_weight, max_weight)
        elif strategy == find_BLS24_curve:
            filename = "bls24_generators/generator_list_%d_%d_%d" % (
                adicity, min_weight, max_weight)
        else:
            filename = "bls48_generators/generator_list_%d_%d_%d" % (
                adicity, min_weight, max_weight)
        if conservative:
            filename += "_conservative"
        if extended:
            filename += "_extended"
        filename += ".csv"
        f = open(filename, 'w')
        writer = csv.writer(f)
        for res in result_list:
            writer.writerow(res)
        f.close()
    else:
        for res in result_list:
            print(res)
    if not silent:
        print("Number of valid candidates: {}".format(len(result_list)))


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
