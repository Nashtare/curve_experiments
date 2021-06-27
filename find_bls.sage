import csv
import sys
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from math import ceil
from itertools import combinations

if sys.version_info[0] == 2: range = xrange

SMALL_PRIMES = (3, 5, 7, 11, 13, 17)

def gcd_small_primes(p):
    return (r for r in SMALL_PRIMES if gcd(p-1, r) == 1)

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

def twoadicity(x, base=0, limit=64):
    return max(i for i in range(base, limit) if ((x-1) % (1<<i) == 0))

def h_weight(x):
    return x.binary().count('1')

# storing parameters of valid BLS12 curves of 2-adicity alpha and Hamming weight <= weight (with log2(p) < 384 or log2(p) < 446)
def find_BLS12_curve(adicity, weight_start = 2, weight_end = 8, conservative = False, wid = 0, processes = 1, verbose=True):
    adicity = adicity // 2
    assert(adicity <= 60)
    L = []
    if conservative:
        limit = 74 # Base field size < 446bits
    else:
        limit = 64 # Base field size < 384bits

    for weight in range(weight_start-1 + wid, weight_end, processes):
        count = 0
        List = list(combinations((i for i in range(1, limit-adicity)), weight))
        output = "Weight %s\n" % (weight+1)
        output += "\tTotal cases: %s\n" % len(List)
        for item in List:
            x = ["0" for k in range(limit)] # to start already at 255 bits for r
            x[0] = "1"
            for i in item:
                x[i] = "1"
            x = Integer(int("0b" + "".join(x),2))
            r = bls_scalar(x)
            if (conservative and r.nbits() > 297) or (not conservative and r.nbits() > 255):
                continue
            if gcd_small_primes(r) == None:
                continue
            if r.is_pseudoprime() :
                w = h_weight(x)
                ad = twoadicity(r)
                p = bls_base(x, r)
                if p.is_pseudoprime():
                    L.append([x, w, ad])
                    count += 1
                p = bls_base(-x, r)
                if p.is_pseudoprime():
                    L.append([-x, w, ad])
                    count += 1
        output += "\tValid cases: %s (%s" % (count, round(1.0 * count / len(List), 4))
        output += " %)\n"
        if verbose:
            print(output)
    return L

########################################################################

def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    conservative = True if "--conservative" in args else False
    sortadicity = True if "--sort_adicity" in args else False
    silent = True if "--silent" in args else False
    save = True if "--save" in args else False
    help = True if "--help" in args else False
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1 or help:
        print("""
Cmd: sage find_bls.sage [--sequential] [--conservative] [--sort_adicity] [--silent] [--save]
                          <min-2adicity> [<min-weight> [<max-weight>]]

Args:
    --sequential        Uses only one process
    --conservative      Uses conservative sizes for BLS security
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
    max_weight = int(args[2]) if len(args) > 2 else 8

    result_list = []
    async_result_list = []

    def collect_result(result):
        async_result_list.append(result)

    if processes == 1 or min_weight == max_weight:
        result_list = find_BLS12_curve(adicity, min_weight, max_weight, conservative)
    else:
        print("Using %d processes." % (processes,))
        pool = Pool(processes=processes)

        for wid in range(processes):
            pool.apply_async(worker, (adicity, min_weight, max_weight, conservative, wid, processes, not silent), callback=collect_result)
    
        pool.close()
        pool.join()
        for i in async_result_list:
            result_list += i

    if len(result_list) > 1:
        if sortadicity:
            result_list.sort(key = lambda x: (-x[2], x[1])) # sorting by 2-adicity first
        else:
            result_list.sort(key = lambda x: (x[1], -x[2])) # sorting by weight first
        result_list.reverse()

    if save:
        filename = "generator_list_%d_%d_%d.csv" % (adicity, min_weight, max_weight)
        f = open(filename, 'w')
        writer = csv.writer(f)
        for res in result_list:
            writer.writerow(res)
        f.close()
    else:
        for res in result_list:
            print(res)
    print("Number of valid candidates: {}".format(len(result_list)))

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
    return find_BLS12_curve(*args)

main()
