# Outputting specific constants for reuse in arkworks ecosystem
# Unformatted means they are not represented in Montgomery form
# https://github.com/arkworks-rs/algebra/tree/master/ff

from util import twoadicity, pretty_element_repr

def params(p):
    F = GF(p)
    output = "\n\n\nParameters for GF(%d) arithmetic construction:\n" % p
    output += "\n\nTWO_ADICITY: %d" % twoadicity(p)
    r = F(2**256)
    r2 = F(2**512)
    s = 2^twoadicity(p)
    t = (p-1) / s
    g = F.multiplicative_generator()
    root_unity = g^t
    output += "\n\nTWO_ADIC ROOT UNITY: %s" % pretty_element_repr(root_unity * r)
    output += "\n\nMODULUS: %s" % pretty_element_repr(p)
    output += "\n\nR: %s" % pretty_element_repr(r)
    output += "\n\nR^2: %s" % pretty_element_repr(r2)
    output += "\n\nMODULUS_MINUS_ONE_DIV_TWO (unformatted): %s" % pretty_element_repr((p-1)/2)
    output += "\n\nT (unformatted): %s" % repr(t)
    output += "\n\nT_MINUS_ONE_DIV_TWO (unformatted): %s" % pretty_element_repr((t-1)/2)
    output += "\n\nGENERATOR: %s" % pretty_element_repr(g * r)
    output += "\n\nMODULUS_BITS: %d" % p.nbits()
    output += "\n\nCAPACITY: %s" % str(p.nbits() - 1)
    output += "\n\nREPR SHAVE BITS: %s" % str(256 - p.nbits())
    output += "\n\nINV: %d" % (-(p^(-1)) % 2^64)
    print(output)
    return

def main():
    args = sys.argv[1:]
    help = True if "--help" in args else False

    if len(args) < 1 or help:
        print("""
Cmd: sage field_params.sage <prime_number>

Args:
    <prime_number>      Modulus of the prime field we want to get constants from
""")
        return

    prime = Integer(str(args[0]))
    params(prime)

main()