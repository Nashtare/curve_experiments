"""
This module aims at providing finite field parameters for Arkworks algebra/ff library.
"""

from dataclasses import field
import sys

from util import twoadicity, pretty_element_repr


def params(p):
    """Return specific field constants for reuse in Arkworks ecosystem

    INPUT:

    - ``p`` -- a prime number to define Fp

    OUTPUT: a string containing all field constants necessary for defining
    a PrimeField of cardinality `p` in the Arkworks ecosystem.

    Unformatted means they are not represented in Montgomery form, see
    https://github.com/arkworks-rs/algebra/tree/master/ff

    """

    F = GF(p)
    field_size_hex = (p.nbits() + 4 - (p.nbits() % 4)) / 4
    output = f"\n\n\nParameters for GF({p}) arithmetic construction:\n"
    output += f"\n\nTWO_ADICITY: {twoadicity(p)}"
    r = F(2)**(field_size_hex * 4)
    r_squared = r ** 2
    s = 2 ** twoadicity(p)
    t = (p-1) / s
    g = F.multiplicative_generator()
    root_unity = g ^ t
    output += f"\n\nTWO_ADIC ROOT UNITY: {pretty_element_repr(root_unity * r)}"
    output += f"\n\nMODULUS: {pretty_element_repr(p, field_size_hex)}"
    output += f"\n\nR: {pretty_element_repr(r, field_size_hex)}"
    output += f"\n\nR^2: {pretty_element_repr(r_squared, field_size_hex)}"
    output += f"\n\nMODULUS_MINUS_ONE_DIV_TWO (unformatted): {pretty_element_repr((p-1)/2, field_size_hex)}"
    output += f"\n\nT (unformatted): {repr(t)}"
    output += f"\n\nT_MINUS_ONE_DIV_TWO (unformatted): {pretty_element_repr((t-1)/2, field_size_hex)}"
    output += f"\n\nGENERATOR: {pretty_element_repr(g * r, field_size_hex)}"
    output += f"\n\nMODULUS_BITS: {p.nbits()}"
    output += f"\n\nCAPACITY: {p.nbits() - 1}"
    output += f"\n\nREPR SHAVE BITS: {2**(p.nbits() - 1).nbits() - p.nbits()}"
    output += f"\n\nINV: {(-(p ^ (-1)) % 2 ^ 64)}"
    print(output)


########################################################################

def main():
    """Main function"""
    args = sys.argv[1:]
    help = "--help" in args

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
