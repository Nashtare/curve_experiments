# curve_experiments

Experimental repo to play with elliptic curves

## Usage

To generate BLS parameters with given 2-adicity for the scalar field:
```bash
sage find_bls.sage 35 1 7 --save --extended --conservative
```
Results are saved by default in the `bls_generators/` folder if `--save` is specified.

---

To find a cycle given a list of BLS generators such that
the scalar field of the BLS is on the cycle:
```bash
sage find_cycle.sage bls_generators/generator_list_32_1_5_conservative.csv
```

## License
[MIT](https://choosealicense.com/licenses/mit/)