# curve_experiments

Experimental repo to play with elliptic curves

## Usage

To generate BLS parameters with given 2-adicity for the scalar field:
```bash
sage find_bls.sage --sort_adicity 35
```

To (at least try) find a cycle given a list of BLS generators such that
the scalar field of the BLS is on the cycle:
```bash
sage find_cycle.sage generator_list_32_2_8.csv
```

## License
[MIT](https://choosealicense.com/licenses/mit/)