#!/bin/bash

declare -a type=("" "--bls24" "--bls48")
declare -a adicity_and_weights=("32 1 5" "64 1 10")
declare -a extended=("" "--extended")
conservative="--conservative"

for i in "${adicity_and_weights[@]}"
do
    for j in "${extended[@]}"
    do
        for k in "${type[@]}"
        do
            str="sage find_bls.sage --silent --save $i $j $k"
            $str
            echo "Finished execution of $str"
        done

        str="sage find_bls.sage --silent --save $i $j --conservative" # only for bls12
        $str
    done
done