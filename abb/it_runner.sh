#!/bin/bash

for i in $(seq 1 9); do
    value=$(bc <<< "scale=1; $i/10")
    echo "Tau: $value"
    Rscript interactive_testing.R $value
done
