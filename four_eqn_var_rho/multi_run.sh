#!/bin/bash
make
for i in {1..20}
do
    for j in {1..20}
    do
        # var1 = $(7+$i*(1.5/20.0))
        var2 = $(("scale=2;5+$i*25/20" | bc))
        # ./four_eqn_model $var1 $var2
    done
done