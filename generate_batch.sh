#!/bin/bash

# phases=(0 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
N=500
for (( index=100; index<N; index++ ))
do
	phase=$index #${phases[$index]}
	sed -e "s/PHPH/${phase}/g"  batchscript_gen_hsw.sh > batchscript_gen_hsw_${phase}.sh
    sbatch batchscript_gen_hsw_${phase}.sh
	rm batchscript_gen_hsw_${phase}.sh
done
