#!/bin/bash

#phases=(0 1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 )
N=1001	
for (( index=1; index<N; index++ ))
do
	phase=$index #${phases[$index]}
	filename=/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B2000G512Z1.1N24000470_b0.345d1.45r40c0.05_seed
	if [[ ! -f "${filename}${phase}.fits" ]]; then
		echo ${filename}${phase}.fits
		sed -e "s/PHPH/${phase}/g"  batchscript_gen_hsw.sh > batchscript_gen_hsw_${phase}.sh
		sbatch batchscript_gen_hsw_${phase}.sh
		rm batchscript_gen_hsw_${phase}.sh
	

fi


done
