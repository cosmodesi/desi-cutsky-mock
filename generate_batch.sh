#!/bin/bash

# phases=(358 1150 1153 1160 1173 1181 1196 1761 1762 1763 1764 1765 1766 1767 1768 1829)

N=4999
for (( index=3001; index<N; index=index+1 ))
do
	# phase=1
	phase=$index
	# phase=${phases[$index]}
	#filename=/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B2000G512Z1.1N24000470_b0.345d1.45r40c0.05_seed
	#if [[ ! -f "${filename}${phase}.fits" ]]; then
	echo ${filename}${phase}.fits
	sed -e "s/PHPH/${phase}/g" -e "s/GALGAL/LRG/g"  batchscript_gen_hsw.sh > batchscript_gen_hsw_${phase}.sh
	sbatch batchscript_gen_hsw_${phase}.sh
	rm batchscript_gen_hsw_${phase}.sh
	

done
