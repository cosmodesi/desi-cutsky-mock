#!/bin/bash

# SNAP=(11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
SNAP=(7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
N=18
# SNAP=(4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
# N=18

path=/global/cfs/cdirs/desi/cosmosim/mocks_SV3/v1/ABACUS/ABACUS_Gcat_FITS/
for (( i=0; i<N; i++ ))
do
    snapshot=${SNAP[$i]}
    echo $snapshot
    # du -ch ${path}/LRG_snap${snapshot}*fits | grep total
    du -ch ${path}/ELG_snap${snapshot}*fits | grep total
    # du -ch ${path}/QSO_snap${snapshot}*fits | grep total
done