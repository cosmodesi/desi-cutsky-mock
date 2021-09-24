#!/usr/bin/env python

# Generate lightcone for DESI mocks
# from HOD catalogues generated by Shadab

#import os
import argparse
#import time
#import tracemalloc
#import numpy as np
#from astropy.io import fits
#import h5py

# from DESILCclass import DESILC
# from UNITLCclass import UNITLC, combine_shells, Paths

from generate_lc import LC, Paths_LC
from foot_nz import FOOT_NZ

def lrg_LC_UNIT(config_file, args):
	input_name = "LRG_snap{}.gcat.sub{}.fits"

	######### Lightcones
	shells_path = "/LC_shells_sv3_nz_radecz_xyz_LRG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	# lc_instance = LC(config_file, args)
	# lc_instance.generate_shells(path_instance, snapshot=None, cutsky=False, nproc=10, Nsubboxes=27)

	foot_nz_instance = FOOT_NZ(config_file, args, galtype="lrg")
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)

	######### CutSky
	shells_path = "/cutsky_shells_sv3_nz_radecz_xyz_LRG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	# lc_instance.generate_shells(path_instance, snapshot=98, cutsky=True, nproc=10)

	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)

def elg_LC_UNIT(config_file, args):
	input_name = "ELG_snap{}.gcat.sub{}.fits"
	
	
	######### Lightcones

	shells_path = "/LC_shells_sv3_nz_radecz_xyz_ELG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	# lc_instance = LC(config_file, args)
	# lc_instance.generate_shells(path_instance, snapshot=None, cutsky=False, nproc=8)

	foot_nz_instance = FOOT_NZ(config_file, args, galtype="elg")
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)
	
	
	# ######### CutSky
	shells_path = "/cutsky_shells_sv3_nz_radecz_xyz_ELG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)
	
	# lc_instance.generate_shells(path_instance, snapshot=93, cutsky=True, nproc=8, Nsubboxes=27)

def lrg_RAN_UNIT(config_file, args):
	input_name = "LRG_snap{}.gcat.sub{}.fits"

	######### Lightcones
	
	lc_instance = LC(config_file, args)

	# foot_nz_instance = FOOT_NZ(config_file, args, galtype="lrg")
	
	######### CutSky
	for i in range(10):
		shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_LRG_h5py/"
		path_instance = Paths_LC(config_file, args, input_name, shells_path)
		
		# lc_instance.generate_shells(path_instance, snapshot=98, cutsky=True, nproc=10, random=i+1)
		foot_nz_instance = FOOT_NZ(config_file, args, galtype="lrg")
		foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
		foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)

def elg_RAN_UNIT(config_file, args):
	input_name = "ELG_snap{}.gcat.sub{}.fits"
	
	
	######### Lightcones

	lc_instance = LC(config_file, args)

	# foot_nz_instance = FOOT_NZ(config_file, args, galtype="elg")
	
	# ######### CutSky
	# for i in range(10):
	i = 9
	shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_ELG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	# lc_instance.generate_shells(path_instance, snapshot=93, cutsky=True, nproc=8, Nsubboxes=27, random=i+1)
	foot_nz_instance = FOOT_NZ(config_file, args, galtype="elg")
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)
	

def lrg_LC_ABACUS(args):
	config_file = "./config_ABACUS_LRG.ini"
	input_name = "LRG_snap{}.gcat.sub{}.fits"

	######### Lightcones
	shells_path = "/LC_shells_sv3_nz_radecz_xyz_LRG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	lc_instance = LC(config_file, args)
	lc_instance.generate_shells(path_instance, snapshot=None, cutsky=False, nproc=20, Nsubboxes=64)

	foot_nz_instance = FOOT_NZ(config_file, args, galtype="lrg")
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)
	
	######### CutSky
	shells_path = "/cutsky_shells_sv3_nz_radecz_xyz_LRG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	lc_instance.generate_shells(path_instance, snapshot=18, cutsky=True, nproc=20, Nsubboxes=64)

	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)

def elg_LC_ABACUS(args):
	config_file = "./config_ABACUS_ELG.ini"
	input_name = "ELG_snap{}.gcat.sub{}.fits"

	######### Lightcones
	# shells_path = "/LC_shells_sv3_nz_radecz_xyz_ELG_h5py/"
	# path_instance = Paths_LC(config_file, args, input_name, shells_path)

	lc_instance = LC(config_file, args)
	# lc_instance.generate_shells(path_instance, snapshot=None, cutsky=False, nproc=15, Nsubboxes=64)

	foot_nz_instance = FOOT_NZ(config_file, args, galtype="elg")
	# foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)

	# ######### CutSky
	# shells_path = "/cutsky_shells_sv3_nz_radecz_xyz_ELG_h5py/"
	# path_instance = Paths_LC(config_file, args, input_name, shells_path)

	# lc_instance.generate_shells(path_instance, snapshot=15, cutsky=True, nproc=15, Nsubboxes=64)
	# foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)

	### Random Cutsky
	# for i in range(10):
	
	i = 9
	shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_ELG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	lc_instance.generate_shells(path_instance, snapshot=15, cutsky=True, nproc=15, Nsubboxes=64, random=i+1)
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=15, fullfootprint=0, todo=0)

def lrg_RAN_ABACUS(args):
	config_file = "./config_ABACUS_LRG.ini"
	input_name = "LRG_snap{}.gcat.sub{}.fits"

	######### Lightcones
	
	lc_instance = LC(config_file, args)

	foot_nz_instance = FOOT_NZ(config_file, args, galtype="lrg")
	
	######### CutSky
	
	# for i in range(10):
	
	i = 8
	shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_LRG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	lc_instance.generate_shells(path_instance, snapshot=18, cutsky=True, nproc=20, random=i+1, Nsubboxes=64)
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)
	
	i = 9
	shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_LRG_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	lc_instance.generate_shells(path_instance, snapshot=18, cutsky=True, nproc=20, random=i+1, Nsubboxes=64)
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)
	
def qso_LC_ABACUS(args):
	config_file = "./config_ABACUS_QSO.ini"
	input_name = "QSO_snap{}.gcat.sub{}.fits"

	######### Lightcones
	shells_path = "/LC_shells_sv3_nz_radecz_xyz_QSO_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	lc_instance = LC(config_file, args)
	lc_instance.generate_shells(path_instance, snapshot=None, cutsky=False, nproc=20, Nsubboxes=64)

	# foot_nz_instance = FOOT_NZ(config_file, args, galtype="qso")
	# foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)
	
	######### CutSky
	shells_path = "/cutsky_shells_sv3_nz_radecz_xyz_QSO_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)

	lc_instance.generate_shells(path_instance, snapshot=11, cutsky=True, nproc=20, Nsubboxes=64)

	# foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)

	######### CutSky
	# for i in range(10):
	i=8
	shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_QSO_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	lc_instance.generate_shells(path_instance, snapshot=11, cutsky=True, nproc=20, random=i+1, Nsubboxes=64)
	
	i=9
	shells_path = f"/ran_{i+1}_cutsky_shells_sv3_nz_radecz_xyz_QSO_h5py/"
	path_instance = Paths_LC(config_file, args, input_name, shells_path)
	
	lc_instance.generate_shells(path_instance, snapshot=11, cutsky=True, nproc=20, random=i+1, Nsubboxes=64)
	
	# foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=2, todo=2)
	# foot_nz_instance.shell(path_instance, nproc=20, fullfootprint=0, todo=0)

def bgs_LC_ABACUS():

	# bgs_file = "/global/cscratch1/sd/avariu/desi/ABACUS_2GPC/BGS/BGS_ANY_LC.fits"

	from foot_nz import apply_footprint
	import fitsio
	import numpy as np
	import glob
	
	files = glob.glob("/global/cscratch1/sd/avariu/desi/ABACUS_2GPC/BGS/*rand*")
	print(len(files))
	for bgs_file in files:
		print(bgs_file)
		fits = fitsio.FITS(bgs_file, "rw")

		print("read ra")
		ra = fits[1].read_column("RA")
		
		print("read dec")
		dec = fits[1].read_column("DEC")

		print("apply sv3 foot")
		footbits_sv3 = apply_footprint(bgs_file, ra, dec, 2)
		print(np.unique(footbits_sv3))

		print("apply Y5 foot")
		footbits_y5 = apply_footprint(bgs_file, ra, dec, 0)
		print(np.unique(footbits_y5))

		print("first sv3 or y5")
		out_arr = np.bitwise_or(footbits_sv3, footbits_y5)
		print(np.unique(out_arr))
		
		print("first foot or downsample")
		out_arr_final = np.bitwise_or(out_arr, np.ones(len(out_arr), dtype=np.int32))
		print(np.unique(out_arr_final))
				
		print("insert column")
		fits[1].insert_column('STATUSGOOD', out_arr_final)
		fits.close()


	


	#['RA', 'DEC', 'Z', 'ZREAL']
	

def main():
	parser = argparse.ArgumentParser()
	# parser.add_argument("config", help="ini file holding configuration", type=str)
	parser.add_argument("--dir_out", type=str, help="output directory (overrides config file)")
	parser.add_argument("--dir_gcat", type=str, help="input directory (same)")
	parser.add_argument("--input_name", type=str, help="name of input catalogs (same)")
	parser.add_argument("--lightcone_path", type=str, help="path of output catalogs (same)")
	args = parser.parse_args()
	# config_file = str(args.config) # config file

	# config_file = "./config_UNIT_LRG.ini"
	# lrg_LC_UNIT(config_file, args)
	# lrg_RAN_UNIT(config_file, args)

	# config_file = "./config_UNIT_ELG.ini"
	# elg_LC_UNIT(config_file, args)
	# elg_RAN_UNIT(config_file, args)

	# qso_LC_ABACUS(args)

	# elg_LC_ABACUS(args)
	
	# lrg_RAN_ABACUS(args)
	# lrg_LC_ABACUS(args)
	bgs_LC_ABACUS()


if __name__ == '__main__':
	main()
