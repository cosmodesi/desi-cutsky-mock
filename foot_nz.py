#!/usr/bin/env python

# Generate lightcone for DESI mocks
# from HOD catalogues generated by Shadab

import sys
import os
import glob
import time
import configparser
import multiprocessing as mp

import numpy as np
from itertools import product
from astropy.table import Table
import desimodel.footprint as foot
import desimodel.io
import h5py


def bits(ask="try"):
	if ask == "LC":              return 0 #(0 0 0 0)
	if ask == "downsample":      return 1 #(0 0 0 1)
	if ask == "entireDESIfoot":  return 2 #(0 0 1 0)
	if ask == "SV3foot":         return 4 #(0 1 0 0)
	if ask == "downsample_main": return 8 #(1 0 0 0)
	sys.exit()

def mask(main=0, nz=0, Y5=0, sv3=0):
    return nz * (2**0) + Y5 * (2**1) + sv3 * (2**2) + main * (2**3)

def get_nz(zz, galtype="test"):
	''' The function where the n(z) is read 
	and the NZ column is computed for the given
	redshifts.
	'''

	if galtype == "LRG":
		filename = "sm_LRG_mycosmo_ev2.1.dat"
		failurerate = 0.

	elif galtype == "LRG_main":
		filename = "LRGzdone_mycosmo_610_nz_DA02.dat"
		failurerate = 0.
	
	elif galtype == "ELG":
		filename = "sm_ELG_mycosmo_ev2.1.dat"
		failurerate = 0.25

	elif galtype == "QSO":
		filename = "sm_QSO_mycosmo_ev2.1.dat"
		failurerate = 0.37
	else:
		print("Unknown galaxy type.")
		sys.exit()
	
	nz_file = f"./nz_files/{filename}"
	
	z, nz = np.loadtxt(nz_file, usecols=(0, 1), unpack=True)
	z_n = z
	
	nz_n = (nz / ( 1 - failurerate )) * 1.0
	

	np.savetxt(f"./nz_files/{filename}_redcor.txt", np.array([z_n, nz_n]).T)

	return np.interp(zz, z_n, nz_n, left=0, right=0)


def downsample_aux(z_cosmo, galtype, ran, n_mean, ask="downsample"):
	""" downsample galaxies following n(z) model specified in galtype"""

	nz = get_nz(z_cosmo, galtype=galtype)

	# downsample
	nz_selected = ran < nz / n_mean
	idx         = np.where(nz_selected)
	print("DOWNSAMPLE: Selected {} out of {} galaxies.".format(len(idx[0]), len(z_cosmo)), flush=True)

	bitval = bits(ask=ask)
	
	newbits = np.zeros(len(z_cosmo), dtype=np.int32)
	newbits[idx] = bitval

	return newbits


def downsample(boxL, galtype, ngalbox, z_cosmo):
	""" downsample galaxies following n(z) model specified in galtype"""

	n_mean = ngalbox/(boxL**3)

	ran         = np.random.rand(len(z_cosmo))

	newbits = downsample_aux(z_cosmo, galtype, ran, n_mean, ask="downsample")
	
	if galtype == "LRG":
		newbits_main = downsample_aux(z_cosmo, "LRG_main", ran, n_mean, ask="downsample_main")
		outbits = np.bitwise_or(newbits, newbits_main)
	
		return outbits, ran

	return newbits, ran


def apply_footprint(ra, dec, footprint_mask):
	""" apply desi footprint """
	
	bitval = 0
	# footprint_mask possibilities
	# 0 - Y5 DESI, 2 - SV3 DESI

	if footprint_mask == 0:
		tiles = desimodel.io.load_tiles()
		point = foot.is_point_in_desi(tiles, ra, dec)
		bitval = bits(ask="entireDESIfoot")
	elif footprint_mask == 2:
		tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
		point = foot.is_point_in_desi(tiles, ra, dec)
		bitval = bits(ask="SV3foot")
	else:
		print("ERROR: Wrong footprint.")
		sys.exit()
	
	idx   = np.where(point)

	print("FOOTPRINT: Selected {} out of {} galaxies.".format(len(idx[0]), len(ra)), flush=True)
	
	newbits = np.zeros(len(ra), dtype=np.int32)
	newbits[idx] = bitval

	return newbits


def determine_return_shellnum(file_, begin_out_shell):
	filename = os.path.basename(file_)
	index_i = filename.find("_shell_")
	index_f = filename.find(".h5py")
	shellnum = filename[index_i + 7: index_f]

	files = glob.glob(begin_out_shell.format("*", "*") + "_shell_" + shellnum + ".h5py")
	print(files)
	if len(files) != 1:
		print("ERROR: return_shellnum", len(files))
		sys.exit()
	elif not os.path.basename(files[0]):
		print("ERROR: return_shellnum", len(files), files[0])
		sys.exit()
	else:
		return int(shellnum) 


def generate_shell(args):
	file_, boxL, galtype, tracer_id, footprint_mask, todo, begin_out_shell, cat_seed = args
	
	shellnum = determine_return_shellnum(file_, begin_out_shell)
	
	print(file_, shellnum)

	unique_seed = tracer_id * 500500 + 250 * cat_seed + shellnum
	print("INFO: UNIQUE SEED:", unique_seed, flush=True)
	np.random.seed(unique_seed)

	
	f = h5py.File(file_, 'r+')
	data = f['galaxy']
	ra = data['RA'][()]
	dec = data['DEC'][()]
	z_cosmo = data['Z_COSMO'][()]
			
	
	foot_bit0 = apply_footprint(ra, dec, 0)
	foot_bit2 = apply_footprint(ra, dec, 2)
	foot_bit = np.bitwise_or(foot_bit0, foot_bit2)
	
	down_bit, ran_arr = downsample(boxL, galtype, f.attrs["NGAL"], z_cosmo)		
	
	out_arr = np.bitwise_or(foot_bit, down_bit)
	out_arr = out_arr.astype(np.int32)
	

	if "STATUS" in data.keys():
		print("STATUS EXISTS")
	
		# status = data["STATUS"][:]
		# new_arr = np.bitwise_or(status, out_arr)
		# new_arr = new_arr.astype(np.int32)
		# print(np.unique(new_arr))
		# data["STATUS"][:] = new_arr[:]

		# data["STATUS"][:] = out_arr[:]
	else:
		f.create_dataset('galaxy/STATUS', data=out_arr,  dtype=np.int32)

	if "RAN_NUM_0_1" in data.keys():
		print("RAN_NUM_0_1 EXISTS")
		# data['RAN_NUM_0_1'][:] = ran_arr[:]
	else:
		f.create_dataset('galaxy/RAN_NUM_0_1', data=ran_arr, dtype=np.float32)		

	f.close()


class FOOT_NZ():
	def __init__(self, config_file, args, galtype=None):
		config     = configparser.ConfigParser()
		config.read(config_file)

		self.boxL       =  config.getint('sim', 'boxL')
		self.zmin       =  config.getfloat('sim', 'zmin')
		self.zmax       =  config.getfloat('sim', 'zmax')
		
		self.galtype = galtype
		
		self.tracer_id = 0

		if galtype == "LRG" or galtype == "LRG_main":
			self.tracer_id = 0
		elif galtype == "ELG":
			self.tracer_id = 1
		elif galtype == "QSO":
			self.tracer_id = 2
		
		print(f"INFO: {self.galtype} with {self.tracer_id} ID")

	
	
	def shell(self, path_instance, cat_seed, nproc=5, footprint_mask=0, todo=1):
		
		infiles = glob.glob(path_instance.shells_out_path + "/*.h5py")

		args = product(infiles, [self.boxL], [self.galtype], [self.tracer_id], [footprint_mask], [todo], [path_instance.begin_out_shell], [cat_seed])

		if nproc > len(infiles):
			nproc = len(infiles)

		pool = mp.Pool(processes=nproc)
	
		pool.map_async(generate_shell, args)
		
		pool.close()
		pool.join()

	def shell_series(self, path_instance, cat_seed, footprint_mask=0, todo=1):
		
		infiles = glob.glob(path_instance.shells_out_path + "/*.h5py")

		for file_ in infiles:
			args = [file_, self.boxL, self.galtype, self.tracer_id, footprint_mask, todo, path_instance.begin_out_shell, cat_seed]
			generate_shell(args)

