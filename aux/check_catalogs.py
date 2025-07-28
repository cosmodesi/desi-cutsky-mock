import matplotlib
# matplotlib.use("Agg")
import glob
import sys
import os
import numpy as np
import matplotlib.pyplot as pt
# from astropy.io import fits, ascii
import scipy.integrate as integrate
import fitsio
import h5py
import pandas as pd

def bits(ask="try"):
    if ask == "LC": return 0  #(0 0 0)
    if ask == "downsample": return 1 #(0 0 1)
    if ask == "entireDESIfoot": return 2 #(0 1 0)
    if ask == "SV3foot": return 4 #(1 0 0)
    sys.exit()

def mask(nz=0, Y5=0, sv3=0):
    return nz * (2**0) + Y5 * (2**1) + sv3 * (2**2)

def E_z(z, Om0, Ol0):
    E = np.sqrt( Om0*(1+z)**3 +Ol0 )
    return E

def compute_dist(z_array):
    c = 299792.458
    dH = c / 100
    # Om0 = 0.292
    # Ol0 = 0.708
    # Om0 = 0.3075
    # Ol0 = 0.6925
    Om0 = 0.315191868
    Ol0 = 1 - Om0
    # Om0 = 0.30
    # Ol0 = 0.70
    
    comov_z = np.zeros(len(z_array))
    for i, z in enumerate(z_array):
        dist = comov_z[i] = dH * integrate.quad(lambda z: 1./E_z(z, Om0, Ol0), 0, z)[0]

    return comov_z

def compute_volume(bins_, surface="SV"):
    if surface == "fulldesi":
        A = 14850.4 * np.pi**2 / (180**2)
    elif surface == "fullsky":
        A = 4 * np.pi
    elif surface == "SVCOMPUTE":
        A = 207.5 * np.pi**2 / (180**2)
    else:
        print('ERROR: Incorrect area')
        exit()

    print(f"Area is= {A}")
    z_val = np.array([(a + b) / 2.0 for a, b in zip(bins_[:-1], bins_[1:])])
    vol = np.array([A*(compute_dist([b])[0]**3 - compute_dist([a])[0]**3) / 3.0 for a, b in zip(bins_[:-1], bins_[1:])])
    
    return np.array(z_val), vol
   
def rdz2xyz(ra, dec, redshift):
    H0 = 100.
    h = H0/100.
    ra_rad = ra * np.pi/180.
    dec_rad = dec * np.pi/180.
    dist = compute_dist(redshift)
    x = dist * np.cos(dec_rad) * np.cos(ra_rad) * h
    y = dist * np.cos(dec_rad) * np.sin(ra_rad) * h
    z = dist * np.sin(dec_rad) * h
    return x, y, z


def get_sorted_df(filename):
    fits_hdu    = fitsio.FITS(filename, "r")    

    z_cosmo = fits_hdu[1]["Z_COSMO"][:]
    status  = fits_hdu[1]["STATUS"][:]

    mask_Y5 = mask(nz=0, Y5=1, sv3=0)
    idx     = np.arange(len(status))
    idx_Y5  = idx[((status & (mask_Y5))==mask_Y5) & (z_cosmo <= 1.65)]

    dict_ = {}

    for col in ["RA", "DEC", "Z_COSMO", "Z", "STATUS", "RAN_NUM_0_1"]:
        print(col)
        col_array  = fits_hdu[1][col][:]
        # dict_[col] = col_array[idx_Y5]
        dict_[col] = np.array(col_array[idx_Y5], dtype=np.float32)
    
    print("Sort")
    df_ = pd.DataFrame(data=dict_)
    df_new = df_.sort_values(by=["RA","DEC","Z_COSMO"])
    fits_hdu.close()

    return df_new


def main():
    # new_cat = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed1_Y5_NGC.fits"
    # old_cat = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed1_NGC.fits"
    
    # new_cat = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B2000G512Z1.1N24000470_b0.345d1.45r40c0.05_seed1.fits"
    # old_cat = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B2000G512Z1.1N24000470_b0.345d1.45r40c0.05_seed1.fits"
    
    # old_cat = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/LightCone/LRG/LRG_LC_AbacusSummit_base_c000_ph000.fits"
    # new_cat = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/AbacusSummit/CubicBox/LRG/LightCone/LRG_LC_AbacusSummit_base_c000_ph000.fits"
    
    # new_cat = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_Y5_NGC.fits"
    # old_cat = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_NGC.fits"
    
    # new_cat = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox_2Gpc_old_code/LRG/z0.800/EZmock_B2000G512Z0.8N8015724_b0.385d4r169c0.3_seed1/seed1.suball_shell_{}.h5py"
    # new_cod = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox/LRG/z0.800/EZmock_B2000G512Z0.8N8015724_b0.385d4r169c0.3_seed1/seed_1_shell_{}.hdf5"
    
    # new_cat = "//global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox_2Gpc_old_code/LRG/z0.800/fits/sort_cutsky_LRG_z0.800_EZmock_B2000G512Z0.8N8015724_b0.385d4r169c0.3_seed1.fits"
    new_cod = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox/LRG/z0.800/sort_cutsky_LRG_z0.800_EZmock_B2000G512Z0.8N8015724_b0.385d4r169c0.3_seed1.fits"
    old_cat = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B2000G512Z0.8N8015724_b0.385d4r169c0.3_seed1.fits"
    

    # for i in range(1, 128 + 1):

    #     f_new = h5py.File(new_cat.format(i), "r")
    #     f_old = h5py.File(new_cod.format(i), "r")
    
    #     diff = f_old["galaxy"]["Z_RSD"][:] - f_new["galaxy"]["Z_RSD"][:]
    #     # print(i, diff, np.sum(np.abs(diff)), len(diff[diff!=0]), len(diff))
    #     print(i, np.sum(np.abs(diff)), len(diff[diff!=0]), len(diff))


    #     f_new.close()
    #     f_old.close()
    mask_Y5 = mask(nz=0, Y5=1, sv3=0)

    print("Cat 1")
    df_1     = get_sorted_df(old_cat)
    # fits_1   = fitsio.FITS(old_cat, "r")
    # status_1 = fits_1[1]["STATUS"][:]
    # idx_1    = np.arange(len(status_1))
    
    # idx_1_Y5  = idx_1[((status_1 & (mask_Y5))==mask_Y5)]

    print("Cat 2")
    df_2 = get_sorted_df(new_cod)
    # fits_2 = fitsio.FITS(new_cod, "r")



    
    for col in ["RA", "DEC", "Z_COSMO", "Z", "STATUS", "RAN_NUM_0_1"]:
        # col_1 = fits_1[1][col][:]
        # col_2 = fits_2[1][col][:]

        diff = df_1[col].values - df_2[col].values
        # diff = col_2 - col_1[idx_1_Y5]
        # print(col, diff, np.sum(np.abs(diff)), len(diff[diff!=0]), len(diff))
        print(col, diff, np.sum(np.abs(diff)), len(diff[diff!=0]), len(diff))

    fits_1.close()
    fits_2.close()

if __name__== '__main__':
    main()
