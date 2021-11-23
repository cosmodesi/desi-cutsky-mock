import matplotlib
# matplotlib.use("Agg")
import glob
import sys
import os
import numpy as np
import matplotlib.pyplot as pt
from astropy.io import fits, ascii
import scipy.integrate as integrate
import fitsio
import h5py

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
      

def plot_ra_dec_z_2foot(filepath, figname="figname", bins=np.linspace(0.5, 2.1, 33), type_="none"):
    
    print(type_)
    fig0, ax0 = pt.subplots(2, 1, sharex=True, figsize=(8,6), gridspec_kw={"height_ratios":[3,1], "hspace":0})
    fig1, ax1 = pt.subplots(1, figsize=(8,6))
    fig2, ax2 = pt.subplots(1, figsize=(8,6))
    fig3, ax3 = pt.subplots(2, 1, sharex=True, figsize=(8,6), gridspec_kw={"height_ratios":[3,1], "hspace":0})

    infiles = glob.glob(filepath + figname + "/*h5py")
    print(len(infiles))

    ra_Y5 = np.empty(0)
    ra_SV = np.empty(0)
    dec_Y5 = np.empty(0)
    dec_SV = np.empty(0)
    z_cosmo_Y5 = np.empty(0)
    z_cosmo_SV = np.empty(0)
    
    n_shells = np.zeros(len(infiles))
    n_means = np.zeros(len(infiles))
    z_s = np.zeros(len(infiles))
    
    for i, file_ in enumerate(infiles):
        f = h5py.File(file_, 'r')
        data = f['galaxy']
        ngal = f.attrs["NGAL"]

        status = data['STATUS'][()]
        ra_tmp = data['RA'][()]
        dec_tmp = data['DEC'][()]
        z_cosmo_tmp = data['Z_COSMO'][()]            
        f.close()

        ### shell number density
        z, vol = compute_volume(np.array([np.min(z_cosmo_tmp), np.max(z_cosmo_tmp)]), surface="fullsky")
        z_s[i] = z

        n_shell = len(z_cosmo_tmp) / vol
        n_shells[i] = n_shell

        n_mean = ngal / (2000* 2000 * 2000)
        n_means[i] = n_mean
        
        ### Masks SV3 Y5
        idx = np.arange(len(status))

        mask_Y5 = mask(nz=1, Y5=1, sv3=0)
        idx_Y5 = idx[(status & (mask_Y5))==mask_Y5]

        mask_SV = mask(nz=1, Y5=0, sv3=1)
        idx_SV = idx[(status & (mask_SV))==mask_SV]

        if len(dec_tmp) != 0:
            z_cosmo_Y5 = np.concatenate((z_cosmo_Y5, z_cosmo_tmp[idx_Y5]))
            ra_Y5 = np.concatenate((ra_Y5, ra_tmp[idx_Y5]))
            dec_Y5 = np.concatenate((dec_Y5, dec_tmp[idx_Y5]))
    
            z_cosmo_SV = np.concatenate((z_cosmo_SV, z_cosmo_tmp[idx_SV]))
            ra_SV = np.concatenate((ra_SV, ra_tmp[idx_SV]))
            dec_SV = np.concatenate((dec_SV, dec_tmp[idx_SV]))


    #####
    ax0[0].scatter(z_s , n_means, label="n box")
    ax0[0].scatter(z_s , n_shells, label="n shells")
    ax0[1].scatter(z_s , n_shells/n_means)
    
    ax0[0].legend()
    fig0.tight_layout()
    fig0.savefig(filepath + "/figures/" + figname + "_shelldens.png")

    
    #####
    values_Y5, bins_, _ = ax2.hist(z_cosmo_Y5, bins=bins, alpha=0.5, color="red", density=False, label="Y5 Footprint")            
    values_SV, bins_, _ = ax2.hist(z_cosmo_SV, bins=bins, alpha=0.5, color="blue", density=False, label="SV3 Footprint")            
    
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(filepath + "/figures/" + figname + "_number_z.png")


    #####
    z, volume_Y5 = compute_volume(bins_, surface="fulldesi")
    z, volume_SV = compute_volume(bins_, surface="SVCOMPUTE")

    ax3[0].plot(z, values_Y5/volume_Y5, label="Y5 output n(z)", marker=".", color="red")
    ax3[0].plot(z, values_SV/volume_SV, label="SV3 output n(z)", marker=".", color="blue")
    
    z_s, nz_s = np.loadtxt("/global/homes/a/avariu/phd/desicodes/generate_survey_mocks/nz_sv3/" + type_.lower() + "_sm_nz_mycosmo_redcor_ev2.1.txt", usecols=(0,1), unpack=True)
   
    ax3[1].set_xlabel("z_cosmo")
    ax3[0].set_ylabel(r"n(z) [Mpc$^{-3}h^{3}$]")
 
    ax3[1].plot(z, nz_s / (values_Y5/volume_Y5), label="input / output Y5", marker=".", color="green")
    ax3[1].plot(z, nz_s / (values_SV/volume_SV), label="input / output SV", marker=".", color="magenta")
    ax3[0].plot(z, nz_s, label="smooth n(z) / (1-failure_rate)", marker="", color="orange", ls="--")
    
    ax3[1].set_ylim([0.91, 1.09])
    ax3[1].axhline(1, color="grey", ls="--")
    
    ax3[0].legend()
    fig3.tight_layout()
    fig3.savefig(filepath + "/figures/" + figname + "_nz_V.png")

    
    # #####
    ax1.scatter(ra_Y5, dec_Y5, s=0.0001, color="red", label="Y5 Footprint")  
    ax1.scatter(ra_SV, dec_SV, s=0.0001, color="blue", label="SV3 Footprint")  
    
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(filepath + "/figures/" + figname + "_RA_DEC.png")

    
def main():
    pt.rcParams.update({'font.size': 12})

    filepath = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/AbacusSummit/CubicBox/"
    # dict_ = {"qso":np.linspace(0.35, 2.1, 71), "lrg":np.linspace(0.01, 1.6, 80), "elg":np.linspace(0.5936708860759494, 1.6, 80-29)}
    dict_ = {"QSO":np.linspace(0.6, 4.5, 79), "LRG":np.linspace(0.01, 1.6, 81), "ELG":np.linspace(0.01, 1.6, 81)}
    

    # for i in range(1, 11):
    #     seed = str(100*i)
    #     key = "QSO"
    #     plot_ra_dec_z_2foot(filepath + f"/{key}/", figname=f"{key}_ran_S{seed}_shells_ph000", bins=dict_[key], type_=key)

    #     key = "ELG"
    #     plot_ra_dec_z_2foot(filepath + f"/{key}/", figname=f"{key}_ran_S{seed}_shells_ph000", bins=dict_[key], type_=key)
        
    #     key = "LRG"
    #     plot_ra_dec_z_2foot(filepath + f"/{key}/", figname=f"{key}_ran_S{seed}_shells_ph000", bins=dict_[key], type_=key)



    # phase = str(int(0)).zfill(3)
    # key = "elg"
    # plot_ra_dec_z_2foot(filepath + "/ELG/z1.100/", figname=f"AbacusSummit_base_c000_ph{phase}", bins=dict_[key], type_=key)
    # exit()
    for i in range(1, 25):
        phase = str(int(i)).zfill(3)
        key = "elg"
        plot_ra_dec_z_2foot(filepath + "/ELG/z1.100/", figname=f"AbacusSummit_base_c000_ph{phase}", bins=dict_[key], type_=key)

    #     key = "lrg"
    #     plot_ra_dec_z_2foot(filepath + "/LRG/z0.800/", figname=f"AbacusSummit_base_c000_ph{phase}", bins=dict_[key], type_=key)


    key = "qso"
    for i in range(13, 25):
        phase = str(int(i)).zfill(3)
        plot_ra_dec_z_2foot(filepath + "/QSO/z1.400/", figname=f"AbacusSummit_base_c000_ph{phase}", bins=dict_[key], type_=key)

    
    

if __name__== '__main__':
    main()
