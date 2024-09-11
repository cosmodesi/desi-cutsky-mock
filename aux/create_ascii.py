import glob
import sys
import numpy as np
import h5py
from astropy.io import fits
import fitsio

def bits(ask="try"):
    if ask == "LC": return 0  #(0 0 0)
    if ask == "downsample": return 1 #(0 0 1)
    if ask == "entireDESIfoot": return 2 #(0 1 0)
    if ask == "SV3foot": return 4 #(1 0 0)
    sys.exit()


def mask(nz=0, Y5=0, sv3=0):
    return nz * (2**0) + Y5 * (2**1) + sv3 * (2**2)


def nz_oneperc(zz, galtype):
    path = "/global/homes/a/avariu/phd/desicodes/generate_survey_mocks/nz_sv3/"
    if galtype == "LRG":
        z, nz = np.loadtxt(path + "sm_LRG_mycosmo_ev2.1.dat", usecols=(0,1), unpack=True)
        failurerate = 0.02
    elif galtype == "ELG":
        z, nz = np.loadtxt(path + "sm_ELG_mycosmo_ev2.1.dat", usecols=(0,1), unpack=True)
        failurerate = 0.25
    elif galtype == "QSO":
        z, nz = np.loadtxt(path + "sm_QSO_mycosmo_ev2.1.dat", usecols=(0,1), unpack=True)
        failurerate = 0.37
    else:
        raise RuntimeError("Unknown galaxy type.")

    nz = (nz / ( 1 - failurerate )) * 1.0
    
    nz_n = nz
    z_n = z
    
    return np.interp(zz, z_n, nz_n, left=0, right=0)


def convert(path, folname, key, redshift):
    files = glob.glob(path + folname + "/*h5py")
    print(len(files))

    counter = 0
    for i, file_ in enumerate(files):
        f = h5py.File(file_, 'r')
        data = f['galaxy']
        ra_tmp      = data['RA'][()]
        status_tmp  = data['STATUS'][()]
        idx = np.arange(len(status_tmp))

        mask_SV = mask(nz=1, Y5=0, sv3=1)
        idx_SV = idx[(status_tmp & (mask_SV))==mask_SV]

        ra_tmp = ra_tmp[idx_SV]

        if len(ra_tmp) == 0:
            continue
        counter = counter + len(ra_tmp)
        f.close()
            
    # ra_     = np.empty(counter, dtype=np.float64)
    # dec_    = np.empty(counter, dtype=np.float64)
    # z_rsd   = np.empty(counter, dtype=np.float64)
    # z_cosmo = np.empty(counter, dtype=np.float64)
    # status  = np.empty(counter, dtype=np.int32)
    
    data_fits = np.zeros(counter, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('NZ', 'f4')])

    index_i = 0
    index_f = 0
    boxL = 2000
    
    for i, file_ in enumerate(files):
        print(f"{i}/{len(files)}")
        f = h5py.File(file_, 'r')
        data = f['galaxy']
        ngalbox    = f.attrs["NGAL"]
        n_mean = ngalbox/(boxL**3)

        ra_tmp      = data['RA'][()]
        dec_tmp     = data['DEC'][()]
        z_rsd_tmp   = data['Z_RSD'][()]
        z_cosmo_tmp = data['Z_COSMO'][()]
        status_tmp  = data['STATUS'][()]
       
        idx = np.arange(len(status_tmp))
        mask_SV = mask(nz=1, Y5=0, sv3=1)
        idx_SV = idx[(status_tmp & (mask_SV))==mask_SV]

        ra_tmp = ra_tmp[idx_SV]
        dec_tmp = dec_tmp[idx_SV]
        z_rsd_tmp = z_rsd_tmp[idx_SV]
        z_cosmo_tmp = z_cosmo_tmp[idx_SV]

        if len(dec_tmp) == 0:
            continue
        index_f = index_i + len(dec_tmp)
        data_fits["RA"][index_i: index_f]      = ra_tmp
        data_fits["DEC"][index_i: index_f]     = dec_tmp
        data_fits["Z"][index_i: index_f]       = z_rsd_tmp
        data_fits["Z_COSMO"][index_i: index_f] = z_cosmo_tmp
        data_fits["NZ"][index_i: index_f]      = nz_oneperc(z_cosmo_tmp, key)
        
        index_i = index_f

        f.close()


    np.savetxt(path + "/fits/cutsky_" + key + "_" + redshift + "_" + folname + ".txt", data_fits)
    # np.savetxt(path + "/fits/" + folname + "_RANDOM_1X.txt", data_fits)

def main():
    
    filepath = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/AbacusSummit/CubicBox/"

    # for i in range(1, 11):
    #     seed = str(100*i)
    #     key = "QSO"      
    #     convert(filepath + f"/{key}/", folname=f"{key}_ran_S{seed}_shells_ph000", key=key, redshift="")

    #     key = "ELG"      
    #     convert(filepath + f"/{key}/", folname=f"{key}_ran_S{seed}_shells_ph000", key=key, redshift="")

    #     key = "LRG"      
    #     convert(filepath + f"/{key}/", folname=f"{key}_ran_S{seed}_shells_ph000", key=key, redshift="")


    for i in range(0, 25):    
        phase = str(int(i)).zfill(3)
        key = "LRG"
        redshift = "z0.800"
        convert(filepath + f"/{key}/{redshift}/", folname=f"AbacusSummit_base_c000_ph{phase}", key=key, redshift=redshift)

    for i in range(0, 25):
        phase = str(int(i)).zfill(3)
        key = "QSO"      
        redshift = "z1.400"
        convert(filepath + f"/{key}/{redshift}/", folname=f"AbacusSummit_base_c000_ph{phase}", key=key, redshift=redshift)

    for i in range(0, 25):
        phase = str(int(i)).zfill(3)
        key = "ELG"
        redshift = "z1.100"
        convert(filepath + f"/{key}/{redshift}/", folname=f"AbacusSummit_base_c000_ph{phase}", key=key, redshift=redshift)


    

    
if __name__== '__main__':
    main()
