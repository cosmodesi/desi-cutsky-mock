import glob
import numpy as np
import h5py
import fitsio
import os

from foot_nz import get_nz, mask
from redshift_error_QSO import sample_redshift_error

def convert(inpath="test", out_file="test", galtype=None, boxL=2000, seed=0, max_seed=0):
    files = glob.glob(inpath + "/*h5py")
    print(len(files))
    counter = 0
    for i, file_ in enumerate(files):
        f = h5py.File(file_, 'r')
        data = f['galaxy']
        ra_tmp      = data['RA'][()]
        status_tmp  = data['STATUS'][()]

        idx = np.arange(len(status_tmp))
        mask_Y5 = mask(main=0, nz=0, Y5=1, sv3=0)
        idx_Y5 = idx[(status_tmp & (mask_Y5))==mask_Y5]


        if len(ra_tmp[idx_Y5]) == 0:
            continue
        counter = counter + len(ra_tmp[idx_Y5])
        f.close()
    
    if galtype == "LRG":
        data_fits = np.zeros(counter, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('NZ_MAIN', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('ID', 'i4')])
    elif galtype == "QSO":
        data_fits = np.zeros(counter, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('Z_ERR_3GAUSS', 'f4'), ('Z_ERR_SIG500', 'f4')])
    else:
        data_fits = np.zeros(counter, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4')])

    index_i = 0
    index_f = 0
    
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
        id_tmp      = data['ID'][()]
        ran_num_0_1_tmp  = data['RAN_NUM_0_1'][()]

        idx = np.arange(len(status_tmp))
        mask_Y5 = mask(main=0, nz=0, Y5=1, sv3=0)
        idx_Y5 = idx[(status_tmp & (mask_Y5))==mask_Y5]

        print(len(ra_tmp), len(dec_tmp[idx_Y5]))

        if len(dec_tmp[idx_Y5]) == 0:
            continue
        index_f = index_i + len(dec_tmp[idx_Y5])
        data_fits["RA"][index_i: index_f]      = ra_tmp[idx_Y5]
        data_fits["DEC"][index_i: index_f]     = dec_tmp[idx_Y5]
        data_fits["Z"][index_i: index_f]       = z_rsd_tmp[idx_Y5]
        data_fits["Z_COSMO"][index_i: index_f] = z_cosmo_tmp[idx_Y5]
        data_fits["STATUS"][index_i: index_f]  = status_tmp[idx_Y5]
        data_fits["NZ"][index_i: index_f]      = get_nz(z_cosmo_tmp[idx_Y5], galtype=galtype)
        data_fits["ID"][index_i: index_f]      = id_tmp[idx_Y5]
        
        if galtype == "LRG":
            data_fits["NZ_MAIN"][index_i: index_f]      = get_nz(z_cosmo_tmp[idx_Y5], galtype="LRG_main")
        elif galtype == "QSO":
            data_fits["Z_ERR_3GAUSS"][index_i: index_f] = sample_redshift_error(z_rsd_tmp[idx_Y5], error_model='3gauss')
            data_fits["Z_ERR_SIG500"][index_i: index_f] = sample_redshift_error(z_rsd_tmp[idx_Y5], error_model='sig500')

        data_fits["RAW_NZ"][index_i: index_f]  = np.ones(len(dec_tmp[idx_Y5])) * n_mean
        data_fits["RAN_NUM_0_1"][index_i: index_f]  = ran_num_0_1_tmp[idx_Y5]

        index_i = index_f

        f.close()

    print(data_fits["RA"][-1], ra_tmp[idx_Y5][-1])    

    out_file = out_file.format(seed) + ".fits"

    hdict = {'SV3_AREA': 207.5, 'Y5_AREA':14850.4}
    fits = fitsio.FITS(out_file+"_tmp", "rw")
    fits.write(data_fits, header=hdict)
    fits.close()

    os.rename(out_file+"_tmp", out_file)
