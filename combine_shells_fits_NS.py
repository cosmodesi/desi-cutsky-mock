import glob
import numpy as np
import h5py
import fitsio
import os

from foot_nz import get_nz, mask
from redshift_error_QSO import sample_redshift_error

def count_aux(status_tmp, ra_tmp):
    idx = np.arange(len(status_tmp))
    mask_Y5 = mask(main=0, nz=0, Y5=1, sv3=0)
        
    range_NGC = (ra_tmp < 303.25) & (ra_tmp > 84)
    range_SGC = ~range_NGC

    idx_Y5_NGC = idx[ ((status_tmp & (mask_Y5))==mask_Y5) & range_NGC ]
    idx_Y5_SGC = idx[ ((status_tmp & (mask_Y5))==mask_Y5) & range_SGC ]

    return idx_Y5_NGC, idx_Y5_SGC

def count_NGC_SGC(files):
    print(len(files))

    counter_NGC = 0
    counter_SGC = 0

    for i, file_ in enumerate(files):
        print(file_)

        f = h5py.File(file_, 'r')
        data = f['galaxy']
        status_tmp  = data['STATUS'][()]
        ra_tmp      = data['RA'][()]

        idx_Y5_NGC, idx_Y5_SGC = count_aux(status_tmp, ra_tmp)

        if len(status_tmp[idx_Y5_NGC]) != 0:
            counter_NGC = counter_NGC + len(status_tmp[idx_Y5_NGC])
        
        if len(status_tmp[idx_Y5_SGC]) != 0:
            counter_SGC = counter_SGC + len(status_tmp[idx_Y5_SGC])
        
        f.close()
    
    return counter_NGC, counter_SGC


def convert(inpath="test", out_file="test", galtype=None, boxL=2000, seed=0, max_seed=0):

    files = glob.glob(inpath + "/*h5py")
    print(len(files))
    
    counter_NGC, counter_SGC = count_NGC_SGC(files)

    if galtype == "LRG":
        data_fits_NGC = np.zeros(counter_NGC, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('NZ_MAIN', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4')])#, ('ID', 'i4')])
        data_fits_SGC = np.zeros(counter_SGC, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('NZ_MAIN', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4')])#, ('ID', 'i4')])
    elif galtype == "QSO":
        data_fits_NGC = np.zeros(counter_NGC, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('Z_ERR_3GAUSS', 'f4'), ('Z_ERR_SIG500', 'f4')])#, ('ID', 'i4')])
        data_fits_SGC = np.zeros(counter_SGC, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('Z_ERR_3GAUSS', 'f4'), ('Z_ERR_SIG500', 'f4')])#, ('ID', 'i4')])
    else:
        data_fits_NGC = np.zeros(counter_NGC, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('ID', 'i4')])
        data_fits_SGC = np.zeros(counter_SGC, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('ID', 'i4')])

    index_i_NGC = 0
    index_f_NGC = 0

    index_i_SGC = 0
    index_f_SGC = 0
    
    
    for i, file_ in enumerate(files):
        print(file_)
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

        ### For randoms
        # id_tmp      = data['ID'][()]
        ran_num_0_1_tmp  = data['RAN_NUM_0_1'][()]

        idx_Y5_NGC, idx_Y5_SGC = count_aux(status_tmp, ra_tmp)

        ### NGC
        print(len(ra_tmp), len(dec_tmp[idx_Y5_NGC]))

        if len(dec_tmp[idx_Y5_NGC]) != 0:

            index_f_NGC = index_i_NGC + len(dec_tmp[idx_Y5_NGC])

            data_fits_NGC["RA"][index_i_NGC: index_f_NGC]      = ra_tmp[idx_Y5_NGC]
            data_fits_NGC["DEC"][index_i_NGC: index_f_NGC]     = dec_tmp[idx_Y5_NGC]
            data_fits_NGC["Z"][index_i_NGC: index_f_NGC]       = z_rsd_tmp[idx_Y5_NGC]
            data_fits_NGC["Z_COSMO"][index_i_NGC: index_f_NGC] = z_cosmo_tmp[idx_Y5_NGC]
            data_fits_NGC["STATUS"][index_i_NGC: index_f_NGC]  = status_tmp[idx_Y5_NGC]
            data_fits_NGC["NZ"][index_i_NGC: index_f_NGC]      = get_nz(z_cosmo_tmp[idx_Y5_NGC], galtype=galtype)

            ### For randoms
            # data_fits_NGC["ID"][index_i_NGC: index_f_NGC]      = id_tmp[idx_Y5_NGC]
            
            if galtype == "LRG":
                data_fits_NGC["NZ_MAIN"][index_i_NGC: index_f_NGC]      = get_nz(z_cosmo_tmp[idx_Y5_NGC], galtype="LRG_main")
            elif galtype == "QSO":
                data_fits_NGC["Z_ERR_3GAUSS"][index_i_NGC: index_f_NGC] = sample_redshift_error(z_rsd_tmp[idx_Y5_NGC], error_model='3gauss')
                data_fits_NGC["Z_ERR_SIG500"][index_i_NGC: index_f_NGC] = sample_redshift_error(z_rsd_tmp[idx_Y5_NGC], error_model='sig500')

            data_fits_NGC["RAW_NZ"][index_i_NGC: index_f_NGC]  = np.ones(len(dec_tmp[idx_Y5_NGC])) * n_mean
            data_fits_NGC["RAN_NUM_0_1"][index_i_NGC: index_f_NGC]  = ran_num_0_1_tmp[idx_Y5_NGC]

            index_i_NGC = index_f_NGC

        if len(dec_tmp[idx_Y5_SGC]) != 0:

            index_f_SGC = index_i_SGC + len(dec_tmp[idx_Y5_SGC])

            data_fits_SGC["RA"][index_i_SGC: index_f_SGC]      = ra_tmp[idx_Y5_SGC]
            data_fits_SGC["DEC"][index_i_SGC: index_f_SGC]     = dec_tmp[idx_Y5_SGC]
            data_fits_SGC["Z"][index_i_SGC: index_f_SGC]       = z_rsd_tmp[idx_Y5_SGC]
            data_fits_SGC["Z_COSMO"][index_i_SGC: index_f_SGC] = z_cosmo_tmp[idx_Y5_SGC]
            data_fits_SGC["STATUS"][index_i_SGC: index_f_SGC]  = status_tmp[idx_Y5_SGC]
            data_fits_SGC["NZ"][index_i_SGC: index_f_SGC]      = get_nz(z_cosmo_tmp[idx_Y5_SGC], galtype=galtype)
            
            ### For randoms
            # data_fits_SGC["ID"][index_i_SGC: index_f_SGC]      = id_tmp[idx_Y5_SGC]
            
            if galtype == "LRG":
                data_fits_SGC["NZ_MAIN"][index_i_SGC: index_f_SGC]      = get_nz(z_cosmo_tmp[idx_Y5_SGC], galtype="LRG_main")
            elif galtype == "QSO":
                data_fits_SGC["Z_ERR_3GAUSS"][index_i_SGC: index_f_SGC] = sample_redshift_error(z_rsd_tmp[idx_Y5_SGC], error_model='3gauss')
                data_fits_SGC["Z_ERR_SIG500"][index_i_SGC: index_f_SGC] = sample_redshift_error(z_rsd_tmp[idx_Y5_SGC], error_model='sig500')

            data_fits_SGC["RAW_NZ"][index_i_SGC: index_f_SGC]  = np.ones(len(dec_tmp[idx_Y5_SGC])) * n_mean
            data_fits_SGC["RAN_NUM_0_1"][index_i_SGC: index_f_SGC]  = ran_num_0_1_tmp[idx_Y5_SGC]

            index_i_SGC = index_f_SGC

        f.close()

    print(data_fits_NGC["RA"][-1], ra_tmp[idx_Y5_NGC][-1])    
    print(data_fits_SGC["RA"][-1], ra_tmp[idx_Y5_SGC][-1])    


    out_file_NGC = out_file.format(seed) + "_NGC.fits"
    out_file_SGC = out_file.format(max_seed - seed + 1) + "_SGC.fits"
    ### For randoms
    # out_file_SGC = out_file.format(max_seed - seed + 100) + "_SGC.fits"
      
    hdict = {'SV3_AREA': 207.5, 'Y5_AREA':14850.4}
    
    fits = fitsio.FITS(out_file_NGC+"_tmp", "rw")
    fits.write(data_fits_NGC, header=hdict)
    fits.close()

    os.rename(out_file_NGC+"_tmp", out_file_NGC)

    fits = fitsio.FITS(out_file_SGC+"_tmp", "rw")
    fits.write(data_fits_SGC, header=hdict)
    fits.close()

    os.rename(out_file_SGC+"_tmp", out_file_SGC)
