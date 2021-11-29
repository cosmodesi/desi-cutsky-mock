import glob
import numpy as np
import h5py
import fitsio

from foot_nz import nz_oneperc


def convert(inpath="test", out_file="test", galtype=None, boxL=2000):
    files = glob.glob(inpath + "/*h5py")
    print(len(files))
    counter = 0
    for i, file_ in enumerate(files):
        f = h5py.File(file_, 'r')
        data = f['galaxy']
        ra_tmp      = data['RA'][()]
        if len(ra_tmp) == 0:
            continue
        counter = counter + len(ra_tmp)
        f.close()
    
    if galtype == "LRG":
        data_fits = np.zeros(counter, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('NZ_MAIN', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4')])
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
        ran_num_0_1_tmp  = data['RAN_NUM_0_1'][()]

        if len(dec_tmp) == 0:
            continue
        index_f = index_i + len(dec_tmp)
        data_fits["RA"][index_i: index_f]      = ra_tmp
        data_fits["DEC"][index_i: index_f]     = dec_tmp
        data_fits["Z"][index_i: index_f]       = z_rsd_tmp
        data_fits["Z_COSMO"][index_i: index_f] = z_cosmo_tmp
        data_fits["STATUS"][index_i: index_f]  = status_tmp
        data_fits["NZ"][index_i: index_f]      = nz_oneperc(z_cosmo_tmp, galtype=galtype)
        if galtype == "LRG":
            data_fits["NZ_MAIN"][index_i: index_f]      = nz_oneperc(z_cosmo_tmp, galtype="LRG_main")

        data_fits["RAW_NZ"][index_i: index_f]  = np.ones(len(dec_tmp)) * n_mean
        data_fits["RAN_NUM_0_1"][index_i: index_f]  = ran_num_0_1_tmp

        index_i = index_f

        f.close()

    print(data_fits["RA"][-1], ra_tmp[-1])    

      
    hdict = {'SV3_AREA': 207.5, 'Y5_AREA':14850.4}
    fits = fitsio.FITS(out_file, "rw")
    fits.write(data_fits, header=hdict)
    fits.close()
