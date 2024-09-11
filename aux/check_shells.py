import h5py
import numpy as np 
import glob
import fitsio
import os


def bits(ask="try"):
    if ask == "LC":              return 0 #(0 0 0 0)
    if ask == "downsample":      return 1 #(0 0 0 1)
    if ask == "entireDESIfoot":  return 2 #(0 0 1 0)
    if ask == "SV3foot":         return 4 #(0 1 0 0)
    if ask == "downsample_main": return 8 #(1 0 0 0)
    if ask == "LRG_NZ":          return 32 # (0 0 0 1 0 0 0 0 0)
    if ask == "LRG_NZ_MAIN":     return 64 # (0 0 1 0 0 0 0 0 0)
    if ask == "ELG_NZ":     return 128 # (0 1 0 0 0 0 0 0 0)
    if ask == "QSO_NZ":     return 256 # (1 0 0 0 0 0 0 0 0)

    sys.exit()

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
    
    nz_file = f"../nz_files/{filename}"
    
    z, nz = np.loadtxt(nz_file, usecols=(0, 1), unpack=True)
    z_n = z
    
    nz_n = (nz / ( 1 - failurerate )) * 1.0
    

    # np.savetxt(f"./nz_files/{filename}_redcor.txt", np.array([z_n, nz_n]).T)

    return np.interp(zz, z_n, nz_n, left=0, right=0)

def mask(main=0, nz=0, Y5=0, sv3=0):
    return nz * (2**0) + Y5 * (2**1) + sv3 * (2**2) + main * (2**3)


def count_aux(status_tmp, ra_tmp):
    idx = np.arange(len(status_tmp))
    mask_Y5 = mask(main=0, nz=0, Y5=1, sv3=0)
    
    range_NGC = (ra_tmp < 303.25) & (ra_tmp > 84)

    
    # idx_Y5 = idx[ ((status_tmp & (mask_Y5))==mask_Y5) & range_NGC]
    idx_Y5 = idx[ ((status_tmp & (mask_Y5))==mask_Y5)]

    return idx_Y5


def count_(files):
    print(len(files))

    counter_ = 0

    for i, file_ in enumerate(files):
        print(file_)

        f = h5py.File(file_, 'r')
        data = f['galaxy']
        status_tmp  = data['STATUS'][()]
        ra_tmp      = data['RA'][()]

        idx_Y5 = count_aux(status_tmp, ra_tmp)

        if len(status_tmp[idx_Y5]) != 0:
            counter_ += len(status_tmp[idx_Y5])
        
        
        f.close()
    
    return counter_


def convert(inpath="test", out_file="test", galtype=None, boxL=2000):
    files = glob.glob(inpath + "/*h5py")
    print(len(files))
    
    counter_ = count_(files)

    print(counter_)
    data_fits = np.zeros(counter_, dtype=[('RA', 'f4'), ('DEC', 'f4'), ('Z_COSMO', 'f4'), ('STATUS', 'i4'), ('NZ', 'f4'), ('NZ_MAIN', 'f4'), ('RAW_NZ', 'f4'), ('RAN_NUM_0_1', 'f4'), ('ONEplusDELTA', 'f4')])

    index_i = 0
    index_f = 0
    
    for i, file_ in enumerate(files):
        print(file_)
        print(f"{i}/{len(files)}")
        f = h5py.File(file_, 'r')
        data = f['galaxy']
        ngalbox    = f.attrs["NGAL"]
        n_mean = ngalbox/(boxL**3)

        ra_tmp      = data['RA'][()]
        dec_tmp     = data['DEC'][()]
        z_cosmo_tmp = data['Z_COSMO'][()]
        status_tmp  = data['STATUS'][()]

        ### For ic
        dens_tmp = 1 + data["DENSITY"][()]
        ### For randoms
        # id_tmp      = data['ID'][()]
        ran_num_0_1_tmp  = data['RAN_NUM_0_1'][()]

        idx_Y5 = count_aux(status_tmp, ra_tmp)

        ### NGC
        print(len(ra_tmp), len(dec_tmp[idx_Y5]))

        if len(dec_tmp[idx_Y5]) != 0:

            index_f += len(dec_tmp[idx_Y5])

            data_fits["RA"][index_i: index_f]      = ra_tmp[idx_Y5]
            data_fits["DEC"][index_i: index_f]     = dec_tmp[idx_Y5]
            data_fits["Z_COSMO"][index_i: index_f] = z_cosmo_tmp[idx_Y5]
            data_fits["STATUS"][index_i: index_f]  = status_tmp[idx_Y5]
            data_fits["NZ"][index_i: index_f]      = get_nz(z_cosmo_tmp[idx_Y5], galtype="LRG")

            ### For ic
            data_fits["ONEplusDELTA"][index_i: index_f]      = dens_tmp[idx_Y5]
            
            ### For randoms            
            data_fits["NZ_MAIN"][index_i: index_f]      = get_nz(z_cosmo_tmp[idx_Y5], galtype="LRG_main")
            data_fits["RAW_NZ"][index_i: index_f]  = np.ones(len(dec_tmp[idx_Y5])) * n_mean
            data_fits["RAN_NUM_0_1"][index_i: index_f]  = ran_num_0_1_tmp[idx_Y5]

            index_i = index_f

       

        f.close()

    print(data_fits["RA"][-1]    , ra_tmp[idx_Y5][-1])    


    ### For randoms
      
    with h5py.File(out_file + "_tmp", 'w') as ff:

        for name in ["RA", "DEC", "Z_COSMO", "NZ", "NZ_MAIN", "RAW_NZ", "RAN_NUM_0_1", "ONEplusDELTA"]:
        # name = "RA"
            ff.create_dataset(name,      data=data_fits[name],    dtype=np.float32)
        
        ff.create_dataset("STATUS",      data=data_fits["STATUS"],    dtype=np.int32)


    # hdict = {'SV3_AREA': 207.5, 'Y5_AREA':14850.4}
    
    # fits = fitsio.FITS(out_file+"_tmp", "rw")
    # fits.write(data_fits, header=hdict)
    # fits.close()

    os.rename(out_file+"_tmp", out_file)



def downsample_aux(ran, nz, raw_nz, ask=):
    """ downsample galaxies following n(z) model specified in galtype"""


    # downsample
    nz_selected = ran < nz / raw_nz
    idx         = np.where(nz_selected)
    print("DOWNSAMPLE: Selected {} out of {} galaxies.".format(len(idx[0]), len(ran)), flush=True)

    bitval = bits(ask=ask)
    
    newbits = np.zeros(len(ran), dtype=np.int32)
    newbits[idx] = bitval

    return newbits



file_ = "/global/cscratch1/sd/avariu/desi/hee_jong_seo/ic_Abacus/ic_dens_N576_AbacusSummit_base_c000_ph000_Y5.hdf5"

# file_ = "/global/cscratch1/sd/avariu/4most/IWG/sweep-020m025-030m020.fits"

inhdf = h5py.File(file_, "r")
print(inhdf.keys())
        
# data = inhdf[1][cols_].read()


print("ran and nz")
zz = inhdf["Z_COSMO"][:]

nz_lrg_main = get_nz(zz, galtype="LRG_main")
nz_lrg = get_nz(zz, galtype="LRG")
np.random.seed(1000)
ran_lrg = np.random.rand(len(zz))

# nz_elg = get_nz(zz, galtype="ELG")
# np.random.seed(2000)
# ran_elg = np.random.rand(len(zz))

# nz_qso = get_nz(zz, galtype="QSO")
# np.random.seed(3000)
# ran_qso = np.random.rand(len(zz))


raw_nz = inhdf["RAW_NZ"][:]


# idx = np.arange(len(zz))
# mask_SV = mask(main=0, nz=0, Y5=0, sv3=1)
# idx_SV = idx[(inhdf["STATUS"][:] & (mask_SV))==mask_SV]

print("bits")

newbits_lrg_main = downsample_aux(ran_lrg, nz_lrg_main, raw_nz)

newbits_lrg = downsample_aux(ran_lrg, nz_lrg, raw_nz)
# newbits_elg = downsample_aux(ran_elg, nz_elg, raw_nz)
# newbits_qso = downsample_aux(ran_qso, nz_qso, raw_nz)


# new_status = np.zeros(len(zz), dtype=np.int32) + int(2)
# new_status = new_status.astype(np.int32)

# print(type(new_status))
print("bits or")

newbits_lrg_main_lrg = np.bitwise_or(newbits_lrg_main, newbits_lrg)
newbits_lrg_main_lrg = newbits_lrg_main_lrg.astype(np.int32)
# newbits_lrg_main_lrg_elg = np.bitwise_or(newbits_lrg_main_lrg, newbits_elg)
# newbits_lrg_main_lrg_elg_qso = np.bitwise_or(newbits_lrg_main_lrg_elg, newbits_qso)

# new_status_to_store = np.bitwise_or(new_status, newbits_lrg_main_lrg)
# new_status_to_store = new_status_to_store.astype(np.int32)

print("write")

file_o = "/global/cscratch1/sd/avariu/desi/hee_jong_seo/ic_Abacus/ic_dens_N576_AbacusSummit_base_c000_ph000_Y5_LRG_ELG_QSO.hdf5"

f = h5py.File(file_o, 'w')
# f.create_dataset('RA', data=inhdf["RA"][:],  dtype=np.float32)
# f.create_dataset('DEC', data=inhdf["DEC"][:],  dtype=np.float32)
# f.create_dataset('Z_COSMO', data=inhdf["Z_COSMO"][:],  dtype=np.float32)
# f.create_dataset('ONEplusDELTA', data=inhdf["ONEplusDELTA"][:],  dtype=np.float32)
# f.create_dataset('RAW_NZ', data=inhdf["RAW_NZ"][:],  dtype=np.float32)

f.create_dataset('STATUS', data=newbits_lrg_main_lrg,  dtype=np.int32)
f.create_dataset('NZ_LRG_MAIN', data=nz_lrg_main,  dtype=np.float32)
f.create_dataset('NZ_LRG', data=nz_lrg,  dtype=np.float32)
# f.create_dataset('NZ_ELG', data=nz_qso,  dtype=np.float32)
# f.create_dataset('NZ_QSO', data=nz_qso,  dtype=np.float32)
f.create_dataset('RAN_LRG', data=ran_lrg,  dtype=np.float32)
# f.create_dataset('RAN_ELG', data=ran_elg,  dtype=np.float32)
# f.create_dataset('RAN_QSO', data=ran_qso,  dtype=np.float32)

f.close()

inhdf.close()

# outfile = "/global/cscratch1/sd/avariu/desi/hee_jong_seo/ic_Abacus/ic_dens_N576_AbacusSummit_base_c000_ph000_Y5.fits"
# path = "/global/cscratch1/sd/avariu/desi/hee_jong_seo/cutsky/AbacusSummit_base_c000_ph000/one/" 

# convert(inpath=path, out_file=outfile, galtype="LRG", boxL=2000)



# f0 = h5py.File(path + "/AbacusSummit_base_c000_ph000/ic_dens_576_SBall_shell_127.h5py", "r")
# f1 = h5py.File(path + "/LRG/AbacusSummit_base_c000_ph000/ic_dens_576_SBall_shell_127.h5py", "r")
# print(f0.keys(), f0["galaxy"].keys(), f1["galaxy"].keys())
# RA0 = f0["galaxy"]["Z_RSD"][:]
# RA1 = f0["galaxy"]["Z_COSMO"][:]
# 
# print(RA0 - RA1)
# print(np.unique((RA0 - RA1)))

# for name in ['DEC', 'DENSITY', 'RA', 'Z_COSMO', 'Z_RSD']:
#     print("\n\n", name)
#     x0 = f0["galaxy"][name][:]
#     x1 = f1["galaxy"][name][:]

#     print(len(x1), len(x0))
#     print(np.unique(x0 - x1))


# f0.close()
# f1.close()