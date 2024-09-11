from astropy.io import fits
import numpy as np
import glob
import h5py 

def count(key="QSO", path="path"):
    snap = np.loadtxt(f"ABACUS_steps_{key}.txt", usecols=(0), unpack=True)
    print(key, snap)
    counters = np.zeros(len(snap))

    for i, snapshot in enumerate(snap):
        files = glob.glob(path + "/" + key + f"*snap{int(snapshot)}.gcat.sub*fits")
        print(len(files))
        counter = 0
        for file_ in files:
            hdul = fits.open(file_)
            data = hdul[1].data
            hdul.close()
            counter += len(data["x"])
        
        print(counter)
        counters[i] = counter
    
    np.savetxt(f"ABACUS_NGAL_{key}.txt", np.array([snap, counters]).T)


def modify_hdf5(key="QSO", path="path"):
    snap, ngal = np.loadtxt(f"ABACUS_NGAL_{key}.txt", usecols=(0,1), unpack=True)

    for i, snapshot in enumerate(snap):
        files = glob.glob(path + f"{key}_snap{int(snapshot)}.gcat.suball._shell_*.h5py")
        print(f"{key}_snap{int(snapshot)}.gcat.suball._shell_*.h5py ", len(files))
        for j, file_ in enumerate(files):
            f = h5py.File(file_, 'r+')
            f.attrs.modify("NGAL", ngal[i])
            f.close()


    

def main():
    # path = "/global/cfs/cdirs/desi/cosmosim/mocks_SV3/v1/ABACUS/ABACUS_Gcat_FITS/"
    path = "/global/cscratch1/sd/avariu/desi/ABACUS_2GPC/"

    key = "QSO"
    # modify_hdf5(key=key, path=path + f"cutsky_shells_sv3_nz_radecz_xyz_{key}_h5py/")
    modify_hdf5(key=key, path=path + f"LC_shells_sv3_nz_radecz_xyz_{key}_h5py/")

    # for i in range (1, 11):
    #     shells_path = f"/ran_{i}_cutsky_shells_sv3_nz_radecz_xyz_{key}_h5py/"
    #     modify_hdf5(key=key, path=path + shells_path)

    # key = "LRG"
    # modify_hdf5(key=key, path=path + f"cutsky_shells_sv3_nz_radecz_xyz_{key}_h5py/")
    # modify_hdf5(key=key, path=path + f"LC_shells_sv3_nz_radecz_xyz_{key}_h5py/")

    # for i in range (1, 11):
    #     shells_path = f"/ran_{i}_cutsky_shells_sv3_nz_radecz_xyz_{key}_h5py/"
    #     modify_hdf5(key=key, path=path + shells_path)

    # key = "ELG"
    # modify_hdf5(key=key, path=path + f"cutsky_shells_sv3_nz_radecz_xyz_{key}_h5py/")
    # modify_hdf5(key=key, path=path + f"LC_shells_sv3_nz_radecz_xyz_{key}_h5py/")
	
    # for i in range (1, 11):
    #     shells_path = f"/ran_{i}_cutsky_shells_sv3_nz_radecz_xyz_{key}_h5py/"
    #     modify_hdf5(key=key, path=path + shells_path)

    
    

if __name__ == '__main__':
	main()
