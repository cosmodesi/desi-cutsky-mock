import numpy as np
from astropy.io import fits
import h5py
import fitsio

def create_random(in_path, random_out_path, tracer, snap, random, partname="ph000.gcat.suball_shell_100.h5py"):
    boxL = 2000

    file_ = in_path + f"/{tracer}_{snap}_{partname}"
    f = h5py.File(file_, 'r')
    length = f.attrs["NGAL"]	
    f.close()
    print(length)


    np.random.seed(int(100*random))
    X = np.random.uniform(low=0, high=boxL, size=length)
    Y = np.random.uniform(low=0, high=boxL, size=length)
    Z = np.random.uniform(low=0, high=boxL, size=length)
    gal_index = np.arange(length)
    counter = 0 
    for i in range(4):
        for j in range(4):
            for k in range(4):
                index_ = i * 4 * 4 + j*4 + k 
                qboxL = boxL/4.
                range_ = (X >= i * qboxL) & (X < (i + 1) * qboxL) & (Y >= j * qboxL) & (Y < (j + 1) * qboxL) & (Z >= k * qboxL) & (Z < (k + 1) * qboxL)
                # print(index_,  i * qboxL, (i + 1) * qboxL, j * qboxL, (j + 1) * qboxL, k * qboxL, (k + 1) * qboxL)
                counter = counter + len(X[range_])
                filename = random_out_path + tracer + "_" + snap + "_SB" + str(index_) + f"_S{100*random}" + "_ph000.fits"
                # print(filename)
                data = np.zeros(len(X[range_]), dtype=[('x','f4'), ('y','f4'), ('z','f4'), ('id', 'i4')])
                fits = fitsio.FITS(filename, 'rw')
                data["x"] = X[range_]
                data["y"] = Y[range_]
                data["z"] = Z[range_]
                data["id"] = gal_index[range_]
                fits.write(data)
                fits.close()
                
    print(length, counter)
        # with h5py.File(path_instance.dir_out + "/ran_box/" + f"Box_RAN{random}_SB{100*subbox}_SN{10000*shellnum}.fits", 'w') as ff:
        #     ff.create_group('galaxy')
        #     ff.create_dataset('galaxy/X', data=data_r["x"], dtype=np.float32)
        #     ff.create_dataset('galaxy/Y', data=data_r["y"], dtype=np.float32)
        #     ff.create_dataset('galaxy/Z', data=data_r["z"], dtype=np.float32)

        # return data_r, preffix, chilow, chiupp



def main():
    main_out_path = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/AbacusSummit/CubicBox/"

    # for i in range(1, 11):
    in_path = main_out_path  + "/ELG/z1.100/AbacusSummit_base_c000_ph000/"
    random_out_path = main_out_path + "/ELG/ran_box/"
    create_random(in_path, random_out_path, "ELGlowDens", "snap16", 10000, partname="ph000.gcat.suball_shell_100.h5py")

        # in_path = main_out_path  + "/QSO/z1.400/AbacusSummit_base_c000_ph000/"
        # random_out_path = main_out_path + "/QSO/ran_box/"
        # create_random(in_path, random_out_path, "QSO", "snap12", i, partname="ph000.gcat.suball._shell_100.h5py")

        # in_path = main_out_path  + "/LRG/z0.800/AbacusSummit_base_c000_ph000/"
        # random_out_path = main_out_path + "/LRG/ran_box/"
        # create_random(in_path, random_out_path, "LRG", "snap20", i, partname="ph000.gcat.suball._shell_100.h5py")

if __name__ == '__main__':
	main()