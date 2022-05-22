import numpy as np
from astropy.io import fits
import h5py
import fitsio
import os

def create_random(in_path, random_out_path, tracer, snap, random, partname="ph000.gcat.suball_shell_100.h5py", length=1, boxL=1):

    # file_ = in_path + f"/{tracer}_{snap}_{partname}"
    # f = h5py.File(file_, 'r')
    # length = f.attrs["NGAL"]	
    # f.close()

    print("NGAL=",length, "BOXsize=", boxL)

    os.mkdir(random_out_path + f"/S{100 * random}/")

    np.random.seed(int(100 * random))
    X = np.random.uniform(low=0, high=boxL, size=length)
    Y = np.random.uniform(low=0, high=boxL, size=length)
    Z = np.random.uniform(low=0, high=boxL, size=length)
    gal_index = np.arange(length)
    counter = 0
    side_n_subbox = 6
    qboxL = boxL / side_n_subbox
    print(qboxL)
    for i in range(side_n_subbox):
        for j in range(side_n_subbox):
            for k in range(side_n_subbox):
                index_ = i * side_n_subbox * side_n_subbox + j * side_n_subbox + k 
                range_ = (X >= i * qboxL) & (X < (i + 1) * qboxL) & (Y >= j * qboxL) & (Y < (j + 1) * qboxL) & (Z >= k * qboxL) & (Z < (k + 1) * qboxL)
                # print(index_,  i * qboxL, (i + 1) * qboxL, j * qboxL, (j + 1) * qboxL, k * qboxL, (k + 1) * qboxL)
                counter = counter + len(X[range_])
                filename = random_out_path + f"/S{100 * random}/" + tracer + "_SB" + str(index_) + f"_S{100 * random}.fits"
                print(filename)
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

def create_random_old(random_out_path, random, length=1, boxL=1):


    print("NGAL=",length, "BOXsize=", boxL)

    # os.mkdir(random_out_path + f"/S{100 * random}/")

    np.random.seed(int(10000 + 100 * random))
    X = np.random.uniform(low=0, high=boxL, size=length)
    Y = np.random.uniform(low=0, high=boxL, size=length)
    Z = np.random.uniform(low=0, high=boxL, size=length)
    gal_index = np.arange(length)
    

    with h5py.File(random_out_path + f"random_x1_BS3000_S{int(10000 + 100 * random)}_elg_box.dat", 'w') as ff:
        ff.create_group('galaxy')
        ff.create_dataset('galaxy/X', data=X, dtype=np.float32)
        ff.create_dataset('galaxy/Y', data=Y, dtype=np.float32)
        ff.create_dataset('galaxy/Z', data=Z, dtype=np.float32)

                

def main():
    main_out_path = "/global/cscratch1/sd/avariu/desi/test_random/box/"
    # in_path = main_out_path  + "/ELG/z1.100/AbacusSummit_base_c000_ph000/"
    # random_out_path = main_out_path
    # create_random(in_path, random_out_path, "LRG", "snap20", 10000, partname="")

    main_out_path = "/global/cscratch1/sd/avariu/desi/FirstGenMocks/EZmock/CubicBox_6Gpc/"
    for i in range(0, 10):
        # in_path = main_out_path  + "/QSO/z1.400/AbacusSummit_base_c000_ph000/"
        # random_out_path = main_out_path + "/ran_box/QSO/"
        # create_random(in_path, random_out_path, "QSO", "snap12", i, partname="ph000.gcat.suball._shell_100.h5py", boxL=6000, length=27432000)

        # in_path = main_out_path  + "/LRG/z0.800/AbacusSummit_base_c000_ph000/"
        # random_out_path = main_out_path + "/ran_box/LRG/"
        # create_random(in_path, random_out_path, "LRG", "snap20", i, partname="ph000.gcat.suball._shell_100.h5py", boxL=6000, length=218160000)

        # in_path = main_out_path  + "/ELG/z1.100/AbacusSummit_base_c000_ph000/"
        # random_out_path = main_out_path + "/ran_box/ELG/"
        # create_random(in_path, random_out_path, "ELG", "snapX", i, partname="ph000.gcat.suball._shell_100.h5py", boxL=6000, length=648010680)

        create_random_old("/global/cscratch1/sd/avariu/desi/small_box_large_box/3gpc_rand/box/", i, length=88448481, boxL=3000)

if __name__ == '__main__':
	main()