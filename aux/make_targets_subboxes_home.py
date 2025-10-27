import numpy as np
import os
from astropy.table import Table
from desitarget.internal import sharedmem
import h5py
import pandas as pd

#targ = 'ELG'



def func(ifile, opath, targ):
        print(targ)
        if not os.path.isdir(opath):
            os.mkdir(opath)
        
        with h5py.File(ifile, 'r') as f:
            keys = list(f.keys())
            print("Datasets:", keys)
            data = {key: f[key][:] for key in keys}

        df = pd.DataFrame(data)
        
        df['x'] += 1000
        df['y'] += 1000
        df['zreal'] += 1000
        
        # Number of subdivisions per dimension
        num_sub_boxes = 4

        # Create bins for x, y, z coordinates
        bins = np.linspace(0, 2000, num_sub_boxes + 1)

        # Digitize the coordinates to determine which sub-box each point falls into
        df['x_bin'] = np.digitize(df['x'], bins) - 1
        df['y_bin'] = np.digitize(df['y'], bins) - 1
        df['z_bin'] = np.digitize(df['zreal'], bins) - 1

        df['sub_box'] = df['x_bin'] * (num_sub_boxes**2) + df['y_bin'] * num_sub_boxes + df['z_bin']
        sub_boxes = [group for _, group in df.groupby('sub_box')]
        print(len(sub_boxes))
        for j,subdf in enumerate(sub_boxes):

            subdf = subdf.drop(columns=['x_bin', 'y_bin', 'z_bin', 'sub_box'])
            t2 = Table.from_pandas(subdf)
            t2.write(os.path.join(opath, f'{targ}_real_space.sub%d.fits.gz' % j), overwrite=True)


func('/global/cfs/projectdirs/desi/mocks/cai/test_HOMe/LRG_mock.h5', '/global/cfs/projectdirs/desi/mocks/cai/test_HOMe/LRG', 'LRG')
func('/global/cfs/projectdirs/desi/mocks/cai/test_HOMe/ELG_mock.h5', '/global/cfs/projectdirs/desi/mocks/cai/test_HOMe/ELG', 'ELG')
#pool = sharedmem.MapReduce(np=25)
#        #with Pool() as pool:#Pool(processes=nproc) as pool:
#inds = np.arange(0,25)
#with pool:
#    res = pool.map(func, inds)

#func(15)
#func(21)
#func(18)
