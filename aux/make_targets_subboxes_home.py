import numpy as np
import os
from astropy.table import Table
from desitarget.internal import sharedmem
import h5py
import pandas as pd

#targ = 'ELG'

def apply_periodic(x, L):
    return (x + 0.5 * L) % L - 0.5 * L

def func(ifile, opath, targ):
        print(targ)
        if not os.path.isdir(opath):
            os.mkdir(opath)
        
        with h5py.File(ifile, 'r') as f:
            keys = list(f.keys())
            print("Datasets:", keys)
            data = {key: f[key][:] for key in keys}

        df = pd.DataFrame(data)
#        print(df.columns)
        
#        print(set(df['class']))

        
#        mask = (df['class'] == b'LRG')
#        print('total', len(df))
#        df = df[mask]
        
        print(targ, len(df), df.columns)
        print(np.min(df['x']), np.max(df['x']))
        print(np.min(df['y']), np.max(df['y']))
        print(np.min(df['z']), np.max(df['z']))
        
        df['x'] += 1000
        df['y'] += 1000
        df['z'] = apply_periodic(df['z'], 2000) + 1000
        
        '''
        lrgs = df.copy()


        if len(lrgs[lrgs['x'] < 0.]) > 0:
            print('X less than 0 is not zero',len(lrgs[lrgs['x'] < 0.]))
            print(lrgs[lrgs['x'] < 0.])

            lrgs['x'][lrgs['x'] < 0.] = 2000 + lrgs['x'][lrgs['x'] < 0.]

        if len(lrgs[lrgs['y'] < 0.]) > 0:
            print('Y less than 0 is not zero',len(lrgs[lrgs['y'] < 0.]))
            print(lrgs[lrgs['y'] < 0.])

            lrgs['y'][lrgs['y'] < 0.] = 2000 + lrgs['y'][lrgs['y'] < 0.]

        if len(lrgs[lrgs['z'] < 0.]) > 0:
            print('Z less than 0 is not zero',len(lrgs[lrgs['z'] < 0.]))
            print(lrgs[lrgs['z'] < 0.])

            lrgs['z'][lrgs['z'] < 0.] = 2000 + lrgs['z'][lrgs['z'] < 0.]
        
        if len(lrgs[lrgs['x'] == 2000.]) > 0:
            print('X equals 2000 is not zero',len(lrgs[lrgs['x'] == 2000.]))
            print(lrgs[lrgs['x'] == 2000.])

            #lrgs['X'][lrgs['X'] == 2000.] = 1999.999

        if len(lrgs[lrgs['y'] == 2000.]) > 0:
            print('Y equals 2000 is not zero',len(lrgs[lrgs['y'] == 2000.]))
            print(lrgs[lrgs['y'] == 2000.])

            #lrgs['Y'][lrgs['Y'] == 2000.] = 1999.999

        if len(lrgs[lrgs['z'] == 2000.]) > 0:
            print('Z equals 2000 is not zero',len(lrgs[lrgs['z'] == 2000.]))
            print(lrgs[lrgs['z'] == 2000.])

            #lrgs['Z'][lrgs['Z'] == 2000.] = 1999.999
        #exit()
        if len(lrgs[lrgs['x'] > 2000.]) > 0:
            print('X greater than 2000 is not zero',len(lrgs[lrgs['x'] > 2000.]))
            print(lrgs[lrgs['x'] > 2000.])

            #lrgs['X'][lrgs['X'] > 2000.] = lrgs['X'][lrgs['X'] > 2000.]-2000

        if len(lrgs[lrgs['y'] > 2000.]) > 0:
            print('Y greater than 2000 is not zero',len(lrgs[lrgs['y'] > 2000.]))
            print(lrgs[lrgs['y'] > 2000.])

            #lrgs['Y'][lrgs['Y'] > 2000.] = lrgs['Y'][lrgs['Y'] > 2000.]-2000

        if len(lrgs[lrgs['z'] > 2000.]) > 0:
            print('Z greater than 2000 is not zero',len(lrgs[lrgs['z'] > 2000.]))
            print(lrgs[lrgs['z'] > 2000.])

            #lrgs['Z'][lrgs['Z'] > 2000.] = lrgs['Z'][lrgs['Z'] > 2000.]-2000


        df = lrgs.copy()
        # Number of subdivisions per dimension
        '''
        num_sub_boxes = 4

        # Create bins for x, y, z coordinates
        bins = np.linspace(0, 2000, num_sub_boxes + 1)

        # Digitize the coordinates to determine which sub-box each point falls into
        df['x_bin'] = np.digitize(df['x'], bins) - 1
        df['y_bin'] = np.digitize(df['y'], bins) - 1
        df['z_bin'] = np.digitize(df['z'], bins) - 1

        print(set(df['x_bin']), set(df['y_bin']),set(df['z_bin']))

        #df.loc[(df['y_bin'] == 4), "y_bin"] = 3


        df['sub_box'] = df['x_bin'] * (num_sub_boxes**2) + df['y_bin'] * num_sub_boxes + df['z_bin']
        sub_boxes = [group for _, group in df.groupby('sub_box')]
        print(len(sub_boxes))
        for j,subdf in enumerate(sub_boxes):

            subdf = subdf.drop(columns=['x_bin', 'y_bin', 'z_bin', 'sub_box'])
            t2 = Table.from_pandas(subdf)
            t2.write(os.path.join(opath, f'{targ}_real_space.sub%d.fits.gz' % j), overwrite=True)


func('/pscratch/sd/a/acarnero/home/LRG_auto_mock.h5', '/pscratch/sd/a/acarnero/home/LRG', 'LRG')
#func('/global/cfs/projectdirs/desi/users/gfavole/aurelio/ELG_LRG_home_Y3_nz_x3.h5', '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/HOMe/Boxes/LRG_period', 'LRG')
#func('/global/cfs/projectdirs/desi/mocks/cai/test_HOMe/ELG_mock.h5', '/global/cfs/projectdirs/desi/mocks/cai/test_HOMe/ELG', 'ELG')
#pool = sharedmem.MapReduce(np=25)
#        #with Pool() as pool:#Pool(processes=nproc) as pool:
#inds = np.arange(0,25)
#with pool:
#    res = pool.map(func, inds)

#func(15)
#func(21)
#func(18)
