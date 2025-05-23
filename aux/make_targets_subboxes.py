import numpy as np
import os
from astropy.table import Table
from desitarget.internal import sharedmem

targ = 'ELG'
#phase = '000'


snapshots = {'LRG':[['0p500','z0.500'],['0p725','z0.725'],['0p950','z0.950']], 'QSO':[['1p400','z1.400']], 'ELG':[['0p950', 'z0.950'], ['1p175', 'z1.175'],['1p475', 'z1.475']]}



def func(i):

#for i in range(0,24):
    phase = str(int(i)).zfill(3)
    for j in range(len(snapshots[targ])):
        opath = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{phase}/Boxes/{targ}/{snapshots[targ][j][1]}'
        if not os.path.isdir(opath):
            os.mkdir(opath)
        ipath = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{phase}/Boxes/{targ}/abacus_HF_{targ}_{snapshots[targ][j][0]}_DR2_v1.0_AbacusSummit_base_c000_ph{phase}_clustering.dat.fits'
        if not os.path.isfile(ipath):
            print(ipath, 'dont exist')
            continue
        print(opath)
        print(ipath)
        lrgs = Table.read(ipath)
        lrgs['X'] += 1000
        lrgs['Y'] += 1000
        lrgs['Z'] += 1000
        lrgs['X_RSD'] += 1000
        lrgs['Y_RSD'] += 1000
        lrgs['Z_RSD'] += 1000

        if len(lrgs[lrgs['X'] == 0.]) > 0:
            print('X equals 0 is not zero',len(lrgs[lrgs['X'] == 0.]))
            print(lrgs[lrgs['X'] == 0.])

            lrgs['X'][lrgs['X'] == 0.] = 0.0001

        if len(lrgs[lrgs['Y'] == 0.]) > 0:
            print('Y equals 0 is not zero',len(lrgs[lrgs['Y'] == 0.]))
            print(lrgs[lrgs['Y'] == 0.])

            lrgs['Y'][lrgs['Y'] == 0.] = 0.0001

        if len(lrgs[lrgs['Z'] == 0.]) > 0:
            print('Z equals 0 is not zero',len(lrgs[lrgs['Z'] == 0.]))
            print(lrgs[lrgs['Z'] == 0.])

            lrgs['Z'][lrgs['Z'] == 0.] = 0.0001
        
        if len(lrgs[lrgs['X'] < 0.]) > 0:
            print('X less than 0 is not zero',len(lrgs[lrgs['X'] < 0.]))
            print(lrgs[lrgs['X'] < 0.])

            lrgs['X'][lrgs['X'] < 0.] = 2000 + lrgs['X'][lrgs['X'] < 0.]

        if len(lrgs[lrgs['Y'] < 0.]) > 0:
            print('Y less than 0 is not zero',len(lrgs[lrgs['Y'] < 0.]))
            print(lrgs[lrgs['Y'] < 0.])

            lrgs['Y'][lrgs['Y'] < 0.] = 2000 + lrgs['Y'][lrgs['Y'] < 0.]

        if len(lrgs[lrgs['Z'] < 0.]) > 0:
            print('Z less than 0 is not zero',len(lrgs[lrgs['Z'] < 0.]))
            print(lrgs[lrgs['Z'] < 0.])

            lrgs['Z'][lrgs['Z'] < 0.] = 2000 + lrgs['Z'][lrgs['Z'] < 0.]
        
        if len(lrgs[lrgs['X'] == 2000.]) > 0:
            print('X equals 2000 is not zero',len(lrgs[lrgs['X'] == 2000.]))
            print(lrgs[lrgs['X'] == 2000.])

            lrgs['X'][lrgs['X'] == 2000.] = 1999.999

        if len(lrgs[lrgs['Y'] == 2000.]) > 0:
            print('Y equals 2000 is not zero',len(lrgs[lrgs['Y'] == 2000.]))
            print(lrgs[lrgs['Y'] == 2000.])

            lrgs['Y'][lrgs['Y'] == 2000.] = 1999.999

        if len(lrgs[lrgs['Z'] == 2000.]) > 0:
            print('Z equals 2000 is not zero',len(lrgs[lrgs['Z'] == 2000.]))
            print(lrgs[lrgs['Z'] == 2000.])

            lrgs['Z'][lrgs['Z'] == 2000.] = 1999.999

        if len(lrgs[lrgs['X'] > 2000.]) > 0:
            print('X greater than 2000 is not zero',len(lrgs[lrgs['X'] > 2000.]))
            print(lrgs[lrgs['X'] > 2000.])

            lrgs['X'][lrgs['X'] > 2000.] = lrgs['X'][lrgs['X'] > 2000.]-2000

        if len(lrgs[lrgs['Y'] > 2000.]) > 0:
            print('Y greater than 2000 is not zero',len(lrgs[lrgs['Y'] > 2000.]))
            print(lrgs[lrgs['Y'] > 2000.])

            lrgs['Y'][lrgs['Y'] > 2000.] = lrgs['Y'][lrgs['Y'] > 2000.]-2000

        if len(lrgs[lrgs['Z'] > 2000.]) > 0:
            print('Z greater than 2000 is not zero',len(lrgs[lrgs['Z'] > 2000.]))
            print(lrgs[lrgs['Z'] > 2000.])

            lrgs['Z'][lrgs['Z'] > 2000.] = lrgs['Z'][lrgs['Z'] > 2000.]-2000




        df = lrgs.to_pandas()

        # Number of subdivisions per dimension
        num_sub_boxes = 4

        # Create bins for x, y, z coordinates
        bins = np.linspace(0, 2000, num_sub_boxes + 1)

        # Digitize the coordinates to determine which sub-box each point falls into
        df['x_bin'] = np.digitize(df['X'], bins) - 1
        df['y_bin'] = np.digitize(df['Y'], bins) - 1
        df['z_bin'] = np.digitize(df['Z'], bins) - 1

        df['sub_box'] = df['x_bin'] * (num_sub_boxes**2) + df['y_bin'] * num_sub_boxes + df['z_bin']
        sub_boxes = [group for _, group in df.groupby('sub_box')]
        print(len(sub_boxes))
        for j,subdf in enumerate(sub_boxes):

            subdf = subdf.drop(columns=['x_bin', 'y_bin', 'z_bin', 'sub_box'])
            t2 = Table.from_pandas(subdf)
            t2.write(os.path.join(opath, f'{targ}_real_space.sub%d.fits.gz' % j), overwrite=True)

'''
pool = sharedmem.MapReduce(np=25)
        #with Pool() as pool:#Pool(processes=nproc) as pool:
inds = np.arange(0,25)
with pool:
    res = pool.map(func, inds)
'''
func(15)
func(21)
func(18)
