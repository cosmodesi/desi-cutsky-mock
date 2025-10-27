import numpy as np
import os
from astropy.table import Table
from desitarget.internal import sharedmem

targ = 'QSO'
#phase = '000'


snapshots = {'LRG':[['0p500','z0.500'],['0p725','z0.725'],['0p950','z0.950']], 'QSO':[['2p000','z2.000'], ['2p500','z2.500'], ['3p000','z3.000']], 'ELG':[['0p950', 'z0.950'], ['1p175', 'z1.175'],['1p475', 'z1.475']]}
#snapshots = {'LRG':[['0p500','z0.500'],['0p725','z0.725'],['0p950','z0.950']], 'QSO':[['1p400','z1.400']], 'ELG':[['0p950', 'z0.950'], ['1p175', 'z1.175'],['1p475', 'z1.475']]}

snaps = ['AbacusSummit_base_c001_ph000', 'AbacusSummit_base_c001_ph004', 'AbacusSummit_base_c002_ph002', 'AbacusSummit_base_c003_ph000', 'AbacusSummit_base_c003_ph004',
         'AbacusSummit_base_c004_ph002', 'AbacusSummit_base_c001_ph001', 'AbacusSummit_base_c001_ph005', 'AbacusSummit_base_c002_ph003', 'AbacusSummit_base_c003_ph001',
         'AbacusSummit_base_c003_ph005', 'AbacusSummit_base_c004_ph003', 'AbacusSummit_base_c001_ph002', 'AbacusSummit_base_c002_ph000', 'AbacusSummit_base_c002_ph004',
         'AbacusSummit_base_c003_ph002', 'AbacusSummit_base_c004_ph000', 'AbacusSummit_base_c004_ph004', 'AbacusSummit_base_c001_ph003', 'AbacusSummit_base_c002_ph001',
         'AbacusSummit_base_c002_ph005', 'AbacusSummit_base_c003_ph003', 'AbacusSummit_base_c004_ph001', 'AbacusSummit_base_c004_ph005']

# /global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/High_dens_Boxes/ELG
def func(i):
    thissnap = snaps[i]
#for i in range(0,24):
    phase = thissnap.split('ph')[-1]
    cosmo = thissnap.split('_')[2]
    for j in range(len(snapshots[targ])):
        opath = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/alternative_cosmologies/{thissnap}/Boxes/{targ}/{snapshots[targ][j][1]}'
        #opath = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{phase}/High_dens_Boxes/{targ}/{snapshots[targ][j][1]}'
        if not os.path.isdir(opath):
            os.mkdir(opath)
        ipath = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/alternative_cosmologies/{thissnap}/Boxes/{targ}/abacus_HF_{targ}_{snapshots[targ][j][0]}_DR2_v1.0_AbacusSummit_base_{cosmo}_ph{phase}_clustering.dat.fits'
        #ipath = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{phase}/High_dens_Boxes/{targ}/abacus_HF_{targ}_{snapshots[targ][j][0]}_DR2_v1.0_AbacusSummit_base_c000_ph{phase}_clustering.dat.fits'
        if not os.path.isfile(ipath):
            print(ipath, 'dont exist')
            continue
        print(opath)
        print(ipath)
        lrgs = Table.read(ipath)
        lrgs['X'] += 1000
        lrgs['Y'] += 1000
        lrgs['Z'] += 1000
        if 'X_RSD' in lrgs.columns:
            lrgs['X_RSD'] += 1000
        if 'Y_RSD' in lrgs.columns:
            lrgs['Y_RSD'] += 1000
        if 'Z_RSD' in lrgs.columns:
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


pool = sharedmem.MapReduce(np=24)
        #with Pool() as pool:#Pool(processes=nproc) as pool:
inds = np.arange(0,24)
with pool:
    res = pool.map(func, inds)

#func(15)
#func(21)
#func(18)
