import numpy as np
import os
from astropy.table import Table
from desitarget.internal import sharedmem
import h5py
import pandas as pd
import numpy as np
from mockfactory import Catalog


def apply_periodic(x, L):
    return (x + 0.5 * L) % L #- 0.5 * L



def catalog_to_rsd(cat, los=np.array([0., 0., 1.])):
    boxsize = cat.header.get('BOXSIZE', 2000.0)
    scalev  = cat.header.get('VELZ2KMS', None)
    if scalev is None:
        raise ValueError("VELZ2KMS not found in catalog.header")

    los = np.asarray(los, dtype=float)
    los /= np.linalg.norm(los)

    pos = np.column_stack([cat['X'],  cat['Y'],  cat['Z']])
    vel = np.column_stack([cat['VX'], cat['VY'], cat['VZ']])

    v_los  = vel @ los + cat['VSMEAR']
    shift  = (v_los / scalev)[:, None] * los
    pos_rsd = pos + shift

    L = boxsize
    pos_rsd = (pos_rsd + L / 2.0) % L - L / 2.0

    return pos_rsd


def func(ifile, opath, targ):
        print(targ)
        if not os.path.isdir(opath):
            os.mkdir(opath)
        if os.path.isfile(os.path.join(opath, f'{targ}_real_space.sub0.fits.gz')):
            return 0
        cat = Catalog.read(ifile)

        test_cat = catalog_to_rsd(cat)
            
        x, y, rsdz = test_cat.T

        cat['RSDZ'] = rsdz
        columns=['HALO_ID', 'ISCENTRAL', 'MASS', 'VSMEAR', 'VX', 'VY', 'VZ', 'X', 'Y', 'Z', 'RSDZ']
        df = pd.DataFrame()
        for col in columns:
            df[col] = cat[col]

        #df = pd.DataFrame(cat)
        
        print(targ, len(df), df.columns)
        df['X'] += 1000
        df['Y'] += 1000
        #testy = apply_periodic(cat['Y'], 2000)
        #print(np.min(testy), np.max(testy))
        #exit()
        df['Z'] += 1000 #apply_periodic(df['z'], 2000) + 1000
        df['RSDZ'] += 1000

        df['X'] = apply_periodic(df['X'],2000)
        df['Y'] = apply_periodic(df['Y'],2000)
        df['Z'] = apply_periodic(df['Z'],2000)
        df['RSDZ'] = apply_periodic(df['RSDZ'],2000)


        print(np.min(df['X']), np.max(df['X']))
        print(np.min(df['Y']), np.max(df['Y']))
        print(np.min(df['Z']), np.max(df['Z']))
        print(np.min(df['RSDZ']), np.max(df['RSDZ']))
        
       
        num_sub_boxes = 4

        # Create bins for x, y, z coordinates
        bins = np.linspace(0, 2000, num_sub_boxes + 1)

        # Digitize the coordinates to determine which sub-box each point falls into
        df['x_bin'] = np.digitize(df['X'], bins) - 1
        df['y_bin'] = np.digitize(df['Y'], bins) - 1
        df['z_bin'] = np.digitize(df['Z'], bins) - 1



        print(set(df['x_bin']), set(df['y_bin']),set(df['z_bin']))


        df['sub_box'] = df['x_bin'] * (num_sub_boxes**2) + df['y_bin'] * num_sub_boxes + df['z_bin']

        

        sub_boxes = [group for _, group in df.groupby('sub_box')]
        for j,subdf in enumerate(sub_boxes):

            subdf = subdf.drop(columns=['x_bin', 'y_bin', 'z_bin', 'sub_box'])
            t2 = Table.from_pandas(subdf)
            t2.write(os.path.join(opath, f'{targ}_real_space.sub%d.fits.gz' % j), overwrite=True)



targ = 'QSO'

#filetosave = open('lrg_list.txt','w')
print('TARG', targ)
for i in range(12, 13):
    print('realization', i)
#    for sn in ['1p850']:  #QSO
    for sn in ['0p950', '1p250', '1p550', '1p850']:  #QSO
#    for sn in ['1p175']:  #ELG
    #for sn in ['0p950', '1p175', '1p475']:  #ELG
    #for sn in ['0p500', '0p725', '0p950']: #LRG
#        print('snapshot', sn)
        #for typ in ['base_B','base']: #LRG
#            print('type fit', typ)
#        for typ in ['base_B','base_B_dv','base','base_dv']: #LRG
        for typ in ['base']: #QSO
            #print('type fit', typ)
        #for typ in ['base_conf_nfwexp']: #ELG
            filename = f'abacus_HF_{targ}_{sn}_DR2_v2.0_AbacusSummit_base_c000_ph{str(int(i)).zfill(3)}_{typ}_clustering.dat.h5'
            input_path = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph{str(int(i)).zfill(3)}/Boxes/{targ}/{filename}'
            output_path = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph{str(int(i)).zfill(3)}/Boxes/{targ}/sn{sn}/{typ}'
            os.makedirs(output_path, exist_ok=True)
#            filetosave.write('%s %s %s\n' %(input_path, output_path, targ))
            print(input_path, output_path)
            func(input_path, output_path, targ)
#filetosave.close()
