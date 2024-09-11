import numpy as np
import os
from astropy.table import Table

opath = '/pscratch/sd/a/acarnero/codes/generate_survey_mocks/third/subboxes/LRG/z0.800'
ipath = '/pscratch/sd/h/hanyuz/Y1HOD/mocks_y1_wponly/AbacusSummit_base_c000_ph000/z0.800/galaxies/LRGs.dat'


lrgs = Table.read(ipath, format='ascii.ecsv')
lrgs['x'] += 1000
lrgs['y'] += 1000
lrgs['z'] += 1000

df = lrgs.to_pandas()

# Number of subdivisions per dimension
num_sub_boxes = 4

# Create bins for x, y, z coordinates
bins = np.linspace(0, 2000, num_sub_boxes + 1)

# Digitize the coordinates to determine which sub-box each point falls into
df['x_bin'] = np.digitize(df['x'], bins) - 1
df['y_bin'] = np.digitize(df['y'], bins) - 1
df['z_bin'] = np.digitize(df['z'], bins) - 1

df['sub_box'] = df['x_bin'] * (num_sub_boxes**2) + df['y_bin'] * num_sub_boxes + df['z_bin']
sub_boxes = [group for _, group in df.groupby('sub_box')]

for j,subdf in enumerate(sub_boxes):

    subdf = subdf.drop(columns=['x_bin', 'y_bin', 'z_bin', 'sub_box'])
    t2 = Table.from_pandas(subdf)
    print('writing ', j)
    t2.write(os.path.join(opath, 'targs_SB%d.fits' % j))


