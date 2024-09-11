import numpy as np
import asdf
import os
from astropy.table import Table
import errno

nmesh = 576
Lbox = 2000. # Mpc/h
cell_size = Lbox/nmesh
bin_edges = np.arange(nmesh+1)*cell_size
bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])-(cell_size/2.)

X,Y,Z = np.meshgrid(bin_centers, bin_centers, bin_centers, indexing='xy')
#print(X)
#print(X.shape)
#print(Y.shape)
#print(Z.shape)
def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
#            print('made %s'%value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


def calculate_unique_id(x, y, z, X, Y):
    """
    Calculate the unique identifier for a node at position (x, y, z)
    in a grid of size XxYxZ.
    """
    return x + y * X + z * (X * Y)

def subDivistion(matr):
    # Sub-box size
    sub_box_size = 144

    # Number of sub-boxes in each dimension
    num_sub_boxes = 4

    # Initialize a list to store the sub-boxes
    sub_boxes = []

    # Iterate over the grid to create sub-boxes
    for i in range(num_sub_boxes):
        for j in range(num_sub_boxes):
            for k in range(num_sub_boxes):
            # Extract the sub-box
                sub_box = matr[
                    i * sub_box_size: (i + 1) * sub_box_size,
                    j * sub_box_size: (j + 1) * sub_box_size,
                    k * sub_box_size: (k + 1) * sub_box_size
                ]
            # Store the sub-box
                sub_boxes.append(sub_box)
    return sub_boxes
X_SB = subDivistion(X)
Y_SB = subDivistion(Y)
Z_SB = subDivistion(Z)


Xs = 576

# Create an array to store the identifiers
identifiers = np.zeros((Xs, Xs, Xs), dtype=np.int64)

# Fill the array with unique identifiers
for x in range(Xs):
    for y in range(Xs):
        for z in range(Xs):
            identifiers[x, y, z] = calculate_unique_id(x, y, z, Xs, Xs)


ID_SB = subDivistion(identifiers)

for i in range(9, 25):
    phase = str(int(i)).zfill(3)
    path = '/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/ic/AbacusSummit_base_c000_ph{PHASE}'.format(PHASE=phase)

    with asdf.open(os.path.join(path, 'ic_dens_N576.asdf'), lazy_load=False) as af:
        dens = af['data']['density']
    print(phase, dens.shape)
    

# Now, `sub_boxes` contains 64 sub-boxes of shape (144, 144, 144)
#print(POS.shape)

    dens_SB = subDivistion(dens)

    pathtosave = '/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/CubicBox/ic/AbacusSummit_base_c000_ph{PHASE}'.format(PHASE=phase)
    test_dir(pathtosave)
    for j in range(64):
        D = dens_SB[j].flatten()
        XX = X_SB[j].flatten()
        YY = Y_SB[j].flatten()
        ZZ = Z_SB[j].flatten()
        IID = ID_SB[j].flatten()
        data_ = {'x':XX, 'y':YY, 'z':ZZ, 'density':D, 'ID_mesh':IID} 

        TT = Table(data_)
        TT.write(os.path.join(pathtosave, 'ic_real_space.sub%d.fits.gz' % j), overwrite=True)
    print(phase, 'done')
#print('subboxes 0 0')
#print('dens', dens_SB[0][0], 'x', X_SB[0][0], 'y', Y_SB[0][0], 'z', Z_SB[0][0])
#print('subboxes 0 -1')

#print(dens_SB[0][-1], X_SB[0][-1], Y_SB[0][-1], Z_SB[0][-1])


