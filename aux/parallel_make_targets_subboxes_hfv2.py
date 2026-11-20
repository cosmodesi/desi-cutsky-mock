import os
import sys
import numpy as np
import pandas as pd
from astropy.table import Table
from mockfactory import Catalog


# ---------------------------------------------------------------------------
# Constants / helpers (unchanged from the original)
# ---------------------------------------------------------------------------

BOXSIZE     = 2000.0
SHIFT       = 1000.0
NUM_SUB_BOX = 4
COLUMNS     = ['HALO_ID', 'ISCENTRAL', 'MASS', 'VSMEAR',
               'VX', 'VY', 'VZ', 'X', 'Y', 'Z']#, 'RSDZ']

LOS = np.array([0., 0., 1.])


def apply_periodic(x, L):
    return (x + 0.5 * L) % L


def catalog_to_rsd(cat, scalev=None):
    
    if scalev is None:
        raise ValueError("VELZ2KMS not found in catalog.header")

    los = LOS / np.linalg.norm(LOS)

    pos = np.column_stack([cat['X'],  cat['Y'],  cat['Z']])
    vel = np.column_stack([cat['VX'], cat['VY'], cat['VZ']])

    v_los  = vel @ los + cat['VSMEAR']
    shift  = np.array(v_los / scalev)[:, None] * los
    pos_rsd = pos + shift

    L = BOXSIZE
    
    # Correct wrapping for [0, L) convention:
    return pos_rsd % L

# ---------------------------------------------------------------------------
# Per-file worker — same logic as your `func`, just tidied up
# ---------------------------------------------------------------------------

def func(ifile, opath, targ):
    print(f"[start] {targ}: {ifile}", flush=True)
    os.makedirs(opath, exist_ok=True)

    marker = os.path.join(opath, f'{targ}_real_space.sub0.fits.gz')
    if os.path.isfile(marker):
        print(f"[skip] {marker}", flush=True)
        return 0

    cat = Catalog.read(ifile)
    scalev = cat.header.get('VELZ2KMS', None)
    print('scalev', scalev)
    print('transforming to df')
    df = pd.DataFrame({c: cat[c] for c in COLUMNS})
    print('applying boundary conditions and putting them in [0-2000]')
    
    for c in ('X', 'Y', 'Z'):
        df[c] = apply_periodic(df[c], BOXSIZE)
    print('start RSD')

    rsd = catalog_to_rsd(df, scalev=scalev)
    
    df['RSDZ'] = rsd[:, 2]

    print(targ, len(df), flush=True)

    bins = np.linspace(0, BOXSIZE, NUM_SUB_BOX + 1)
    b = NUM_SUB_BOX
    df['sub_box'] = (
        (np.digitize(df['X'], bins) - 1) * (b * b)
        + (np.digitize(df['Y'], bins) - 1) * b
        + (np.digitize(df['Z'], bins) - 1)
    )
    print("sub_boxes:", sorted(set(df['sub_box'])), flush=True)

    # sort=True (default) -> same file numbering as your original code,
    # which enumerated sorted groups and wrote sub0.fits.gz, sub1.fits.gz, ...
    for j, (_, subdf) in enumerate(df.groupby('sub_box')):
        subdf = subdf.drop(columns=['sub_box'])
        out_path = os.path.join(opath, f'{targ}_real_space.sub{j}.fits.gz')
        Table.from_pandas(subdf).write(out_path, overwrite=True)

    print(f"[done] {ifile}", flush=True)
    return 0


# ---------------------------------------------------------------------------
# Build the 110-task flat list, then dispatch to a process pool
# ---------------------------------------------------------------------------

TARG     = 'QSO'
PHASES   = range(25)                                  # 22 realizations
N_WORKERS = int(os.environ.get('N_WORKERS', '16'))       # tune me

FILE_TEMPLATES = [
    'abacus_HF_LRG_0p500_DR2_v2.0_AbacusSummit_base_c000_ph{ph}_base_B_clustering.dat.h5',
#    'abacus_HF_QSO_0p950_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z0.8_1.1_clustering.dat.h5',
#    'abacus_HF_QSO_1p250_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z1.1_1.4_clustering.dat.h5',
#    'abacus_HF_QSO_1p475_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z0.8_2.1_clustering.dat.h5',
#    'abacus_HF_QSO_1p550_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z1.4_1.7_clustering.dat.h5',
#    'abacus_HF_QSO_1p850_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z1.7_2.1_clustering.dat.h5',
]
#BASE = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph{ph}/Boxes/LRG'
BASE = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph{ph}/Boxes/QSO_lorentzian'

tmpl = 'abacus_HF_QSO_1p475_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z0.8_2.1_clustering.dat.h5'

def build_tasks_old():
    tasks = []
    for i in PHASES:
        ph = str(i).zfill(3)
        base = BASE.format(ph=ph)
        for j, tmpl in enumerate(FILE_TEMPLATES):
            tasks.append((os.path.join(base, tmpl.format(ph=ph)),
                          os.path.join(base, str(j)),
                          TARG))
    return tasks

def build_tasks():
    tasks = []
    for i in PHASES:
        ph = str(i).zfill(3)
        base = BASE.format(ph=ph)
        #for j, tmpl in enumerate(FILE_TEMPLATES):
        tasks.append((os.path.join(base, tmpl.format(ph=ph)),
                          os.path.join(base, '0_z0.8_2.1'),
                          TARG))
    return tasks


def _worker(task):
    """Top-level wrapper so the call target is picklable on `spawn`."""
    return func(*task)


def main():
    tasks = build_tasks()
    print(f"Total tasks: {len(tasks)} | workers: {N_WORKERS}", flush=True)

    # Pre-create dirs once on the main process to avoid a thundering herd
    # of mkdir() calls hitting the Lustre/CFS metadata server.
    for _, opath, _ in tasks:
        os.makedirs(opath, exist_ok=True)

    from concurrent.futures import ProcessPoolExecutor, as_completed

    n_ok = n_err = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(_worker, t): t for t in tasks}
        for fut in as_completed(futures):
            ifile, opath, targ = futures[fut]
            try:
                fut.result()
                n_ok += 1
                print(f"[{n_ok + n_err}/{len(tasks)}] ok  {ifile}", flush=True)
            except Exception as e:
                n_err += 1
                print(f"[{n_ok + n_err}/{len(tasks)}] FAIL {ifile}: {e}",
                      flush=True)

    print(f"All done. {n_ok} ok, {n_err} failed.", flush=True)


if __name__ == '__main__':
    main()

