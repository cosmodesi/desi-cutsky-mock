import os
import sys
import numpy as np
import pandas as pd
from astropy.table import Table


# ---------------------------------------------------------------------------
# Constants / helpers (unchanged from the original)
# ---------------------------------------------------------------------------

BOXSIZE     = 2000.0
NUM_SUB_BOX = 4

LOS = np.array([0., 0., 1.])






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

    cat = Table.read(ifile)
    df = cat.to_pandas()

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

TARG     = 'IC'
PHASES   = range(1, 25)                                  # 22 realizations
N_WORKERS = int(os.environ.get('N_WORKERS', '16'))       # tune me

FILE_TEMPLATES = [
    'atgrid_density_catalog_corners_phase{ph}.fits',
#    'abacus_HF_QSO_0p950_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z0.8_1.1_clustering.dat.h5',
#    'abacus_HF_QSO_1p250_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z1.1_1.4_clustering.dat.h5',
#    'abacus_HF_QSO_1p475_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z0.8_2.1_clustering.dat.h5',
#    'abacus_HF_QSO_1p550_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z1.4_1.7_clustering.dat.h5',
#    'abacus_HF_QSO_1p850_QSO_fit_dv_AbacusSummit_base_c000_ph{ph}_base_lorentzian_dv_cat_fixed_z1.7_2.1_clustering.dat.h5',
]
base = '/global/cfs/cdirs/desi/users/acarnero/abacus_ic'
#BASE = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph{ph}/Boxes/QSO_lorentzian'

tmpl = 'atgrid_density_catalog_corners_phase{ph}.fits'

def build_tasks():
    tasks = []
    for i in PHASES:
        ph = str(i).zfill(3)
        
        #for j, tmpl in enumerate(FILE_TEMPLATES):
        tasks.append((os.path.join(base, tmpl.format(ph=ph)),
                          os.path.join(base, str(i)),
                          TARG))
    return tasks

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

