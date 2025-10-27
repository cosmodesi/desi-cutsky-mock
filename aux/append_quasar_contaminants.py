import numpy as np
import os
import sys
from astropy.io import fits

import fitsio
from astropy.table import Table,unique,join,vstack
import LSS.common_tools as common
from LSS.globals import main

# Script to append quasar contaminants to the high-fidelity mocks
# Starts by identifying contaminants in the full file of the data
# and then adds them to the mock

append_stars = True
append_unclassified = True
append_galaxies = False


basedir = '/global/cfs/cdirs/desi/survey/catalogs/'
survey = 'DA2'
data = 'LSS'
verspec = 'loa-v1'
version = 'v1.1'
mock_path = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400/forclustering/cutsky_abacusHF_DR2_QSO_z1p400_zcut_0p8to3p5_clustering.dat.fits'
path_out = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400/forclustering/'
name_out = 'cutsky_abacusHF_DR2_QSO_z1p400_zcut_0p8to3p5_clustering_append_stars_unclassified.dat.fits'
tp = 'QSO'

indir = os.path.join(basedir, survey, data, verspec, 'LSScats', version)


full = Table(fitsio.read(os.path.join(indir, tp + '_full_HPmapcut.dat.fits')))
## the selection of valid samples
sel_obs = full['ZWARN'] != 999999
sel_obs&= full['ZWARN']*0 == 0

selection = sel_obs

gz = common.goodz_infull(tp[:3], full)
selection_gz = selection & gz

emline = fits.open(os.path.join(basedir, survey, data, verspec, 'emlin_catalog.fits'))[1].data

ss = np.searchsorted(sorted(emline['TARGETID']), full['TARGETID'])

argss = np.argsort(emline['TARGETID'])

oii_flux = emline['OII_FLUX'][argss][ss]
oii_flux_ivar = emline['OII_FLUX_IVAR'][argss][ss]

oiii_flux = emline['OIII_FLUX'][argss][ss]
oiii_flux_ivar = emline['OIII_FLUX_IVAR'][argss][ss]

dchi_cut = 30
o2c_cut = 0.9
oiii_cut = 5
selgal = ((selection  & ~gz & (full['Z_RR'] > 0.01)  & (full['Z_RR'] < 1.625)) & 
          ((full['DELTACHI2'] > dchi_cut) | (np.log10(oii_flux * oii_flux_ivar**0.5) > o2c_cut - 0.2 * np.log10(full['DELTACHI2']))
           | (oiii_flux * oiii_flux_ivar**0.5 > oiii_cut))
          )  #...now 40,1.2,5...earlier 30,0.9,5

selstar = (selection & ~gz & (full['Z_RR'] < 0.01))


print('Fraction of stars',np.sum(selstar) / np.sum(selection))
print('Fraction of galaxies',np.sum(selgal) / np.sum(selection))
print('Fraction of qsos',np.sum(selection&gz)/np.sum(selection))
print('Fraction of junk',1-( np.sum(selstar) + np.sum(selgal) + np.sum(selection&gz) )/np.sum(selection))
flavour = 'abacus'
sim_data = Table.read(mock_path)#fits.open(mock_path)[1].data

if append_stars:
    stars = full[selstar]
    if flavour == 'abacus':
        '''
            name = 'RA'; format = 'E',             name = 'DEC'; format = 'E',           name = 'TRUEZ'; format = 'E',
            name = 'STATUS'; format = 'J',            name = 'RAW_NZ'; format = 'E',            name = 'RAN_NUM_0_1'; format = 'E',
            name = 'NZ'; format = 'E',             name = 'Z'; format = 'E',             name = 'HALO_ID'; format = 'K',
            name = 'HALO_MASS'; format = 'E',             name = 'IS_CENTRAL'; format = 'L',             name = 'Z_ERR_3GAUSS'; format = 'E',
            name = 'Z_ERR_SIG500'; format = 'E',             name = 'WEIGHT'; format = 'D',             name = 'DESI_TARGET'; format = 'K',
            name = 'PRIORITY_INIT'; format = 'K',             name = 'PRIORITY'; format = 'K',             name = 'NUMOBS_MORE'; format = 'K',
            name = 'NUMOBS_INIT'; format = 'K',             name = 'BGS_TARGET'; format = 'K',             name = 'TARGETID'; format = 'K',
            name = 'MWS_TARGET'; format = 'K',             name = 'SUBPRIORITY'; format = 'D',             name = 'BRICKNAME'; format = '8A',
            name = 'OBSCONDITIONS'; format = 'K',            name = 'SCND_TARGET'; format = 'K',            name = 'ZWARN'; format = 'K',             )
        '''

        star_arr = np.empty(len(stars), dtype=sim_data.dtype)
        star_arr['RA'] = stars['RA'].astype(sim_data.dtype['RA'])
        star_arr['DEC'] = stars['DEC'].astype(sim_data.dtype['DEC'])
        
        star_out = Table(star_arr)
    elif flavour == 'uchuu':
        star_out = Table(np.array([-9 * np.ones(len(stars)).astype('int'), -9 * np.ones(len(stars)).astype('int'), stars['RA'], stars['DEC'], np.zeros_like(stars['RA']), np.zeros_like(stars['DEC']),
                                   np.array(['STAR'] * len(stars)), np.zeros_like(stars['RA']), np.zeros_like(stars['RA'])]).T, dtype = sim_data.dtype) #names=('GALAXYID',
        #'PID','RA','DEC','Z_NORSD','Z','TRACER_TYPE',
        #'NX','WEIGHT_FKP'))
    else:
        raise Exception('should be uchuu or abacus')
    sim_data = vstack((Table(sim_data), star_out))
print(len(sim_data), len(star_out))
if append_unclassified:
    unclassified = full[selection & ~selstar & ~selgal & ~(selection &gz)]
    if flavour == 'abacus':
        unclassified_arr = np.empty(len(unclassified), dtype=sim_data.dtype)
        unclassified_arr['RA'] = unclassified['RA'].astype(sim_data.dtype['RA'])
        unclassified_arr['DEC'] = unclassified['DEC'].astype(sim_data.dtype['DEC'])
        unclassified_out = Table(unclassified_arr)
    elif flavour == 'uchuu':
        unclassified_out = Table(np.array([-9 * np.ones(len(unclassified)).astype('int'),
                                       -9 * np.ones(len(unclassified)).astype('int'),
                                       unclassified['RA'],
                                       unclassified['DEC'],
                                       np.zeros_like(unclassified['RA']),
                                       np.zeros_like(unclassified['DEC']),
                                       np.array(['UNKNOWN'] * len(unclassified)),
                                       np.zeros_like(unclassified['RA']),
                                       np.zeros_like(unclassified['RA'])]).T,dtype=sim_data.dtype) #names=('GALAXYID',
    #'PID','RA','DEC','Z_NORSD','Z','TRACER_TYPE',
        #'NX','WEIGHT_FKP'))
    else:
        raise Exception('should be uchuu or abacus')

    sim_data = vstack((Table(sim_data), unclassified_out))
print(len(sim_data), len(unclassified_out))

if append_galaxies:
    galaxies = full[selgal]
    galaxies_out = Table(np.array([-9 * np.ones(len(galaxies)).astype('int'),
                                   -9 * np.ones(len(galaxies)).astype('int'),
                                   galaxies['RA'],
                                   galaxies['DEC'],
                                   np.zeros_like(galaxies['RA']),
                                   np.zeros_like(galaxies['DEC']),
                                   np.array(['GALAXY'] * len(galaxies)),
                                   np.zeros_like(galaxies['RA']),
                                   np.zeros_like(galaxies['RA'])]).T,dtype=sim_data.dtype) #names=('GALAXYID',
    #'PID','RA','DEC','Z_NORSD','Z','TRACER_TYPE',
        #'NX','WEIGHT_FKP'))
    sim_data = vstack((Table(sim_data), galaxies_out))

sim_data.write(os.path.join(path_out, name_out))
