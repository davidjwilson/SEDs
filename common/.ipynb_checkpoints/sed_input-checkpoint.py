"""
@verison: 4

@author: David Wilson

@date 20200512

Organises all of the Mega-MUSCLES data and produces the SEDs. New 2020 workflow orgainising data by source rather than star.

"""

import instruments
import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from scipy.io.idl import readsav
from astropy.table import Table, vstack
from astropy.io import ascii
import astropy.units as u
import make_mm_sed as sed
import prepare_cos
import prepare_stis
import prepare_model
import prepare_xmm
import prepare_euv
import prepare_chandra
from craftroom import resample
from scipy.interpolate import interp1d
import make_sed_files
from shutil import copyfile

"""File structure"""

path = '/media/david/5tb_storage1/muscles/' #data are in the external harddrive
sources = ['cos','stis', 'lya','phoenix', 'xmm', 'chandra', 'apec', 'euv']


#'2MASS-J23062928-0502285' leaving out Trappist-1
stars = ['L-980-5',
        'GJ674', 
        'GJ676A',
        'GJ649',
        'GJ699',
        'GJ163',
        'GJ849',
        'GJ1132',
        'LHS-2686',
        'GJ729',
        'GJ15A']

stars = ['GJ674']
airglow =  [1207, 1222, 1300, 1310, 1353, 1356]
cos_gratings = ['G130M', 'G160M']
stis_gratings = ['G140M','E140M','G140L', 'G230L', 'G230LB', 'G430L']

#lya_ecsvs = glob.glob('{}*ecsv'.format(lya_path))
#star_params = Table.read(star_params_path)
#print(lya_ecsvs)
###################
version = 1
###################

def make_repo(star, path, version):
    """
    Makes directories to store the produced files
    """
    repo = '{}sed_repo/{}/'.format(path, star)
    component_repo = '{}components_v{}/'.format(repo, version)
    if os.path.exists(repo) == False: #makes the parent directory then puts another directory in it for the components
        os.mkdir(repo)
    if os.path.exists(component_repo) == False:
        os.mkdir(component_repo)
    return repo, component_repo

def sort_components(star, path, sources, component_repo):
    """
    Moves all components to individual star directories
    """
    for source in sources:
        spath = '{}{}_hlsp/'.format(path, source)
        starfiles = glob.glob('{}*{}*'.format(spath, star.lower()))
        for sf in starfiles:
            fname = os.path.split(sf)[1]
            copyfile(sf, '{}{}'.format(component_repo, fname))

      #  print(spec)
              #  filename = os.path.split(spec)[1]
               # copyfile(spec, 'nuv_collection/spectra/{}'.format(filename))

def which_xray(component_repo):
    xscope = 'none'
    if len(glob.glob('{}*xmm*'.format(component_repo))) > 0:
           xscope = 'xmm'
    if len(glob.glob('{}*cxo*'.format(component_repo))) > 0:
           xscope = 'cxo'
    return xscope 

for star in stars:
    print(star)
    repo, component_repo = make_repo(star, path, version)
    sort_components(star, path, sources, component_repo)
    
   #COS
    
    sed_table, instrument_list = sed.build_cos_fuv(component_repo, airglow)
    
    #STIS and Lya
    
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, component_repo, airglow[0:2], instrument_list, airglow[2:], norm=False)
    
    #PHOENIX
    
    sed_table, instrument_list = sed.add_phoenix_and_g430l(sed_table, component_repo, instrument_list, scale=False)
    
    #X-ray
    
    sed_table, instrument_list, euv_gap = sed.add_xray_spectrum(sed_table, component_repo, instrument_list, which_xray(component_repo), add_apec = True, find_gap=True)
    
    #EUV
    euv_name = 'euv-scaling'
    sed_table, instrument_list = sed.add_euv(sed_table, component_repo, instrument_list, euv_gap, euv_name)
    
    sed_table.sort(['WAVELENGTH'])
    lim = np.mean(sed_table['FLUX'][(sed_table['WAVELENGTH'] > 2e5) & (sed_table['WAVELENGTH'] < 3e5)])
 
    savdat = Table([sed_table['WAVELENGTH']*u.AA, sed_table['FLUX']*u.erg/u.s/u.cm**2/u.AA, sed_table['ERROR']*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
    ascii.write(savdat, '{}/basic_seds/{}_basic_v1.ecsv'.format(path, star), format='ecsv', overwrite=True)

    plt.figure(star, figsize=(7, 5))
    plt.plot(sed_table['WAVELENGTH'], sed_table['FLUX'], label=star, rasterized=True)
   # plt.plot(sed_table['WAVELENGTH'], sed_table['ERROR'], rasterized=True)
    plt.ylim(lim)
    plt.xlim(5, 3e5)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Wavelength (\AA)')
    plt.ylabel('Flux erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)')
    plt.legend(loc=1)
    plt.tight_layout()
    #plt.savefig('plots/first_seds/{}_v{}_sed.png'.format(star, version))
    plt.show()
    plt.close()
    
    

    

 
print('Done')

"""    

    #savdat = Table([w_full*u.AA, f_full*u.erg/u.s/u.cm**2/u.AA, e_full*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
    #ascii.write(savdat, 'uv_first_pass/{}_hst+opt_v1.ecsv'.format(star), format='ecsv', overwrite=True)

"""

