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
import make_fits
import instruments

"""File structure"""

# path = '/media/david/5tb_storage1/muscles/' #data are in the external harddrive
# path = '/media/david/1tb_storage1/emergency_data/mega_muscles/' #backup hd
path = '/media/david/2tb_ext_hd/hddata/mega_muscles/' #new hd
sources = ['cos','stis', 'lya','phoenix', 'xmm', 'chandra', 'apec', 'euv']


# stars= []# leaving out Trappist-1
# stars = ['2MASS-J23062928-0502285',
#         'L-980-5',
#         'GJ674', 
#         'GJ676A',
#         'GJ649',
#         'GJ699',
#         'GJ163',
#         'GJ849',
#         'GJ1132',
#         'LHS-2686',
#         'GJ729',
#         'GJ15A']
# stars = ['L-980-5']
stars = ['GJ15A']
# stars = ['GJ699']
# stars = ['L-980-5']#'GJ676A']
# stars = ['GJ1132']
stars = ['GJ15A', 'GJ163', 'GJ699', 'GJ849', 'LHS-2686']

airglow =  [1207, 1222, 1300, 1310, 1353, 1356]
cos_gratings = ['G130M', 'G160M']
stis_gratings = ['G140M','E140M','G140L', 'G230L', 'G230LB', 'G430L']

#lya_ecsvs = glob.glob('{}*ecsv'.format(lya_path))

#print(lya_ecsvs)
###################
version = 2
###################

def make_repo(star, path, version):
    """
    Makes directories to store the produced files
    """
    repo = '{}hlsp/{}/'.format(path, star)
    component_repo = '{}components_v1/'.format(repo)
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

def table_name(star): 
    """
    Turns the input name into the name in the stellar params table
    """
    if star[:2] == 'GJ':
        star = star.replace(' ', '')
    elif star == 'Trappist-1':
        star = '2MASS-J23062928-0502285'
    elif star[0] == 'L':
        star = star.replace(' ', '-')
    return star
    
star_params_path = '/home/david/work/muscles/SEDs/optical/stellar_parameters.csv'
star_params = Table.read(star_params_path)
targets = np.array([table_name(star) for star in star_params['Target']])
# print(targets)
    
for star in stars:
    print(star)
    repo, component_repo = make_repo(star, path, version)
    print(component_repo)
#     sort_components(star, path, sources, component_repo)

    row = star_params[np.where(targets == star)[0][0]]
    
   #COS
    
    sed_table, instrument_list = sed.add_cos(component_repo, airglow, remove_negs=True)
    
    #STIS and Lya
    
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, component_repo, airglow[0:2], instrument_list, airglow[2:], norm=False, remove_negs=True)
    
    #PHOENIX
    

    
    sed_table, instrument_list = sed.add_phoenix_and_g430l(sed_table, component_repo, instrument_list, scale=False)
#     sed_table, instrument_list= sed.add_phx_spectrum(sed_table, component_repo, instrument_list)
    
    #X-ray
    
    sed_table, instrument_list, euv_gap = sed.add_xray_spectrum(sed_table, component_repo, instrument_list, which_xray(component_repo), add_apec = True, find_gap=True)
    
    #EUV
    euv_name = 'euv-scaling'
    if len(glob.glob('{}*dem*'.format(component_repo))) > 1:
        euv_name = 'dem'
    
    sed_table, instrument_list = sed.add_euv(sed_table, component_repo, instrument_list, euv_gap, euv_name)
    
    sed_table.sort(['WAVELENGTH'])
#     print(sed_table.meta)
    #bolometric flux
    sed_table = sed.add_bolometric_flux(sed_table, component_repo, row)
    
    print(sed_table['BOLOFLUX'][100])
    
    #final meta keys
    sed_table.meta['WAVEMIN'] = min(sed_table['WAVELENGTH'])
    sed_table.meta['WAVEMAX'] = max(sed_table['WAVELENGTH'])
    sed_table.meta['FLUXMIN'] = min(sed_table['FLUX'])
    sed_table.meta['FLUXMAX'] = max(sed_table['FLUX'])

    
#     np.save('test_to_fits/ti_instlist', instrument_list)
#     sed_table.write('test_to_fits/t1_table_test.ecsv', overwrite=True)
    make_fits.make_mm_fits(component_repo, sed_table, instrument_list, version,sed_type='adapt-var')
    
#     sed_table_1A = sed.sed_to_const_res(sed_table)
#     sed_table_1A.meta['WAVEMIN'] = min(sed_table_1A['WAVELENGTH'])
#     sed_table_1A.meta['WAVEMAX'] = max(sed_table_1A['WAVELENGTH'])
#     sed_table_1A.meta['FLUXMIN'] = min(sed_table_1A['FLUX'])
#     sed_table_1A.meta['FLUXMAX'] = max(sed_table_1A['FLUX'])
#     print(min(sed_table_1A['WAVELENGTH']), max(sed_table_1A['WAVELENGTH']),min(sed_table_1A['FLUX']), max(sed_table_1A['FLUX']))
#     make_fits.make_mm_fits(component_repo, sed_table_1A, instrument_list, version,sed_type='const')
    
#     print (sed_table.dtype.names)
    
#     print(sed_table[1300])
#     print(sed_table.meta)
#     print(instrument_list)
 
#     savdat = Table([sed_table['WAVELENGTH']*u.AA, sed_table['FLUX']*u.erg/u.s/u.cm**2/u.AA, sed_table['ERROR']*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
#     ascii.write(savdat, '{}/basic_seds/{}_basic_v1.ecsv'.format(path, star), format='ecsv', overwrite=True)

    plt.figure(star, figsize=(7, 5))
    plt.plot(sed_table['WAVELENGTH'], sed_table['FLUX'], label=star)
    plt.plot(sed_table['WAVELENGTH'], sed_table['ERROR'])
#     plt.plot(sed_table_1A['WAVELENGTH'], sed_table_1A['FLUX'], label=star)
#     plt.plot(sed_table_1A['WAVELENGTH'], sed_table_1A['ERROR'])

#     plt.ylim(lim)
#     plt.ylim(1e-17, 1e-13)
#     plt.xlim(5, 3e5)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Wavelength (\AA)')
    plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)')
#     plt.legend(loc=1)
    plt.tight_layout()
    #plt.savefig('plots/first_seds/{}_v{}_sed.png'.format(star, version))
    plt.show()
    plt.close()
    
    

    

print('Done')

"""    

    #savdat = Table([w_full*u.AA, f_full*u.erg/u.s/u.cm**2/u.AA, e_full*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
    #ascii.write(savdat, 'uv_first_pass/{}_hst+opt_v1.ecsv'.format(star), format='ecsv', overwrite=True)

"""

