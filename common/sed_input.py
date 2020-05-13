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
        if source == 'stis':
            spath = '{}{}_hlsp/{}/'.format(path, source, star)
        else:
            spath = '{}{}_hlsp/'.format(path, source)
        starfiles = glob.glob('{}*{}*'.format(spath, star.lower()))
        for sf in starfiles:
            fname = os.path.split(sf)[1]
            copyfile(sf, '{}{}'.format(component_repo, fname))

      #  print(spec)
              #  filename = os.path.split(spec)[1]
               # copyfile(spec, 'nuv_collection/spectra/{}'.format(filename))



for star in stars[0:2]:
    repo, component_repo = make_repo(star, path, version)
   # sort_components(star, path, sources, component_repo)
    
   #COS
    
    sed_table, instrument_list = sed.build_cos_fuv(component_repo, airglow)
    
    #STIS and Lya
    
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, component_repo, airglow[0:2], instrument_list, airglow[2:])
    
    #PHOENIX
    
    sed_table, instrument_list = sed.add_phoenix_and_g430l(sed_table, component_repo, instrument_list)
    
    
    args = np.argsort(sed_table['WAVELENGTH'])
    plt.plot(sed_table['WAVELENGTH'][args], sed_table['FLUX'][args])
    plt.show()
    
"""    
    
    cos_savs = glob.glob('{}{}/*.sav'.format(cos_path, star))
    #print(cos_savs)
    stis_ecsvs = glob.glob('{}{}/*.ecsv'.format(stis_path, star))

    w_full = np.array([], dtype=float)
    f_full = np.array([], dtype=float)
    e_full = np.array([], dtype=float)
    
    w1 = 1160
    
    #COS
    
    sed_table, instrument_list = sed.build_cos_fuv(star, '{}cos_hlsp'.format(path), airglow)
    
    #Lya + Stis
    
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, component_repo, airglow[0:2], instrument_list, airglow[2:])
    
    lw0, lw1 = 1000000, 0 #placeholder so L-980-5 works
    for lya in lya_ecsvs:
        data = Table.read(lya)
      #  print(data.meta['TARGNAME'])
        if data.meta['TARGNAME'] == star:
            w, f = data['WAVELENGTH'], data['FLUX']hlsp_muscles_hst_cos_gj1132_g130m_v1_component-spec.fits
            plt.step(w, f, where='mid', label=r'Ly $\alpha$')
            lw0, lw1 = w[0], w[-1]
            w_full = np.concatenate((w_full, np.array(w)))
            f_full = np.concatenate((f_full, np.array(f)))
            e_full = np.concatenate((e_full, np.zeros(len(w))))
            
    #STIS
            
    for grating in stis_gratings:
        for spec in stis_ecsvs:
            data= Table.read(spec)
            if data.meta['GRATING'] == grating:
                w, f, e = data['WAVELENGTH'], data['FLUX'], data['ERROR']
                if grating == 'G140M':
                    mask = (w > 1207) & (w < lw0) | (w > lw1) & (w < 1222) 
               # elif grating == 'G140L':
                #    mask = (w > 1300) & (w < 1310) | (w > 1353) & (w < 1356) #| (w > w1)
                 #   w1 = w[-1]
                #elif grating == 'E140M':
                 #   mask =  (w > w1) & (w < lw0) #| (w > lw1)
                  #  w1 = w[-1]
                else:
                  #  mask = (w > w1)
                    mask = (w > w[0]-10)
                    w1 = w[-1]
                w, f, e = w[mask], f[mask], e[mask]
                smooth = 2
               # if grating =='G430L':
                #    f = convolve(np.array(f),Box1DKernel(smooth))
                 #   e = convolve(np.array(e),Box1DKernel(smooth))/(smooth**0.5)
                plt.step(w, f, where='mid', label=grating)
                w_full = np.concatenate((w_full, np.array(w)))
                f_full = np.concatenate((f_full, np.array(f)))
                e_full = np.concatenate((e_full, np.array(e)))
                
    #PHOENIX        
    
    opath = glob.glob(phoenix_path+star+'*')
 #   print(opath)
    pdata = Table.read(opath[0])
    w, f, e  = pdata['WAVELENGTH'], pdata['FLUX']*pdata.meta['NORMFAC'], np.zeros(len(pdata['WAVELENGTH']))
    mask = w > max(w_full)
    w, f, e = w[mask], f[mask], e[mask]
    plt.plot(w, f, ls='--')
    w_full = np.concatenate((w_full, np.array(w)))
    f_full = np.concatenate((f_full, np.array(f)))
    e_full = np.concatenate((e_full, np.array(e)))
    
    #EUV
     
    #dem_file = glob.glob('{}{}_v*.fits')
   # if len(dem_file) > 1:
     #   print('More than one dem file for this star')
    #if len(dem_file > 0):
     #   euv_name = 'dem'
      #  prepare_euv.make_euv(dem_file, dem_path, euv_inputs=euv_inputs)
       # prepare_model.make_model_spectrum(, version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, model_name=euv_name)
  #  if len(dem_file) = 0:
    euv_files = glob.glob('{}hlsp_muscles_model_euv-scaling_{}_na_v1_component-spec.ecsv'.format(euv_repo, star.lower()))
    if len(euv_files) == 1:
        data = Table.read(euv_files[0])
        w, f = data['WAVELENGTH'], data['FLUX']
        mask = w < min(w_full)
        w, f = w[mask], f[mask]
        plt.step(w, f, where='mid', label='EUV scaling')
        w_full = np.concatenate((w_full, np.array(w)))
        f_full = np.concatenate((f_full, np.array(f)))
        e_full = np.concatenate((e_full, np.zeros(len(w))))
    
                
    plt.legend(loc= 1)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1e-17)
    plt.show()
    args = np.argsort(w_full)
    w_full, f_full, e_full = w_full[args], f_full[args], e_full[args] 
    #savdat = Table([w_full*u.AA, f_full*u.erg/u.s/u.cm**2/u.AA, e_full*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
    #ascii.write(savdat, 'uv_first_pass/{}_hst+opt_v1.ecsv'.format(star), format='ecsv', overwrite=True)

"""

