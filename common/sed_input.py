"""
@verison: 3

@author: David Wilson

@date 20200504

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

"""File structure"""

path = '/media/david/5tb_storage1/muscles/' #data are in the external harddrive
phoenix_path = path + 'phoenix_models/interpolated_models/'
stis_path = path + 'stis_ecsvs/'
cos_path = path + 'cos_savs/'
stis_x1ds = path + 'stis_x1ds/'
cos_x1ds = path + 'cos_x1ds'
chandra_path = path + 'xray/chandra/'
xmm_path = path + 'xray/xmm/'
dem_path = path + 'dem_models'
star_params_path = path + 'mega_muscles_stellar_parameters.csv'
lya_path = path + 'lya_hlsp/' 
euv_repo = path + 'euv_repo/'

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

airmask =  [1207, 1222, 1300, 1310, 1353, 1356]
cos_gratings = ['G130M', 'G160M']
stis_gratings = ['G140M','E140M','G140L', 'G230L', 'G230LB', 'G430L']

lya_ecsvs = glob.glob('{}*ecsv'.format(lya_path))
star_params = Table.read(star_params_path)
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


for star in stars[1:2]:
    repo, component_repo = make_repo(star, path, version)
    
    print(star)
    cos_savs = glob.glob('{}{}/*.sav'.format(cos_path, star))
    #print(cos_savs)
    stis_ecsvs = glob.glob('{}{}/*.ecsv'.format(stis_path, star))

    w_full = np.array([], dtype=float)
    f_full = np.array([], dtype=float)
    e_full = np.array([], dtype=float)
    
    w1 = 1160
    
    #COS
    
    for grating in cos_gratings:
        for sav in cos_savs:
            data = readsav(sav)
            if str(data['grating'][0]) == "b'{}'".format(grating) :
                w, f, e = data['wave'], data['flux'], data['err']
                mask = (w < airmask[0]) | (w > airmask[1]) & (w < airmask[2]) | (w > airmask[3]) & (w < airmask[4]) | (w > airmask[5]) 
                w, f, e = w[mask], f[mask], e[mask]
              #  mask = (w > w1)
            #    w, f, e = w[mask], f[mask], e[mask]
                plt.step(w,f, where='mid', label=grating)
                w1 = w[-1]
                w_full = np.concatenate((w_full, np.array(w)))
                f_full = np.concatenate((f_full, np.array(f)))
                e_full = np.concatenate((e_full, np.array(e)))
    
    #Lya
    
    lw0, lw1 = 1000000, 0 #placeholder so L-980-5 works
    for lya in lya_ecsvs:
        data = Table.read(lya)
      #  print(data.meta['TARGNAME'])
        if data.meta['TARGNAME'] == star:
            w, f = data['WAVELENGTH'], data['FLUX']
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
     
    dem_file = glob.glob('{}{}_v*.fits')
    if len(dem_file) > 1:
        print('More than one dem file for this star')
    if len(dem_file > 0):
        euv_name = 'dem'
        prepare_euv.make_euv(dem_file, dem_path, euv_inputs=euv_inputs)
        prepare_model.make_model_spectrum(, version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, model_name=euv_name)
    
    
                
    plt.legend(loc= 1)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1e-17)
    plt.show()
    args = np.argsort(w_full)
    w_full, f_full, e_full = w_full[args], f_full[args], e_full[args] 
    #savdat = Table([w_full*u.AA, f_full*u.erg/u.s/u.cm**2/u.AA, e_full*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
    #ascii.write(savdat, 'uv_first_pass/{}_hst+opt_v1.ecsv'.format(star), format='ecsv', overwrite=True)



