"""
@verison: 1

@author: David Wilson

@date 20190805

Using this to wite the main function for make_mm_sed.
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
import prepare_lya

#paths = some way of storing all the paths to the different spectra 
"""
input paths = dict('COS_readsav', 'COS_x1d', 'STIS_FUV', lya_model)
airglow (lya)
"""


def make_sed(input_paths, savepath, version, lya_range, other_airglow, save_components = False):
    airglow = lya_range+other_airglow
    #COS FUV 
    component_repo = savepath+'components/' #directory where the component spectra are saved
    prepare_cos.make_cos_spectrum(input_paths['COS_readsav'], version,  input_paths['COS_x1d'], savepath = component_repo, plot=False, save_ecsv=save_components, save_fits=save_components)
    sed_table, instrument_list = sed.build_cos_fuv(component_repo, airglow)
    
    #STIS FUV and Lya
    gratings = prepare_stis.make_stis_spectum(input_paths['STIS_FUV'], version, savepath = component_repo, save_ecsv=save_components, return_gratings=True, save_fits = save_components)
    print(gratings)
    if 'G140L' in gratings:
        stis_normfac = sed.find_stis_normfac(component_repo, airglow, 'fuv')
       # prepare_stis.make_stis_spectum(input_paths['STIS_FUV'], version, savepath = component_repo, save_ecsv=save_components, return_gratings=True, save_fits = save_components, normfac=stis_normfac, sx1=False)
    else:
        stis_normfac= 1.0
    prepare_lya.make_lya_spectrum(input_paths['lya_model'], version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, normfac=stis_normfac)
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, component_repo, lya_range, instrument_list, other_airglow)
    if 'G140L' not in gratings:
        print('adding polynomial fills')
        sed_table, instrument_list = sed.fill_cos_airglow(sed_table, other_airglow, instrument_list)
    if 'G230L' in gratings:
        nuv_normfac = sed.find_stis_normfac(component_repo, airglow, 'nuv')
        sed_table, instrument_list = sed.add_stis_nuv(sed_table, component_repo, instrument_list)
       # print (nuv_normfac)
    
    #works to here
    #NUV- STIS or COS
  #  if 'G230L' in gratings:
        
    
    
    
    sed_table.sort(['WAVELENGTH'])
                                              
    return sed_table, instrument_list
        

    
    
def gj_674_test():
    """
    Testing each stage with gj674
    """
    star = 'gj_674' #as it appears in the filepath
    star_up = 'GJ_674'
    path = '/home/david/work/muscles/SEDs/'+star+'/'
    muscles_path = '/home/david/work/muscles/MegaMUSCLES/'+star_up+'/'
    input_paths = dict(COS_readsav = path+'COS/', COS_x1d = muscles_path+'HST/COS/',STIS_FUV = muscles_path+'HST/STIS/', 
                       lya_model = path + 'lya/GJ674_intrinsic_LyA_profile.txt')
    lya_range = [1207, 1225] #lyman alpha region to remove
    other_airglow = [1304, 1304.5, 1355, 1356] #oi airglow to remove
    save_path = path + 'test_files/'
    version = 1
    sed_table, instrument_list = make_sed(input_paths, save_path, version, lya_range, other_airglow, save_components=True)
    
    #print(sed_table.sort('WAVELENGTH'))
    plt.figure(star+'_test')
    plt.step(sed_table['WAVELENGTH'], sed_table['FLUX'], where='mid')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    
def trappist_1_test():
    """
    Testing each stage with trappist-1
    """
    star = 'trappist-1' #as it appears in the filepath
    star_up = 'TRAPPIST-1'
    path = '/home/david/work/muscles/SEDs/'+star+'/'
    muscles_path = '/home/david/work/muscles/MegaMUSCLES/'+star_up+'/'
    input_paths = dict(COS_readsav = path+'COS/', COS_x1d = '/home/david/work/muscles/trappist-1/hst/data/',
                       STIS_FUV = '/home/david/work/muscles/trappist-1/hst/g140m_cals/picked_trace_extracts/' , 
                       lya_model = path + 'lya/Trappist-1_lya_simple.txt')
    lya_range = [1207, 1225] #lyman alpha region to remove
    other_airglow =  [1273.9, 1287.3, 1301, 1307]  #oi airglow to remove
    save_path = path + 'test_files/'
    version = 1
    sed_table, instrument_list = make_sed(input_paths, save_path, version, lya_range, other_airglow)
    
    #print(sed_table.sort('WAVELENGTH'))
    plt.figure(star+'_test')
    plt.step(sed_table['WAVELENGTH'], sed_table['FLUX'], where='mid')
    plt.show()
    
    
    
gj_674_test()
#trappist_1_test()
                                              
    #filepaths = {'xmm':'/xmm/GJ674.fits',
     #        'cos_g130m':'/COS/GJ674_COS130M_Mm1_NOSCL_03apr18.sav',
      #       'lya':'/lya/GJ674_intrinsic_LyA_profile.txt',
       ##      'stis_g140m':'/STIS/GJ674_G140M_coadd.ecsv',
            # 'stis_g140l':'/STIS/GJ674_G140L_noflare_x1d.fits',
         ##    'stis_g230l':'/STIS/GJ674_G230L_x1d.fits',
             #'stis_g430l':'/STIS/odlm21010_sx1.fits',
           #  'phoenix':'/PHOENIX/lte03400-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
            # 'phoenix_wave' :'/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits' }
    

