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
import prepare_model
import prepare_xmm
import prepare_euv

#paths = some way of storing all the paths to the different spectra 
"""
input paths = dict('COS_readsav', 'COS_x1d', 'STIS_FUV', lya_model)
airglow (lya)
"""


def make_sed(input_paths, savepath, version, lya_range, other_airglow, save_components = False, star_params={}, phoenix_repo = '/home/david/work/muscles/phoenix/model_repo/', phoenix_wave = '/home/david/work/muscles/SEDs/common/wavegrid_hires.fits', phoenix_cut = 4000, do_phoenix=True, euv_inputs={}):
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
    prepare_model.make_model_spectrum(input_paths['lya_model'], version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, normfac=stis_normfac)
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, component_repo, lya_range, instrument_list, other_airglow)
    if 'G140L' not in gratings:
        print('adding polynomial fills')
        sed_table, instrument_list = sed.fill_cos_airglow(sed_table, other_airglow, instrument_list)
    
    #NUV- COS or STIS?
    if 'G230L' in gratings:
        nuv_normfac = sed.find_stis_normfac(component_repo, airglow, 'nuv')
        sed_table, instrument_list = sed.add_stis_nuv(sed_table, component_repo, instrument_list)
    else:
        gap_edges = prepare_cos.make_cos_nuv(input_paths['COS_x1d'], version, savepath = component_repo, plot=False, save_ecsv=save_components, save_fits=save_components, find_gap=True)
        #print (gap)
        sed_table, instrument_list = sed.add_cos_nuv(sed_table, component_repo, instrument_list, gap_edges)
       # print (nuv_normfac)
    #OPTICAL
    
    sed_table, instrument_list = sed.add_stis_optical(sed_table, component_repo, instrument_list)
    
    if do_phoenix: #phoenix takes ages so I'm adding the option to turn it off for testing purposes
        if star_params == {}:
            star_params = prepare_model.load_star_params(input_paths['STAR_PARAMS'], FeH=0.0, aM=0.0)
        if len(os.listdir(input_paths['PHOENIX'])) == 0:
            prepare_model.make_phoenix_spectrum(phoenix_wave,input_paths['PHOENIX'], phoenix_repo, star_params, save_ecsv=True, plot=False)
        prepare_model.make_model_spectrum(input_paths['PHOENIX']+os.listdir(input_paths['PHOENIX'])[0], version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, model_name='PHX')
        phoenix_normfac = sed.phoenix_norm(component_repo, cut=phoenix_cut)
        sed_table, instrument_list = sed.add_phx_spectrum(sed_table, component_repo, instrument_list)
                     
    #xray- xmm/chandra +model
    if 'XMM_path' in input_paths:
        prepare_xmm.make_xmm_spectra(input_paths['XMM_path'], component_repo, sed_table.meta, version, apec_repo=input_paths['APEC'], save_ecsv=save_components, save_fits=save_components)
        scope = 'xmm'
    prepare_model.make_model_spectrum(input_paths['APEC']+os.listdir(input_paths['APEC'])[0], version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, model_name='apec')
    
    sed_table, instrument_list, euv_gap = sed.add_xray_spectrum(sed_table, component_repo, instrument_list, scope, add_apec = True, find_gap=True)
    
    #EUV
    dem_path = ''
    euv_name = 'euv-scaling'
    if 'DEM_path' in input_paths:
        dem_path = input_paths['DEM_path']
        euv_name = 'dem'
    prepare_euv.make_euv(input_paths['EUV'], dem_path, euv_inputs=euv_inputs)
    prepare_model.make_model_spectrum(input_paths['EUV']+os.listdir(input_paths['EUV'])[0], version, sed_table ,savepath = component_repo, save_ecsv=save_components, save_fits=save_components, model_name=euv_name)
    
    sed_table, instrument_list = sed.add_euv(sed_table, component_repo, instrument_list, euv_gap, euv_name)
        
    
    sed_table.sort(['WAVELENGTH'])
                                              
    return sed_table, instrument_list
        

    
    
def gj_674_test():
    """
    Testing each stage with gj674
    """
    star = 'gj_674' #as it appears in the filepath
    star_up = 'GJ_674'
    star_params = {'Teff':3400, 'logg':4.5, 'FeH':0.0 , 'aM':0.0 }
    path = '/home/david/work/muscles/SEDs/'+star+'/'
    muscles_path = '/home/david/work/muscles/MegaMUSCLES/'+star_up+'/'
    input_paths = dict(COS_readsav = path+'COS/', 
                       COS_x1d = muscles_path+'HST/COS/',
                       STIS_FUV = muscles_path+'HST/STIS/', 
                       lya_model = path + 'lya/GJ674_intrinsic_LyA_profile.txt', 
                       PHOENIX= path+'phoenix_repo/',
                       XMM_path = path+'xmm/GJ674.fits',
                       APEC = path+'apec/',
                       EUV = path+'euv_repo/'
                       )
    lya_range = [1207, 1225] #lyman alpha region to remove
    other_airglow = [1304, 1304.5, 1355, 1356] #oi airglow to remove
    save_path = path + 'test_files/'
    version = 1
    euv_inputs = dict(lya=2.9*2.06e-12, distance=4.54 )
    sed_table, instrument_list = make_sed(input_paths, save_path, version, lya_range, other_airglow, save_components=True, star_params=star_params, do_phoenix=False, euv_inputs = euv_inputs)
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
    input_paths = dict(COS_readsav = path+'COS/', 
                       COS_x1d = '/home/david/work/muscles/trappist-1/hst/data/',
                       STIS_FUV = '/home/david/work/muscles/trappist-1/hst/stis_collection/' , 
                       lya_model = path + 'lya/Trappist-1_lya_simple.txt',
                       PHOENIX= path+'phoenix_repo/',
                       XMM_path = path+'xmm/Trappist-1.fits',
                       APEC = path+'apec/',
                       EUV = path+'euv_repo/',
                       DEM_path = path + 'dem/trappist-1_dem_spectra.fits'
                      )
    lya_range = [1207, 1225] #lyman alpha region to remove
    other_airglow =  [1273.9, 1287.3, 1301, 1307]  #oi airglow to remove
    save_path = path + 'test_files/'
    version = 1
    star_params = {'Teff':2560, 'logg':5.0, 'FeH':0.0 , 'aM':0.0 }
    sed_table, instrument_list = make_sed(input_paths, save_path, version, lya_range, other_airglow, save_components=True, star_params=star_params, do_phoenix=False)
    
    #print(sed_table.sort('WAVELENGTH'))
    plt.figure(star+'_test')
    plt.step(sed_table['WAVELENGTH'], sed_table['FLUX'], where='mid')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    
    
#gj_674_test()
trappist_1_test()
                                              
    #filepaths = {'xmm':'/xmm/GJ674.fits',
     #        'cos_g130m':'/COS/GJ674_COS130M_Mm1_NOSCL_03apr18.sav',
      #       'lya':'/lya/GJ674_intrinsic_LyA_profile.txt',
       ##      'stis_g140m':'/STIS/GJ674_G140M_coadd.ecsv',
            # 'stis_g140l':'/STIS/GJ674_G140L_noflare_x1d.fits',
         ##    'stis_g230l':'/STIS/GJ674_G230L_x1d.fits',
             #'stis_g430l':'/STIS/odlm21010_sx1.fits',
           #  'phoenix':'/PHOENIX/lte03400-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
            # 'phoenix_wave' :'/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits' }
    

