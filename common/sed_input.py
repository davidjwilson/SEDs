"""
@verison: 1

@author: David Wilson

@date 20190805

Using this to wite the main function for make_mm_sed.
"""

import make_mm_sed as sed
import prepare_cos
import prepare_stis

#paths = some way of storing all the paths to the different spectra 
"""
input paths = dict('COS_readsav', 'COS_x1d', 'STIS_FUV')
airglow (lya)
"""


def make_sed(datapaths, savepath, version, lya_range, other_airglow):
    airglow = lya_range+other_airglow
    #COS FUV 
    component_repo = savepath+'componets' #directory where the component spectra are saved
    prepare_cos.make_cos_spectrum(input_paths['COS_readsav'], version,  input_paths['COS_x1d'], savepath = component_repo, plot=False, save_ecsv=True, save_fits=True)
    sed_table, instrument_list = sed.build_cos_fuv(component_repo, airglow)
    
    #STIS FUV
    gratings = make_stis_spectum(input_paths['STIS_FUV'], version, savepath = componet_repo, save_ecsv=True, return_gratings=True, save_fit = True) 
    if 'G140L' in gratings:
        stis_normfac = sed.find_stis_norm(component_repo, airglow)
        
        

    

