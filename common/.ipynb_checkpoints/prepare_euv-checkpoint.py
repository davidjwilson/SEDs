import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import scipy.interpolate as interpolate
import os
import astropy.io.fits as fits

"""

@verison: 3

@author: David Wilson

@date 20200504

Makes txt files for the EUV portion of a spectrum that can then be fed into prepare_model. Either uses a DEM model or the Linsky+14 relationships

v3 added ability to put star name in files

Calculating the EUV fluxes using the relationships of Linsky + 14 (https://ui.adsabs.harvard.edu/abs/2014ApJ...780...61L/abstract)

log[F(delta lambda) /F(lya)]=

10–20 nm (stars) 	 	−0.491 	 
20–30 nm (stars) 	 	−0.548 	 
30–40 nm (stars) 		−0.602 	 
40–50 nm (models) 	  	  	−2.294+0.258 log[f (Lyα)]
50–60 nm (models) 	  	  	−2.098+0.572 log[f (Lyα)]
60–70 nm (models) 	  	  	−1.920+0.240 log[f (Lyα)]
70–80 nm (models) 	  	  	−1.894+0.518 log[f (Lyα)]
80–91.2 nm (models) 	  	  	−1.811+0.764 log[f (Lyα)]
91.2–117 nm (models) 	  	  	−1.004+0.065 log[f (Lyα)]

"""

def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff0 = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    diff1 = np.concatenate((diff, np.array([diff[-1]]))) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff0/2.
    w1 = w + diff1/2.
    return w0, w1

def euv_estimator(euv_inputs, save_path, star, save=True):
    
    lya, distance = euv_inputs['lya'], euv_inputs['distance']
        
    distance_conversion = ((1*u.au.to(u.m))/(distance*u.pc.to(u.m)))**2

    lya_1au = lya / distance_conversion

    w1 = np.array([100,200,300,400,500,600,700,800,912], dtype=float) #A
    w2 = np.array([200,300,400,500,600,700,800,912,1170], dtype=float)
    bandwidth = w2-w1

    a = np.array([-0.491,-0.548,-0.602,-2.294,-2.098,-1.920,-1.894,-1.811,-1.004], dtype=float)
    b = np.array([ 0.,    0.,    0.,    0.258, 0.572, 0.240, 0.518, 0.764, 0.065], dtype=float)

    #log(f/lya) = a + b *log(lya)
    f = a + b*np.log10(lya_1au)
    

   # print('Total EUV=',np.sum(f))
    f = (lya_1au * 10**f)/bandwidth

    f *= distance_conversion

    #extrapolate onto 1A grid
    wavelength = np.arange((w1[0])+0.5, (w2[-1])+0.5, 1.0)
    flux = interpolate.interp1d(np.mean([w1, w2], axis=0), f, kind='nearest', bounds_error=False, fill_value='extrapolate')(wavelength)
 

    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    if star != '':
        star = '{}_'.format(star)
    filename = '{}l14euv_lya{}_d{}.txt'.format(star, euv_inputs['lya'], euv_inputs['distance'])
    savedat = Table([wavelength, flux], names=['WAVELENGTH', 'FLUX'])
    ascii.write(savedat, save_path+filename, overwrite=True)
    
def dem_to_1A(w,f):
    """
    Converts a DEM model at 5A resolution to 1A resolution
    """
    w1 = np.arange(w[0], w[-1], 1.)
    f1 = interp1d(w, f, fill_value='extrapolate', kind='nearest')(w1)
    return w1, f1, 
    
def make_dem(dem_path, save_path):
    """
    Extracts the wavelength and flux from a DEM and saves it as a text file
    """
    data = fits.getdata(dem_path, 1)
    wavelength, flux = data['Wavelength'], data['Flux_density']
    w0, w1 = wavelength_edges(wavelength)
    #flux = bin_flux/(w1-w0) #convert from bin-intergrated flux to flux -not required for new-generation seds
    name = 'dem.txt' 
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    savedat = Table([wavelength, flux], names=['WAVELENGTH', 'FLUX'])
    ascii.write(savedat, save_path+name, overwrite=True)


def make_euv(savepath, star = '', dem_path = '', euv_inputs = {}):
    """
    Main fuction. Uses a DEM if available, if not makes the EUV spectrum from euveuv_txt_files_inputes (lya flux, distance)
    """
    if dem_path == '':
        euv_estimator(euv_inputs, savepath, star)
    else:
        make_dem(dem_path, savepath)
        
        
        
#star = 'GJ674'
#lya = 2.9*2.06e-12 #erg /s/cm2 for GJ674, 2.9 is flux scaling factor 
#distance = 4.54
#star = 'GJ176'
#lya = 3.9e-13
#distance = 9.3 
#star = 'HD_97658'
#lya = 9.1e-13
#distance = 21.1
#print(lya)
#star = 'TRAPPIST-1'
#distance = 12.1
#lya = 2.12e-15