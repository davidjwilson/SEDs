import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u


"""
@author: David Wilson

@version: 1 

@date :20190805

Turns AYs Lyman alpha fits into a standard MUSCLES file.
"""

def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff/2.
    w1 = w + diff/2.
    return w0, w1
    
def make_lya_metadata(new_data, normfac):
    """
    Makes the metadata for the ecsv files- eventually will be converted into the fits file
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    exptimes = []
    start_times = []
    end_times = []
    dates = []
    data = readsav(sav)
    for name in data['files']:
        x = x1dpath+str(name)[2:-1]+'_x1d.fits'
        hdr0 = fits.getheader(x,0)
        hdr1 = fits.getheader(x,1)
        exptimes.append(hdr1['EXPTIME'])
        start_times.append(hdr1['EXPSTART'])
        end_times.append(hdr1['EXPEND'])
        dates.append(hdr1['DATE-OBS'])
    muscles_name = 'Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanet Host Stars'
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','APERTURE','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN','EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['','',hdr0['OPT_ELEM'],'','','','','',muscles_name,'MUSCLES','David J. Wilson','','',min(dates),min(start_times),max(end_times),sum(exptimes),'SUM', 
                min(exptimes), max(exptimes), np.median(exptimes),1.0,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = hdr0[name]
        else:
            metadata[name] = filler
    return metadata


def make_lya_data(lya_path):
    """
    Makes the lya data array, assuming an input .txt file with WAVELENGTH and flux
    """
    