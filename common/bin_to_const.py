import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import interpolate
import math as mt
from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler
from astropy.nddata import StdDevUncertainty
from astropy.units import cds
cds.enable()

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

def spectrum_to_const_res(spectrum, res=1):
    """
    Rebins an MUSCLES-style component spectrum to a wavelength grid with a bin size of res, default = 1A
    Expecting a MM spectrum with 
    'WAVELENGTH',
    'WAVELENGTH0',
    'WAVELENGTH1',
    'FLUX',
    'ERROR',
    'EXPTIME',
    'DQ',
    'EXPSTART',
    'EXPEND'    
    or
    'WAVELENGTH',
    'WAVELENGTH0',
    'WAVELENGTH1',
    'FLUX'
    if a model
    """
    fluxcon = FluxConservingResampler(extrapolation_treatment='zero_fill')
    start, end= mt.ceil(spectrum['WAVELENGTH'][0]), mt.floor(spectrum['WAVELENGTH'][-1])
    new_wavelength = np.arange(start,end+res, res)
    
    if 'ERROR' in spectrum.dtype.names:                  
        print('binning an observed spectrum')
        input_spec = Spectrum1D(spectral_axis=(spectrum['WAVELENGTH'].value)*u.AA, flux=spectrum['FLUX']*u.Unit('erg cm-2 s-1 AA-1') , uncertainty= StdDevUncertainty(spectrum['ERROR']))
        new_spec_fluxcon = fluxcon(input_spec, new_wavelength*u.AA)
        new_wavelength = (new_spec_fluxcon.spectral_axis.value)
        new_flux = (new_spec_fluxcon.flux.value)
        new_error = (1/(new_spec_fluxcon.uncertainty.array**0.5))
        
        new_flux = np.nan_to_num(new_flux, 0)
        new_error = np.nan_to_num(new_error, 0)
        
        new_w0, new_w1 = wavelength_edges(new_wavelength)

        new_exptime = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPTIME'])(new_wavelength)

        #dq - interploate, then look for unusual values and correct them, summing if the values to either side are different.

        new_dq = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['DQ'], kind='previous')(new_wavelength)
        new_dq = new_dq.astype(int)

        #expstart - minumum expstart in each bin
        startups = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPSTART'], kind='next')(new_wavelength)
        startdowns = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPSTART'], kind='previous')(new_wavelength)
        new_expstart = np.min([startups, startdowns], axis=0)

        #expends - maximum expend in each bin
        endups = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPEND'], kind='next')(new_wavelength)
        enddowns = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPEND'], kind='previous')(new_wavelength)
        new_expend = np.max([endups, enddowns], axis=0)

        names = spectrum.dtype.names
        # print(names)
        new_spectrum = Table([new_wavelength*u.AA, new_w0*u.AA, new_w1*u.AA, new_flux*u.erg/u.s/u.cm**2/u.AA, new_error*u.erg/u.s/u.cm**2/u.AA, new_exptime*u.s, 
                               new_dq,new_expstart*cds.MJD, new_expend*cds.MJD], names=names, meta= spectrum.meta)
        
    else: #if it's a phoenix model it breaks the data down into chunks and makes a spectrum out to 1.3e5 A, where the native resolution drops below 1A
        if spectrum.meta['INSTRUME'] == 'PHX':
            print('binning a PHOENIX model')
            mask = spectrum['WAVELENGTH'] < 1.3e5
            w1, f1 = spectrum['WAVELENGTH'][mask], spectrum['FLUX'][mask]
            start = w1[0]
            step = 1000
            new_wavelength = np.array([], dtype=float)
            new_flux = np.array([], dtype=float)
            while start <= w1[-1]:
                mask = (w1 >= start-1) & (w1 < start+step+1) #plus one fixes nans
                w2, f2 = w1[mask], f1[mask]
                neww = np.arange(w2[0]+1, w2[-1]-1, 1) #but you have to take it off again
                input_spec = Spectrum1D(spectral_axis=(w2.value)*u.AA, flux=f2*u.Unit('erg cm-2 s-1 AA-1') )
                new_spec_fluxcon = fluxcon(input_spec, neww*u.AA)
                new_wavelength = np.concatenate((new_wavelength, new_spec_fluxcon.spectral_axis.value))
                new_flux = np.concatenate((new_flux, new_spec_fluxcon.flux.value))
                start += step
            new_flux = np.nan_to_num(new_flux, 0)
                
        else:
            print('binning a non-PHOENIX model')
            input_spec = Spectrum1D(spectral_axis=(spectrum['WAVELENGTH'].value)*u.AA, flux=spectrum['FLUX']*u.Unit('erg cm-2 s-1 AA-1'))
            new_spec_fluxcon = fluxcon(input_spec, new_wavelength*u.AA)
            new_wavelength = new_spec_fluxcon.spectral_axis.value
            new_flux = new_spec_fluxcon.flux.value
            new_flux = np.nan_to_num(new_flux, 0)
            
        new_w0, new_w1 = wavelength_edges(new_wavelength)
        
        names = spectrum.dtype.names
        new_spectrum = Table([new_wavelength*u.AA, new_w0*u.AA, new_w1*u.AA, new_flux*u.erg/u.s/u.cm**2/u.AA], names=names, meta= spectrum.meta)
        
    return new_spectrum


        
        
    