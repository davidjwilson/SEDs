import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import interpolate
from craftroom import resample
# import stistools
from astropy.convolution import convolve, Box1DKernel
from astropy.modeling import models, fitting
from scipy.io.idl import readsav
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian1DKernel

import math as mt



def sed_to_const_res(sed_table, res=1):
    """
    Rebins an SED to a wavelength grid with a bin size of res, default = 1A
    """
    
        ('WAVELENGTH',
     'WAVELENGTH0',
     'WAVELENGTH1',
     'FLUX',
     'ERROR',
     'EXPTIME',
     'DQ',
     'EXPSTART',
     'EXPEND',
     'INSTRUMENT',
     'NORMFAC',
     'BOLOFLUX',
     'BOLOERR')
    
    #wavelength 
    start, end= mt.ceil(w[0]), mt.floor(w[-1])
    new_wavelength = np.arange(start,end+res, res)
    new_w0 = new_wavelength - (0.5 * res)
    new_w1 = new_wavelength + (0.5 * res)
    
    #flux and error
    new_wavelength, new_flux, new_error = resample.bintogrid(sed_table['WAVELENGTH'], sed_table['FLUX'], newxnew_wavelength, unc= sed_table['ERROR'])
    
    #exptime - linear extrapolation is similar to averaged to bin widths
    new_exptime = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPTIME'])(new_wavelength)
    
    #dq - interploate, then look for unusual values and correct them, summing if the values to either side are different.
    dqs = np.unique(sed_table['DQ'])
    new_dq = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['DQ'])(new_wavelength)
    for i in range(new_dq):
        if new_dq[i] not in dqs:
            new_dq[i] == new_dq[i-1] + new_dq[i+1]
    
    #expstart - minumum expstart in each bin
    startups = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPSTART'], kind='next')(new_wavelength)
    startdowns = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPSTART'], kind='previous')(new_wavelength)
    new_expstart = np.min([startups, startdowns], axis=0)
    
    #expends - maximum expend in each bin
    endups = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPEND'], kind='next')(new_wavelength)
    enddowns = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPEND'], kind='previous')(new_wavelength)
    new_expend = np.max([endups, enddowns], axis=0)
    
    #instrument - as dqs
    insts = np.unique(sed_table['DQ'])
    new_instrument = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['INSTRUMENT'])(new_wavelength)
    for i in range(new_instrument):
        if new_instrument[i] not in insts:
            new_instrument[i] == new_instrument[i-1] + new_instrument[i+1]
    #combine loops with dq later
    
    #normfac - linear extrapolation
    new_normfac = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['NORMFAC'])(new_wavelength)
    
    #boloflux -use original boloflux for consitency
    bolo_int = sed_table.meta['BOLOFLUX']*(u.erg/u.s/u.cm**2)
    new_boloflux = (new_flux/bolo_int).value
    new_boloerr = (new_error/bolo_int).value
    
    