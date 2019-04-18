import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from scipy.io.idl import readsav
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u

"""
v1 20190418
Script to make Mega-Muscles spectra. Will develop over time.

objectives for this version
output: a file with everything that MUSCLES did i.e.
'WAVELENGTH' 
'WAVELENGTH0' = wavelength bin edges
'WAVELENGTH1'
'FLUX'
'ERROR'
'EXPTIME' exptime of that bit of spectrum
'DQ',
'EXPSTART' exposure dates, mjd
'EXPEND',
'INSTRUMENT', bitwise- could we do this with words or is it what hlsp requires?
'NORMFAC',
'BOLOFLUX' =flux normalised to bolometric flux
'BOLOERR'




"""


#totals = [w_full, f_full, e_full, dq_full, n_full,i_full]
#additions =[w, f, e, dq, n, instrument]

def add_spec(totals, additions):
    """running combinations of wavelength, flux, error, data quality flags, normalisation and instrument"""
    totals[i] =  np.concatenate((totals, additions), axis=1)
    return totals

def sort_totals(totals):
    """sorts totals array by wavelength"""
    arr1inds = totals[0].argsort()
    for t in totals:
        t = t[arr1inds]
    return totals
    
def totals_setup():
    """set up arrays to keep running totals of data"""
    names = ['WAVELENGTH','WAVELENGTH0', 'WAVELENGTH1','FLUX','ERROR','EXPTIME',
             'DQ','EXPSTART','EXPEND','INSTRUMENT','NORMFAC','BOLOFLUX' ,'BOLOERR']
    totals = []
    i = 0
    while i < len(names):
        totals.append([])
        i += 1
    totals = np.array(totals, dtype = float)
    return names, totals
    
    
def plot_spectrum(w, f):
    """plots a spectrum section"""
    plt.step(w,f)

def get_data(file_path, file_format, data_type):
    """reads spectrum, gets information. This is the complicated one.
    
    Information I need to provide for each spectrum:
    -where to get the data
    -how to extract it
    -what range to plot
    -what scaling to use
    
    Information I need to extract:
    -wavelength
    -wavelength edges (for some)
    -flux
    -error
    -dq
    -exptime
    -expstart
    -expend
    -instrument
    
    Information I need to caculate:
    -wavelength edges (for some)
    -boloflux
    -boloerror
    
    """
    if filetype = "sav":
        
        