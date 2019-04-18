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

def euv_estimator(star, lya, distance, save=False, plot=False):
    """
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
    distance_conversion = ((1*u.au.to(u.m))/(distance*u.pc.to(u.m)))**2

    lya_1au = lya / distance_conversion

    w1 = np.array([100,200,300,400,500,600,700,800,912], dtype=float) #A
    w2 = np.array([200,300,400,500,600,700,800,912,1170], dtype=float)
    bandwidth = w2-w1

    a = np.array([-0.491,-0.548,-0.602,-2.294,-2.098,-1.920,-1.894,-1.811,-1.004], dtype=float)
    b = np.array([ 0.,    0.,    0.,    0.258, 0.572, 0.240, 0.518, 0.764, 0.065], dtype=float)

    #log(f/lya) = a + b *log(lya)
    f = a + b*np.log10(lya_1au)

    print('Total EUV=',np.sum(f))
    f = (lya_1au * 10**f)/bandwidth

    f *= distance_conversion

    #extrapolate onto 1A grid
    wav = np.arange((w1[0])+0.5, (w2[-1])+0.5, 1.0)
    flux = []
    for w1i, w2i, fi in zip(w1, w2,f):
        for wi in wav:
            if wi > w1i and wi < w2i :
                flux.append(fi)

    if save == True:
        data = Table([wav*u.AA,  flux*u.erg/u.s/u.cm**2], names=['WAVELENGTH', 'FLUX'])
        ascii.write(data, star+'_1Aeuv_estimate.ecsv', delimiter=',', format='ecsv', overwrite=True)

    if plot == True:
        plt.figure(star+'_EUV', figsize=(8,6))
        plt.subplots_adjust(top=0.99, right=0.99)
        plt.plot(wav, flux)
        plt.xlabel('Wavelength (\AA)', size=20)
        plt.ylabel('Flux (erg s$^{-1}$\AA$^{-1}$cm$^{-2}$)', size=20)
        plt.yscale('log')
        plt.show()

    return(wav, flux) 

def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff = np.concatenate((diff[0], diff)) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff/2.
    w1 = w + diff/2.
    return w0, w1

def read_idl(filepath):
    """
    Extracts data from KF's IDL files
    """
    data = readsav(filepath)
    w, f, e, exptime = data['wave'], data['flux'], data['err'], data['exptime']
    dq = np.zeros(len(w)) #consult with kf on this
    w0, w1 = wavelength_edges(w)


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

