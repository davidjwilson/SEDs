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

from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler
from astropy.nddata import StdDevUncertainty


def spectrum_to_const_res(sed_table, res=1):
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
    """
    
    
    