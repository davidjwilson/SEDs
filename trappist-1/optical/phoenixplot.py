import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from scipy.interpolate import interp1d
from scipy.io.idl import readsav
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian1DKernel

def smear(w,f, R, w_sample=2):
    '''
    Smears a model spectrum with a gaussian kernel to the given resolution, R.
    Adapeted from https://github.com/spacetelescope/pysynphot/issues/78

    Parameters
    -----------

    w,f:  spectrum to smear

    R: int
        The resolution (dL/L) to smear to

    w_sample: int
        Oversampling factor for smoothing

    Returns
    -----------

    sp: PySynphot Source Spectrum
        The smeared spectrum
    '''

    # Save original wavelength grid and units
    w_grid = w
    

    # Generate logarithmic wavelength grid for smoothing
    w_logmin = np.log10(np.nanmin(w_grid))
    w_logmax = np.log10(np.nanmax(w_grid))
    n_w = np.size(w_grid)*w_sample
    w_log = np.logspace(w_logmin, w_logmax, num=n_w)

    # Find stddev of Gaussian kernel for smoothing
    R_grid = (w_log[1:-1]+w_log[0:-2])/(w_log[1:-1]-w_log[0:-2])/2
    sigma = np.median(R_grid)/R
    if sigma < 1:
        sigma = 1

    # Interpolate on logarithmic grid
    f_log = np.interp(w_log, w_grid, f)

    # Smooth convolving with Gaussian kernel
    gauss = Gaussian1DKernel(stddev=sigma)
    f_conv = convolve_fft(f_log, gauss)

    # Interpolate back on original wavelength grid
    f_sm = np.interp(w_grid, w_log, f_conv)

    # Write smoothed spectrum back into Spectrum object
    return w_grid, f_sm

def error_cut(w, f, e, bin_width = 30): #cut region before a rolling 30pt mean SN > 1
    sn = np.array([np.mean(f[i:i+bin_width]/e[i:i+bin_width]) for i in range(len(w[:-bin_width]))])
    start = w[:-bin_width][np.where(sn > 1)[0][0]]
    mask = (w > start) & (f > 0)
    return w[mask], f[mask], e[mask]


specs = glob.glob('interpolated_models/*ecsv')

sdata = fits.getdata('../combined/odlm41010_sx1.fits',1)[0]
wo, fo, eo = sdata['WAVELENGTH'], sdata['FLUX'], sdata['ERROR']
wo2, fo2, eo2 = error_cut(wo, fo, eo)

irpath = '/home/david/work/muscles/SEDs/trappist-1/ir_data/'

g19data = glob.glob(irpath+'PS_Gaia*')
g19data

gw, gf, ge = np.loadtxt(g19data[0], unpack=True)
gw *=10000
gphot = np.genfromtxt(g19data[1], dtype=None, delimiter=',', names=True, encoding=None)
gpn, gpw, gpf, gpe = gphot['Band'], gphot['Wavelength'], gphot['Flux'], gphot['Error'] 
gpw*=10000


plt.figure()
plt.plot(gw, gf)
plt.errorbar(gpw, gpf, yerr=gpe, marker='o', ls='none')

#plt.xlim(gw[0], gw[-1])
#print(gw[0])

for sp in specs:
    data = Table.read(sp)
    #print (data.meta)
    w, f = data['WAVELENGTH'], data['FLUX']*data.meta['NORMFAC']
    w, f = smear(w, f, 4000) 
    plt.plot(w, f, label='{}'.format(data.meta['TEFF']))
plt.xscale('log')
plt.yscale('log')

plt.plot(wo2, fo2)

plt.show()
#plt.xlim(gw[0], gw[-1])
#plt.ylim(1e-16, 2.5e-14)
