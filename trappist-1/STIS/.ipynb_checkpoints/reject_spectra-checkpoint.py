import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
from pylab import cm  
from matplotlib.colors import LogNorm
from astropy.convolution import convolve, Box1DKernel


def no_zero_errors(flux, error):
    """
    Corrects instances where negative flux measurements have very small errors
    """
    e_new = error
    for i in range(len(error)):
        if flux[i] < 0.0 and error[i] < 0.1*abs(flux[i]):
            e_new[i] = abs(flux[i])
    return e_new

y2 = []

def on_key(event):
    global y2
    if event.key == 'w':        
        y2.append(0.0)
        plt.close()
    if event.key == 'e':
        y2.append(1.0)
        plt.close()
    
#path = '/home/david/work/muscles/trappist-1/hst/g140m_cals/all_obs/'
path = 'new_x1ds/'
x1ds = glob.glob(path +'*x1d.fits')

roots = []
dates = []
#o_roots = ['od3v02010', 'od3v03010', 'od3v01020', 'od3v01010']

smooth =2
for i, x1d in enumerate(x1ds):
    print(i)
    hdul = fits.open(x1d)
    rootname = hdul[0].header['ROOTNAME']
    date = hdul[0].header['TEXPSTRT']
    roots.append(rootname)
    dates.append(date)
    data = hdul[1].data[0]
    fig = plt.figure(rootname, figsize=(14,6))
    wavelength, flux, error = data['WAVELENGTH'], data['FLUX'], data['ERROR'] 
    error = no_zero_errors(flux, error)
    flux = convolve(flux,Box1DKernel(smooth))
    error = convolve(error,Box1DKernel(smooth))/(smooth**0.5)
    
    
    plt.step(wavelength, flux, where='mid')
    plt.step(wavelength, error, where='mid')
    plt.tight_layout()
    plt.axhline(0, c='0.5', ls ='--')
    plt.axvline(1215.44, ls ='--', c='r')
    plt.ylim(-1e-15, 1.4e-14)
    plt.xlim(1210, 1220)

    #plt.xlim(300, 600)
    #plt.ylim(350,650)
    #plt.grid()
    cid = fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()
  
  
    
savedat = Table([roots, dates, y2], names=['ROOTNAME', 'DATE', 'KEEP'])
ascii.write(savedat, 'clean_spectra.ecsv', format='ecsv', overwrite=True)