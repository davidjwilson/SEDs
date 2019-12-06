import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
from pylab import cm  
from matplotlib.colors import LogNorm


    
path = 'raw_data/'
flts = glob.glob(path +'*flt.fits')


for flt in flts:
    rootname = fits.getheader(flt,0)['ROOTNAME']
    data = fits.getdata(flt,1)
    fig = plt.figure(rootname, figsize=(13,13))
    plt.imshow(data, cmap=cm.gray_r, norm = LogNorm())
    plt.tight_layout()
    plt.savefig('trace_plots/'+rootname+'_trace.png')
    plt.close()