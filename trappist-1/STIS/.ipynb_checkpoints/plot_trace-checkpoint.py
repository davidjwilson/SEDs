import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
from pylab import cm  
from matplotlib.colors import LogNorm


    
path = '/home/david/work/muscles/trappist-1/hst/g140m_cals/all_obs/'
xpath = path +'../picked_trace_extracts/'
flts = glob.glob(path +'*flt.fits')
x1ds = glob.glob(xpath+'*x1d.fits')
i = 1
for flt in flts:
    print(i)
    rootname = fits.getheader(flt,0)['ROOTNAME']
    data = fits.getdata(flt,1)
    fig = plt.figure(rootname, figsize=(10,10))
    plt.imshow(data, cmap=cm.gray_r, norm = LogNorm())
    x1d_name = xpath+rootname+'_new_x1d.fits'
    if x1d_name in x1ds:
        ys = fits.getdata(x1d_name,1)[0]['EXTRLOCY'][::-1]
        plt.plot(np.arange(len(ys)),ys)
    #plt.axhline(500, c='r', ls ='--')
    #plt.xlim(300, 600)
    #plt.ylim(350,650)
    
    plt.show()
    plt.close()
    i+=1

