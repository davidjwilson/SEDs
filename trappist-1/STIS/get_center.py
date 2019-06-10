import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
from pylab import cm  
from matplotlib.colors import LogNorm

y1, y2, x1 = [], [], []

def on_key(event):
    global y1, y2, x1
    if event.key == 'w':
        y1.append(event.ydata)
        y2.append(0.0)
        x1.append(event.xdata)
        print('%.3f' %event.xdata, '%.3f' %event.ydata)
        plt.close()
    if event.key == 'e':
        y2.append(event.ydata)
        y1.append(0.0)
        x1.append(event.xdata)
        print('%.3f' %event.xdata, '%.3f' %event.ydata)
        plt.close()
    
path = '/home/david/work/muscles/trappist-1/hst/g140m_cals/all_obs/'
flts = glob.glob(path +'*flt.fits')

roots = []
#o_roots = ['od3v02010', 'od3v03010', 'od3v01020', 'od3v01010']
i = 1
for flt in flts[0:5]:
    print(i)
    rootname = fits.getheader(flt,0)['ROOTNAME']
    #if rootname in o_roots:
    roots.append(rootname)
    data = fits.getdata(flt,1)
    fig = plt.figure(rootname, figsize=(13,13))
    plt.imshow(data, cmap=cm.gray_r, norm = LogNorm())
    plt.tight_layout()
    #plt.axhline(500, c='r', ls ='--')
    #plt.xlim(300, 600)
    #plt.ylim(350,650)
    #plt.grid()
    cid = fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()
    #plt.close()
    i+=1

    
#savedat = Table([roots, x1, y1, y2], names=['ROOTNAME', 'XCOL', 'A2CENTER', 'UNCERTAIN'])
#ascii.write(savedat, 'new_trace_centers.ecsv', format='ecsv', overwrite=True)