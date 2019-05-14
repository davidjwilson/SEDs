import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii


x1 = []

def on_key(event):
  global x1, x2, first
  if event.key == 'w':
    x1.append(event.xdata)
    print('%.3f' %event.xdata)
    plt.close()
    
path = '/home/david/work/muscles/trappist-1/hst/g140m_cals/all_obs/'
flts = glob.glob(path +'*flt.fits')

roots = []

i = 1
for flt in flts:
    print(i)
    rootname = fits.getheader(flt,0)['ROOTNAME']
    roots.append(rootname)
    data = fits.getdata(flt,1)
    fig = plt.figure('ROOTNAME', figsize=(10,10))
    plt.imshow(data)
    plt.axhline(500, c='r', ls ='--')
    plt.xlim(300, 600)
    plt.ylim(350,650)
    cid = fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()
    #plt.close()
    i+=1

    
savedat = Table([roots, x1], names=['ROOTNAME', 'A2CENTER'])
ascii.write(savedat, 'trace_centers.ecsv', format='ecsv', overwrite=True)