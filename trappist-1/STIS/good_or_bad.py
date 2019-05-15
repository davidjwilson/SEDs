import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
#isolate bad spectra

path = '/home/david/work/muscles/trappist-1/hst/g140m_cals/all_obs/'
flts = glob.glob(path +'*flt.fits')

good = []
bad = []
ugly = []


i = 1
for flt in flts:
    print(i)
    rootname = fits.getheader(flt,0)['ROOTNAME']
    data = fits.getdata(flt,1)
    fig = plt.figure('ROOTNAME', figsize=(10,10))
    plt.imshow(data)
    plt.axhline(500, c='r', ls ='--')
    plt.xlim(300, 600)
    plt.ylim(350,650)
    plt.show()
    choice = True
    while choice:
        qual = input('Good (g), bad (b) or ugly (u):')
        if qual == 'g':
            good.append(rootname)
            choice = False
           # plt.close()
        elif qual == 'b':
            bad.append(rootname)
            choice = False
            #plt.close()
        elif qual == 'u':
            ugly.append(rootname)
            choice = False
            #plt.close()
            

  #  plt.show()
    #plt.close()
    i+=1
for a, n in zip([good, bad, ugly], ['good', 'bad', 'ugly']):
    if len(a) > 0:
        ascii.write(Table([a]), n+'_roots.csv', format='csv', overwrite=True)    
#savedat = Table([good, bad, ugly], names=['GOOD', 'BAD', 'UGLY'])
#ascii.write(savedat, 'flt_quality.ecsv', format='ecsv', overwrite=True)