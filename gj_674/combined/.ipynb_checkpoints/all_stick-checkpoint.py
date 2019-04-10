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
Script to help stick Mega-Muscules COS and STIS spectra together.
"""
star='GJ_674'

#with masks
plt.figure(star+'_combined',figsize = (13, 7))
plt.subplots_adjust(top = 0.95, right = 0.99, left = 0.07, bottom = 0.11)


#COS
#G130M
data = readsav('../COS/GJ674_COS130M_Mm1_NOSCL_03apr18.sav')
#mask = data['Flux'] > 0
#lymask = (data['Wave'] <1207)|(data['Wave'] >1225) #mask out lya
glowmask = (data['Wave'] <1207)|(data['Wave'] >1225)&(data['Wave'] <1304)|(data['Wave'] >1304.5)&(data['Wave'] <1355)|(data['Wave'] >1356)   #mask out lya, airglow

plt.step(data['Wave'][glowmask], data['Flux'][glowmask])
#plt.step(data['Wave'], data['Flux'])
cos_end = data['Wave'][-1]


#STIS
#G140M
stis_scale = 1.8 #scaled to COS data
data = Table.read('../STIS/GJ674_G140M_coadd.ecsv')
#mask = data['FLUX'] > 0
lyinc = (data['WAVELENGTH'] >1207)&(data['WAVELENGTH'] <1225) #include just lya
plt.step(data['WAVELENGTH'][lyinc], data['FLUX'][lyinc]*stis_scale)
#plt.step(data['WAVELENGTH'], data['FLUX'])

#G140L
data = fits.getdata('../STIS/GJ674_G140L_noflare_x1d.fits', 1)[0]
#mask = data['FLUX'] > 0
mask = (data['WAVELENGTH']>=cos_end) #only need the bit not covered by COS nb here no difference between > and >=, could change for other stars

plt.step(data['WAVELENGTH'][mask], data['FLUX'][mask]*stis_scale)
#plt.step(data['WAVELENGTH'], data['FLUX'])
g140L_end = data['WAVELENGTH'][-1]


#G230L
nuv_scale = 0.86
data = fits.getdata('../STIS/GJ674_G230L_x1d.fits')[0]
#clip_st, clip_end = 30,-6 
clip_end = -6 #don't need to clip the start any more
#mask = data['FLUX'][clipLst:clip_end] > 0
mask = data['WAVELENGTH'][:clip_end]>g140L_end
plt.step(data['WAVELENGTH'][:clip_end][mask], data['FLUX'][:clip_end][mask]*nuv_scale)
#plt.step(data['WAVELENGTH'][clip_st:clip_end], data['FLUX'][clip_st:clip_end])
g230L_end = data['WAVELENGTH'][clip_end] 

#G430L nb normalise to photometry, PHOENIX
scale = 1.08
ccd = '../STIS/odlm21010_sx1.fits'
data = fits.getdata(ccd)[0]
clip_end = -1 #don't need to clip the start any more
#mask = data['FLUX'][clipLst:clip_end] > 0
mask = data['WAVELENGTH'][:clip_end]>g230L_end
#mask = data['WAVELENGTH'][:clip_end]>3500.
plt.step(data['WAVELENGTH'][:clip_end][mask], data['FLUX'][:clip_end][mask]*scale)
w_end = data['WAVELENGTH'][:clip_end][mask][-1]

#Phoenix
data = Table.read('../photometry/scaled_03400-4.50-0.0_phoenix_gj674.ecsv')
mask = data['WAVELENGTH'] > w_end
mw, mf = data['WAVELENGTH'][mask], data['FLUX'][mask]
#mw, mf = data['WAVELENGTH'], data['FLUX']
plt.step(mw, mf, zorder=-10, )
#plt.plot(mw, mf, zorder=-10, c='0.5', ls='--')
#xmm
xdt = fits.getdata('../xmm/GJ674.fits',1)
plt.step(xdt['wave'], xdt['cflux'])


plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wavelength (\AA)', size=20)
plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)', size=20)
plt.axhline(0, ls='--', c='k')




plt.show()
