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
w_full = np.array([], dtype=float)
f_full = np.array([], dtype=float)
e_full = np.array([], dtype=float)
n_full = np.array([], dtype=float)

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
cos_start, cos_end = data['Wave'][0], data['Wave'][-1]
w_full = np.concatenate((w_full, data['Wave'][glowmask]))
f_full = np.concatenate((f_full, data['Flux'][glowmask]))
e_full = np.concatenate((e_full, data['Err'][glowmask]))
n_full = np.concatenate((n_full, np.full(len(data['Wave'][glowmask]), 1.0)))

#lya
g140m_scale = 2.9 #1.8 #scaled to COS data
lya = Table.read('../lya/GJ674_intrinsic_LyA_profile.txt', format='ascii')
plt.plot(lya['WAVELENGTH'], lya['FLUX']*g140m_scale)
#lyast, lyaed = lya['WAVELENGTH'][0], lya['WAVELENGTH'][-1]

w_full = np.concatenate((w_full, lya['WAVELENGTH']))
f_full = np.concatenate((f_full, lya['FLUX']*g140m_scale))
e_full = np.concatenate((e_full, np.zeros(len(lya['WAVELENGTH']))))
n_full = np.concatenate((n_full, np.full(len(lya['WAVELENGTH']),g140m_scale)))
             
#STIS
#G140M
data = Table.read('../STIS/GJ674_G140M_coadd.ecsv')
#mask = data['FLUX'] > 0
lyinc = (data['WAVELENGTH'] >1207)&(data['WAVELENGTH'] <lya['WAVELENGTH'][0])|(data['WAVELENGTH'] >lya['WAVELENGTH'][-1])&(data['WAVELENGTH'] <1225) #include just bits not covered by cos or lya reconstruction
#lyinc = (data['WAVELENGTH'] >1207)&(data['WAVELENGTH'] <1225)
plt.step(data['WAVELENGTH'][lyinc], data['FLUX'][lyinc]*g140m_scale)
#plt.step(data['WAVELENGTH'], data['FLUX'])

w_full = np.concatenate((w_full, data['WAVELENGTH'][lyinc]))
f_full = np.concatenate((f_full, data['FLUX'][lyinc]*g140m_scale))
e_full = np.concatenate((e_full, data['ERROR'][lyinc]*g140m_scale))
n_full = np.concatenate((n_full, np.full(len(data['WAVELENGTH'][lyinc]), g140m_scale)))


#G140L
g140l_scale = 3.5
data = fits.getdata('../STIS/GJ674_G140L_noflare_x1d.fits', 1)[0]
#mask = data['FLUX'] > 0
mask = (data['WAVELENGTH'] > 1304)&(data['WAVELENGTH'] <1304.5)|(data['WAVELENGTH'] > 1355)&(data['WAVELENGTH'] < 1356)|(data['WAVELENGTH']>=cos_end) #only need the bit not covered by COS nb here no difference between > and >=, could change for other stars

plt.step(data['WAVELENGTH'][mask], data['FLUX'][mask]*g140l_scale)
#plt.step(data['WAVELENGTH'], data['FLUX'])
g140L_end = data['WAVELENGTH'][-1]


w_full = np.concatenate((w_full, data['WAVELENGTH'][mask]))
f_full = np.concatenate((f_full, data['FLUX'][mask]*g140l_scale))
e_full = np.concatenate((e_full, data['ERROR'][mask]*g140l_scale))
n_full = np.concatenate((n_full, np.full(len(data['WAVELENGTH'][mask]), g140l_scale)))



#G230L
nuv_scale = 1.43
data = fits.getdata('../STIS/GJ674_G230L_x1d.fits')[0]
#clip_st, clip_end = 30,-6 
clip_end = -6 #don't need to clip the start any more
#mask = data['FLUX'][clipLst:clip_end] > 0
mask = data['WAVELENGTH'][:clip_end]>g140L_end
plt.step(data['WAVELENGTH'][:clip_end][mask], data['FLUX'][:clip_end][mask]*nuv_scale)
#plt.step(data['WAVELENGTH'][clip_st:clip_end], data['FLUX'][clip_st:clip_end])
g230L_end = data['WAVELENGTH'][clip_end] 

w_full = np.concatenate((w_full, data['WAVELENGTH'][:clip_end][mask]))
f_full = np.concatenate((f_full, data['FLUX'][:clip_end][mask]*nuv_scale))
e_full = np.concatenate((e_full, data['ERROR'][:clip_end][mask]*nuv_scale))
n_full = np.concatenate((n_full, np.full(len(data['WAVELENGTH'][:clip_end][mask]), nuv_scale)))

#G430L nb normalise to photometry, PHOENIX
ccd_scale = 1.11
ccd = '../STIS/odlm21010_sx1.fits'
data = fits.getdata(ccd)[0]
clip_end = -1 #don't need to clip the start any more
#mask = data['FLUX'][clipLst:clip_end] > 0
mask = data['WAVELENGTH'][:clip_end]>g230L_end
#mask = data['WAVELENGTH'][:clip_end]>3500.
plt.step(data['WAVELENGTH'][:clip_end][mask], data['FLUX'][:clip_end][mask]*ccd_scale)
w_end = data['WAVELENGTH'][:clip_end][mask][-1]


w_full = np.concatenate((w_full, data['WAVELENGTH'][:clip_end][mask]))
f_full = np.concatenate((f_full, data['FLUX'][:clip_end][mask]*ccd_scale))
e_full = np.concatenate((e_full, data['ERROR'][:clip_end][mask]*ccd_scale))
n_full = np.concatenate((n_full, np.full(len(data['WAVELENGTH'][:clip_end][mask]), ccd_scale)))


#Phoenix
data = Table.read('../photometry/scaled_03400-4.50-0.0_phoenix_gj674.ecsv')
mask = data['WAVELENGTH'] > w_end
mw, mf = data['WAVELENGTH'][mask], data['FLUX'][mask]
#mw, mf = data['WAVELENGTH'], data['FLUX']
plt.step(mw, mf, zorder=-10, )
#plt.plot(mw, mf, zorder=-10, c='0.5', ls='--')

w_full = np.concatenate((w_full, mw))
f_full = np.concatenate((f_full, mf))
e_full = np.concatenate((e_full, np.full(len(mw),0.0)))
n_full = np.concatenate((n_full, np.full(len(mw), 3.19143165e-26)))


#xmm
xdt = fits.getdata('../xmm/GJ674.fits',1)
plt.step(xdt['wave'], xdt['cflux'])
xray_end = xdt['wave'][-1]

w_full = np.concatenate((w_full, xdt['wave']))
f_full = np.concatenate((f_full, xdt['cflux']))
e_full = np.concatenate((e_full, xdt['cflux_err']))
n_full = np.concatenate((n_full, np.full(len(xdt['wave']), 1.0)))

#euv estimates
euv = Table.read('GJ674_1Aeuv_estimate.ecsv')
mask = (euv['WAVELENGTH'] < cos_start)
plt.step(euv['WAVELENGTH'][mask], euv['FLUX'][mask])
euv_start = euv['WAVELENGTH'][0]

w_full = np.concatenate((w_full, euv['WAVELENGTH']))
f_full = np.concatenate((f_full, euv['FLUX']))
e_full = np.concatenate((e_full, np.full(len(euv['WAVELENGTH']), 0.0)))
n_full = np.concatenate((n_full, np.full(len(euv['WAVELENGTH']), 1.0)))

#xray model
xdt = fits.getdata('../xmm/GJ674.fits',2)
mask = (xdt['wave'] > xray_end) & (xdt['wave'] < euv_start)
plt.step(xdt['wave'][mask], xdt['flux'][mask])

w_full = np.concatenate((w_full, xdt['wave'][mask]))
f_full = np.concatenate((f_full, xdt['flux'][mask]))
e_full = np.concatenate((e_full, np.full(len(xdt['wave'][mask]), 0.0)))
n_full = np.concatenate((n_full, np.full(len(xdt['wave'][mask]), 1.0)))

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wavelength (\AA)', size=20)
plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)', size=20)
plt.axhline(0, ls='--', c='k')

#sort complete arrays by wavelength

arr1inds = w_full.argsort()
w_full = w_full[arr1inds]
f_full = f_full[arr1inds]
e_full = e_full[arr1inds]
n_full = n_full[arr1inds]

data = Table([w_full*u.AA, f_full*u.erg/u.cm**2/u.s/u.AA, e_full*u.erg/u.cm**2/u.s/u.AA, n_full], names = ['WAVELENGTH', 'FLUX', 'ERROR', 'NORMFAC'] )
ascii.write(data, 'gj674_data+models_v1.ecsv', delimiter=',', format='ecsv', overwrite=True)

plt.show()
