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
Script to stick Trappist-1 data together
"""

def add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, w, f, e, dq, n, instrument):
    """running combinations of wavelength, flux, error, data quality flags, normalisation and instrument"""
    w_full = np.concatenate((w_full, w))
    f_full = np.concatenate((f_full, f))
    e_full = np.concatenate((e_full, e))
    dq_full = np.concatenate((dq_full, dq))
    n_full = np.concatenate((n_full, np.full(len(w),n)))
    i_full = np.concatenate((i_full, np.full(len(w), instrument)))
    return w_full, f_full, e_full, dq_full, n_full, i_full
                            

star = 'TRAPPIST-1'

w_full = np.array([], dtype=float)
f_full = np.array([], dtype=float)
e_full = np.array([], dtype=float)
dq_full = np.array([], dtype=float)
n_full = np.array([], dtype=float)
i_full = np.array([], dtype=str)

plt.figure(star+'_combined',figsize = (13, 7))
plt.subplots_adjust(top = 0.95, right = 0.99, left = 0.07, bottom = 0.11)

#COSG130M
instrument = 'COS_G130M'
data = readsav('../COS/TRAPPIST1_G130M_Mm1_NOSCL_10dec2018.sav')
glowmask = (data['Wave'] <1207)|(data['Wave'] >1225)&(data['Wave'] <1301)|(data['Wave'] >1307.)&(data['Wave'] <1355)|(data['Wave'] >1356)
plt.step(data['wave'][glowmask], data['flux'][glowmask])
end_130 = data['wave'][-1]
w_full, f_full, e_full, dq_full, n_full, i_full = add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, 
                                                   data['wave'][glowmask], data['flux'][glowmask], 
                                                   data['Err'][glowmask], np.zeros(len(data['flux'][glowmask]), dtype=int), 1.0, instrument)


#COSG160M
instrument = 'COS_G160M'
data = readsav('../COS/TRAPPIST1_G160M_3orb_Mm1_NOSCL_09dec2018.sav')
mask = (data['wave'] > end_130)
plt.step(data['wave'][mask], data['flux'][mask])
#plt.step(data['wave'][mask], data['err'][mask])
end_160 = data['wave'][mask][-1]
w_full, f_full, e_full, dq_full, n_full, i_full = add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, 
                                                   data['wave'][mask], data['flux'][mask], 
                                                   data['Err'][mask], np.zeros(len(data['flux'][mask]), dtype=int), 1.0, instrument)
                                                   
                                                   

#STIS G140M
instrument = 'STIS_G140M'
stis_scale = 1.0 #scaled to COS data
data = Table.read('../STIS/archival_TRAPPIST-1_G140M_mean.ecsv')
lyinc = (data['WAVELENGTH'] >1207)&(data['WAVELENGTH'] <1225) #include just lya
w, f = data['WAVELENGTH'][lyinc], data['FLUX'][lyinc]*stis_scale
plt.step(w, f)
w_full, f_full, e_full, dq_full, n_full, i_full = add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, 
                                                   w,f, 
                                                   np.zeros(len(w)), np.zeros(len(w), dtype=int), stis_scale, instrument)
                                                   
#COSG230L
instrument = 'COS_G230L'
specpath = '../../../trappist-1/hst/data/ldlm42010_x1dsum.fits'
cnw = np.array([], dtype=float)
cnf = np.array([], dtype=float)
cne = np.array([], dtype=float)
cndq = np.array([], dtype=int)
for dt in fits.getdata(specpath,1)[0:2]:
    cnw= np.concatenate((cnw, dt['WAVELENGTH']))
    cnf = np.concatenate((cnf, dt['FLUX']))
    cne = np.concatenate((cne, dt['ERROR']))
    cndq = np.concatenate((cndq, dt['DQ']))
mask = (cnw > end_160)
cnw, cnf, cne, cndq = cnw[mask], cnf[mask], cne[mask], cndq[mask] 
plt.step(cnw, cnf)
w_full, f_full, e_full, dq_full, n_full, i_full = add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, 
                                                   cnw,cnf, cne, cndq, 1.0, instrument)
nuv_end = cnw[-1]

#ccd
instrument = 'STIS_G430L'
stis_opt = '../../../trappist-1/hst/data/odlm41010_sx1.fits'
data = fits.getdata(stis_opt,1)[0]
w, f, e, dq = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
mask = (w > nuv_end)
w, f, e, dq = w[mask], f[mask], e[mask], dq[mask]
plt.step(w,f)
w_full, f_full, e_full, dq_full, n_full, i_full = add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, 
                                                   cnw,cnf, cne, cndq, 1.0, instrument)
w_end = w[-1]


#phoenix
m_scale = 4.71774681e-28
instrument = 'PHOENIX'
data = Table.read('../optical/scaled_02560-5.00-0.0_phoenix_trappist-1.ecsv')
mask = data['WAVELENGTH'] > w_end
mw, mf = data['WAVELENGTH'][mask], data['FLUX'][mask]
#mw, mf = data['WAVELENGTH'], data['FLUX']
plt.step(mw, mf, zorder=-10, )
#plt.plot(mw, mf, zorder=-10, c='0.5', ls='--')
w_full, f_full, e_full, dq_full, n_full, i_full = add_spec(w_full, f_full, e_full, dq_full, n_full,i_full, 
                                                   mw,mf, np.zeros(len(mw)),np.zeros(len(mw), dtype=int), m_scale, instrument)

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Wavelength (\AA)', size=20)
plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)', size=20)
plt.axhline(0, ls='--', c='k')


#sort complete arrays by wavelength

arr1inds = w_full.argsort()
w_full = w_full[arr1inds]
f_full = f_full[arr1inds]
e_full = e_full[arr1inds]
dq_full = e_full[arr1inds]
n_full = n_full[arr1inds]
i_full = i_full[arr1inds]

data = Table([w_full*u.AA, f_full*u.erg/u.cm**2/u.s/u.AA, e_full*u.erg/u.cm**2/u.s/u.AA, dq_full, n_full, i_full], names = ['WAVELENGTH', 'FLUX', 'ERROR', 'DQ', 'NORMFAC', 'INSTRUMENT'] )
ascii.write(data, star+'_hst+phoenix_v1.ecsv', delimiter=',', format='ecsv', overwrite=True)

plt.show()