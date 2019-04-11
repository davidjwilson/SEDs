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
    n_full = np.concatenate((n_full, n))
    i_full = np.concatenate((i_full, instrument))
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
instrument = 'COSG130M'
data = readsav('../COS/TRAPPIST1_G130M_Mm1_NOSCL_10dec2018.sav')
glowmask = (data['Wave'] <1207)|(data['Wave'] >1225)&(data['Wave'] <1301)|(data['Wave'] >1307.)&(data['Wave'] <1355)|(data['Wave'] >1356)
plt.step(data['wave'][glowmask], data['flux'][glowmask])
end_130 = data['wave'][-1]

#COSG160M
instrument = 'COSG130M'
data = readsav('../COS/TRAPPIST1_G160M_3orb_Mm1_NOSCL_09dec2018.sav')
mask = (data['wave'] > end_130)
plt.step(data['wave'][mask], data['flux'][mask])


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

#data = Table([w_full*u.AA, f_full*u.erg/u.cm**2/u.s/u.AA, e_full*u.erg/u.cm**2/u.s/u.AA, n_full], names = ['WAVELENGTH', 'FLUX', 'ERROR', 'NORMFAC'] )
#ascii.write(data, star+'_data+phoenix_v1.ecsv', delimiter=',', format='ecsv', overwrite=True)

plt.show()