import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from scipy.io.idl import readsav
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.time import Time
import spectra_combine as sc


names, totals = sc.totals_setup()
star = 'gj_674' #as it appears in the filepath
path = '/home/david/work/muscles/SEDs/'+star

plt.figure(star+'_combined',figsize = (13, 7))
plt.subplots_adjust(top = 0.95, right = 0.99, left = 0.07, bottom = 0.11)

#file_locations
filepaths = {'xmm':'xmm/GJ674.fits',
             'cos_g130m':'/COS/GJ674_COS130M_Mm1_NOSCL_03apr18.sav',
             'lya':'/lya/GJ674_intrinsic_LyA_profile.txt',
             'stis_g140m':'/STIS/GJ674_G140M_coadd.ecsv',
             'stis_g140l':'/STIS/GJ674_G140L_noflare_x1d.fits',
             'stis_g230L':'/STIS/GJ674_G230L_x1d.fits',
             'stis_g430L':'/STIS/odlm21010_sx1.fits',
             'phoenix':'/photometry/scaled_03400-4.50-0.0_phoenix_gj674.ecsv'}

cosdata = sc.read_idl(path+filepaths['cos_g130m'])
w = cosdata['w']
glowmask = (w <1207)|(w >1225)&(w <1304)|(w >1304.5)&(w <1355)|(w >1356)   #mask out lya, airglow
w, f, cos_additions = sc.masked_additions(cosdata, glowmask)
sc.plot_spectrum(w,f)
totals = sc.add_spec(totals, cos_additions)





#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Wavelength (\AA)', size=20)
plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)', size=20)
plt.axhline(0, ls='--', c='k')

totals = sc.sort_totals(totals)

plt.show()