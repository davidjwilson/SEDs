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
star = 'trappist-1' #as it appears in the filepath
path = '/home/david/work/muscles/SEDs/'+star

plt.figure(star+'_combined',figsize = (13, 7))
plt.subplots_adjust(top = 0.95, right = 0.99, left = 0.07, bottom = 0.11)

#file_locations
filepaths = {'xmm':'/xmm/Trappist-1.fits',
             'cos_g130m':'/COS/TRAPPIST1_G130M_Mm1_NOSCL_10dec2018.sav',
             'cos_g160m':'/COS/TRAPPIST1_G160M_3orb_Mm1_NOSCL_09dec2018.sav'
             'lya':'/lya/GJ674_intrinsic_LyA_profile.txt',
             'cos_g230l': '/COS/ldlm42010_x1dsum.fits'
             'stis_g140m':'/STIS/TRAPPIST-1_G140M_mean.ecsv',
             'stis_g430l':'/STIS/odlm41010_sx1.fits',
             'phoenix':'/optical/unscaled_02560-5.00-0.0_phoenix_trappist-1.ecsv',
             'phoenix_wave' :'/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits' }


#COS
#G130M
g130m_data = sc.read_idl(path+filepaths['cos_g130m'])
w = g130m_data['w']
lya_edges = [1207, 1225] #lyman alpha region to remove
airglow = [1301, 1307, 1355, 1356] #oi airglow to remove
glowmask = (w <lya_edges[0])|(w > lya_edges[1])&(w <airglow[0])|(w >airglow[1])&(w <airglow[2])|(w >airglow[3])   #mask out lya, airglow
g130m_start, g130m_end, totals = sc.make_section(totals, g130m_data, mask =glowmask)

#G160M
g160m_data = sc.read_idl(path+filepaths['cos_g160m'])
w = g160m_data['w']
mask = (w > g130m_end)
g160m_start, g160m_end, totals = sc.make_section(totals, g160m_data, mask =mask)

#lya 
g140m_normfac = 1.0 #to do - calculate scales in spectra_combine
#lya_data = sc.read_lyamod(path+filepaths['lya'])
#lya_start, lya_end,totals = sc.make_section(totals, lya_data, normfac=g140m_normfac)

#G140M
g140m_data = sc.read_ecsv(path+filepaths['stis_g140m'])
w = g140m_data['w']
mask = (w >lya_edges[0])&(w <lya_start)|(w >lya_end)&(w < lya_edges[1]) #include just bits not co
g140m_start, g140m_end, totals = sc.make_section(totals, g140m_data, mask =mask, normfac=g140m_normfac)

#G230L

#G430L 
g430l_data = sc.read_stis_ccd(path+filepaths['stis_g430l'])
w = g430l_data['w']
g430l_clip = [0,-1]
mask = (w > g230l_end)
g430l_start, g430l_end, totals = sc.make_section(totals, g430l_data, mask =mask, clip = g430l_clip)

#phoenix
phoenix_normfac = 3.19143165e-26
phx_data = sc.read_phoenix(path+filepaths['phoenix'], path+filepaths['phoenix_wave'])
w = phx_data['w']
mask = (w > g430l_end)
phx_start, phx_end, totals = sc.make_section(totals, phx_data, mask=mask, normfac=phoenix_normfac)

#xmm and apec
xmm_data, apec_data = sc.read_xmm(path+filepaths['xmm'])
xmm_start, xmm_end, totals = sc.make_section(totals, xmm_data)
w = apec_data['w']
mask = (w > xmm_end)
apec_start, apec_end, totals = sc.make_section(totals, apec_data, mask=mask)

#euv
lya_flux = 2.06e-12
lya_flux *= g140m_normfac
distance = 4.54 #pc
euv_data = sc.make_euv(lya_flux, distance)
w = euv_data['w']
mask = (w > apec_end) & (w < g130m_start)
euv_start, euv_end, totals = sc.make_section(totals, euv_data, mask=mask)

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Wavelength (\AA)', size=20)
plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)', size=20)
plt.axhline(0, ls='--', c='k')

totals = sc.sort_totals(totals)

#sc.save_to_ecsv(totals, names, star, 'v2')

plt.show()