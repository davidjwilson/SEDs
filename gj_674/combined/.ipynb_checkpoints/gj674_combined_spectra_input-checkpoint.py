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



totals = sc.totals_setup()
star = 'gj_674' #as it appears in the filepath
path = '/home/david/work/muscles/SEDs/'+star

plt.figure(star+'_combined',figsize = (13, 7))
plt.subplots_adjust(top = 0.95, right = 0.99, left = 0.07, bottom = 0.11)

#file_locations
filepaths = {'xmm':'/xmm/GJ674.fits',
             'cos_g130m':'/COS/GJ674_COS130M_Mm1_NOSCL_03apr18.sav',
             'lya':'/lya/GJ674_intrinsic_LyA_profile.txt',
             'stis_g140m':'/STIS/GJ674_G140M_coadd.ecsv',
             'stis_g140l':'/STIS/GJ674_G140L_noflare_x1d.fits',
             'stis_g230l':'/STIS/GJ674_G230L_x1d.fits',
             'stis_g430l':'/STIS/odlm21010_sx1.fits',
             'phoenix':'/PHOENIX/lte03400-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits',
             'phoenix_wave' :'/PHOENIX/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits' }


#COS
#G130M
g130m_data = sc.read_idl(path+filepaths['cos_g130m'])
w = g130m_data['w']
lya_edges = [1207, 1225] #lyman alpha region to remove
airglow = [1304, 1304.5, 1355, 1356] #oi airglow to remove
glowmask = (w <lya_edges[0])|(w > lya_edges[1])&(w <airglow[0])|(w >airglow[1])&(w <airglow[2])|(w >airglow[3])   #mask out lya, airglow
g130m_start, g130m_end, totals = sc.make_section(totals, g130m_data, mask =glowmask)

#lya 
g140m_normfac = 2.9 #to do - calculate scales in spectra_combine
lya_data = sc.read_lyamod(path+filepaths['lya'])
lya_start, lya_end,totals = sc.make_section(totals, lya_data, normfac=g140m_normfac)

#G140M
g140m_data = sc.read_ecsv(path+filepaths['stis_g140m'])
w = g140m_data['w']
mask = (w >lya_edges[0])&(w <lya_start)|(w >lya_end)&(w < lya_edges[1]) #include just bits not co
g140m_start, g140m_end, totals = sc.make_section(totals, g140m_data, mask =mask, normfac=g140m_normfac)

#G140L
#includes all data blue-wards of the g130m and fills in the airglow gaps
g140l_normfac = 3.5
g140l_data = sc.read_stis_x1d(path+filepaths['stis_g140l'])
w = g140l_data['w']
mask = (w > airglow[0])&(w < airglow[1])|(w > airglow[2])&(w < airglow[3])|(w >g130m_end) #only need the bit not covered by cos, plus airglow filler
g140l_start, g140l_end, totals = sc.make_section(totals, g140l_data, mask =mask, normfac=g140l_normfac)

#G230L
g230l_normfac = 1.43
g230l_data = sc.read_stis_x1d(path+filepaths['stis_g230l'])
w = g230l_data['w']
g230l_clip = [0,-6]
mask = (w > g140l_end)
g230l_start, g230l_end, totals = sc.make_section(totals, g230l_data, mask =mask, normfac=g230l_normfac, clip = g230l_clip)

#G430L 
# normalised to PHOENIX
g430l_normfac = 1.11
g430l_data = sc.read_stis_ccd(path+filepaths['stis_g430l'])
w = g430l_data['w']
g430l_clip = [0,-1]
mask = (w > g230l_end)
g430l_start, g430l_end, totals = sc.make_section(totals, g430l_data, mask =mask, normfac=g430l_normfac, clip = g430l_clip)

#phoenix
phoenix_normfac = 3.19143165e-26
phx_data = sc.read_phoenix(path+filepaths['phoenix'], path+filepaths['phoenix_wave'])
phw = phx_data['w']
phf = phx_data['f']
mask = (phw > g430l_end)
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

#bolometric flux
teff = 3400*u.K
bolometric_flux, bf, be = sc.calculate_bolometric_flux(teff, totals, phw, phf, phoenix_normfac)


sc.save_to_ecsv(totals, bf, be, star, 'v3')

plt.show()