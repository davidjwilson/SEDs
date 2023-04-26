"""
@verison: 6

@author: David Wilson

@date 20230223

The big one. Draft here, will spin off to modules as required. 
v6 updating to use fits files, not ecsv

"""

import instruments
import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from scipy.io.idl import readsav
from astropy.table import Table, vstack
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from astropy.modeling import models, fitting
from craftroom import resample
from scipy.interpolate import interpolate
from astropy.convolution import convolve, Box1DKernel, convolve_fft, Gaussian1DKernel
from astropy.units import cds
from scipy.io import readsav
from scipy.optimize import leastsq
from scipy.signal import argrelmax
from scipy.integrate import quad
import math as mt
from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler
from astropy.nddata import StdDevUncertainty
import remove_negatives as negs
import bin_to_const as bin1A
cds.enable()


def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff0 = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    diff1 = np.concatenate((diff, np.array([diff[-1]]))) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff0/2.
    w1 = w + diff1/2.
    return w0, w1

def mask_maker(x, pairs, include=True):
    """
    creates a mask for a spectrum that excudes between pairs from an array
    """
    b = pairs[::2]
    r = pairs[1::2]
    C = np.zeros_like(x,dtype='bool')
    for i in range(len(b)):
        C |= (x>b[i])&(x<r[i])
    if include:
        return ~C
    else:
        return C

def hst_instrument_column(table, hdr):
    """
    Builds an instrument column and adds it to data. For HST data.
    """
    telescope, instrument, grating = hdr['TELESCOP'], hdr['INSTRUME'], hdr['GRATING']
    if instrument == 'STIS':
        instrument = 'sts'
    inst_string = '%s_%s_%s' % (telescope.lower(), instrument.lower(), grating.lower())
    inst_code = instruments.getinsti(inst_string)
    inst_array = np.full(len(table['WAVELENGTH']), inst_code, dtype=int)
    table['INSTRUMENT'] = inst_array
    return inst_code, table

def normfac_column(table, hdr):
    """
    Adds a normfac column to data
    """
    norm_array = np.full(len(table['WAVELENGTH']), hdr['NORMFAC'], dtype =float)
    table['NORMFAC'] = norm_array
    return table

def add_cos(cospath, airglow, remove_negs=False, to_1A=False):
    """
    cospath is a path to where the output from prepare_cos are stored. 
    Airglow is a list of airglow regions to mask out (inculding the Lyman alpha). Defined by visual inspection of each spectrum.
    """
    
    instrument_list = [] #starts a running count of all instruments
    g130m_path = glob.glob(cospath+'*g130m*.fits')
    g160m_path = glob.glob(cospath+'*g160m*.fits')
    g230l_path = glob.glob(cospath+'*cos*g230l*.fits')
    
    if len(g130m_path) == 1:
    
        g130m = Table(fits.getdata(g130m_path[0], 1))
        hdr = fits.getheader(g130m_path[0], 0)
        if remove_negs:
            print('removing negatives from {}'.format(g130m_path[0]))
            g130m = negs.make_clean_spectrum(g130m)
        if to_1A:
            print('binning {}'.format(g130m_path[0]))
            g130m = bin1A.spectrum_to_const_res(g130m)
        
        instrument_code, g130m = hst_instrument_column(g130m, hdr)
        g130m = normfac_column(g130m, hdr)
        instrument_list.append(instrument_code)
        airglow_mask = mask_maker(g130m ['WAVELENGTH'], airglow)
        sed_table = g130m[airglow_mask] #start the full SED table
        sed_table.meta = dict(hdr)
    else: 
        sed_table = {'Message':'nothing here yet, this star uses E140M'}
    
    if len(g160m_path) > 0:
        g160m = Table(fits.getdata(g160m_path[0], 1))
        hdr = fits.getheader(g160m_path[0], 0)
        if remove_negs:
            print('removing negatives from {}'.format(g160m_path[0]))
            g160m = negs.make_clean_spectrum(g160m)
        if to_1A:
            print('binning {}'.format(g160m_path[0]))
            g160m = bin1A.spectrum_to_const_res(g160m)
       
        instrument_code, g160m = hst_instrument_column(g160m,  hdr)
        instrument_list.append(instrument_code)
        g160m = normfac_column(g160m, hdr)
        g160m = g160m[g160m['WAVELENGTH'] > sed_table['WAVELENGTH'][-1]] #cut off everything covered by g130m
        sed_table = vstack([sed_table, g160m], metadata_conflicts = 'silent')
      
    if len(g230l_path) > 0:
        g230l = Table(fits.getdata(g230l_path[0], 1))
        hdr = fits.getheader(g230l_path[0], 0)
        if remove_negs:
            print('removing negatives from {}'.format(g230l_path[0]))
            g230l = negs.make_clean_spectrum(g230l)
        if to_1A:
            print('binning {}'.format(g230l_path[0]))
            g230l = bin1A.spectrum_to_const_res(g230l)
        
        gap_edges = [2085.0, 2777.0]
        gap_mask = mask_maker(g230l['WAVELENGTH'], gap_edges)
        g230l = g230l[gap_mask] 
        instrument_code, g230l = hst_instrument_column(g230l, hdr)
        instrument_list.append(instrument_code)
        g230l = normfac_column(g230l, hdr)
        g230l = g230l[g230l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        sed_table = vstack([sed_table, g230l], metadata_conflicts = 'silent')
        sed_table, instrument_list = fill_cos_airglow(sed_table, gap_edges, instrument_list, hdr, nuv=True)

        
    
    return sed_table, instrument_list #sed table is the main thing.

def fill_cos_airglow(sed_table, airglow, instrument_list, hdr, nuv = False):
    """
    Fills in the gaps in cos airglow if stis spectra are unavailable. Fits to specta 5A on either side. If nuv =True then it instead fills the gap in the NUV spectrum, which requires different treatment
    """
    if nuv:
        b, r = airglow[0], airglow[1]
        gap_w = np.arange(b, r+1, 1)
        w, f = sed_table['WAVELENGTH'], sed_table['FLUX'] 
        mask = (w > 1700) & (w < 2790) | (w > 2805) & (w < 3150) #cut to nuv range and remove mg ii
        w, f = w[mask], f[mask]
        gap_f = np.polyval((np.polyfit(w,f,1)), gap_w)
#         print(np.mean(gap_f))
    else:
        b = airglow[::2]
        r = airglow[1::2]
        gap_w = np.array([], dtype=float)
        gap_f = np.array([], dtype=float)
        for i in range(len(b)):
            mask = (sed_table['WAVELENGTH'] > b[i] - 5) & (sed_table['WAVELENGTH'] < r[i] + 5)
            wi = np.arange(b[i], r[i], 1.0)
            gap_w = np.concatenate((gap_w, wi))
            fi = np.polyval((np.polyfit(sed_table['WAVELENGTH'][mask], sed_table['FLUX'][mask], 2)), wi)
            gap_f = np.concatenate((gap_f, fi))
    w0, w1 = wavelength_edges(gap_w)
    fill_table = Table([gap_w*u.AA, w0*u.AA, w1*u.AA, gap_f*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'WAVELENGTH0', 'WAVELENGTH1','FLUX'], meta={'NORMFAC': 1.0})
    instrument_code, fill_table = fill_model(fill_table, 'mod_gap_fill-', hdr)
    sed_table = vstack([sed_table, fill_table], metadata_conflicts = 'silent')
    instrument_list.append(instrument_code)
    return sed_table, instrument_list
    
def update_norm(fits_file, normfac):
    """
    Updates the normalisation factors in stis ecsv and fits files
    """
    h = fits.open(fits_file)
    h[0].header['NORMFAC'] = normfac
    h.writeto(fits_file, overwrite=True)
    h.close()

def fill_model(table, model_name, hdr): 
    """
    Fills out the missing columns from a model ecsv
    """
    table_length = len(table['WAVELENGTH'])
    fill_zeros = np.zeros(table_length)
    extra_names = ['EXPTIME','DQ','EXPSTART','EXPEND']
    for i in range(len(extra_names)):
        table[extra_names[i]] = fill_zeros#*units[i]
    inst_code = instruments.getinsti(model_name)
    inst_array = np.full(table_length, inst_code, dtype=int)
    table['INSTRUMENT'] = inst_array    
    norm_array = np.full(len(table['WAVELENGTH']), hdr['NORMFAC'], dtype =float)
    table['NORMFAC'] = norm_array
    return inst_code, table

def add_stis_and_lya(sed_table, component_repo, lya_range, instrument_list, other_airglow, norm=False, error_cut=True, optical = False, remove_negs=False, to_1A=False):
    """
    Add the stis fuv spectra and lya model to the sed
    """
    stis_gratings = ['E140M', 'G140M','G140L', 'G230L', 'G230LB']
    if optical:
        stis_gratings.append('G430L') #usually want to add the optical spectrum with the phoenix model, but retaining the option here
   # g140l_path = glob.glob(component_repo+'*g140l*.ecsv')
   # g140m_path = glob.glob(component_repo+'*g140m*.ecsv')
    lya = dict(WAVELENGTH = [10000, 0]) #filler for the star that doesn't have a lya measurement 
    lya_path = glob.glob(component_repo+'*lya*.fits')
    if len(lya_path) == 1:
        lya = Table(fits.getdata(lya_path[0], 1))
        hdr = fits.getheader(lya_path[0], 0)
        if to_1A:
            print('binning {}'.format(lya_path[0]))
            lya = bin1A.spectrum_to_const_res(lya)
        instrument_code, lya = fill_model(lya, 'mod_lya_young', hdr)
        instrument_list.append(instrument_code)
        lya = normfac_column(lya, hdr)        
    normfac = 1.0
    uses_e140m = False #if nether present fill in COS airglow with a polynomial
    used_g140l = False
    for grating in stis_gratings:
        specpath = glob.glob('{}*{}_v*.fits'.format(component_repo, grating.lower()))  
        if len(specpath) == 1:
            data= Table(fits.getdata(specpath[0], 1))
            hdr = fits.getheader(specpath[0], 0)
            if hdr['INSTRUME'] =='STIS':
                if remove_negs:
                    print('removing negatives from {}'.format(specpath))
                    data = negs.make_clean_spectrum(data)
                if to_1A:
                    print('binning {}'.format(specpath))
                    data = bin1A.spectrum_to_const_res(data)
                instrument_code, data = hst_instrument_column(data,  hdr)
                instrument_list.append(instrument_code)
                if grating != 'E140M':
                    if norm:
                        normfac = find_normfac(sed_table, specpath[0], np.concatenate((lya_range, other_airglow)), normfac) 
                    update_norm(specpath[0], normfac)
                if grating == 'E140M':
                    uses_e140m = True
                    mask = (data['WAVELENGTH'] > 1160) & (data['WAVELENGTH'] < lya['WAVELENGTH'][0]) | (data['WAVELENGTH'] > lya['WAVELENGTH'][-1]) 

                elif grating == 'G140M':
                    mask = (data['WAVELENGTH'] > lya_range[0]) & (data['WAVELENGTH'] < lya['WAVELENGTH'][0]) | (data['WAVELENGTH'] > lya['WAVELENGTH'][-1]) & (data['WAVELENGTH'] < lya_range[1])
                elif grating == 'G140L':
                    used_g140l = True
                    mask = mask_maker(data['WAVELENGTH'], other_airglow, include=False) #fill in airglow gaps
                    mask |= (data['WAVELENGTH'] > max(sed_table['WAVELENGTH']))
                    
                elif grating == 'G430L':
                    if error_cut: #cut region before a rolling 30pt mean SN > 1
                        bin_width = 30
                        w, f, e = data['WAVELENGTH'], data['FLUX'], data['ERROR']
                        sn = np.array([np.mean(f[i:i+bin_width]/e[i:i+bin_width]) for i in range(len(w[:-bin_width]))])
                        start = w[:-bin_width][np.where(sn > 1)[0][0]]
                        mask = (w > start) & (f > 0)
                        data = data[mask]
                    mask = (data['WAVELENGTH'] > max(sed_table['WAVELENGTH']))
                else:
                    mask = (data['WAVELENGTH'] > max(sed_table['WAVELENGTH']))
                data = data[mask]
                if norm:
                    data['FLUX'] = data['FLUX'] * normfac
                    data['ERROR'] = data['ERROR'] * normfac
                data = normfac_column(data, hdr)
                if grating == 'E140M':
                    sed_table = data
                    sed_table.meta = dict(hdr)
                else:                
                    sed_table = vstack([sed_table, data], metadata_conflicts = 'silent')
    if len(lya_path) == 1:    #lya needs to be added after e140m  
        print('adding lya')
        sed_table = vstack([sed_table, lya], metadata_conflicts = 'silent')
    sed_table.sort(['WAVELENGTH'])

    if  uses_e140m == False and used_g140l == False:
        print('filling COS airglow with polynomials')
        sed_table, instrument_list = fill_cos_airglow(sed_table, other_airglow, instrument_list, hdr)

                
    return sed_table, instrument_list





    
def residuals(scale, f, mf):
    return f - mf/scale

        
def add_phoenix_and_g430l(sed_table, component_repo, instrument_list, error_cut=True, scale=False, remove_negs=False, to_1A=False):
    """
    Adds both the phoenix model and the g430l spectrum, triming the g430l spectrum by and error cut and filling in any gap with the phoenix model. 
    """
    phx_path = glob.glob(component_repo+'*phx*.fits')
    g430l_path = glob.glob(component_repo+'*g430l*.fits')
    if len(phx_path) == 1 and len(g430l_path) == 1:
        phx = Table(fits.getdata(phx_path[0], 1))
        hdr = fits.getheader(phx_path[0], 0)
        if to_1A:
            print('binning {}'.format(phx_path[0]))
            phx = bin1A.spectrum_to_const_res(phx)
        instrument_code, phx = fill_model(phx, 'mod_phx_-----', hdr)
        instrument_list.append(instrument_code)
        phx = normfac_column(phx, hdr)
       # print(phx.meta['NORMFAC'])
        phx['FLUX'] *= hdr['NORMFAC']
        phx['ERROR'] *= hdr['NORMFAC']
        
        g430l = Table(fits.getdata(g430l_path[0], 1))
        hdr = fits.getheader(g430l_path[0], 0)
        if error_cut: #cut region before a rolling 30pt mean SN > 1
            bin_width = 30
            w, f, e = g430l['WAVELENGTH'], g430l['FLUX'], g430l['ERROR']
            sn = np.array([np.mean(f[i:i+bin_width]/e[i:i+bin_width]) for i in range(len(w[:-bin_width]))])
            start = w[:-bin_width][np.where(sn > 1)[0][0]]
            mask = (w > start)
            g430l = g430l[mask]
        
        if remove_negs:
            print('removing negatives from {}'.format(g430l_path[0]))
            g430l = negs.make_clean_spectrum(g430l)
        if to_1A:
            print('binning {}'.format(g430l_path[0]))
            g430l = bin1A.spectrum_to_const_res(g430l)
        instrument_code, g430l = hst_instrument_column(g430l,  hdr)
        instrument_list.append(instrument_code)
        
       
     
        if scale: #scale g430l spectrum to the phoenix spectrum
            mask = (phx['WAVELENGTH'] >= g430l['WAVELENGTH'][0]) & (phx['WAVELENGTH'] <= g430l['WAVELENGTH'][-1]) 
            mw, mf = phx['WAVELENGTH'][mask], phx['FLUX'][mask]
            mw, mf = smear(mw, mf, 1000)
            pfr = interpolate.interp1d(mw, mf, fill_value='extrapolate')(g430l['WAVELENGTH'])
            normfac = leastsq(residuals, 1., args=(g430l['FLUX'], pfr))[0]
            g430l['FLUX'] *= normfac
            g430l['ERROR'] *= normfac
            update_norm(g430l_path[0], normfac[0])
        
        g430l = normfac_column(g430l, hdr)
        g430l = g430l[g430l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        # print('optical gap', max(sed_table['WAVELENGTH']), min(g430l['WAVELENGTH']))
        phx = phx[(phx['WAVELENGTH'] > max(sed_table['WAVELENGTH'])) & (phx['WAVELENGTH'] < min(g430l['WAVELENGTH'])) | (phx['WAVELENGTH'] > max(g430l['WAVELENGTH']))]

        sed_table = vstack([sed_table, g430l], metadata_conflicts = 'silent')
        sed_table = vstack([sed_table, phx], metadata_conflicts = 'silent')
        
    return sed_table, instrument_list
        
    
    
def add_xray_spectrum(sed_table, component_repo, instrument_list, scope, add_apec = True, find_gap=True, to_1A=False):
    """
    Adds either a Chandra or and XMM spectrum and an APEC model. Can also return the gap that the EUV/DEM will fit into.
    """
    if scope == 'xmm':
        instrument_name = 'xmm_epc_multi'
    if scope == 'cxo':
        instrument_name = 'cxo_acs_-----'
    cos_start = min(sed_table['WAVELENGTH']) #save the lowest wavelength on the table before we add anything to it
    xray_path = glob.glob(component_repo+'*'+scope+'*.fits')
    xray_end = 0
    if len(xray_path) > 0:
        xray = Table(fits.getdata(xray_path[0], 1))
        hdr = fits.getheader(xray_path[0], 0)
        if to_1A:
            print('binning {}'.format(xray_path[0]))
            xray = bin1A.spectrum_to_const_res(xray)
        error = xray['ERROR'] #retain error and exptime 
        exptime = xray['EXPTIME']
        instrument_code, xray = fill_model(xray, instrument_name, hdr)
        xray['ERROR'] = error
        xray['EXPTIME'] = exptime
        instrument_list.append(instrument_code)
        #xray = normfac_column(xray)
        xray_end = max(xray['WAVELENGTH'])
        sed_table = vstack([sed_table, xray], metadata_conflicts = 'silent')
    if add_apec:
        apec_path = glob.glob(component_repo+'*apec*.fits')
        if len(apec_path) > 0:
#             print(apec_path)
            apec = Table(fits.getdata(apec_path[0], 1))
            hdr = fits.getheader(apec_path[0], 0)
            if to_1A:
                print('binning {}'.format(apec_path[0]))
                apec = bin1A.spectrum_to_const_res(apec)
            instrument_code, apec = fill_model(apec, 'mod_apc_-----', hdr)
            instrument_list.append(instrument_code)
            apec = normfac_column(apec, hdr)
            apec = apec[apec['WAVELENGTH'] > xray_end]
            xray_end = max(apec['WAVELENGTH'])
            sed_table = vstack([sed_table, apec], metadata_conflicts = 'silent')
    if find_gap:
        return sed_table, instrument_list, [xray_end, cos_start]
    else:
        return sed_table, instrument_list
    
def add_euv(sed_table, component_repo, instrument_list, euv_gap, euv_type, to_1A=False):
    """
    Add the euv portion of the spectrum, either a Linsky_14 estmation or a DEM.
    """
    instrument_name = 'mod_euv_young'
    if euv_type == 'dem':
        instrument_name = 'mod_dem_-----'
    euv_path = glob.glob(component_repo+'*'+euv_type+'*.fits')
    if len(euv_path) > 0:
        euv = Table(fits.getdata(euv_path[0], 1))
        hdr = fits.getheader(euv_path[0], 0)
        if to_1A:
            print('binning {}'.format(euv_path[0]))
            euv = bin1A.spectrum_to_const_res(euv)
        instrument_code, euv = fill_model(euv, instrument_name, hdr)
        instrument_list.append(instrument_code)
        euv = normfac_column(euv, hdr)
        euv = euv[(euv['WAVELENGTH'] > euv_gap[0]) & (euv['WAVELENGTH'] < euv_gap[1])]
        sed_table = vstack([sed_table, euv], metadata_conflicts = 'silent')
    return sed_table, instrument_list






def add_bolometric_flux(sed_table, component_repo, star_params):
    """
    Creates and adds the bolometric flux column to the sed
    """
    # phx = Table(fits.getdata(glob.glob(component_repo+'*phx*fits')[0], 1))
    # print(sed_table['WAVELENGTH'][100])
    bolo_int = np.trapz(sed_table['FLUX'], sed_table['WAVELENGTH'])*(u.erg/u.s/u.cm**2)
    # print(bolo_int)
#     bolo_int = bolo_integral(sed_table,phx,star_params['Teff'])*(u.erg/u.s/u.cm**2)
    boloflux = (sed_table['FLUX']/bolo_int).value
    # print(sed_table['FLUX'][100])
    # print(boloflux)
    boloerr = (sed_table['ERROR']/bolo_int).value
    sed_table['BOLOFLUX'] = boloflux*(1/u.AA)
    sed_table['BOLOERR'] = boloerr*(1/u.AA)
    sed_table.meta['BOLOFLUX'] = bolo_int.value
    return sed_table

def sed_to_const_res(sed_table, res=1, start_cut=0, end_cut = 1e5):
    """
    Rebins an SED to a wavelength grid with a bin size of res, default = 1A
    """
    
    #wavelength 
    start, end= mt.ceil(sed_table['WAVELENGTH'][0]), mt.floor(sed_table['WAVELENGTH'][-1])
    if start < start_cut: #cut MM SEDs down to 1e5 A
        start = start_cut
    if end > end_cut:
        end = end_cut 
    
    new_wavelength = np.arange(start,end+res, res)
    new_w0 = new_wavelength - (0.5 * res)
    new_w1 = new_wavelength + (0.5 * res)
    
    #flux and error
    w, f, e= np.array(sed_table['WAVELENGTH']), np.array(sed_table['FLUX']), np.array(sed_table['ERROR'])
    
    #add an error array to models so that bintogrid works. Take it away later on
    model_instruments = []
    for i in range(len(w)):
        if e[i] == 0.0:
            model_instruments.append(sed_table['INSTRUMENT'][i])
            e[i]  = 0.1*f[i]
    print('length of new w:', len(new_wavelength))
#     cut = 5700 #spectutils struggles with large numbers, do top of spectrum separatly.
    mask = sed_table['INSTRUMENT'] != 131072 #cut of the phoenix spectrum and do it in craftroom
    cut = w[~mask][0] #where to cut the new wavelength grid
#     mask = ( w < cut)
#     print('here', e[mask])
    
    input_spec = Spectrum1D(spectral_axis=w[mask]*u.AA, flux=f[mask]*u.Unit('erg cm-2 s-1 AA-1') , uncertainty= StdDevUncertainty(e[mask]))
    fluxcon = FluxConservingResampler()
    new_spec_fluxcon = fluxcon(input_spec, new_wavelength[new_wavelength < cut]*u.AA)
                    
    new_wavelength, new_flux = resample.bintogrid(w[~mask], f[~mask], newx=new_wavelength[new_wavelength >= cut])
    
    new_error = np.concatenate(((1/new_spec_fluxcon.uncertainty.array**0.5), np.zeros(len(new_wavelength))))
    new_wavelength = np.concatenate((new_spec_fluxcon.spectral_axis.value, new_wavelength))
    new_flux = np.concatenate((new_spec_fluxcon.flux.value, new_flux))
    print(len(new_wavelength))

#     #flux
#     new_wavelength, new_flux = resample.bintogrid(sed_table['WAVELENGTH'], sed_table['FLUX'], newx=new_wavelength)
    
#     #error
#     new_error =  interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['ERROR'])(new_wavelength)# this doesn't work
    
    #exptime - linear extrapolation is similar to averaged to bin widths
    new_exptime = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPTIME'])(new_wavelength)
    
    #dq - interploate, then look for unusual values and correct them, summing if the values to either side are different.
   
    new_dq = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['DQ'], kind='previous')(new_wavelength)
    new_dq = new_dq.astype(int)
    
    #expstart - minumum expstart in each bin
    startups = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPSTART'], kind='next')(new_wavelength)
    startdowns = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPSTART'], kind='previous')(new_wavelength)
    new_expstart = np.min([startups, startdowns], axis=0)
    
    #expends - maximum expend in each bin
    endups = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPEND'], kind='next')(new_wavelength)
    enddowns = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['EXPEND'], kind='previous')(new_wavelength)
    new_expend = np.max([endups, enddowns], axis=0)
    
    #instrument - as dqs
    new_instrument = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['INSTRUMENT'], kind='previous')(new_wavelength)
    new_instrument = new_instrument.astype(int)
    
    #error, dq and instrument loop
    for i in range(len(new_wavelength)):
        if i < len(new_wavelength)-1:
            if new_dq[i] != new_dq[i+1]:
                new_dq[i] = new_dq[i] + new_dq[i+1]
            if new_instrument[i] != new_instrument[i+1]:
                new_instrument[i] = new_instrument[i] + new_instrument[i+1]
        if new_instrument[i] in instruments.getmodelcodes():
            new_error[i] = 0.0
#         if i > 0 and i < len(new_wavelength):
#             if new_instrument[i-1] in instruments.getmodelcodes() and if new_instrument[i+1] in instruments.getmodelcodes():
#                 new_error[i] = 0 #fudge for when two instruments next to each other
        if np.isnan(new_flux[i]) == True:
            print('yes')
            new_flux[i] = 0.0
        if np.isnan(new_error[i]) == True:
            new_error[i] = 0.0

    
    #normfac - linear extrapolation
    new_normfac = interpolate.interp1d(sed_table['WAVELENGTH'], sed_table['NORMFAC'])(new_wavelength)
    
    #boloflux -use original boloflux for consitency
    bolo_int = sed_table.meta['BOLOFLUX']*(u.erg/u.s/u.cm**2)
    new_boloflux = (new_flux/bolo_int).value
    new_boloerr = (new_error/bolo_int).value
    

    names = sed_table.dtype.names
    new_sed_table = Table([new_wavelength*u.AA, new_w0*u.AA, new_w1*u.AA, new_flux*u.erg/u.s/u.cm**2/u.AA, new_error*u.erg/u.s/u.cm**2/u.AA, new_exptime*u.s, 
                           new_dq,new_expstart*cds.MJD, new_expend*cds.MJD, new_instrument, new_normfac, new_boloflux*(1/u.AA), new_boloerr*(1/u.AA)], names=names, meta= sed_table.meta)
    return new_sed_table
          

