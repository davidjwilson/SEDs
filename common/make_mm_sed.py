"""
@verison: 3

@author: David Wilson

@date 20202811

The big one. Draft here, will spin off to modules as required. 

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
from scipy.interpolate import interp1d
from astropy.convolution import convolve, Box1DKernel, convolve_fft, Gaussian1DKernel
import prepare_cos 
import prepare_stis
from astropy.units import cds
from scipy.io import readsav
from scipy.optimize import leastsq
from scipy.signal import argrelmax
from scipy.integrate import quad
cds.enable()


def smear(w,f, R, w_sample=1):
    '''
    Smears a model spectrum with a gaussian kernel to the given resolution, R.
    Adapeted from https://github.com/spacetelescope/pysynphot/issues/78

    Parameters
    -----------

    w,f:  spectrum to smear

    R: int
        The resolution (dL/L) to smear to

    w_sample: int
        Oversampling factor for smoothing
photometry
    Returns
    -----------

    sp: PySynphot Source Spectrum
        The smeared spectrum
    '''

    # Save original wavelength grid and units
    w_grid = w
    

    # Generate logarithmic wavelength grid for smoothing
    w_logmin = np.log10(np.nanmin(w_grid))
    w_logmax = np.log10(np.nanmax(w_grid))
    n_w = np.size(w_grid)*w_sample
    w_log = np.logspace(w_logmin, w_logmax, num=n_w)

    # Find stddev of Gaussian kernel for smoothing
    R_grid = (w_log[1:-1]+w_log[0:-2])/(w_log[1:-1]-w_log[0:-2])/2
    sigma = np.median(R_grid)/R
    if sigma < 1:
        sigma = 1

    # Interpolate on logarithmic grid
    f_log = np.interp(w_log, w_grid, f)

    # Smooth convolving with Gaussian kernel
    gauss = Gaussian1DKernel(stddev=sigma)
    f_conv = convolve_fft(f_log, gauss)

    # Interpolate back on original wavelength grid
    f_sm = np.interp(w_grid, w_log, f_conv)

    # Write smoothed spectrum back into Spectrum object
    return w_grid, f_sm


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

def hst_instrument_column(table):
    """
    Builds an instrument column and adds it to data. For HST data.
    """
    telescope, instrument, grating = table.meta['TELESCOP'], table.meta['INSTRUME'], table.meta['GRATING']
    if instrument == 'STIS':
        instrument = 'sts'
    inst_string = '%s_%s_%s' % (telescope.lower(), instrument.lower(), grating.lower())
    inst_code = instruments.getinsti(inst_string)
    inst_array = np.full(len(table['WAVELENGTH']), inst_code, dtype=int)
    table['INSTRUMENT'] = inst_array
    return inst_code, table

def normfac_column(table):
    """
    Adds a normfac column to data
    """
    norm_array = np.full(len(table['WAVELENGTH']), table.meta['NORMFAC'], dtype =float)
    table['NORMFAC'] = norm_array
    return table

def add_cos(cospath, airglow):
    """
    cospath is a path to where the output from prepare_cos are stored. 
    Airglow is a list of airglow regions to mask out (inculding the Lyman alpha). Defined by visual inspection of each spectrum.
    """
    
    instrument_list = [] #starts a running count of all instruments
    g130m_path = glob.glob(cospath+'*g130m*.ecsv')
    g160m_path = glob.glob(cospath+'*g160m*.ecsv')
    g230l_path = glob.glob(cospath+'*cos*g230l*.ecsv')
    
    if len(g130m_path) == 1:
    
        g130m = Table.read(g130m_path[0])
        instrument_code, g130m = hst_instrument_column(g130m)
        g130m = normfac_column(g130m)
        instrument_list.append(instrument_code)
        airglow_mask = mask_maker(g130m ['WAVELENGTH'], airglow)
        sed_table = g130m[airglow_mask] #start the full SED table
    else: 
        sed_table = {'Message':'nothing here yet, this star uses E140M'}
    
    if len(g160m_path) > 0:
        g160m = Table.read(g160m_path[0])
        instrument_code, g160m = hst_instrument_column(g160m)
        instrument_list.append(instrument_code)
        g130m = normfac_column(g130m)
        g160m = g160m[g160m['WAVELENGTH'] > sed_table['WAVELENGTH'][-1]] #cut off everything covered by g130m
        sed_table = vstack([sed_table, g160m], metadata_conflicts = 'silent')
      
    if len(g230l_path) > 0:
        g230l = Table.read(g230l_path[0])
        gap_edges = [2085.0, 2777.0]
        gap_mask = mask_maker(g230l['WAVELENGTH'], gap_edges)
        g230l = g230l[gap_mask] 
        instrument_code, g230l = hst_instrument_column(g230l)
        instrument_list.append(instrument_code)
        g230l = normfac_column(g230l)
        g230l = g230l[g230l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        sed_table = vstack([sed_table, g230l], metadata_conflicts = 'silent')
        sed_table, instrument_list = fill_cos_airglow(sed_table, gap_edges, instrument_list, nuv=True)

        
    
    return sed_table, instrument_list #sed table is the main thing.

def fill_cos_airglow(sed_table, airglow, instrument_list, nuv = False):
    """
    Fills in the gaps in cos airglow if stis spectra are unavailable. Fits to specta 5A on either side. If nuv =True then it instead fills the gap in the NUV spectrum, which requires different treatment
    """
    if nuv:
        b, r = airglow[0], airglow[1]
        gap_w = np.arange(b, r, 1)
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
    fill_table = Table([gap_w*u.AA, gap_f*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX'], meta={'NORMFAC': 1.0})
    instrument_code, fill_table = fill_model(fill_table, 'mod_gap_fill-')
    sed_table = vstack([sed_table, fill_table], metadata_conflicts = 'silent')
    instrument_list.append(instrument_code)
    return sed_table, instrument_list
    
def update_norm(ecsv_file, fits_file, normfac):
    """
    Updates the normalisation factors in stis ecsv and fits files
    """
    t = Table.read(ecsv_file)
    t.meta['NORMFAC'] = normfac
    t.write(ecsv_file, overwrite=True)
    h = fits.open(fits_file)
    h[0].header['NORMFAC'] = normfac
    h.writeto(fits_file, overwrite=True)
    h.close()
    


def fill_model(table, model_name): 
    """
    Fills out the missing columns from a model ecsv
    """
    table_length = len(table['WAVELENGTH'])
    fill_zeros = np.zeros(table_length)
    extra_names = ['ERROR','EXPTIME','DQ','EXPSTART','EXPEND']
    #units = [*u.erg/u.s/u.cm**2/u.AA, *u.s,None, *cds.MJD, *cds.MJD]
    for i in range(len(extra_names)):
        table[extra_names[i]] = fill_zeros#*units[i]
    #[table[name]=fill_zeros*unit for name, unit in zip(extra_names, units)]
    
    inst_code = instruments.getinsti(model_name)
    inst_array = np.full(table_length, inst_code, dtype=int)
    table['INSTRUMENT'] = inst_array 
    
    norm_array = np.full(len(table['WAVELENGTH']), table.meta['NORMFAC'], dtype =float)
    table['NORMFAC'] = norm_array
    return inst_code, table

    #'ERROR':e_new*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq_new,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD keep incase units are an issue
    

def add_stis_and_lya(sed_table, component_repo, lya_range, instrument_list, other_airglow, norm=True, error_cut=True):
    """
    Add the stis fuv spectra and lya model to the sed
    """
    stis_gratings = ['E140M', 'G140M','G140L', 'G230L', 'G230LB', 'G430L']
   # g140l_path = glob.glob(component_repo+'*g140l*.ecsv')
   # g140m_path = glob.glob(component_repo+'*g140m*.ecsv')
    lya = dict(WAVELENGTH = [10000, 0]) #filler for the star that doesn't have a lya measurement 
    lya_path = glob.glob(component_repo+'*lya*.ecsv')
    if len(lya_path) == 1:
        lya = Table.read(lya_path[0])
        instrument_code, lya = fill_model(lya, 'mod_lya_young')
        instrument_list.append(instrument_code)
        lya = normfac_column(lya)
        
    normfac = 1.0
    for grating in stis_gratings:
        specpath = glob.glob('{}*{}_v*.ecsv'.format(component_repo, grating.lower()))
        
        if len(specpath) == 1:
           
                 
            data= Table.read(specpath[0])
            if data.meta['INSTRUME'] =='STIS':
                instrument_code, data = hst_instrument_column(data)
                instrument_list.append(instrument_code)
                if grating != 'E140M':
                    if norm:
                        normfac = find_normfac(sed_table, specpath[0], np.concatenate((lya_range, other_airglow)), normfac) #normalise to COS spectra then sequentially 
                    update_norm(specpath[0], '{}.fits'.format(specpath[0][:-5]), normfac)
                if grating == 'E140M':
                    mask = (data['WAVELENGTH'] > 1160) & (data['WAVELENGTH'] < lya['WAVELENGTH'][0]) | (data['WAVELENGTH'] > lya['WAVELENGTH'][-1]) 

                elif grating == 'G140M':
                    mask = (data['WAVELENGTH'] > lya_range[0]) & (data['WAVELENGTH'] < lya['WAVELENGTH'][0]) | (data['WAVELENGTH'] > lya['WAVELENGTH'][-1]) & (data['WAVELENGTH'] < lya_range[1])
                elif grating == 'G140L':
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
                if grating == 'E140M':
                    sed_table = data
                else:                
                    sed_table = vstack([sed_table, data], metadata_conflicts = 'silent')
                if len(lya_path) == 1:    #lya needs to be added after e140m  
                    sed_table = vstack([sed_table, lya], metadata_conflicts = 'silent')
                sed_table.sort(['WAVELENGTH'])

                
    return sed_table, instrument_list





    
def residuals(scale, f, mf):
    return f - mf/scale
    
def phoenix_norm(component_repo, star_params, plot=False): 
    """
    find the normalisation factor between the phoenix model and the stis ccd (ccd_path)
    """
    norm = Table.read(glob.glob(component_repo+'*phx*.ecsv')[0])
    radius, distance = star_params['radius'], star_params['distance']
    normfac = ((radius.to(u.cm)/distance.to(u.cm))**2).value
    print('PHOENIX NORMFAC =', normfac)
    update_norm(glob.glob(component_repo+'*phx*.ecsv')[0], glob.glob(component_repo+'*phx*.fits')[0], normfac)

    if plot:
        plt.figure(star+'_scaled')
        plt.plot(w_phx, f_phx*normfac)
        #plt.step(w,f, where='mid')
        plt.step(w1, f1, where='mid')
        plt.xlabel('Wavelength (\AA)', size=20)
        plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)', size=20)
        plt.xlim(2000, 6000)
        plt.yscale('log')
        plt.axvline(cut, c='r', ls='--')
        plt.tight_layout()
        plt.show()
    return normfac

def add_stis_optical(sed_table, component_repo, instrument_list):
    """
    Adds the G430L spectrum
    """
    g430l_path = glob.glob(component_repo+'*g430l*.ecsv')
    if len(g430l_path) > 0:
        g430l = Table.read(g430l_path[0])
        instrument_code, g430l = hst_instrument_column(g430l)
        instrument_list.append(instrument_code)
        g430l = normfac_column(g430l)
        g430l = g430l[g430l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        sed_table = vstack([sed_table, g430l], metadata_conflicts = 'silent')
    return sed_table, instrument_list

def add_phx_spectrum(sed_table, component_repo, instrument_list):
    """
    Adds the scaled phoenix spectrum
    """
    phx_path = glob.glob(component_repo+'*phx*.ecsv')
    if len(phx_path) > 0:
        phx = Table.read(phx_path[0])
        instrument_code, phx = fill_model(phx, 'mod_phx_-----')
        instrument_list.append(instrument_code)
        phx = normfac_column(phx)
        phx = phx[phx['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        phx['FLUX'] *= phx.meta['NORMFAC']
        sed_table = vstack([sed_table, phx], metadata_conflicts = 'silent')
    return sed_table, instrument_list
        
def add_phoenix_and_g430l(sed_table, component_repo, instrument_list, error_cut=True, scale=True):
    """
    Adds both the phoenix model and the g430l spectrum, triming the g430l spectrum by and error cut and filling in any gap with the phoenix model. 
    """
    phx_path = glob.glob(component_repo+'*phx*.ecsv')
    g430l_path = glob.glob(component_repo+'*g430l*.ecsv')
    if len(phx_path) == 1 and len(g430l_path) == 1:
        phx = Table.read(phx_path[0])
        instrument_code, phx = fill_model(phx, 'mod_phx_-----')
        instrument_list.append(instrument_code)
        phx = normfac_column(phx)
       # print(phx.meta['NORMFAC'])
        phx['FLUX'] *= phx.meta['NORMFAC']
        
        g430l = Table.read(g430l_path[0])
        instrument_code, g430l = hst_instrument_column(g430l)
        instrument_list.append(instrument_code)
        
        if error_cut: #cut region before a rolling 30pt mean SN > 1
            bin_width = 30
            w, f, e = g430l['WAVELENGTH'], g430l['FLUX'], g430l['ERROR']
            sn = np.array([np.mean(f[i:i+bin_width]/e[i:i+bin_width]) for i in range(len(w[:-bin_width]))])
            start = w[:-bin_width][np.where(sn > 1)[0][0]]
            mask = (w > start) & (f > 0)
            g430l = g430l[mask]
        
        if scale: #scale g430l spectrum to the phoenix spectrum
            mask = (phx['WAVELENGTH'] >= g430l['WAVELENGTH'][0]) & (phx['WAVELENGTH'] <= g430l['WAVELENGTH'][-1]) 
            mw, mf = phx['WAVELENGTH'][mask], phx['FLUX'][mask]
            mw, mf = smear(mw, mf, 1000)
            pfr = interp1d(mw, mf, fill_value='extrapolate')(g430l['WAVELENGTH'])
            normfac = leastsq(residuals, 1., args=(g430l['FLUX'], pfr))[0]
            g430l['FLUX'] *= normfac
            g430l['ERROR'] *= normfac
            update_norm(g430l_path[0], '{}.fits'.format(g430l_path[0][:-5]), normfac[0])
        
        g430l = normfac_column(g430l)
        g430l = g430l[g430l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        phx = phx[(phx['WAVELENGTH'] > max(sed_table['WAVELENGTH'])) & (phx['WAVELENGTH'] < min(g430l['WAVELENGTH'])) | (phx['WAVELENGTH'] > max(g430l['WAVELENGTH']))]

        sed_table = vstack([sed_table, g430l], metadata_conflicts = 'silent')
        sed_table = vstack([sed_table, phx], metadata_conflicts = 'silent')
        
    return sed_table, instrument_list
        
    
    
def add_xray_spectrum(sed_table, component_repo, instrument_list, scope, add_apec = True, find_gap=True):
    """
    Adds either a Chandra or and XMM spectrum and an APEC model. Can also return the gap that the EUV/DEM will fit into.
    """
    if scope == 'xmm':
        instrument_name = 'xmm_epc_multi'
    if scope == 'cxo':
        instrument_name = 'cxo_acs_-----'
    cos_start = min(sed_table['WAVELENGTH']) #save the lowest wavelength on the table before we add anything to it
    xray_path = glob.glob(component_repo+'*'+scope+'*.ecsv')
    xray_end = 0
    if len(xray_path) > 0:
        xray = Table.read(xray_path[0])
        error = xray['ERROR']
        instrument_code, xray = fill_model(xray, instrument_name)
        xray['ERROR'] = error
        instrument_list.append(instrument_code)
        #xray = normfac_column(xray)
        xray_end = max(xray['WAVELENGTH'])
        sed_table = vstack([sed_table, xray], metadata_conflicts = 'silent')
    if add_apec:
        apec_path = glob.glob(component_repo+'*apec*.ecsv')
        if len(apec_path) > 0:
#             print(apec_path)
            apec = Table.read(apec_path[0])
            instrument_code, apec = fill_model(apec, 'mod_apc_-----')
            instrument_list.append(instrument_code)
            apec = normfac_column(apec)
            apec = apec[apec['WAVELENGTH'] > xray_end]
            xray_end = max(apec['WAVELENGTH'])
            sed_table = vstack([sed_table, apec], metadata_conflicts = 'silent')
    if find_gap:
        return sed_table, instrument_list, [xray_end, cos_start]
    else:
        return sed_table, instrument_list
    
def add_euv(sed_table, component_repo, instrument_list, euv_gap, euv_type):
    """
    Add the euv portion of the spectrum, either a Linsky_14 estmation or a DEM.
    """
    instrument_name = 'mod_euv_young'
    if euv_type == 'dem':
        instrument_name = 'mod_dem_-----'
    euv_path = glob.glob(component_repo+'*'+euv_type+'*.ecsv')
    if len(euv_path) > 0:
        euv = Table.read(euv_path[0])
        instrument_code, euv = fill_model(euv, instrument_name)
        instrument_list.append(instrument_code)
        euv = normfac_column(euv)
        euv = euv[(euv['WAVELENGTH'] > euv_gap[0]) & (euv['WAVELENGTH'] < euv_gap[1])]
        sed_table = vstack([sed_table, euv], metadata_conflicts = 'silent')
    return sed_table, instrument_list

def blackbody_fit(phx, Teff):
    """Return a function that is a blackbody fit to the phoenix spectrum for the star. The fit is to the unnormalized
    phoenix spectrum, so the fit function values must be multiplied by the appropriate normalization factor to match
    the normalized spectrum. From PLs code"""

    
    # recursively identify relative maxima until there are fewer than N points
    N = 10000
    keep = np.arange(len(phx))
    while len(keep) > N:
        temp, = argrelmax(phx['FLUX'][keep])
        keep = keep[temp]

    efac = const.h * const.c / const.k_B / (Teff * u.K)
    efac  = efac.to(u.angstrom).value
    w = phx['WAVELENGTH']
    w = w[keep]
    planck_shape = 1.0/w**5/(np.exp(efac/w) - 1)
    y = phx['FLUX'][keep]

    Sfy = np.sum(planck_shape * y)
    Sff = np.sum(planck_shape**2)

    norm = Sfy/Sff

    return lambda w: norm/w**5/(np.exp(efac/w) - 1)

def bolo_integral(pan,phx,teff,uplim=np.inf, tail=False):
    """
    Calculates the bolometric integral flux of the SED by adding a blackbody fit to the end of the sed
    """
    fit_unnormed = blackbody_fit(phx, teff)
    normfac = pan[-1]['NORMFAC'] 
    #Ibody = flux_integral(pan)[0]
    Ibody = np.trapz(pan['FLUX'], pan['WAVELENGTH'])
    if tail:
        Itail = normfac*quad(fit_unnormed, pan['WAVELENGTH'][-1], uplim)[0]
        I = Ibody + Itail
    else:
        I = Ibody
#     print(I)

    return I


def add_bolometric_flux(sed_table, component_repo, star_params):
    """
    Creates and adds the bolometric flux column to the sed
    """
    phx = Table.read(glob.glob(component_repo+'*phx*ecsv')[0])
    bolo_int = bolo_integral(sed_table,phx,star_params['Teff'])
    sed_table['BOLOFLUX'] = (sed_table['FLUX']/bolo_int)*1/u.A
    sed_table['BOLOERR'] = (sed_table['ERROR']/bolo_int)*1/u.A
    #boloerr = np.zeros(len(sed_table['ERROR']))
    #for i in range(len(boloerr)):
     #   if sed_table['ERROR'][i] > 0.0:
      #      boloerr[i] = sed_table['ERROR'][i]/boloflux
    return sed_table