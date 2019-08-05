"""
@verison: 1

@author: David Wilson

@date 20190805

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
from astropy.convolution import convolve, Box1DKernel
import prepare_cos 
import prepare_stis
from astropy.units import cds
from scipy.io import readsav
cds.enable()


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

def build_cos_fuv(cospath, airglow):
    """
    cospath is a path to where the output from prepare_cos are stored. 
    Airglow is a list of airglow regions to mask out (inculding the Lyman alpha). Defined by visual inspection of each spectrum.
    """
    
    instrument_list = [] #starts a running count of all instruments
    g130m_path = glob.glob(cospath+'*g130m*.ecsv')
    g160m_path = glob.glob(cospath+'*g160m*.ecsv')
    g130m = Table.read(g130m_path[0])
    instrument_code, g130m = hst_instrument_column(g130m)
    g130m = normfac_column(g130m)
    instrument_list.append(instrument_code)
    airglow_mask = mask_maker(g130m ['WAVELENGTH'], airglow)
    sed_table = g130m[airglow_mask] #start the full SED table
    
    if len(g160m_path) > 0:
        g160m = Table.read(g160m_path[0])
        instrument_code, g160m = hst_instrument_column(g160m)
        instrument_list.append(instrument_code)
        g130m = normfac_column(g130m)
        g160m = g160m[g160m['WAVELENGTH'] > sed_table['WAVELENGTH'][-1]] #cut off everything covered by g130m
        sed_table = vstack([sed_table, g160m], metadata_conflicts = 'silent')
    
    return sed_table, instrument_list #sed table is the main thing.

def find_stis_normfac(component_repo, airglow):
    """
    Finds the normaliastion factor between the COS FUV data and the STIS G140L spectrum, if present
    """
    g140l = Table.read(glob.glob(component_repo+'*g140l*.ecsv')[0])
    g130m = Table.read(glob.glob(component_repo+'*g130m*.ecsv')[0])
    cw, cf, cdq = g130m['WAVELENGTH'], g130m['FLUX'], g130m['DQ']
    sw, sf, sdq = g140l['WAVELENGTH'], g140l['FLUX'], g140l['DQ']
    c_mask = (cw >= sw[0]) & (cw <=sw[-1]) & (cdq == 0)
    s_mask = (sw >= cw[0]) & (sw <=cw[-1]) & (sdq == 0)#mask to same same waveband and cut dq flags
    cw, cf, sw, sf = cw[c_mask], cf[c_mask], sw[s_mask], sf[s_mask]
    cw1, cf1 = resample.bintogrid(cw, cf, newx=sw) #rebin to stis wavelength grid
    stis_airglow_mask = mask_maker(sw, airglow)
    cos_airglow_mask = mask_maker(cw1, airglow)
    sw, sf, cw1, cf1 = sw[stis_airglow_mask], sf[stis_airglow_mask], cw1[cos_airglow_mask], cf1[cos_airglow_mask] #remove airglow
    c_int = np.trapz(cf1,cw1)
    s_int =  np.trapz(sf,sw)
    normfac = c_int/s_int
    print('STIS normfac = ', normfac)
    return normfac

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
    
    inst_code = instruments.getinsti('mod_lya_young')
    inst_array = np.full(table_length, inst_code, dtype=int)
    table['INSTRUMENT'] = inst_array 
    
    norm_array = np.full(len(table['WAVELENGTH']), table.meta['NORMFAC'], dtype =float)
    table['NORMFAC'] = norm_array
    return inst_code, table

    #'ERROR':e_new*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq_new,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD keep incase units are an issue
    

def add_stis_and_lya(sed_table, component_repo, lya_range, instrument_list, other_airglow):
    """
    Add the stis spectra and lya modelto the sed
    """
    g140l_path = glob.glob(component_repo+'*g140l*.ecsv')
    g140m_path = glob.glob(component_repo+'*g140m*.ecsv')
    lya = Table.read(glob.glob(component_repo+'*lya*.ecsv')[0])
    instrument_code, lya = fill_model(lya, 'mod_lya_young')
    instrument_list.append(instrument_code)
    
    if len(g140m_path) > 0:
        g140m = Table.read(g140m_path[0])
        instrument_code, g140m = hst_instrument_column(g140m)
        instrument_list.append(instrument_code)
        g140m = normfac_column(g140m)
        g140m_mask = (g140m['WAVELENGTH'] > lya_range[0]) & (g140m['WAVELENGTH'] < lya['WAVELENGTH'][0]) | (g140m['WAVELENGTH'] > lya['WAVELENGTH'][0]) & (g140m['WAVELENGTH'] > lya_range[0])
        g140m = g140m[g140m_mask]
        g140m['FLUX'] = g140m['FLUX'] * g140m.meta['NORMFAC']
        
        sed_table = vstack([sed_table, g140m], metadata_conflicts = 'silent')
    
    if len(g140l_path) > 0:
        g140l = Table.read(g140l_path[0])
        instrument_code, g140l = hst_instrument_column(g140l)
        instrument_list.append(instrument_code)
        g140l = normfac_column(g140l)
        g140l_mask = mask_maker(g140l['WAVELENGTH'], other_airglow, include=False) #fill in airglow gaps
        g140l_mask |= (g140l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])) #add in everything beyond COS range
        g140l = g140l[g140l_mask]
        sed_table = vstack([sed_table, g140l], metadata_conflicts = 'silent')
        
    return sed_table, instrument_list

    
                                                               
        
        
        

                     
    