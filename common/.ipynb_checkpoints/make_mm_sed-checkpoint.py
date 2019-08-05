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
    g130m_path = glob.glob(cospath+'*g130m*.ecsv')[0]
    g160m_path = glob.glob(cospath+'*g160m*.ecsv')[0]
    g130m = Table.read(g130m_path)
    instrument_code, g130m = hst_instrument_column(g130m)
    g130m = normfac_column(g130m)
    instrument_list.append(instrument_code)
    airglow_mask = mask_maker(g130m ['WAVELENGTH'], airglow)
    sed_table = g130m[airglow_mask] #start the full SED table
    
    if len(g160m_path) > 0:
        g160m = Table.read(g160m_path)
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
    g140l = Table.read(glob.glob(componet_repo+'*g140l*.ecsv')[0])
    g130m = Table.read(glob.glob(componet_repo+'*g130m*.ecsv')[0])
    cw, cf, cdq = g130m['WAVELENGTH'], g130m['FLUX'], g130m['DQ']
    sw, sf, sdq = g140l['WAVELENGTH'], g140l['FLUX'], g140l['DQ']
    c_mask = (cw >= sw[0]) & (cw <=sw[-1]) & (cdq == 0)
    s_mask = (sw >= cw[0]) & (sw <=cw[-1]) & (sdq == 0)#mask to same same waveband and cut dq flags
    cw, cf, sw, sf = cw[c_mask], cf[c_mask], sw[s_mask], sf[s_mask]
    cw1, cf1 = resample.bintogrid(cw, cf, newx=sw) #rebin to stis wavelength grid
    airglow_mask = mask_maker(sw, airglow)
    sw, sf, cf1 = sw[airglow_mask], sf[airglow_mask], cw[airglow_mask] #remove airglow
    c_int = np.trapz(cf1,sw)
    s_int =  np.trapz(sf,sw)
    normfac = c_int/s_int
    return normfac

def add_stis(sed_table, component_repo):
    """
    Add the stis_spectrum to the sed
    """