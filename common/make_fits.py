import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.units import cds

"""
Making hslp standard fits files for Mega-Muscles
"""

def data_header(hdr):
    """
    Makes the header for the table extension. 
    
    Keywords to add from MUSCLES:

    TDESC1  = 'midpoint of the wavelength bin'                                      
    TDESC2  = 'left/blue edge of the wavelength bin'                                
    TDESC3  = 'right/red edge of the wavelength bin'                                
    TDESC4  = 'average flux over the bin'                                           
    TDESC5  = 'error on the flux'                                                   
    TDESC6  = 'cumulative exposure time for the bin'                                
    TDESC7  = 'data quality flags (HST data only)'                                  
    TDESC8  = 'modified julian date of start of first exposure'                     
    TDESC9  = 'modified julian date of end of last exposure'                        
    TDESC10 = 'bitmask identifying the source instrument(s). See "instlgnd" & extension for a legend.'                                             
    TDESC11 = 'normalization factor applied to the source spectrum'                 
    TDESC12 = 'flux density normalized by the bolometric flux'                      
    TDESC13 = 'error on bolometrically-normalized flux density'

    """
    new_keywords = ('TDESC1','TDESC2','TDESC3','TDESC4', 'TDESC5',
               'TDESC6', 'TDESC7', 'TDESC8', 'TDESC9', 'TDESC10',
               'TDESC11','TDESC12','TDESC13') 
    new_values =  ('midpoint of the wavelength bin',
                   'left/blue edge of the wavelength bin',
                   'right/red edge of the wavelength bin',
                   'average flux over the bin',
                   'error on the flux',
                   'cumulative exposure time for the bin',
                   'data quality flags (HST data only)',
                   'modified julian date of start of first exposure',
                   'modified julian date of end of last exposure',
                   'bitmask identifying the source instrument(s). See "instlgnd" & extension for a legend.',
                   'normalization factor applied to the source spectrum',
                   'flux density normalized by the bolometric flux',
                   'error on bolometrically-normalized flux density')
    for i, n, v in zip(range(len(new_keywords)), new_keywords, new_values):
        hdr.insert(i+8, (new_keywords[i], new_values[i]))
    return hdr


def make_data_ext(table):
    """
    The table extension, takes an astroy table 
    """
    hdu = fits.table_to_hdu(table)
    hdu.header = data_header(hdu.header)
    return hdu

def make_primary_ext():
    """
    Make the primary header
    """
    hdr = fits.Header()
    primary_hdu = fits.PrimaryHDU(header=hdr)
    return primary_hdu

