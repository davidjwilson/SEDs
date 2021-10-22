import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.units import cds
import instruments
cds.enable()
"""
Making hslp standard fits files for Mega-Muscles
"""
def make_primary_header(hdr, sed_table, instrument_list):
    meta = sed_table.meta
    telescopes = [instruments.getinststr(inst)[0:3] for inst in instrument_list]
    instrus = [instruments.getinststr(inst)[4:7] for inst in instrument_list]
    gratings = [instruments.getinststr(inst)[8:] for inst in instrument_list]
    hdr.append(('TELESCOP', 'MULTI'))
    hdr.append(('INSTRUME', 'MULTI'))
    hdr.append(('GRATING', 'MULTI'))
    for i in range(len(telescopes)):
        hdr.append(('TELESC{:02.0f}'.format(i), telescopes[i].upper()))
        hdr.append(('INSTRU{:02.0f}'.format(i), instrus[i].upper()))
        hdr.append(('GRATIN{:02.0f}'.format(i), gratings[i].upper()))
    extra_keys =  ['TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','FLUXMIN',
                  'FLUXMAX','FLUXUNIT', 'BOLOFLUX', 'LNZ_NORM','LNZ_GAM']
    for key in extra_keys[:-2]: #temporary while I don't have the last ones in there.
        hdr.append((key,meta[key]))
    return hdr



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


def make_instrument_extension(sed_table, instlist, version):
    """
    Make the instrument extension for a MEGA-MUSCLES file, containing translations of the instrument codes and the list of component file names
    
    """
    telescopes = [instruments.getinststr(inst)[0:3].upper() for inst in instlist]
    instrus = [instruments.getinststr(inst)[4:7].upper() for inst in instlist]
    gratings = [instruments.getinststr(inst)[8:].replace('-----','na').upper() for inst in instlist]
    target = sed_table.meta['TARGNAME']
    filenames = ['hlsp_muscles_{}_{}_{}_{}_v{}_component-spec.fits'.format(tel, inst, target, grating, version).lower() for tel, inst, grating in zip(telescopes, instrus,gratings)]
    data = Table([instlist, telescopes, instrus, gratings , filenames], names = ['BITVALUE','TELESCOPE','INSTRUMENT','GRATING','HLSP_FILE'])
    hdu = fits.table_to_hdu(data)
    
    comment = 'This extension is a legend for the integer identifiers in the instrument column of the '\
              'spectrum extension. Instruments are identified by bitwise flags so that any combination of '\
              'instruments contributing to the data within a spectral element can be identified together. '\
              'For example, if instruments 4 and 16 (100 and 10000 in binary) both contribute to the data '\
              'in a bin, then that bin will have the value 20, or 10100 in binary, to signify that both '\
              'instruments 4 and 16 have contributed. This is identical to the handling of bitwise data '\
              'quality flags. Note that polynomial fits for filling spectral gaps were not saved as separate spectra.'\

    
    
    hdu.header.append(('COMMENT',comment))
    return hdu
    
    


def make_data_ext(sed_table):
    """
    The table extension, takes an astropy table 
    """
    hdu = fits.table_to_hdu(Table(dict(sed_table)))
    hdu.header = data_header(hdu.header)
    hdu.name = 'SPECTRUM' #Monday me see if this works.
    return hdu


def make_primary_ext(sed_table, instrument_list):
    """
    Make the primary header
    """
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header = make_primary_header(primary_hdu.header, sed_table, instrument_list)
    return primary_hdu

def make_mm_fits(savepath, sed_table, instrument_list, version,sed_type='var'):
    """
    Saves an SED as a Mega-MUSCLES fits file
    """
    primary_hdu = make_primary_ext(sed_table, instrument_list)
    data_ext = make_data_ext(sed_table)
    inst_ext = make_instrument_extension(sed_table, instrument_list, version)
    star = sed_table.meta['TARGNAME'].lower()
    if star == '2mass-j23062928-0502285':
        star = 'trappist-1'
    file_name = 'hlsp_muscles_multi_multi_{}_broadband_v{}_{}-res-sed'.format(star, version, sed_type)
    hdul = fits.HDUList([primary_hdu, data_ext,inst_ext])
    hdul.writeto('{}{}.fits'.format(savepath,file_name), overwrite=True)
    print('sed saved as {}'.format(file_name))