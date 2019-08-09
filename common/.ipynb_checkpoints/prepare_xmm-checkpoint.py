import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.units import cds
cds.enable()


"""
@verison: 1

@author: David Wilson

@date 20190809

Turns XMM data in to HLSP format ecsv and fits files, ready to be added to an SED
"""

def nan_clean(array):
    """
    replace nans in arrays with zeros
    """
    for i in range(len(array)):
        if np.isnan(array[i]) == True:
            array[i] = 0.0
    return array

def apec_to_ecsv(data, hdr0, save_path):
    """
    save the apec model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    wavelength, flux = data['Wave'], data['Flux']
    target = hdr0['TARGET'].replace(' ','_')
    savedat = Table([wavelength, flux], names=['WAVELENGTH', 'FLUX'])
    name = target+'apec.txt' 
    ascii.write(savedat, save_path+name, overwrite=True)

def build_xmm_data(data, hdr0):
    w, bins, f, e = data['Wave'], data['bin_width'], data['CFlux'], data['CFlux_err']
    w0, w1 = w - (bins/2), w+(bins/2)
    exptime = np.full(len(w), hdr0['pn_DURATION'])
    start = np.full(len(w), (Time(hdr0['pn_DATE-OBS']).mjd))
    end = np.full(len(w), (Time(hdr0['pn_DATE-END']).mjd))
    dq = np.zeros(len(w), dtype=int)
    f, e = nan_clean(f), nan_clean(e)
    new_data = {'WAVELENGTH':w*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':e*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
    return new_data

def build_xmm_metadata(hdr0, new_data, sed_meta):
    """
    Makes the metadata for the xmm data table
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    start, end, exptime = np.min(new_data['EXPSTART'][new_data['EXPSTART']>0]).value, np.max(new_data['EXPEND']).value, np.max(new_data['EXPTIME']).value
    meta_names =['TELESCOP','INSTRUME','GRATING','DETECTOR','DETECT00','DETECT01','DETECT02','FILTER','FILTER00',
                 'FILTER01','FILTER02','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD',
                 'PR_INV_L','PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN',
                 'EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['XMM','EPIC','NA','MULTI','PN','MOS1','MOS2','MULTI', hdr0['pn_FILTER'],
                 hdr0['mos1_FILTER'],hdr0['mos2_FILTER'],'sed','sed','sed','sed','sed','sed','sed',
                 'sed','sed',hdr0['pn_DATE-OBS'], start, end, exptime, 'MEAN', 
                 exptime, 1.0, min(wavelength), max(wavelength), 'ang', 'vac', specres, waveres,np.min(flux[np.isnan(flux)==False]),
                 np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']  
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == 'sed':
            metadata[name] = sed_meta[name]
        else:
            metadata[name] = filler
    return metadata
    
def save_to_ecsv(data, metadata, save_path, version):
    """
    save the new model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    file_name = make_component_filename(metadata, version)
    savedat = Table(data, meta=metadata)
    savedat.write(save_path+file_name+'.ecsv', overwrite=True, format='ascii.ecsv')
    print('Spectrum saved as '+file_name+'.ecsv')

    
def save_to_fits(data, metadata, dataset_hdu, savepath, version):
    """
    Saves to a MUSCLES-standard fits file
    """
    if os.path.exists(savepath) == False:
        os.mkdir(savepath)
    file_name = make_component_filename(metadata, version)
    hdr = fits.Header(metadata)
    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdu = fits.table_to_hdu(Table(data))
    descriptions =['midpoint of the wavelength bin', 'left/blue edge of the wavelength bin','right/red edge of the wavelength bin','average flux over the bin',
                'error on the flux','cumulative exposure time for the bin','data quality flags (HST data only)','modified julian date of start of first exposure', 
                'modified julian date of end of last exposure']
    hdu.header.insert(8, ('EXTNAME','SPECTRUM'))
    hdu.header.insert(9, ('EXTNO',2))
    [hdu.header.insert(i[0]+10, ('TDESC%s' %(i[0]), i[1])) for i in enumerate(descriptions)]
    hdul = fits.HDUList([primary_hdu, hdu, dataset_hdu])
    hdul.writeto(savepath+file_name+'.fits', overwrite=True)
    print('Spectrum saved as '+file_name+'.fits')    
    
    
def make_component_filename(metadata, version):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def make_dataset_extension(hdr):
    """
    Makes a fits extension containg a list of rootnames and dataset IDs used in the spectrum.
    """
    description_text = 'This extension contains a list of observation IDs (DATASET_ID used for consistency with HST data) that can be used to locate the data in the XMM archives. XMM data all come from only a single observation (unlike the HST observations), but this extension is retained in place of a keyword for consistency with the HST files. '
    
    rootnames = [hdr['OBS_ID'],]
    datasets = ['']
    dataset_table = Table([rootnames,datasets], names=[('ROOTNAME'),('DATASET_ID')])
    hdu = fits.table_to_hdu(dataset_table)
    hdu.header.insert(8, ('EXTNAME','SRCSPECS'))
    hdu.header.insert(9, ('EXTNO',3))
    hdu.header['COMMENT'] = description_text
    return hdu
    
    
def make_xmm_spectra(xmm_path, savepath, sed_meta, version, apec_repo='', make_apec=True, save_ecsv=False, save_fits=False):
    hdul = fits.open(xmm_path)
    hdr0 = hdul[0].header
    data = hdul[1].data
    apec_data = hdul[2].data
    hdul.close
    data = build_xmm_data(data, hdr0)
    metadata = build_xmm_metadata(hdr0, data, sed_meta)
    if make_apec:
        apec_to_ecsv(apec_data, hdr0, apec_repo)
    if save_ecsv:
        save_to_ecsv(data, metadata, savepath, version)
    if save_fits:
        data_set_hdu = make_dataset_extension(hdr0)
        save_to_fits(data, metadata, data_set_hdu, savepath, version)
    
"""
                                             
TARGET  = 'GJ 674  '                                                            
OBS_ID  = '0810210301'                                                          
HIERARCH pn_FILTER = 'Medium  '                                                 
HIERARCH pn_MODE = 'IMAGING '                                                   
HIERARCH pn_SUBMODE = 'PrimeSmallWindow'                                        
HIERARCH pn_DATE-OBS = '2018-04-03T06:37:28'                                    
HIERARCH pn_DATE-END = '2018-04-03T14:39:08'                                    
HIERARCH pn_EXP_ID = '0810210301001'                                            
HIERARCH pn_DURATION = 29317.0                                                  
HIERARCH mos1_FILTER = 'Medium  '                                               
HIERARCH mos1_MODE = 'IMAGING '                                                 
HIERARCH mos1_SUBMODE = 'PrimePartialW2'                                        
HIERARCH mos1_DATE-OBS = '2018-04-03T06:31:59'                                  
HIERARCH mos1_DATE-END = '2018-04-03T14:36:15'                                  
HIERARCH mos1_EXP_ID = '0810210301002'                                          
HIERARCH mos2_FILTER = 'Medium  '                                               
HIERARCH mos2_MODE = 'IMAGING '                                                 
HIERARCH mos2_SUBMODE = 'PrimePartialW3'                                        
HIERARCH mos2_DATE-OBS = '2018-04-03T06:32:29'                                  
HIERARCH mos2_DATE-END = '2018-04-03T14:36:20'                                  
HIERARCH mos2_EXP_ID = '0810210301003'                                          
HIERARCH Instrument = 'EPIC    '                                                





TELESCOP= 'XMM     '                                                            
INSTRUME= 'EPIC    '                                                            
GRATING = 'NA      '                                                            
DETECTOR= 'MULTI   '                                                            
DETECT00= 'PN      '                                                            
DETECT01= 'MOS1    '                                                            
DETECT02= 'MOS2    '                                                            
FILTER  = 'MULTI   '                                                            
FILTER00= 'Thick   '                                                            
FILTER01= 'Thick   '                                                            
FILTER02= 'Medium  '                                                            
TARGNAME= 'GJ832   '                                                            
RA_TARG =           323.391564                                                  
DEC_TARG=           -49.009005                                                  
PROPOSID=                13650                                                  
HLSPNAME= 'Measurements of the Ultraviolet Spectral Characteristics of &'       
CONTINUE  'Low-mass Exoplanet Host Stars'                                       
HLSPACRN= 'MUSCLES '                                                            
HLSPLEAD= 'R. O. Parke Loyd'                                                    
PR_INV_L= 'France  '                                                            
PR_INV_F= 'Kevin   '                                                            
DATE-OBS= '2014-10-11T03:08:08.000'                                             
EXPSTART=    56941.13064814815                                                  
EXPEND  =    56941.26559027778                                                  
EXPTIME =    10633.33333333333                                                  
EXPDEFN = 'MEAN    '                                                            
NORMFAC =                  1.0 / normalization factor used by MUSCLES           
WAVEMIN =    7.424886226654053                                                  
WAVEMAX =    56.66274642944336                                                  
WAVEUNIT= 'ang     '                                                            
AIRORVAC= 'vac     '                                                            
SPECRES =    22.45938587188721                                                  
WAVERES =   0.5318784713745117                                                  
FLUXMIN = 1.61376280359525E-16                                                  
FLUXMAX = 3.44657261545477E-15                                                  
FLUXUNIT= 'erg/s/cm2/ang'   
"""
