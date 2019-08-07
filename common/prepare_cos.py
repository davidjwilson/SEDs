import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from scipy.interpolate import interpolate
from astropy.units import cds
from scipy.io import readsav
cds.enable()

"""
@author: David Wilson

version 1 20190723

Turns KF's COS files into Muscles Stanard fits files.

"""
def no_zero_errors(flux, error):
    """
    Corrects instances where negative flux measurements have very small errors
    """
    e_new = error
    for i in range(len(error)):
        if flux[i] < 0.0 and error[i] < 0.1*abs(flux[i]):
            e_new[i] = abs(flux[i])
    return e_new

def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff/2.
    w1 = w + diff/2.
    return w0, w1

def nan_clean(array):
    """
    replace nans in arrays with zeros
    """
    for i in range(len(array)):
        if np.isnan(array[i]) == True:
            array[i] = 0.0
    return array

def make_cos_mjds(w_new, wave_array, x1d_array, x1dpath):
    start = []
    end = []
    for w, name in zip(wave_array,x1d_array):
        x1d = x1dpath+str(name)[2:-1]+'_x1d.fits'
        hdr = fits.getheader(x1d,1)
        starti = np.full(len(w), hdr['EXPSTART'])
        starti = interpolate.interp1d(w, starti, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
        endi = np.full(len(w), hdr['EXPEND'])
        endi = interpolate.interp1d(w, endi, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
        start.append(starti)
        end.append(endi)
    start = np.min(np.ma.masked_array(start, mask=[np.array(start) == 0.]), axis=0).filled(0)
    end = np.max(np.array(end), axis=0)
    return start, end


def make_cos_data(sav, x1dpath, correct_error=True):
    """
    Turns a readsav file into the data array. Needs the file and path to the x1d files making that file
    """
    data = readsav(sav)
    w_new, f_new, e_new, exptime = np.array(data['wave']), np.array(data['flux']), np.array(data['err']), np.array(data['exptime']) 
    dq_new = []
    
    for i in range(len(data['dqin'])):
        dqi =  interpolate.interp1d(data['wavein'][i], data['dqin'][i], kind='nearest',bounds_error=False, fill_value=0.)(w_new)
        dq_new.append(dqi)



    dq_new = np.array(dq_new, dtype=int)
    dq_new = [(np.sum(np.unique(dq_new[:,i]))) for i in range(len(dq_new[0]))] 

    start, end = make_cos_mjds(w_new, data['wavein'], data['files'], x1dpath)
    
    if correct_error:    
            e_new = no_zero_errors(f_new, e_new)
    f_new, e_new = nan_clean(f_new), nan_clean(e_new)
    w0, w1 = wavelength_edges(w_new)
    new_data = {'WAVELENGTH':w_new*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f_new*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':e_new*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq_new,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
   

    return new_data

def make_nuv_data(g230l_path, correct_error=True):
    """
    Turns a cos g230l file into a data array
    """
    hdul = fits.open(g230l_path) 
    data = fits.getdata(filepath,1)[0:2] #not using reflected order
    hdr0 = fits.getheader(filepath,0)
    hdr1 = fits.getheader(filepath,1)
    w = np.array([], dtype=float)
    f = np.array([], dtype=float)
    e = np.array([], dtype=float)
    dq = np.array([], dtype=int)
    if fillgap == True:
        gap_w, gap_f = nuv_fill(data)
        for dt in data:
            mask = (dt['WAVELENGTH'] < gap_w[0]) | (dt['WAVELENGTH'] > gap_w[-1])
            w= np.concatenate((w, dt['WAVELENGTH'][mask]))
            f = np.concatenate((f, dt['FLUX'][mask]))
            e = np.concatenate((e, dt['ERROR'][mask]))
            dq = np.concatenate((dq, dt['DQ'][mask]))
        gap_w0, gap_w1 = wavelength_edges(gap_w)
        gap_code = inst.getinsti('oth_---_other')
        gap_collection = dict_builder(gap_w0, gap_w1, gap_w, gap_f, np.zeros(len(gap_w)), np.zeros(len(gap_w)),np.zeros(len(gap_w)), 0., 0., gap_code)        
    else:
        for dt in data:
            w= np.concatenate((w, dt['WAVELENGTH']))
            f = np.concatenate((f, dt['FLUX']))
            e = np.concatenate((e, dt['ERROR']))
            dq = np.concatenate((dq, dt['DQ']))
        gap_collection = {}
    w0, w1 = wavelength_edges(w)
    exptime = np.full(len(w), hdr1['EXPTIME'])
    expstart, expend = hdr1['EXPSTART'], hdr1['EXPEND']
    instrument_name = hdr0['TELESCOP']+'_cos_'+hdr0['OPT_ELEM']   
    instrument_code = inst.getinsti(instrument_name.lower())
    cos_nuv_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    return cos_nuv_collection, gap_collection

def make_cos_metadata(sav, new_data, x1dpath):
    """
    Makes the metadata for the ecsv files- eventually will be converted into the fits file
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    exptimes = []
    start_times = []
    end_times = []
    dates = []
    data = readsav(sav)
    for name in data['files']:
        x = x1dpath+str(name)[2:-1]+'_x1d.fits'
        hdr0 = fits.getheader(x,0)
        hdr1 = fits.getheader(x,1)
        exptimes.append(hdr1['EXPTIME'])
        start_times.append(hdr1['EXPSTART'])
        end_times.append(hdr1['EXPEND'])
        dates.append(hdr1['DATE-OBS'])
    muscles_name = 'Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanet Host Stars'
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','APERTURE','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN','EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['','',hdr0['OPT_ELEM'],'','','','','',muscles_name,'MUSCLES','David J. Wilson','','',min(dates),min(start_times),max(end_times),sum(exptimes),'SUM', 
                min(exptimes), max(exptimes), np.median(exptimes),1.0,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = hdr0[name]
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
    

def plot_spectrum(data, metadata):
    plt.figure('%s_%s' % (metadata['TARGNAME'],metadata['GRATING']))
    plt.step(data['WAVELENGTH'], data['FLUX'], where='mid')
    plt.step(data['WAVELENGTH'], data['ERROR'], where='mid')
    plt.xlabel('Wavelength (\AA)', size=20)
    plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)', size=20)
    plt.tight_layout()
    plt.show()

def make_component_filename(metadata, version):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def make_dataset_extension(sav,x1dpath):
    """
    Makes a fits extension containg a list of rootnames and dataset IDs used in the spectrum.
    """
    description_text = 'This extension contains a list of HST rootnames (9 character string in HST files downloaded from MAST) and dataset IDs of the exposures used to create this spectrum file. The dataset IDs can be used to directly locate the observations through the MAST HST data archive search interface. Multiple identifiers indicate the spectra were coadded. In some cases the spectra have been rextracted and will be different from those availbale from MAST.' 
    
    rootnames = []
    datasets = []
    data = readsav(sav)
    for name in data['files']:
        x = x1dpath+str(name)[2:-1]+'_x1d.fits'
        hdr = fits.getheader(x,0)
        rootnames.append(hdr['ROOTNAME'],)
        datasets.append(hdr['ASN_ID'])
    dataset_table = Table([rootnames,datasets], names=[('ROOTNAME'),('DATASET_ID')])
    hdu = fits.table_to_hdu(dataset_table)
    hdu.header.insert(8, ('EXTNAME','SRCSPECS'))
    hdu.header.insert(9, ('EXTNO',3))
    hdu.header['COMMENT'] = description_text
    return hdu

def make_cos_spectrum(savpath, version, x1dpath,savepath = '', plot=False, save_ecsv=False, save_fits=False):
    """
    Main function. Take one of Kevin France's coadded x1d files (they come added by grating) and make it into a muscles fits file
    """
    savs = glob.glob(savpath+'*.sav')
    for sav in savs:
        data = make_cos_data(sav, x1dpath)
        metadata = make_cos_metadata(sav, data, x1dpath)
        if plot:
            plot_spectrum(data, metadata)
        if save_ecsv:
            save_to_ecsv(data, metadata, savepath, version)
        if save_fits:
            data_set_hdu = make_dataset_extension(sav,x1dpath)
            save_to_fits(data, metadata, data_set_hdu, savepath, version)
        
def make_cos_nuv():
    """
    Makes and saves an nuv spectrum
    """
    

def test():
    """
    testing with Trappist-1
    """
    savpath = ''
    x1dpath = '/home/david/work/muscles/trappist-1/hst/data/'
    version = 1
    savepath = 'test_files/'
    make_cos_spectrum(savpath, version,  x1dpath, savepath = savepath, plot=True, save_ecsv=True, save_fits=True)
    
#test()