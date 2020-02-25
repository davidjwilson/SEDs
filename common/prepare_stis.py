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
from scipy.io.idl import readsav



cds.enable()

"""
@author: David Wilson

version 1 20190717


Finds all STIS x1d files, groups them by grating, coadds them and saves to file with required metadata

"""
def coadd_flux(f_array, e_array):
    """
    Returns a variance-weighted coadd with standard error of the weighted mean (variance weights, scale corrected).
    f_array and e_arrays are collections of flux and error arrays, which should have the same lenth and wavelength scale
    """
    weights = 1 / (e_array**2)
    flux = np.average(f_array, axis =0, weights = weights)
    var = 1 / np.sum(weights, axis=0)
    rcs = np.sum((((flux - f_array)**2) * weights), axis=0) / (len(f_array)-1) #reduced chi-squared
    error = (var * rcs)**0.5
    return flux, error

def no_zero_errors(flux, error):
    """
    Corrects instances where negative flux measurements have very small errors
    """
    e_new = error
    for i in range(len(error)):
        if flux[i] < 0.0 and error[i] < 0.1*abs(flux[i]):
            e_new[i] = abs(flux[i])
    return e_new

def build_wavelength(x1ds):
    """
    builds a wavelength array covering all wavelength ranges in x1d (different cenwaves)
    """
    starts = []
    ends = []
    diffs = []
    for x in x1ds:
        w = fits.getdata(x, 1)[0]['WAVELENGTH']
        starts.append(w[0])
        ends.append(w[-1])
        diffs.append(np.max(np.diff(w)))
    w_new = np.arange(min(starts),max(ends), max(diffs))
    return w_new


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

 
def get_ayres_e140m(x1ds): 
    """ 
    reads T Ayres combined e140m files and bulds the data arrays from them. Hacky for now, will replace with my own routines when the new e140m calibrations are available.
    """
    target = fits.getheader(x1ds[0])['TARGNAME']
    savpath = '/home/david/work/muscles/SEDs/common/ayres_e140m/{}_E140M_coadd.sav'.format(target)
    data = readsav(savpath)
    w_new, f_new, e_new, dq_new = data['wave'], data['flux'], data['photerr'], data['epsilon']
    return w_new, f_new, e_new, dq_new
    
def combine_x1ds(x1ds, correct_error=True):
    """
    coadds a collection of x1d fluxes and adds columns for exposure time detials. Input is a list of paths to x1d files with the same grating. Also works for sx1 files

    """
    if len(x1ds) > 1:
        f_new = []
        e_new = []
        dq_new = []
        exptime = []
        start = []
        end = []
        w_new = build_wavelength(x1ds)
        for x in x1ds:
            data_extension = 1
            if x[-8:-5] == 'sx1':
                data_extension = 0
            data = fits.getdata(x,data_extension)[0]
            hdr = fits.getheader(x,0)
            fi = interpolate.interp1d(data['WAVELENGTH'], data['FLUX'], bounds_error=False, fill_value=0.)(w_new)
            ei = interpolate.interp1d(data['WAVELENGTH'], data['ERROR'], bounds_error=False, fill_value=0.)(w_new)
            dqi =  interpolate.interp1d(data['WAVELENGTH'], data['DQ'], kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            expi = np.full(len(data['WAVELENGTH']), hdr['TEXPTIME'])
            expi = interpolate.interp1d(data['WAVELENGTH'], expi, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            starti = np.full(len(data['WAVELENGTH']), hdr['TEXPSTRT'])
            starti = interpolate.interp1d(data['WAVELENGTH'], starti, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            endi = np.full(len(data['WAVELENGTH']), hdr['TEXPEND'])
            endi = interpolate.interp1d(data['WAVELENGTH'], endi, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            
            if correct_error:    
                ei = no_zero_errors(fi, ei)
            f_new.append(fi)
            e_new.append(ei)
            dq_new.append(dqi)
            exptime.append(expi)
            start.append(starti)
            end.append(endi)
            

        f_new, e_new = coadd_flux(np.array(f_new), np.array(e_new))
        dq_new = np.array(dq_new, dtype=int)
        dq_new = [(np.sum(np.unique(dq_new[:,i]))) for i in range(len(dq_new[0]))]
        exptime = np.sum(np.array(exptime), axis=0)
        start = np.min(np.ma.masked_array(start, mask=[np.array(start) == 0.]), axis=0)
        end = np.max(np.array(end), axis=0)
    
    else: #in the case where there's only one available spectrum
        data_extension = 1
        if x1ds[0][-8:-5] == 'sx1': #modified 1 off for t1 spectrum, must improve later
            data_extension = 0
          #  data = Table.read('/home/david/work/muscles/SEDs/trappist-1/optical/t1_g430m_edit.ecsv')
        #else:
        data = fits.getdata(x1ds[0],data_extension)[0]
        hdr = fits.getheader(x1ds[0],0)
        w_new, f_new, e_new, dq_new = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
        exptime, start, end = np.full(len(data['WAVELENGTH']), hdr['TEXPTIME']), np.full(len(data['WAVELENGTH']), hdr['TEXPSTRT']), np.full(len(data['WAVELENGTH']), hdr['TEXPEND'])
        if correct_error:    
                e_new = no_zero_errors(f_new, e_new)
    if fits.getheader(x1ds[0])['OPT_ELEM'] == 'E140M':
        print('yes')
        w_new, f_new, e_new, dq_new = get_ayres_e140m(x1ds)
        exptime, start, end = np.full(len(w_new), 0),np.full(len(w_new), 0), np.full(len(w_new), 0)
    f_new, e_new = nan_clean(f_new), nan_clean(e_new)
    w0, w1 = wavelength_edges(w_new)
    new_data = {'WAVELENGTH':w_new*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f_new*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':e_new*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq_new,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
   

    return new_data


 
def make_metadata(x1ds, new_data, normfac):
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
    for x in x1ds:      
        hdr = fits.getheader(x,0)
        exptimes.append(hdr['TEXPTIME'])
        start_times.append(hdr['TEXPSTRT'])
        end_times.append(hdr['TEXPEND'])
        dates.append(hdr['TDATEOBS'])
    muscles_name = 'Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanet Host Stars'
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','APERTURE','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN','EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['','',hdr['OPT_ELEM'],'','','','','',muscles_name,'MUSCLES','David J. Wilson','','',min(dates),min(start_times),max(end_times),sum(exptimes),'SUM', 
                min(exptimes), max(exptimes), np.median(exptimes),normfac,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = hdr[name]
        else:
            metadata[name] = filler
    return metadata
    
def setup_list(x1ds):
    """
    Takes all x1ds in selection and sorts them by instrument setup
    """
    gratings = []
    x1ds_by_setup = []
    for x in x1ds:
        hdr = fits.getheader(x,0)
        grating = hdr['OPT_ELEM']
      
        gratings.append(grating)
    setups = np.unique(gratings, axis=0)
    for s in setups:
        collection = []
        for i in range(len(x1ds)):
            if gratings[i] == s:
                collection.append(x1ds[i])
        if len(collection) > 0:
            x1ds_by_setup.append(collection)
    return setups, x1ds_by_setup

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

def stis_clean(x1ds):
    """
    checks that all x1d files are stis spectra
    """
    stis_x1ds =[]
    for x in x1ds:
        if fits.getheader(x,0)['INSTRUME'] == 'STIS':
            stis_x1ds.append(x)
    return stis_x1ds

def make_component_filename(metadata, version):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def make_dataset_extension(x1ds):
    """
    Makes a fits extension containg a list of rootnames and dataset IDs used in the spectrum.
    """
    description_text = 'This extension contains a list of HST rootnames (9 character string in HST files downloaded from MAST) and dataset IDs of the exposures used to create this spectrum file. The dataset IDs can be used to directly locate the observations through the MAST HST data archive search interface. Multiple identifiers indicate the spectra were coadded. In some cases the spectra have been rextracted and will be different from those availbale from MAST.' 
    
    rootnames = []
    datasets = []
    for x in x1ds:
        hdr = fits.getheader(x)
        rootnames.append(hdr['ROOTNAME'],)
        datasets.append(hdr['ASN_ID'])
    dataset_table = Table([rootnames,datasets], names=[('ROOTNAME'),('DATASET_ID')])
    hdu = fits.table_to_hdu(dataset_table)
    hdu.header.insert(8, ('EXTNAME','SRCSPECS'))
    hdu.header.insert(9, ('EXTNO',3))
    hdu.header['COMMENT'] = description_text
    return hdu

    
def make_stis_spectum(x1dpath, version,savepath = '', plot=False, save_ecsv=False, save_fits=False, return_data=False, return_gratings = False, normfac=1.0, sx1 = True):
    """
    main function
    """
    all_x1ds = glob.glob(x1dpath+'*x1d.fits')
    stis_x1ds = stis_clean(all_x1ds) #get rid of any not-stis x1ds
    if sx1:
        all_sx1 = glob.glob(x1dpath+'*sx1.fits') #do the ccd as well
        stis_x1ds = np.concatenate((stis_x1ds, all_sx1))
    
    if len(stis_x1ds) > 0:
        gratings, x1ds_by_grating = setup_list(stis_x1ds)
        for x1ds in x1ds_by_grating:
            data = combine_x1ds(x1ds)
           # data = [wavelength*u.AA, flux*u.erg/u.s/u.cm**2/u.AA, error*u.erg/u.s/u.cm**2/u.AA, dq]
            metadata = make_metadata(x1ds, data, normfac)
            if plot:
                plot_spectrum(data, metadata)
            if save_ecsv:
                save_to_ecsv(data, metadata, savepath, version)
            if save_fits:
                data_set_hdu = make_dataset_extension(x1ds)
                save_to_fits(data, metadata, data_set_hdu, savepath, version)
    if return_data:
        return data
    if return_gratings:
        return gratings


def test():
    """
    testing with GJ 699
    """
    x1dpath = '/home/david/work/muscles/MegaMUSCLES/GJ_699/HST/STIS/'
    version = 1
    savepath = 'test_files/'
    make_stis_spectum(x1dpath, version, savepath = savepath, plot=True, save_ecsv=True, save_fits=True)
    
#test()