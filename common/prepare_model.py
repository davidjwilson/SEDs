import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import griddata, interp1d
from urllib.request import urlretrieve, urlopen

"""
@author: David Wilson

@version: 2 

@date :20190808

Turns models into standard MUSCLES file. Added all of prepare_phoenix in here as well to keep things cleaner
"""

def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff/2.
    w1 = w + diff/2.
    return w0, w1
    
def get_model_data(model_path):
    """
    Makes the model data array, assuming an input .txt file with WAVELENGTH and FLUX
    """
    data = Table.read(model_path, format = 'ascii')
    w, f = data['WAVELENGTH'], data['FLUX']
    w0, w1 = wavelength_edges(w)
    new_data = {'WAVELENGTH':w*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f*u.erg/u.s/u.cm**2/u.AA}
    return new_data
    
def make_model_metadata(new_data, normfac, sed_metadata, model_name):
    """
    Makes the metadata for the lya file -sed metadata is from the total sed so far (should just be COS?)
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['Model',model_name,'NA','','','','','','','','','',normfac,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = sed_metadata[name]
        else:
            metadata[name] = filler
    return metadata

def make_component_filename(metadata, version):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

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

    
def model_save_to_fits(data, metadata, savepath, version):
    """
    Saves to a MUSCLES-standard fits file for models
    """
    if os.path.exists(savepath) == False:
        os.mkdir(savepath)
    file_name = make_component_filename(metadata, version)
    hdr = fits.Header(metadata)
    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdu = fits.table_to_hdu(Table(data))
    descriptions =['midpoint of the wavelength bin', 'left/blue edge of the wavelength bin','right/red edge of the wavelength bin','average flux over the bin']
    hdu.header.insert(8, ('EXTNAME','SPECTRUM'))
    hdu.header.insert(9, ('EXTNO',2))
    [hdu.header.insert(i[0]+10, ('TDESC%s' %(i[0]), i[1])) for i in enumerate(descriptions)]
    hdul = fits.HDUList([primary_hdu, hdu])
    hdul.writeto(savepath+file_name+'.fits', overwrite=True)
    print('Spectrum saved as '+file_name+'.fits')

def make_model_spectrum(model_path, version, sed_data ,savepath = '', save_ecsv=False, save_fits=False, normfac=1.0, model_name='LYA-RECONSTRUCTION'):
    """
    Main function.
    """
    data = get_model_data(model_path)
    metadata = make_model_metadata(data,  normfac, sed_data.meta, model_name)
    if save_ecsv:
        save_to_ecsv(data, metadata, savepath, version)
    if save_fits:
        model_save_to_fits(data, metadata, savepath, version)
        
        
"""
PHOENIX
-----------------------------------------------------------------------------

"""
        
def phxurl(Teff, logg=4.5, FeH=0.0, aM=0.0, repo='ftp'):
    """
    Constructs the URL for the phoenix spectrum file for a star with effective
    temperature Teff, log surface gravity logg, metalicity FeH, and alpha
    elemnt abundance aM.

    Does not check that the URL is actually valid, and digits beyond the
    precision of the numbers used in the path will be truncated.
    """
    phoenixbaseurl = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/'
    zstr = '{:+4.1f}'.format(FeH)
    if FeH == 0.0: zstr = '-' + zstr[1:]
        
    #20190628 alphas for fe = -4.0 are missing, don't know why but adding these lines to account
   # if FeH == -4.0:
    #    aM = 0.0
    
    astr = '.Alpha={:+5.2f}'.format(aM) if aM != 0.0 else ''
    name = ('lte{T:05.0f}-{g:4.2f}{z}{a}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
            ''.format(T=Teff, g=logg, z=zstr, a=astr))

    if repo == 'ftp':
        folder = 'Z' + zstr + astr + '/'
        return phoenixbaseurl + folder + name
    else:
        return os.path.join(repo, name)

def fetchphxfile(Teff, logg, FeH, aM, repo):
    loc, ftp = [phxurl(Teff, logg, FeH, aM, repo=r) for r in [repo, 'ftp']]
    urlretrieve(ftp, loc)

def make_dicts(param_list):
    """
    makes array of dictionaries with parameters to download
    """
    param_dicts = []
    for teff in param_list[0]:
        for logg in param_list[1]:
            for feh in param_list[2]:
                for aM in param_list[3]:
                    param_dict = {'Teff':teff, 'logg':logg, 'FeH': feh, 'aM': aM}
                    if param_dict not in param_dicts:
                        param_dicts.append(param_dict)
    return param_dicts

def make_param_list(star_params, grids):
    """
    makes a list of required atmospheic parameters to be retreived, also records which params need interpolation
    """
    params_to_interp = []
    param_names = ['Teff', 'logg', 'FeH', 'aM']
    param_list = []
    for param, grid, name in zip([star_params['Teff'],star_params['logg'] ,star_params['FeH'], star_params['aM']], grids, param_names):
        if param in grid:
            param_list.append([param, param])
        else:
            idx = np.searchsorted(grid, param)
            param_list.append([grid[idx-1],grid[idx]])
            params_to_interp.append(name)
    return param_list, params_to_interp

def get_grids():
    """
    arrays storing the available phoenix spectra
    """
    phxTgrid = np.hstack([np.arange(2300,7000,100),
                   np.arange(7000,12001,200)])
    phxggrid = np.arange(0.0, 6.1, 0.5)
    phxZgrid = np.hstack([np.arange(-4.0, -2.0, 1.0),
                       np.arange(-2.0, 1.1, 0.5)])
    phxagrid = np.arange(-0.2, 1.3, 0.2)
    return phxTgrid,phxggrid,phxZgrid, phxagrid
    
def get_models(repo,param_dicts):
    """
    downloads phoenix spectra. Returns "spectra" param_dicts but with the model flux added to each dictionary
    """
    if os.path.exists(repo) == False:
        os.mkdir(repo)
    spectra = []
    for params in param_dicts:
        Teff, logg, FeH, aM = params['Teff'], params['logg'], params['FeH'], params['aM']
        file_path = phxurl(Teff, logg, FeH, aM, repo=repo)
        if os.path.exists(file_path) == False: #only download if we need to
            fetchphxfile(Teff, logg, FeH, aM, repo=repo)
        params.update({'flux':fits.getdata(file_path)})
        spectra.append(params)
    return spectra
    
def interp_flux(spectra, params_to_interp, star_params):
    """
    build the new spectrum
    """
    out_vals = [star_params[p] for p in params_to_interp]
    in_vals = [[s[p] for p in params_to_interp] for s in spectra]
    fluxes = [s['flux'] for s in spectra]
    if len(params_to_interp) == 1:
        in_vals = [s[params_to_interp[0]] for s in spectra]
        new_flux = interp1d(in_vals, fluxes, axis=0, fill_value='extrapolate')(star_params[params_to_interp[0]])
    else:
        out_vals = [star_params[p] for p in params_to_interp]
        in_vals = [[s[p] for p in params_to_interp] for s in spectra]
        new_flux = griddata(in_vals, fluxes, out_vals)[0]
    return new_flux
    
def get_wavelength(wave_file):
    """
    Load the wavelenth array from file
    """
    return fits.getdata(wave_file)
    
def phx_save_to_ecsv(wavelength, flux, save_path, star_params):
    """
    save the new model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    savedat = Table([wavelength*u.AA, flux], names=['WAVELENGTH', 'FLUX'])
    name = '%s_%s_%s_%s_phoenix_interpolated.txt' % (star_params['Teff'], star_params['logg'], star_params['FeH'], star_params['aM'])
    ascii.write(savedat, save_path+name, overwrite=True)
    
def plot_spectrum(wavelength, flux, star):
    plt.figure(star)
    plt.plot(wavelength, flux)
    plt.xlabel('Wavelength (\AA)', size=20)
    plt.ylabel('Flux (adu)', size=20)
    plt.tight_layout()
    plt.show()

def get_existing_model(star_params, repo):
    """
    Get the flux if there's already a good phoenix model
    """
    if os.path.exists(repo) == False:
        os.mkdir(repo)
    Teff, logg, FeH, aM = star_params['Teff'], star_params['logg'], star_params['FeH'], star_params['aM']
    file_path = phxurl(Teff, logg, FeH, aM, repo=repo)
    if os.path.exists(file_path) == False: #only download if we need to
        fetchphxfile(Teff, logg, FeH, aM, repo=repo)
    flux = fits.getdata(file_path)
    return flux
        
    
def make_phoenix_spectrum(wave_file, save_path, repo, star_params, save_ecsv=False, plot=False):
    """
    Main Function. Takes a list of stellar parameters and makes a phoenix spectrum out of. Save_path is where you want the final spectrum to go, repo is where downloaded phoenix files go. wave_file is where the wavelength array is 
    """
    tgrid, ggrig,fgrid, agrid = get_grids()
    param_list, params_to_interp = make_param_list(star_params, [tgrid, ggrig,fgrid, agrid])
    if len(params_to_interp) == 0: #i.e. if there's an existing model
        print('phoenix model available')
        flux = get_existing_model(star_params, repo)
    else:
        param_dicts = make_dicts(param_list)
        spectra = get_models(repo,param_dicts)
        flux = interp_flux(spectra, params_to_interp, star_params)
    wavelength = get_wavelength(wave_file)
    if save_ecsv == True:
        phx_save_to_ecsv(wavelength, flux, save_path, star_params)
    if plot == True:
        plot_spectrum(wavelength, flux, star)
    #return wavelength, flux

def load_star_params(star_table_path, FeH=0.0, aM=0.0):
    """
    Load one of SP's tables. Does not include metalicty and FeH, will proabaly need my own table.
    """
    data = Table.read(star_table_path, format ='ascii')[0]
    mass, radius = data['Mass__0'] * u.M_sun, data['Radius__0'] * u.Rsun
    g_star = (const.G * mass) / radius**2
    logg = np.log10((g_star.to(u.cm/u.s**2)).value)
    star_params = {'Teff':data['Teff__0'], 'logg': logg, 'FeH': FeH, 'aM': aM}
    return star_params

def residuals(scale, f, mf):
    return f - mf/scale

def phoenix_norm(star, w_phx, f_phx, ccd_path, plot=False, cut=2000): 
    """
    find the normalisation factor between the phoenix model and the stis ccd (ccd_path)
    """
    data = fits.getdata(ccd_path)[0]
    w, f, dq = data['WAVELENGTH'], data['FLUX'], data['DQ']
    #mask =  (dq ==0)
    mask = (w > cut) & (f > 0) & (dq ==0)
    w1, f1 = w[mask], f[mask]
    norm_mask = (w_phx >= w1[0]) & (w_phx <= w1[-1])
    phx_flux = interp1d(w_phx[norm_mask], f_phx[norm_mask], fill_value='extrapolate')(w1)
    scale, flag = leastsq(residuals, 1., args=(f1, phx_flux))
    normfac = 1/scale[0]

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
    