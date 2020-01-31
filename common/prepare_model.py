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
import lzma

"""
@author: David Wilson

@version: 23 

@date :20191212

Turns models into standard MUSCLES file. Added all of prepare_phoenix in here as well to keep things cleaner. Adapted to use Lyon models.
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
    print(save_path)
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
        
def phxurl(Teff, logg=4.5, repo='ftp'):
    """
    Constructs the URL for the phoenix spectrum file for a star with effective
    temperature Teff, log surface gravity logg. FeH and aM are fixed at 0.0 as the lyon database does not have other options

    Does not check that the URL is actually valid, and digits beyond the
    precision of the numbers used in the path will be truncated.
    """
    phoenixbaseurl = 'https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/SPECTRA/'
  
    name = 'lte{T:05.1f}-{g:3.1f}-0.0a+0.0.BT-Settl.spec.7.xz'.format(T=Teff/100.0, g=logg)
    print(name)

    if repo == 'ftp':
        return phoenixbaseurl + name
    else:
        return os.path.join(repo, name)

def fetchphxfile(Teff, logg, FeH, aM, repo, source = 'lyon'):
    #keeping a source option in here in case I want to use the gottigen versions again
    if source == 'lyon':
        loc, ftp = [phxurl(Teff, logg, repo=r) for r in [repo, 'ftp']]
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
    makes a list of required atmospheic parameters to be retreived, also records which params need interpolation. Fixing FeH and aM = 0.0 for now.
    """
    params_to_interp = []
    param_names = ['Teff', 'logg', 'FeH', 'aM']
    param_list = []
    for param, grid, name in zip([star_params['Teff'],star_params['logg'] ,0.0, 0.0], grids, param_names):
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
    phxTgrid = np.arange(1200, 7001, 100)
    phxTgrid = np.hstack([np.arange(2300,7000,100),
                   np.arange(7000,12001,200)])
    phxggrid = np.arange(0.0, 6.1, 0.5)
    phxZgrid = np.array([0.0])
    phxagrid = np.array([0.0])
    return phxTgrid,phxggrid,phxZgrid, phxagrid

def extract_spectrum(filepath):
    """
    Extracts the spectrum from the Lyon files, which is non-trivial. Adapts code by JSP. 
    """
    DF = -8
    if filepath[-2:] == 'xz':
        nameout = filepath[:-3]
        with lzma.open(filepath) as f, open(filepath, 'wb') as fout: #https://stackoverflow.com/a/33718185
            file_content = f.read()
            fout.write(file_content)
    else:
        nameout = filepath
    phoenixR = ascii.read(nameout,format="fixed_width_no_header",col_starts=(0,14),col_ends=(12,25),delimiter=" ",names=('Wave','Spec'))
    ph1, jj = np.unique(np.array(phoenixR['Wave']),return_index=True)
    phoenix = np.zeros((len(ph1),2))
    for kk in range(len(jj)):
        phoenix[kk,1] = np.float64(phoenixR['Spec'][jj[kk]].replace("D","E"))
    phoenix[:,0] = ph1
    xraw = phoenix[:,0]
    yraw = np.power(10.,phoenix[:,1] + DF)
    mask = (xraw <= 100000) #original w is way too big
    return xraw[mask], yraw[mask]

    
def get_models(repo,param_dicts):
    """
    downloads phoenix spectra. Returns "spectra" param_dicts but with the model flux added to each dictionary
    """
    if os.path.exists(repo) == False:
        os.mkdir(repo)
    spectra = []
    for params in param_dicts:
        Teff, logg, FeH, aM = params['Teff'], params['logg'], params['FeH'], params['aM']
        file_path = phxurl(Teff, logg, repo=repo)
        if os.path.exists(file_path) == False and os.path.exists(file_path[:-3]) == False: #only download if we need to
            fetchphxfile(Teff, logg, FeH, aM, repo=repo)
        if os.path.exists(file_path[:-3]) == True:
            file_path = file_path[:-3]
        wavelength, flux =  extract_spectrum(file_path)
        params.update({'wavelength':wavelength})
        params.update({'flux':flux})
        spectra.append(params)
    return spectra
    
def interp_flux(spectra, params_to_interp, star_params):
    """
    build the new spectrum, interpolation each phoenix model onto the shortest wavelength array then interploating to the correct parameters
    """
    out_vals = [star_params[p] for p in params_to_interp]
    in_vals = [[s[p] for p in params_to_interp] for s in spectra]
    wavelengths = [s['wavelength'] for s in spectra]
    nwave = np.min([len(w) for w in wavelengths])
    for w in wavelengths: 
        if len(w) == nwave:
            wavelength = w
    fluxes = []
    for s in spectra:
        if len(s['flux']) == nwave:
            fluxes.append(s['flux'])
        else:
            fi = interp1d(s['wavelength'], s['flux'], fill_value='extrapolate')(wavelength)
            fluxes.append(fi)
        
    if len(params_to_interp) == 1:
        in_vals = [s[params_to_interp[0]] for s in spectra]
        new_flux = interp1d(in_vals, fluxes, axis=0, fill_value='extrapolate')(star_params[params_to_interp[0]])
    else:
        out_vals = [star_params[p] for p in params_to_interp]
        in_vals = [[s[p] for p in params_to_interp] for s in spectra]
        new_flux = griddata(in_vals, fluxes, out_vals)[0]
    return wavelength, new_flux

    
def save_phoenix(wavelength, flux, save_path):
    """
    save the new model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    savedat = Table([wavelength*u.AA, flux], names=['WAVELENGTH', 'FLUX'])
    ascii.write(savedat, save_path+'phoenix_interpolated.ecsv', format='ecsv', overwrite=True)
    
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
    wavelength, flux = extract_spectrum(filepath)
    return wavelength, flux
        
    
def make_phoenix_spectrum(save_path, repo, star_params, save_ecsv=False, plot=False):
    """
    Main array. Takes a list of stellar parameters and makes a phoenix spectrum out of. Save_path is where you want the final spectrum to go, repo is where downloaded phoenix files go. wave_file is where the wavelength array is 
    """
    tgrid, ggrig,fgrid, agrid = get_grids()
    param_list, params_to_interp = make_param_list(star_params, [tgrid, ggrig,fgrid, agrid])
    if len(params_to_interp) == 0: #i.e. if there's an existing model
        print('phoenix model available')
        wavelength, flux = get_existing_model(star_params, repo)
    else:
        param_dicts = make_dicts(param_list)
        spectra = get_models(repo,param_dicts)
        wavelength, flux = interp_flux(spectra, params_to_interp, star_params) 
    if save_ecsv:
        save_phoenix(wavelength, flux, save_path)
    if plot == True:
        plot_spectrum(wavelength, flux, star)
    return wavelength, flux

def distance_scale(radius, distance, flux):
    """
    Scales the phoenix model using the distance to the star
    """
    scale = (radius.to(u.cm)/distance.to(u.cm))**2
    return flux * scale

def load_star_params(star_table_path, FeH=0.0, aM=0.0):
    """
    Load one of SP's tables. Does not include metalicty and FeH, will proabaly need my own table.
    """
    data = Table.read(star_table_path, format ='ascii')[0]
    mass, radius = data['Mass__0'] * u.M_sun, data['Radius__0'] * u.Rsun
    g_star = (const.G * mass) / radius**2
    logg = np.log10((g_star.to(u.cm/u.s**2)).value)
    star_params = {'Teff':data['Teff__0'], 'logg': logg, 'FeH': FeH, 'aM': aM, 'radius':data['Radius__0'] * u.Rsun, 'distance':data['Distance__0']*u.pc }
    return star_params
