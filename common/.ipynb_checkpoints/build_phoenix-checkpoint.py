import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import griddata, interp1d
from urllib.request import urlretrieve, urlopen

"""
@author David Wilson

version 2 20190709

Script to retreive phoenix models and interpolate them onto the correct values. "phxurl" and "fetchphxfile" adaped from Parke Loyds scripts
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
    
def save_to_ecsv(star,wavelength, flux, save_path):
    """
    save the new model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    savedat = Table([wavelength*u.AA, flux], names=['WAVELENGTH', 'FLUX'])
    ascii.write(savedat, save_path+star+'_phoenix_interpolated.ecsv', overwrite=True)
    
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
        
    
def make_phoenix_spectrum(star,wave_file, save_path, repo, star_params, save_ecsv=False, plot=False):
    """
    Main array. Takes a list of stellar parameters and makes a phoenix spectrum out of. Save_path is where you want the final spectrum to go, repo is where downloaded phoenix files go. wave_file is where the wavelength array is 
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
        save_to_ecsv(star, wavelength, flux, save_path)
    if plot == True:
        plot_spectrum(wavelength, flux, star)
    return wavelength, flux

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


def test():
    star = 'GJ_699'
    star_params = load_star_params(star+'_ParamStats.txt', FeH=-0.4)   
    path = '/home/david/work/muscles/phoenix/file_dump/'
    wave_file = path+ 'phoenix/wavegrid_hires.fits'
    repo = path + star+'_repo/'
    save_path = path + star+'_output/'
    #star_params = {'Teff':3300, 'logg': 4.00, 'FeH': -2.0, 'aM': 0.0}
    w_phx, f_phx = make_phoenix_spectrum(star,wave_file, save_path, repo, star_params, save_ecsv=False, plot=False)
    ccd_path = 'g430l/odlm24010_sx1.fits'
    normfac = phoenix_norm(star, w_phx, f_phx, ccd_path, plot=True)
    print(normfac)
    #facs = []
    #cuts = np.arange(3000,4000, 100)
    #for cut in cuts:
     #   normfac = phoenix_norm(cut,star, w_phx, f_phx, ccd_path, plot=False)
      #  facs.append(normfac)
    #plt.plot(cuts, facs/np.mean(facs))
    #plt.show()

#test()