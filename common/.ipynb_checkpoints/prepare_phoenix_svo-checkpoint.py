#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import griddata, interp1d
import sys
import glob
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian1DKernel



__author__ ='David Wilson, Parke Loyd'
__version__=6.01
__date__=20220921


"""
BT Settl models are now availabe on the SVO, and are much easier to work with.
Updating to optionally return an error estimate, based on effective temperature
"""


def make_filepath(Teff, logg=4.5, repo='ftp'):
    """
    Constructs the filepath for a phoenix spectrum file for a star with effective
    temperature Teff, log surface gravity logg. 
    """
   
    name = 'lte{T:05.1f}-{g:3.1f}-0.0a+0.0.BT-Settl.spec.7.dat'.format(T=Teff/100.0, g=logg)
    #print(name)
    
    return os.path.join(repo, name)

def make_dicts(param_list):
    """
    makes array of dictionaries with parameters to load
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
    makes a list of required atmospheric parameters to be retreived, also records which params need interpolation. Fixing FeH and aM = 0.0 for now.
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
  #  print (param_list)
    return param_list, params_to_interp

def get_grids():
    """
    arrays storing the available phoenix spectra. NOTE: check svo is the same
    """
    phxTgrid = np.arange(1200, 7001, 100)
    phxTgrid = np.hstack([np.arange(2300,7000,100),
                   np.arange(7000,12001,200)])
    phxggrid = np.arange(0.0, 6.1, 0.5)
    phxZgrid = np.array([0.0])
    phxagrid = np.array([0.0])
    return phxTgrid,phxggrid,phxZgrid, phxagrid

def get_models(repo,param_dicts):
    """
    Returns "spectra" param_dicts but with the model flux added to each dictionary
    """
    spectra = []
    for params in param_dicts:
        Teff, logg, FeH, aM = params['Teff'], params['logg'], params['FeH'], params['aM']
        filepath = make_filepath(Teff, logg, repo=repo) 
        wavelength, flux =  extract_spectrum(filepath)
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
           # print(len(s['wavelength']), s['wavelength'][0], s['wavelength'][-1])
            fluxes.append(s['flux'])
        else:
           # print(len(s['wavelength']), s['wavelength'][0], s['wavelength'][-1])
            fi = interp1d(s['wavelength'], s['flux'], fill_value='extrapolate')(wavelength)
            fluxes.append(fi)
        
    if len(params_to_interp) == 1:
        in_vals = [s[params_to_interp[0]] for s in spectra]
        new_flux = interp1d(in_vals, fluxes, axis=0, fill_value='extrapolate')(star_params[params_to_interp[0]])
    else:
        out_vals = [star_params[p] for p in params_to_interp]
        in_vals = [[s[p] for p in params_to_interp] for s in spectra]
     #   print(in_vals)
      #  print(out_vals)
       # print(len(fluxes))
        new_flux = griddata(in_vals, fluxes, out_vals)[0]
    return wavelength, new_flux

    
def save_to_ecsv(star,wavelength, flux, error, save_path, star_params, normfac, make_error):
    """
  #  save the new model to an ecsv file
  #  """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    if make_error:
        metadata = {'OBJECT':star, 'TEFF':star_params['Teff'], 'TEFF_e':star_params['Teff_e'], 'LOGG':star_params['logg'], 'NORMFAC':normfac}
        savedat = Table([wavelength*u.AA, flux, error], names=['WAVELENGTH', 'FLUX', 'ERROR'], meta=metadata)
    else:
        metadata = {'OBJECT':star, 'TEFF':star_params['Teff'], 'LOGG':star_params['logg'], 'NORMFAC':normfac}
        savedat = Table([wavelength*u.AA, flux], names=['WAVELENGTH', 'FLUX'], meta=metadata)
    star = star.replace(' ', '')
    #ascii.write(savedat, save_path+star+'_phoenix_interpolated.ecsv', overwrite=True, format='ecsv')
    savedat.write(save_path+star+'_phoenix_interpolated.ecsv', overwrite=True, format='ascii.ecsv')
    
def plot_spectrum(wavelength, flux, star, normfac):
    plt.figure(star, figsize=(5, 5))
    plt.subplot(211)
    plt.plot(wavelength, flux, label = 'Flux at stellar surface')
    plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=1, frameon=True)
    plt.subplot(212)
    plt.plot(wavelength, flux*normfac, label = 'Flux at Earth')
    plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
    plt.xlabel('Wavelength (\AA)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=1, frameon=True)
    plt.tight_layout()
    plt.show()

def get_existing_model(star_params, repo):
    """
    Get the flux if there's already a good phoenix model
    """
    Teff, logg, FeH, aM = star_params['Teff'], star_params['logg'], star_params['FeH'], star_params['aM']
    file_path = make_filepath(Teff, logg, FeH, aM, repo=repo)
    wavelength, flux = extract_spectrum(filepath)
    return wavelength, flux
  
def air_to_vac(w_air, flux, flux_interp = False):
    """
    Converts the air wavelengths to vaccum wavelengths via the formular from https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    """
    if w_air[0] == 0.0: #correct a divide by zero problem by adding a very small number to the first wavelength element
        w_air[0] += 0.01 * w_air[1]
        print(w_air[0])
    print(w_air[0])   
    s = 1e4/w_air
    n = 1. + 0.00008336624212083 + (0.02408926869968 / (130.1065924522 - s**2)) + (0.0001599740894897 / (38.92568793293 - s**2))
    w_vac = w_air * n
    if flux_interp: #interpolate flux back onto old wavelength grid
        flux = interp1d(w_vac, flux, fill_value='extrapolate')(w_air)
        w_vac = w_air
    return w_vac, flux
    
    
def extract_spectrum(filepath):
    """
    Open and extract a svo txt file. So much easier than before!
    """
    try:
        w_raw, f_raw = np.loadtxt(filepath, unpack=True)
    except:
        print ('model {} not in repo'.format(os.path.split(filepath)[1]))
        sys.exit(1)
    return w_raw, f_raw

def build_spectrum(repo, star_params, tgrid, ggrig,fgrid, agrid):
    """
    Returns the phoenix spectrum for the given parameters, or interpolates a new one.
    """
    param_list, params_to_interp = make_param_list(star_params, [tgrid, ggrig,fgrid, agrid])
    if len(params_to_interp) == 0: #i.e. if there's an existing model
        print('phoenix model available')
        wavelength, flux = get_existing_model(star_params, repo)
    else:
        param_dicts = make_dicts(param_list)
        spectra = get_models(repo,param_dicts)
        wavelength, flux = interp_flux(spectra, params_to_interp, star_params)
    wavelength, flux = wavelength[wavelength >= 501.0], flux[wavelength >= 501.0] #spectrum does funny things at lambda < 501
    return wavelength, flux
    
def make_phoenix_spectrum(star, save_path, repo, star_params, save_ecsv=False, plot=False, to_vac=False, make_error=False):
    """
    Main array. Takes a list of stellar parameters and makes a phoenix spectrum out of. Save_path is where you want the final spectrum to go, repo is where downloaded phoenix files go. wave_file is where the wavelength array is 
    """
    tgrid, ggrig,fgrid, agrid = get_grids()
    wavelength, flux = build_spectrum(repo,star_params, tgrid, ggrig,fgrid, agrid)
    error = np.zeros(len(flux))
    normfac = find_normfac(star_params['Radius'], star_params['Distance'])
    if make_error:
        teff, teff_e = star_params['Teff'], star_params['Teff_e']
        star_params['Teff'] = teff+teff_e
        wave_up, flux_up = build_spectrum(repo,star_params, tgrid, ggrig,fgrid, agrid)
        star_params['Teff'] = teff-teff_e
        wave_down, flux_down = build_spectrum(repo, star_params, tgrid, ggrig,fgrid, agrid)
        star_params['Teff'] = teff
        if not np.array_equal(wave_up, wavelength): #check if wavelength grids are the same. Should be, but you might get unlucky and hit a grid change
            flux_up = interp1d(wave_up, flux_up,fill_value='extrapolate')(wavelength)
        if not np.array_equal(wave_down, wavelength): 
            flux_down = interp1d(wave_down, flux_down,fill_value='extrapolate')(wavelength)
        error = np.mean([abs(flux_up-flux), abs(flux-flux_down)], axis=0)
    if to_vac:
        wavelength, flux = air_to_vac(wavelength, flux)
    if save_ecsv:
        save_to_ecsv(star, wavelength, flux, error, save_path, star_params, normfac, make_error)
    if plot:
        plot_spectrum(wavelength, flux, star, normfac)
    # print(make_error)
    if make_error == True and plot==True:
        plt.figure()
        plt.plot(wavelength, flux, zorder=10)
        plt.plot(wavelength, flux_up)
        plt.plot(wavelength, flux_down)
        plt.yscale('log')
        plt.xscale('log')
        plt.show()
    if make_error == True:
        return wavelength, flux, error
    else:
        return wavelength, flux


def find_normfac(radius, distance):
    """
    finds the scaling factor for the spectrum
    """
    return ((radius.to(u.cm)/distance.to(u.cm))**2).value
 

def smear(w,f, R, w_sample=1):
    '''
    Smears a model spectrum with a gaussian kernel to the given resolution, R.
    Adapeted from https://github.com/spacetelescope/pysynphot/issues/78

    Parameters
    -----------

    w,f:  spectrum to smear

    R: int
        The resolution (dL/L) to smear to

    w_sample: int
        Oversampling factor for smoothing

    Returns
    -----------

    sp: PySynphot Source Spectrum
        The smeared spectrum
    '''

    # Save original wavelength grid and units
    w_grid = w
    

    # Generate logarithmic wavelength grid for smoothing
    w_logmin = np.log10(np.nanmin(w_grid))
    w_logmax = np.log10(np.nanmax(w_grid))
    n_w = np.size(w_grid)*w_sample
    w_log = np.logspace(w_logmin, w_logmax, num=n_w)

    # Find stddev of Gaussian kernel for smoothing
    R_grid = (w_log[1:-1]+w_log[0:-2])/(w_log[1:-1]-w_log[0:-2])/2
    sigma = np.median(R_grid)/R
    if sigma < 1:
        sigma = 1

    # Interpolate on logarithmic grid
    f_log = np.interp(w_log, w_grid, f)

    # Smooth convolving with Gaussian kernel
    gauss = Gaussian1DKernel(stddev=sigma)
    f_conv = convolve_fft(f_log, gauss)

    # Interpolate back on original wavelength grid
    f_sm = np.interp(w_grid, w_log, f_conv)

    # Write smoothed spectrum back into Spectrum object
    return w_grid, f_sm

def test():
    star = 'Trappist-1_test' 
    repo = '/media/david/5tb_storage1/btsettl_test/t1_test/' #where the files to be interpolated are
    save_path = 'test_output/' #where you want the ecsv files to be saved
    star_params = {'Teff': 2628, 'logg': 5.21, 'FeH':0.0, 'aM':0.0, 'Radius':1.16*u.R_jup, 'Distance':12.43*u.pc}
    w_phx, f_phx = make_phoenix_spectrum(star, save_path, repo, star_params, save_ecsv=True, plot=True)

def test_load():
    path = 'test_output/'
    spectra = glob.glob('{}*.ecsv'.format(path))
    if len(spectra) > 0:
        for spectrum in spectra:
            data = Table.read(spectrum)
            plot_spectrum(data['WAVELENGTH'], data['FLUX'], data.meta['OBJECT'],data.meta['NORMFAC'])
    else:
        print('No ecsv files in path')
        
    
def test2():
    G = const.G
    M = const.M_sun.to(u.kg)
    R = const.R_sun.to(u.m)
    mass = 1.22

    teff = 5675
    teff_e = 75
    radius = 1.38
    distance = (1000/7.8288000)
    save_path = 'models/'
    star = 'NGTS_10'
    g = ((G*mass*M)/(radius*R)**2).to(u.cm/u.s**2)
    print(np.log10(g.value))
    repo = '/media/david/2tb_ext_hd/hddata/mega_muscles/data-vacuum/'

    star_params = {'Teff': teff, 'logg': np.log10(g.value), 'FeH': 0.00, 'aM': 0, 'Radius':radius*u.R_sun, 'Distance':distance*u.pc, 'Teff_e':teff_e}
    pw, pf, pe = make_phoenix_spectrum(star, save_path, repo, star_params, save_ecsv=False, plot=False, make_error=True)
    normfac = ((radius*R)/(distance*u.pc.to(u.m)))**2

    


# test2()
# test_load() 