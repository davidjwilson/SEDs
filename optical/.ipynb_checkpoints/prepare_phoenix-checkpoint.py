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
import lzma
import requests

"""
@author David Wilson

version 3 20191211

Script to retreive phoenix models and interpolate them onto the correct values. "phxurl" and "fetchphxfile" adaped from Parke Loyds scripts. Adapetd further to use the lyon models
"""

def phxurl(Teff, logg=4.5, repo='ftp'):
    """
    Constructs the URL for the phoenix spectrum file for a star with effective
    temperature Teff, log surface gravity logg. FeH and aM are fixed at 0.0 as the lyon database does not have other options

    Does not check that the URL is actually valid, and digits beyond the
    precision of the numbers used in the path will be truncated.
    """
    phoenixbaseurl = 'https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/SPECTRA/'
  
    name = 'lte{T:05.1f}-{g:3.1f}-0.0a+0.0.BT-Settl.spec.7'.format(T=Teff/100.0, g=logg)
    print(name)

    if repo == 'ftp':
        return phoenixbaseurl + name
    else:
        return os.path.join(repo, name)

def fetchphxfile(Teff, logg, FeH, aM, repo, source = 'lyon'):
    #keeping a source option in here in case I want to use the gottigen versions again
    if source == 'lyon':
        loc, ftp = [phxurl(Teff, logg, repo=r) for r in [repo, 'ftp']]
    #urlretrieve(ftp, loc)
    r = requests.get(ftp)
    with open(loc,'wb') as f:
        f.write(r.content)

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
    print (param_list)
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

def unzip_file(filepath):
    """
    Moving the lzma work to a separate directory
    """
    nameout = filepath[:-3]
    with lzma.open(filepath) as f, open(filepath, 'wb') as fout: #https://stackoverflow.com/a/33718185
        file_content = f.read()
        fout.write(file_content)
    return nameout

def extract_spectrum(filepath):
    """
    Extracts the spectrum from the Lyon files, which is non-trivial. Adapts code by JSP. 
    """
    #nameout = unzip_file(filepath)
    wavemin, wavemax, DF = 1000, 1000000, 8
    
    try:
        phoenixR = ascii.read(filepath,format="fixed_width_no_header",col_starts=(0,14),col_ends=(12,25),delimiter=" ",names=('Wave','Spec'))
        ph1, jj = np.unique(np.array(phoenixR['Wave']),return_index=True)
        phoenix = np.zeros((len(ph1),2))
        for kk in range(len(jj)):
            phoenix[kk,1] = np.float64(phoenixR['Spec'][jj[kk]].replace("D","E"))
        phoenix[:,0] = ph1
    except:
            try:
                print("File Badly Formatted --- trying again...")
                phoenixR = ascii.read(filepath,format="fixed_width_no_header",col_starts=(0,13),col_ends=(12,24),delimiter=" ",names=('Wave','Spec'))
                ph1, jj = np.unique(np.array(phoenixR['Wave']),return_index=True)
                phoenix = np.zeros((len(ph1),2))
                for kk in range(len(jj)):
                    phoenix[kk,1] = np.float64(phoenixR['Spec'][jj[kk]].replace("D","E"))
                phoenix[:,0] = ph1
            except:
                print("... and again ... ")
                phoenixR = ascii.read(filepath,format="no_header",delimiter=" ")
                temp = np.zeros(len(phoenixR['col1']))
                for kk in range(len(temp)):
                    temp[kk] = np.float64(phoenixR['col1'][kk].replace("D","E"))
                ph1, jj = np.unique(temp,return_index=True)
                phoenix = np.zeros((len(ph1),2))
                for kk in range(len(jj)):
                    phoenix[kk,0] = np.float64(phoenixR['col1'][jj[kk]].replace("D","E"))
                    phoenix[kk,1] = np.float64(phoenixR['col2'][jj[kk]].replace("D","E"))
    
 
    ind = np.where( (phoenix[:,0] <= wavemax) & (phoenix[:,0] >= wavemin))[0]  
    xraw = phoenix[ind,0]
    yraw = np.power(10.,phoenix[ind,1] + DF)
    return xraw, yraw

    
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
        if os.path.exists(file_path) == False: #only download if we need to
            fetchphxfile(Teff, logg, FeH, aM, repo=repo)
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
            print(len(s['wavelength']), s['wavelength'][0], s['wavelength'][-1])
            fluxes.append(s['flux'])
        else:
            print(len(s['wavelength']), s['wavelength'][0], s['wavelength'][-1])
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

    
def save_to_ecsv(star,wavelength, flux, save_path, star_params):
    """
  #  save the new model to an ecsv file
  #  """
    normfac = find_normfac(star_params['Radius'], star_params['Distance'] )
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    metadata = {'TEFF':star_params['Teff'], 'LOGG':star_params['logg'], 'NORMFAC':normfac}
    savedat = Table([wavelength*u.AA, flux], names=['WAVELENGTH', 'FLUX'], meta=metadata)
    star = star.replace(' ', '')
    #ascii.write(savedat, save_path+star+'_phoenix_interpolated.ecsv', overwrite=True, format='ecsv')
    savedat.write(save_path+star+'_phoenix_interpolated.ecsv', overwrite=True, format='ascii.ecsv')
    
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
        
    
def make_phoenix_spectrum(star, save_path, repo, star_params, save_ecsv=False, plot=False):
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
        print(param_dicts)
        spectra = get_models(repo,param_dicts)
        wavelength, flux = interp_flux(spectra, params_to_interp, star_params) 
    if save_ecsv:
        save_to_ecsv(star, wavelength, flux, save_path, star_params)
    if plot == True:
        plot_spectrum(wavelength, flux, star)
    return wavelength, flux

def distance_scale(radius, distance, flux):
    """
    Scales the phoenix model using the distance to the star
    """
    scale = (radius.to(u.cm)/distance.to(u.cm))**2
    return flux * scale

def find_normfac(radius, distance):
    """
    finds the scaling factor for the spectrum
    """
    return (radius.to(u.cm)/distance.to(u.cm))**2
 


def test():
    star = 'Trappist-1_test'
    #star_params = load_star_params(star+'_ParamStats.txt', FeH=-0.4)   
    path = '/home/david/work/muscles/phoenix/file_dump/'
    #wave_file = path+ 'phoenix/wavegrid_hires.fits'
    repo = path + star+'_repo/'
    save_path = path + star+'_output/'
    #star_params = {'Teff':3300, 'logg': 4.00, 'FeH': -2.0, 'aM': 0.0}
    star_params = {'Teff': 2628, 'logg': 5.21, 'FeH': 0.00, 'aM': 0, 'radius':1.16*u.R_jup, 'distance':12.43*u.pc}
    w_phx, f_phx = make_phoenix_spectrum(star, save_path, repo, star_params, save_ecsv=False, plot=True)
   # ccd_path = 'g430l/odlm24010_sx1.fits'
    #normfac = phoenix_norm(star, w_phx, f_phx, ccd_path, plot=True)
    print(normfac)
    #facs = []
    #cuts = np.arange(3000,4000, 100)
    #for cut in cuts:
     #   normfac = phoenix_norm(cut,star, w_phx, f_phx, ccd_path, plot=False)
      #  facs.append(normfac)
    #plt.plot(cuts, facs/np.mean(facs))
    #plt.show()

#test()