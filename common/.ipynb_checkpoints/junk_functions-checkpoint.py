    
def phoenix_norm(component_repo, star_params, plot=False): 
    """
    find the normalisation factor between the phoenix model and the stis ccd (ccd_path)
    """
    norm = Table.read(glob.glob(component_repo+'*phx*.ecsv')[0])
    radius, distance = star_params['radius'], star_params['distance']
    normfac = ((radius.to(u.cm)/distance.to(u.cm))**2).value
    print('PHOENIX NORMFAC =', normfac)
    update_norm(glob.glob(component_repo+'*phx*.ecsv')[0], glob.glob(component_repo+'*phx*.fits')[0], normfac)

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

def add_stis_optical(sed_table, component_repo, instrument_list):
    """
    Adds the G430L spectrum
    """
    g430l_path = glob.glob(component_repo+'*g430l*.ecsv')
    if len(g430l_path) > 0:
        g430l = Table.read(g430l_path[0])
        instrument_code, g430l = hst_instrument_column(g430l)
        instrument_list.append(instrument_code)
        g430l = normfac_column(g430l)
        g430l = g430l[g430l['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        sed_table = vstack([sed_table, g430l], metadata_conflicts = 'silent')
    return sed_table, instrument_list

def add_phx_spectrum(sed_table, component_repo, instrument_list):
    """
    Adds the scaled phoenix spectrum
    """
    phx_path = glob.glob(component_repo+'*phx*.ecsv')
    if len(phx_path) > 0:
        phx = Table.read(phx_path[0])
        instrument_code, phx = fill_model(phx, 'mod_phx_-----')
        instrument_list.append(instrument_code)
        phx = normfac_column(phx)
        phx = phx[phx['WAVELENGTH'] > max(sed_table['WAVELENGTH'])]
        phx['FLUX'] *= phx.meta['NORMFAC']
        sed_table = vstack([sed_table, phx], metadata_conflicts = 'silent')
    return sed_table, instrument_list
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
photometry
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
def blackbody_fit(phx, Teff):
    """Return a function that is a blackbody fit to the phoenix spectrum for the star. The fit is to the unnormalized
    phoenix spectrum, so the fit function values must be multiplied by the appropriate normalization factor to match
    the normalized spectrum. From PLs code"""

    
    # recursively identify relative maxima until there are fewer than N points
    N = 10000
    keep = np.arange(len(phx))
    while len(keep) > N:
        temp, = argrelmax(phx['FLUX'][keep])
        keep = keep[temp]

    efac = const.h * const.c / const.k_B / (Teff * u.K)
    efac  = efac.to(u.angstrom).value
    w = phx['WAVELENGTH']
    w = w[keep]
    planck_shape = 1.0/w**5/(np.exp(efac/w) - 1)
    y = phx['FLUX'][keep]

    Sfy = np.sum(planck_shape * y)
    Sff = np.sum(planck_shape**2)

    norm = Sfy/Sff

    return lambda w: norm/w**5/(np.exp(efac/w) - 1)
def bolo_integral(pan,phx,teff,uplim=np.inf, tail=False):
    """
    Calculates the bolometric integral flux of the SED by adding a blackbody fit to the end of the sed
    """
#     fit_unnormed = blackbody_fit(phx, teff)
    normfac = pan[-1]['NORMFAC'] 
    #Ibody = flux_integral(pan)[0]
    Ibody = np.trapz(pan['FLUX'], pan['WAVELENGTH'])
#     if tail:
#         Itail = normfac*quad(fit_unnormed, pan['WAVELENGTH'][-1], uplim)[0]
#         I = Ibody + Itail
    else:
        I = Ibody
#     print(I)

    return I