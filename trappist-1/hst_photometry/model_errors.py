import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.modeling import models, fitting

mg2 = [2796.352, 2803.53]
os = [6.08e-01, 3.03e-01]

def nuv_bb_model(amp_in, scale_in):
    mg_mod = models.Gaussian1D(amplitude=amp_in*u.erg/u.s/u.cm**2/u.AA, mean = mg2[0]*u.AA, stddev=0.5*u.AA, fixed = dict(mean=True, stddev=True)) + models.Gaussian1D(amplitude=amp_in*u.erg/u.s/u.cm**2/u.AA, mean = mg2[1]*u.AA, stddev=0.5*u.AA, fixed = dict(mean=True, stddev=True))
    
    def tiedamp(mod):
        amplitude = mod.amplitude_0.value/(os[0]/os[1])
        return amplitude
    mg_mod.amplitude_1.tied = tiedamp
        
    bb_mod = models.BlackBody(temperature=8500*u.K, scale=scale_in*u.erg/u.s/u.cm**2/u.AA/u.sr, fixed=dict(temperature=True, scale=True))
    

    
    # nuv_mod = mg_mod+bb_mod
    mod_w = np.arange(1990, 3501, 0.1)*u.AA


    bb_mod = bb_mod(mod_w).value
    mg_mod = mg_mod(mod_w).value
    mod_w = mod_w.value

    mod_spec = bb_mod + mg_mod

    return mod_w, mod_spec

ntries = 10000
amp_in = np.array([5.7, 4.9, 4.4])*1e-16
amp_in_e = np.array([4.8, 5.5, 4.7])*1e-17
scale_in = np.array([2.0, 2.0, 1.8])*1e-25
scale_in_e = np.array([1.2, 1.3, 1.2])*1e-26

# 5.714344493884374e-16 1.9967318694715743e-25 [4.77795165e-17 1.24331445e-26]
# 4.921099180036549e-16 2.0491163652767218e-25 [5.52970285e-17 1.30596556e-26]
# 4.369600366965082e-16 1.84750622443636e-25 [4.70821212e-17 1.17079044e-26]

fluxes = []

for i in range(len(amp_in)):
    amps = np.random.normal(amp_in[i], amp_in_e[i], ntries)
    scales = np.random.normal(scale_in[i], scale_in_e[i], ntries)
    n = 0

    while n < ntries:
        w, f = nuv_bb_model(amps[n], scales[n])
        fluxes.append(f)
        print(i, n)
        n +=1

    fluxmean, fluxstd = np.mean(fluxes, axis=0), np.std(fluxes, axis=0)

    fig, ax = plt.subplots(nrows=2)
    ax[0].plot(w, fluxmean)
    ax[0].plot(w, fluxstd, alpha=0.5)
    ax[1].step(w, fluxmean/fluxstd, where='mid')
    ax[0].set_ylabel('Flux (10$^{-18}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
    ax[1].set_ylabel('SNR')
    ax[1].set_xlabel('Wavelength (\AA)')
    ax[0].set_title('Epoch {}'.format(i+1))
    fig.tight_layout()

    savdat = Table([w*u.AA, fluxmean*u.erg/u.s/u.cm**2/u.AA, fluxstd*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
    savdat.write('model_spectra/epoch{}_nuv_mod.ecsv'.format(i+1), format='ascii.ecsv', overwrite=True)



plt.show()












