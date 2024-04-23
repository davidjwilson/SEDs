import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.modeling import models, fitting

mg2 = [2796.352, 2803.53]
os = [6.08e-01, 3.03e-01]

def nuv_bb_model(amp_in, scale_in):
    mg_mod = models.Gaussian1D(amplitude=amp_in*u.erg/u.s/u.cm**2/u.AA, mean = mg2[0]*u.AA, stddev=0.5*u.AA, fixed = dict(mean=True, stddev=True)) + models.Gaussian1D(amplitude=amp_in/(os[0]/os[1])*u.erg/u.s/u.cm**2/u.AA, mean = mg2[1]*u.AA, stddev=0.5*u.AA, fixed = dict(mean=True, stddev=True))
    
    # def tiedamp(mod):
    #     amplitude = mod.amplitude_0.value/(os[0]/os[1])
    #     return amplitude
    # mg_mod.amplitude_1.tied = tiedamp
        
    bb_mod = models.BlackBody(temperature=8500*u.K, scale=scale_in*u.erg/u.s/u.cm**2/u.AA/u.sr, fixed=dict(temperature=True, scale=True))
    

    
    # nuv_mod = mg_mod+bb_mod
    mod_w = np.arange(1100, 3501, 0.1)*u.AA


    bb_mod = bb_mod(mod_w).value
    mg_mod = mg_mod(mod_w).value
    mod_w = mod_w.value

    mod_spec = bb_mod + mg_mod

    return mod_w, mod_spec

ntries = 10000


# amp_in = np.array([5.7, 4.9, 4.4])*1e-16
# amp_in_e = np.array([4.8, 5.5, 4.7])*1e-17
# scale_in = np.array([2.0, 2.0, 1.8])*1e-25
# scale_in_e = np.array([1.2, 1.3, 1.2])*1e-26

#with fixed line ratios
amp_in = np.array([7.638781293208803e-16, 6.568414968276414e-16, 5.841014739084463e-16])
amp_in_e = np.array([6.397343721256946e-17, 7.397898101536548e-17, 6.304554133610627e-17])
scale_in = np.array([1.9917254457999382e-25, 2.0459168790762651e-25, 1.8442737618238748e-25])
scale_in_e = np.array([1.2466942704625804e-26, 1.309791110479485e-26, 1.1731741923033736e-26])

# 5.714344493884374e-16 1.9967318694715743e-25 [4.77795165e-17 1.24331445e-26]
# 4.921099180036549e-16 2.0491163652767218e-25 [5.52970285e-17 1.30596556e-26]
# 4.369600366965082e-16 1.84750622443636e-25 [4.70821212e-17 1.17079044e-26]

#7.638781293208803e-16 1.9917254457999382e-25 [6.39734372e-17 1.24669427e-26]
#6.568414968276414e-16 2.0459168790762651e-25 [7.39789810e-17 1.30979111e-26]
#5.841014739084463e-16 1.8442737618238748e-25 [6.30455413e-17 1.17317419e-26]

fluxes = []

for i in range(len(amp_in)):
    if i == 2:
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

        w, f = nuv_bb_model(amp_in[i], scale_in[i])
        e = f * (fluxstd/fluxmean)
    
        savdat = Table([w*u.AA, f*u.erg/u.s/u.cm**2/u.AA, e*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
        savdat.write('model_spectra/epoch{}_nuv_mod.ecsv'.format(i+1), format='ascii.ecsv', overwrite=True)
    


plt.show()












