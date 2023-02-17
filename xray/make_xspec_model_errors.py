import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from  xspec import *

#trappist-1


kts = np.array([0.2, 0.4])
ktes = np.array([0, 0])
norms = np.array([0.35, 0.07])*1e-4
normes = np.array([0.2,0.02])*1e-4



nsteps = 10000
i = 0
kt1s = np.random.normal(kts[0], ktes[0], nsteps)
kt2s = np.random.normal(kts[1], ktes[1], nsteps)
norm1s = np.random.normal(norms[0], normes[0], nsteps)
norm2s = np.random.normal(norms[1], normes[1], nsteps)
# print(len(norm2s[norm2s < 0]))

fluxes = []

plt.figure('models')

while i < nsteps:
    print(i)
    if kt1s[i] > 0 and kt2s[i] >0 and norm1s[i] > 0 and norm2s[i] > 0: #some small amount go below 0
        mod = Model('phabs*(apec+apec)', setPars={1:1e-3, 2:kt1s[i], 3:0.4, 5:norm1s[i], 6:kt2s[i], 7:0.4, 9:norm2s[i]})
        Plot.xAxis = "angstrom"
        Plot.perHz = False
        Plot.area=True
        fluxnum = mod.flux[0]
        AllModels.setEnergies("0.1 2.5 2400")
        Plot("model")
        xVals = Plot.x()
        yVals = Plot.model()
        wx = xVals*u.AA
        fx  = (yVals * (u.photon/u.s/u.cm**2/u.AA)).to(u.erg/u.s/u.cm**2/u.AA, equivalencies=u.spectral_density(wx))
        plt.plot(wx, fx, alpha=0.01, c='C0', rasterized=True)
        fluxes.append(fx.value)
    i +=1
plt.xlim(5, 50)
plt.yscale('log')
# plt.show()

# print(len(wx))

plt.figure('final')

flux = np.mean(fluxes, axis=0)
flux_err = np.std(fluxes, axis=0)
# print(len(flux))
plt.plot(wx, flux)
plt.plot(wx, flux+flux_err, c='C0', alpha=0.5)
plt.plot(wx, flux-flux_err, c='C0', alpha=0.5)
plt.yscale('log')

savdat = Table([wx[::-1], flux[::-1]*u.erg/u.s/u.cm**2/u.AA, flux_err[::-1]*u.erg/u.s/u.cm**2/u.AA], names=['WAVELENGTH', 'FLUX', 'ERROR'])
ascii.write(savdat, 't1_test.ecsv', format='ecsv', overwrite=True)

plt.show()