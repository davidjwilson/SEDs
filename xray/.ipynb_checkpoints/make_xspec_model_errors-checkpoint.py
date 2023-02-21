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

ntries = 10000

xtab = Table.read('mega-muscles_x-ray_params.csv')
hlsppath = 'hlsp/'

for x in xtab[6:]:
    n = 0
    fluxes = []
    while n < ntries:
        print(x['Star'], 'Iteration:', n)
        if x['kt3'] != 0.0:
            print(x['Star'], '3temps')
            if x['fxeu'] == 0.0:
                fx = x['fx']
            else:
                fx = np.random.normal(x['fx'], np.mean([x['fxeu'], x['fxel']]))
                if fx < 0.0:
                    fx = x['fx']
            if x['kt1eu'] == 0.0:
                kt1 = x['kt1']
            else:
                kt1 = np.random.normal(x['kt1'], np.mean([x['kt1eu'], x['kt1el']]))
                if kt1 < 0.008:
                    kt1 = 0.008
            if x['kt2eu'] == 0.0:
                kt2 = x['kt2']
            else:
                kt2 = np.random.normal(x['kt2'], np.mean([x['kt2eu'], x['kt2el']]))
                if kt2 < 0.008:
                    kt2 = 0.008
            if x['kt3eu'] == 0.0:
                kt3 = x['kt3']
            else:
                kt3 = np.random.normal(x['kt3'], np.mean([x['kt3eu'], x['kt3el']]))
                if kt3 < 0.008:
                    kt3 = 0.008
            abd, nh = x['abd'], x['nh'] * 0.01
            




            # fx, kt1, kt2, kt3, abd, nh = x['fx'], x['kt1'], x['kt2'], x['kt3'], x['abd'], x['nh']*0.01

            mod = Model('(apec+apec+apec)*phabs', setPars={1:kt1, 2:abd,
                                                          5:kt2, 6:abd,
                                                          9:kt3, 10:abd,
                                                          13:nh})
            Plot.xAxis = "angstrom"
            Plot.perHz = False
            Plot.area=True
            AllModels.setEnergies(".3 10. 1000")
            flux = AllModels.calcFlux(".3 10")
            fluxnum = mod.flux[0]
            norm = (fx*1e-14)/fluxnum
            mod.setPars({4:norm})
            mod.setPars({8:norm})
            mod.setPars({12:norm}) 

        elif x['kt2'] != 0.0:
            print(x['Star'], '2temps')
            if x['fxeu'] == 0.0:
                fx = x['fx']
            else:
                fx = np.random.normal(x['fx'], np.mean([x['fxeu'], x['fxel']]))
                if fx < 0.0:
                    fx = x['fx']
            if x['kt1eu'] == 0.0:
                kt1 = x['kt1']
            else:
                kt1 = np.random.normal(x['kt1'], np.mean([x['kt1eu'], x['kt1el']]))
                if kt1 < 0.008:
                    kt1 = 0.008
            if x['kt2eu'] == 0.0:
                kt2 = x['kt2']
            else:
                kt2 = np.random.normal(x['kt2'], np.mean([x['kt2eu'], x['kt2el']]))
                if kt1 < 0.008:
                    kt1 = 0.008
            abd, nh = x['abd'], x['nh'] * 0.01
            # fx, kt1, kt2, abd, nh = x['fx'], x['kt1'], x['kt2'], x['abd'], x['nh']*0.01

            mod = Model('(apec+apec)*phabs', setPars={1:kt1, 2:abd,
                                                          5:kt2, 6:abd,
                                                          9:nh})
            Plot.xAxis = "angstrom"
            Plot.perHz = False
            Plot.area=True
            AllModels.setEnergies(".3 10. 1000")
            flux = AllModels.calcFlux(".3 10")
            fluxnum = mod.flux[0]
            norm = (fx*1e-14)/fluxnum
            mod.setPars({4:norm})
            mod.setPars({8:norm})

        else:
            print(x['Star'], '1temp')
            if x['fxeu'] == 0.0:
                fx = x['fx']
            else:
                fx = np.random.normal(x['fx'], np.mean([x['fxeu'], x['fxel']]))
                if fx < 0.0:
                    fx = x['fx']
            if x['kt1eu'] == 0.0:
                kt1 = x['kt1']
            else:
                kt1 = np.random.normal(x['kt1'], np.mean([x['kt1eu'], x['kt1el']]))
                if kt1 < 0.008:
                    kt1 = 0.008
            abd, nh = x['abd'], x['nh'] * 0.01
            # fx, kt1, abd, nh = x['fx'], x['kt1'], x['abd'], x['nh']*0.01

            mod = Model('(apec)*phabs', setPars={1:kt1, 2:abd,
                                                          5:nh})
            Plot.xAxis = "angstrom"
            Plot.perHz = False
            Plot.area=True
            AllModels.setEnergies(".3 10. 1000")
            flux = AllModels.calcFlux(".3 10")
            fluxnum = mod.flux[0]
            norm = (fx*1e-14)/fluxnum
            mod.setPars({4:norm})

        if x['tel'] == 'xmm':
            AllModels.setEnergies("0.102 2.5 2401")
        else:
            AllModels.setEnergies("0.1 2.5 2400")
        # AllModels.setEnergies("0.1 2.5 2400")
        Plot("model")
        xVals = Plot.x()
        yVals = Plot.model()
        wx = xVals*u.AA
        fx  = (yVals * (u.photon/u.s/u.cm**2/u.AA)).to(u.erg/u.s/u.cm**2/u.AA, equivalencies=u.spectral_density(wx))
        # plt.plot(wx, fx, c='C0', alpha=0.2)
        # plt.xlim(5, 50)
        # plt.yscale('log')
        
        fluxes.append(fx[::-1])
        
        n += 1

    starx = glob.glob('{}*apec*{}*.fits'.format(hlsppath, x['Star'].lower()))
    if len(starx) > 0:
        print('found the hlsp')
        data = fits.getdata(starx[0], 1)
        # plt.plot(data['WAVELENGTH'], data['FLUX'], c='C1', ls='--')
        
        fluxes = np.array(fluxes)
        fluxmean, fluxstd = np.mean(fluxes, axis=0), np.std(fluxes, axis=0)
        errper = fluxstd/fluxmean
        flux_err = data['FLUX'] * errper
        # plt.plot(data['WAVELENGTH'], flux_err, c='C2')
        
        # print('SN', np.median(1/errper))
        
        savdat = Table((data['WAVELENGTH']*u.AA, data['FLUX']*u.erg/u.s/u.cm**2/u.AA , flux_err*u.erg/u.s/u.cm**2/u.AA), 
                       names=['WAVELENGTH', 'FLUX', 'ERROR'])
        ascii.write(savdat, 'new_apec_specs/{}_apec_errs.ecsv'.format(x['Star']), format='ecsv', overwrite=True)

print('DONE')
    # plt.show()
