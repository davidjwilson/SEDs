import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u

"""
Calculating the EUV fluxes using the relationships of Linsky + 14 (https://ui.adsabs.harvard.edu/abs/2014ApJ...780...61L/abstract)

log[F(delta lambda) /F(lya)]=

10–20 nm (stars) 	 	−0.491 	 
20–30 nm (stars) 	 	−0.548 	 
30–40 nm (stars) 		−0.602 	 
40–50 nm (models) 	  	  	−2.294+0.258 log[f (Lyα)]
50–60 nm (models) 	  	  	−2.098+0.572 log[f (Lyα)]
60–70 nm (models) 	  	  	−1.920+0.240 log[f (Lyα)]
70–80 nm (models) 	  	  	−1.894+0.518 log[f (Lyα)]
80–91.2 nm (models) 	  	  	−1.811+0.764 log[f (Lyα)]
91.2–117 nm (models) 	  	  	−1.004+0.065 log[f (Lyα)]

"""

def euv_estimator(lya, distance, star='', save=False, plot=False):
        
    distance_conversion = ((1*u.au.to(u.m))/(distance*u.pc.to(u.m)))**2

    lya_1au = lya / distance_conversion

    w1 = np.array([100,200,300,400,500,600,700,800,912], dtype=float) #A
    w2 = np.array([200,300,400,500,600,700,800,912,1170], dtype=float)
    bandwidth = w2-w1

    a = np.array([-0.491,-0.548,-0.602,-2.294,-2.098,-1.920,-1.894,-1.811,-1.004], dtype=float)
    b = np.array([ 0.,    0.,    0.,    0.258, 0.572, 0.240, 0.518, 0.764, 0.065], dtype=float)

    #log(f/lya) = a + b *log(lya)
    f = a + b*np.log10(lya_1au)

   # print('Total EUV=',np.sum(f))
    f = (lya_1au * 10**f)/bandwidth

    f *= distance_conversion

    #extrapolate onto 1A grid
    wav = np.arange((w1[0])+0.5, (w2[-1])+0.5, 1.0)
    flux = []
    for w1i, w2i, fi in zip(w1, w2,f):
        for wi in wav:
            if wi > w1i and wi < w2i :
                flux.append(fi)

    if save == True:
        data = Table([wav*u.AA,  flux*u.erg/u.s/u.cm**2], names=['WAVELENGTH', 'FLUX'])
        ascii.write(data, star+'_1Aeuv_estimate.ecsv', delimiter=',', format='ecsv', overwrite=True)

    if plot == True:
        plt.figure(star+'_EUV', figsize=(8,6))
        plt.subplots_adjust(top=0.99, right=0.99)
        plt.plot(wav, flux)
        plt.xlabel('Wavelength (\AA)', size=20)
        plt.ylabel('Flux (erg s$^{-1}$\AA$^{-1}$cm$^{-2}$)', size=20)
        plt.yscale('log')
        plt.show()
    
    return(wav, np.array(flux)) 

        
#star = 'GJ674'
#lya = 2.9*2.06e-12 #erg /s/cm2 for GJ674, 2.9 is flux scaling factor 
#distance = 4.54
#star = 'GJ176'
#lya = 3.9e-13
#distance = 9.3 
#star = 'HD_97658'
#lya = 9.1e-13
#distance = 21.1
#print(lya)
#star = 'TRAPPIST-1'
#distance = 12.1
#lya = 8e-15
#euv_estimator(lya,  distance, star=star, plot=True, save=True)