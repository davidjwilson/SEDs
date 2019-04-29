import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel

plt.figure()

#COSG130M
instrument = 'COS_G130M'
data = readsav('../COS/TRAPPIST1_G130M_Mm1_NOSCL_10dec2018.sav')
glowmask = (data['Wave'] <1207)|(data['Wave'] >1225)&(data['Wave'] <1301)|(data['Wave'] >1307.)&(data['Wave'] <1355)|(data['Wave'] >1356)
f = convolve(data['flux'][glowmask],Box1DKernel(5))
#plt.step(data['wave'][mask], data['flux'][mask])
plt.step(data['wave'][glowmask],f)
end_130 = data['wave'][-1]


#COSG160M
instrument = 'COS_G160M'
data = readsav('../COS/TRAPPIST1_G160M_3orb_Mm1_NOSCL_09dec2018.sav')
mask = (data['wave'] > end_130)
f = convolve(data['flux'][mask],Box1DKernel(5))
#plt.step(data['wave'][mask], data['flux'][mask])
plt.step(data['wave'][mask],f)

#plt.step(data['wave'][mask], data['err'][mask])
end_160 = data['wave'][mask][-1]

instrument = 'COS_G230L'
specpath = '../../../trappist-1/hst/data/ldlm42010_x1dsum.fits'
cnw = np.array([], dtype=float)
cnf = np.array([], dtype=float)
cne = np.array([], dtype=float)
cndq = np.array([], dtype=int)
ndata = fits.getdata(specpath,1)[0:2]
for dt in ndata:
    cnw= np.concatenate((cnw, dt['WAVELENGTH']))
    cnf = np.concatenate((cnf, dt['FLUX']))
    cne = np.concatenate((cne, dt['ERROR']))
    cndq = np.concatenate((cndq, dt['DQ']))
mask = (cnw > end_160)
cnw, cnf, cne, cndq = cnw[mask], cnf[mask], cne[mask], cndq[mask] 
plt.step(cnw, cnf)

plt.show()

           