import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
import astropy.constants as const

plt.figure()

#COSG130M
instrument = 'COS_G130M'
data = readsav('../COS/TRAPPIST1_G130M_Mm1_NOSCL_10dec2018.sav')
mask = (data['Wave'] <1207)|(data['Wave'] >1225)&(data['Wave'] <1301)|(data['Wave'] >1307.)&(data['Wave'] <1355)|(data['Wave'] >1356)
f = convolve(data['flux'][mask],Box1DKernel(5))
e = convolve(data['err'][mask],Box1DKernel(5))/(5**0.5)
#plt.step(data['wave'][mask], data['flux'][mask])
plt.step(data['wave'][mask],f, where='mid', c='C0')
plt.step(data['wave'][mask],e, where='mid', c='C1')

end_130 = data['wave'][-1]


#COSG160M
instrument = 'COS_G160M'
data = readsav('../COS/TRAPPIST1_G160M_3orb_Mm1_NOSCL_09dec2018.sav')
mask = (data['wave'] > end_130)
f = convolve(data['flux'][mask],Box1DKernel(5))
e = convolve(data['err'][mask],Box1DKernel(5))/(5**0.5)
#plt.step(data['wave'][mask], data['flux'][mask])
plt.step(data['wave'][mask],f, where='mid', c='C0')
plt.step(data['wave'][mask],e, where='mid', c='C1')

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
plt.step(cnw, cnf, where='mid', c='C0')

plt.axhline(0, ls='--', c='k')

lines = [1561, 1657.00751,1657.37863,1657.90661,1486,1666.153,1407.382,1670.787, 1533.45  ,1393.755, 1402.770]
line_names = ['Ci', 'Ci','Ci','Ci','Oiii','Oiv','Alii','Siii', 'Siiv','Siiv']

dv = -56300 #Bourrier+17a 
c = const.c.value
dshift = (1.0+(dv/c))
[plt.axvline(line*dshift, ls='--',c='r') for line in lines]
#plt.ylim(5e-16,-3e-16)

plt.show()

           