import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy import interpolate
from astropy.convolution import convolve, Box1DKernel,convolve_fft, Gaussian1DKernel
from astropy.modeling import models, fitting

data = Table.read('../combined/TRAPPIST-1_hst+phoenix_v1.ecsv')



#all_lines = [ciii_1175, feii_1200, siiii_1206, nv_1240, ci_1288+feii_1292, cii_1335, feii_1527, civ_1550, mgii_2800]

#making spectrum for plot
w, f, e  = data['WAVELENGTH'], data['FLUX'], data['ERROR']
lya_edges = [1207, 1225]
glowmask = (w <lya_edges[0])|(w > lya_edges[1])
w, f, e = w[glowmask], f[glowmask], e[glowmask]




lines = [1561, 1657.00751,1657.37863,1657.90661,1486,1666.153,1407.382,1670.787, 1533.45  ,1393.755, 1402.770]
line_names = ['Ci', 'Ci','Ci','Ci','Oiii','Oiv','Alii','Siii', 'Siiv','Siiv']

"""
line_edges = [[1174.6,1176.1],[1198.4,1200.95],[1205.9,1206.8],[1238.2,1238.9],[1242.3,1242.75],[1288.5,1288.8],[1333.9, 1334.55],[1334.9,1335.9],[1526.8, 1527.2],[1547.65, 1548.3],[1550.25,1550.75],[2794.,2797]]
b_bg = [[1173., 1174.],[1197., 1198.],[1205.5, 1205.75],[1237.5, 1238.],[1241.5, 1242.], [1287.5, 1288.],[1332., 1333.],[1332., 1333.],[1525., 1526.],[1546.,1547.],[1549.,1550],[2791,2793.]]

r_bg = [[1176.25, 1176.5],[1201., 1203.],[1204., 1205.],[1239.5, 1240.],[1243., 1244.],[1289., 1290.],[1334.6,1334.8],[1334.6,1334.8],[1527.5,1528.5],[1548.5,1549.],[1551,1552],[2798., 2799.]]

line_width = np.diff(line_edges)
b_width = np.diff(b_bg)
b_scale = line_width/b_width

r_width = np.diff(r_bg)
r_scale = line_width/b_width


#all_lines[4] = ci_1288


#all_lines[-1] =[2795.528]
"""
dv = -56300 #Bourrier+17a 
c = const.c.value
dshift = (1.0+(dv/c))
#print(dshift)

xlims = []

gauss_int  = []
gauss_e = []
#line_b = []
#line_r = []

#fitter = fitting.SLSQPLSQFitter()
fitter = fitting.LevMarLSQFitter()

plt.figure(figsize=(18,15))
plt.subplots_adjust(hspace=0.1, wspace=0.1,top=0.99, right=0.99, left = 0.05, bottom =0.05)
for i in range(len(lines)):
    pl = lines[i]
    xlims.append([4, 4])
    mask = (w > pl-7) & (w < pl+7)
    w1, f1, e1  = w[mask], f[mask], e[mask]
    plt.subplot(3,4,i+1)
    #if i < 11:
    f1 = convolve(f1,Box1DKernel(5))
    e1 = convolve(e1,Box1DKernel(5)) / (5**0.5)
    f1 *= 1e16
    e1 *= 1e16
    plt.step(w1, f1, where='mid')
    mask2 = (w1 > pl-xlims[i][0]) & (w1 < pl+xlims[i][1])
    top = max(f1[mask2])
    bottom = min(f1[mask2])
    #nmask = (lines > pl-5) & (lines < pl+5)
    plt.ylim(bottom*1.1, top*1.6)
    plt.xlim(pl - xlims[i][0], pl + xlims[i][1])
    #if i ==3:
     #   plt.ylabel('Flux ($10^{-16}$ erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)')
    #if i == 7:
     #   plt.xlabel('Wavelength (\AA)')
    #[plt.annotate(name,(line*dshift, top*1.05), xytext=(line*dshift, top*1.3), horizontalalignment='center') for name, line in zip(plot_name[nmask], lines[nmask])]
    #[plt.annotate('',(line*dshift, top*1.05), xytext=(line*dshift, top*1.25),arrowprops=dict(arrowstyle='-'), horizontalalignment='center') for line in all_lines[i]]
    plt.axhline(0, c='k', ls='--')
    guess_sigma = 0.05
    #if i == 11:
     #   guess_sigma = 0.1
    gg_init = models.Gaussian1D(1.0, pl*dshift, guess_sigma, bounds = {'mean': [pl*dshift-1,  pl*dshift+1], 'amplitude':[0, None] })
    #for line in all_lines[i][1:]:
     #   gg_init = gg_init + models.Gaussian1D(1.0, line*dshift, guess_sigma, bounds = {'mean': [all_lines[i][0]*dshift-1,  all_lines[i][0]*dshift+1], 'amplitude':[0, None]})
  #  if i !=2:
    gg_init = gg_init+ models.Const1D(-1)
    gg_fit = fitter(gg_init, w1, f1, maxiter=1000)
    plt.plot(w1, gg_fit(w1))
    
  #  l_flux = np.trapz(gg_fit(w1)[(w1 > line_edges[i][0]) & (w1 < line_edges[i][1])], w1[(w1 > line_edges[i][0]) & (w1 < line_edges[i][1])])
  #  gauss_int.append(l_flux)
  #  b_bg_flux = np.trapz(f1[(w1 > b_bg[i][0]) & (w1 < b_bg[i][1])], w1[(w1 > b_bg[i][0]) & (w1 < b_bg[i][1])])*b_scale[i][0]
  #  r_bg_flux = np.trapz(f1[(w1 > r_bg[i][0]) & (w1 < r_bg[i][1])], w1[(w1 > r_bg[i][0]) & (w1 < r_bg[i][1])])*r_scale[i][0]
  #  gauss_e.append((abs(b_bg_flux)+abs(r_bg_flux))/2.)
    
    #print(plot_name[i],lines[i],l_flux,(abs(b_bg_flux)+abs(r_bg_flux))/2. )
    
   # gauss_int.append(np.trapz(gg_fit(w1), w1))
  #  if i == 0:
   #     lines_model = gg_fit
   # else:
    #    lines_model = lines_model + gg_fit
#for a,b, c, d in zip(plot_name,lines,gauss_int,gauss_e):
 #   print(a,'&',b,'&',c*100,'$\pm$',d*100,'\\\\')

#print(gauss_e)
#savdat = Table([plot_name, lines, gauss_int, gauss_e, np.array(line_edges)[:,0], np.array(line_edges)[:,1]], names = ['Species', 'lambda', 'int_flux', 'int_error', 'blue_edge', 'red_edge'])
#ascii.write(savdat, 'int_flux_with_const_table.ecsv', format='ecsv', overwrite=True)
    
plt.show()
        