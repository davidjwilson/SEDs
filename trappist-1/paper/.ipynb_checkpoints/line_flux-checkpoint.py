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

#name line(A) flux(e-17 erg s-1 cm-2) error (e-17 erg s-1 cm-2) 
#species = np.array(['ciii', 'feii', 'siiii','nv', 'nv', 'ci', 'feii', 'cii', 'cii', 'feii', 'civ', 'civ', 'mgii', 'mgii'] )#
#lines = np.array([1175.5, 1200, 1206.53,1238.821, 1242.804,1288.8, 1292.5,1334.532, 1335.708 ,1527.3 ,1548.202,1550.774 ,2796.35 ,2803.53])
#fluxes = np.array([6.0, 27.7, 4.35,2.31,1.3,-1,3.20,5.37,2.2, 10.8,4.39, 71.0, 33])
#flux_errors = np.array([0.2, 1.1, 0.7, 0.88, 0.1, 0.0, 0.51, 0.80, 0.3, 0.90, 0.67, 6.5,6.0])

#ciii_1175 = [1174.93,1175.26,1175.59,1175.71,1175.99,1176.37]
ciii_1175 = [1174.93,1175.26,1175.59,1175.71,1175.99]
#feii_1200 =[1197.434,1197.499,1198.051,1198.194,1198.368,1198.398,1198.660,1198.924,1199.231,1199.681,1200.889,1201.415,1201.506,1201.548,1201.651,1201.960,1202.452,1202.580]
feii_1200 =[1198.924,1199.231,1199.681,1200.889]
siiii_1206 = [1206.51, 1206.53]
nv_1240 = [1238.821, 1242.804]
#ci_1288 = [1287.609,1287.809,1288.038,1288.040,1288.422,1288.711,1288.918,1289.893,1289.975]
ci_1288 = [1288.711,1288.918]
#feii_1292 = [1290.187,1290.353,1291.168,1291.274,1291.579,1291.824,1291.862,1292.405]
feii_1292 = [1292.405]
cii_1335 = [1334.532, 1335.708]
feii_1527 = [1527.2397] 
#siii_1526 = [1526.72]
civ_1550 = [1548.202,1550.774]
mgii_2800 = [2795.528, 2802.704]


#plotlines = [1175, 1200, 1206.5, 1240, 1290, 1335, 1527, 1549, 2800]
plot_name = np.array(['C\,{\sc iii}', 'Fe\,{\sc ii}', 'Si\,{\sc iii}', 'N\,{\sc v}', 'N\,{\sc v}', 'C\,{\sc i}', 'C\,{\sc ii}', 'C\,{\sc ii}', 'Fe\,{\sc ii}','C\,{\sc iv}', 'C\,{\sc iv}', 'Mg\,{\sc ii}'])  

lines = np.array([1175.5, 1200, 1206.53,1238.821, 1242.804,1288.8,1334.532, 1335.708 ,1527.3 ,1548.202,1550.774 ,2796.35 ])

xlims = [[2.9,2.9], [3.9,3.9],[1.9,0.9],[2.9,3.9],[2.9,3.9],[2.9,2.9],[2.9,2.9],[2.9,2.9],[1.9,1.9],[3.9,3.9],[3.9,3.9],[6.9,6.9]]


#all_lines = [ciii_1175, feii_1200, siiii_1206, nv_1240, ci_1288+feii_1292, cii_1335, feii_1527, civ_1550, mgii_2800]

#making spectrum for plot
w, f, e  = data['WAVELENGTH'], data['FLUX'], data['ERROR']
lya_edges = [1207, 1225]
glowmask = (w <lya_edges[0])|(w > lya_edges[1])
w, f, e = w[glowmask], f[glowmask], e[glowmask]

#note these ones are tweaked for the fit
plotlines = [1175, 1200, 1206.5, 1239, 1243, 1289, 1334, 1336, 1527, 1548, 1550, 2795]


ciii_1175 = [1174.93,1175.26,1175.59,1175.71,1175.99]
feii_1200 =[1199.231,1199.681,1200.889]
siiii_1206 = [1206.51, 1206.53]
nv_1240 = [1238.821, 1242.804]
ci_1288 = [1288.918]
feii_1292 = [1292.405]
cii_1335 = [1334.532, 1335.708]
feii_1527 = [1527.2397] 
civ_1550 = [1548.202,1550.774]
mgii_2800 = [2795.528]


all_lines = [ciii_1175, feii_1200, siiii_1206, [nv_1240[0]], [nv_1240[1]], ci_1288, [cii_1335[0]],[cii_1335[1]], feii_1527, [civ_1550[0]],[civ_1550[1]], mgii_2800]

line_edges = [[1174.6,1176.1],[1198.4,1200.95],[1205.9,1206.8],[1238.2,1239.9],[1242.3,1242.75],[1288.5,1288.8],[1333.9, 1334.55],[1334.9,1335.9],[1526.8, 1527.2],[1547.65, 1548.3],[1550.25,1550.75],[2794.,2797]]
b_bg = [[1173., 1174.],[1197., 1198.],[1205.5, 1205.75],[1237.5, 1238.],[1241.5, 1242.], [1287.5, 1288.],[1332., 1333.],[1332., 1333.],[1525., 1526.],[1546.,1547.],[1549.,1550],[2791,2793.]]

r_bg = [[1176.25, 1176.5],[1201., 1203.],[1204., 1205.],[1239.5, 1240.],[1243., 1244.],[1289., 1290.],[1334.6,1334.8],[1334.6,1334.8],[1527.5,1528.5],[1548.5,1549.],[1551,1552],[2798., 2799.]]

line_width = np.diff(line_edges)
b_width = np.diff(b_bg)
b_scale = line_width/b_width

r_width = np.diff(r_bg)
r_scale = line_width/b_width


#all_lines[4] = ci_1288


#all_lines[-1] =[2795.528]

dv = -56300 #Bourrier+17a 
c = const.c.value
dshift = (1.0+(dv/c))
#print(dshift)



gauss_int  = []
gauss_e = []

fitter = fitting.SLSQPLSQFitter()


plt.figure(figsize=(18,15))
plt.subplots_adjust(hspace=0.1, wspace=0.1,top=0.99, right=0.99, left = 0.05, bottom =0.05)
for i in range(len(plotlines)):
    pl = plotlines[i]
    mask = (w > pl-7) & (w < pl+7)
    w1, f1, e1  = w[mask], f[mask], e[mask]
    plt.subplot(3,4,i+1)
    if i < 11:
        f1 = convolve(f1,Box1DKernel(5))
        e1 = convolve(e1,Box1DKernel(5))
    f1 *= 1e16
    e1 *= 1e16
    plt.step(w1, f1)
    mask2 = (w1 > pl-xlims[i][0]) & (w1 < pl+xlims[i][1])
    top = max(f1[mask2])
    bottom = min(f1[mask2])
    nmask = (lines > pl-5) & (lines < pl+5)
    plt.ylim(bottom*1.1, top*1.6)
    plt.xlim(pl - xlims[i][0], pl + xlims[i][1])
    #if i ==3:
     #   plt.ylabel('Flux ($10^{-16}$ erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)')
    #if i == 7:
     #   plt.xlabel('Wavelength (\AA)')
    [plt.annotate(name,(line*dshift, top*1.05), xytext=(line*dshift, top*1.3), horizontalalignment='center') for name, line in zip(plot_name[nmask], lines[nmask])]
    [plt.annotate('',(line*dshift, top*1.05), xytext=(line*dshift, top*1.25),arrowprops=dict(arrowstyle='-'), horizontalalignment='center') for line in all_lines[i]]
    plt.axhline(0, c='k', ls='--')
    guess_sigma = 0.05
    if i == 11:
        guess_sigma = 0.1
    gg_init = models.Gaussian1D(1.0, all_lines[i][0]*dshift, guess_sigma)
    for line in all_lines[i][1:]:
        gg_init = gg_init + models.Gaussian1D(1.0, line*dshift, guess_sigma)
    gg_fit = fitter(gg_init, w1, f1)
    plt.plot(w1, gg_fit(w1))
    
    l_flux = np.trapz(gg_fit(w1)[(w1 > line_edges[i][0]) & (w1 < line_edges[i][1])], w1[(w1 > line_edges[i][0]) & (w1 < line_edges[i][1])])
    gauss_int.append(l_flux)
    b_bg_flux = np.trapz(f1[(w1 > b_bg[i][0]) & (w1 < b_bg[i][1])], w1[(w1 > b_bg[i][0]) & (w1 < b_bg[i][1])])*b_scale[i]
    r_bg_flux = np.trapz(f1[(w1 > r_bg[i][0]) & (w1 < r_bg[i][1])], w1[(w1 > r_bg[i][0]) & (w1 < r_bg[i][1])])*r_scale[i]
    gauss_e.append((abs(b_bg_flux)+abs(r_bg_flux))/2.)
    
    #print(plot_name[i],lines[i],l_flux,(abs(b_bg_flux)+abs(r_bg_flux))/2. )
    
   # gauss_int.append(np.trapz(gg_fit(w1), w1))
  #  if i == 0:
   #     lines_model = gg_fit
   # else:
    #    lines_model = lines_model + gg_fit
for a,b, c, d in zip(plot_name,lines,gauss_int,gauss_e):
    print(a,b,c,d)

plt.show()
        