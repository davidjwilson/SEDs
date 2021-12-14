#************************************************************************************
# This code calculates habitable zone 'fluxes' using the expression given in the 
# Kopparapu et al.(2014) paper. The corresponding output file is 'HZs.dat'. 
# It also generates a file 'HZ_coefficients.dat' that gives the coefficients for 
# the analytical expression.
#
# Ravi kumar Kopparapu April 19 2014
#
# Translated to python code by John Armstrong (jcarmstrong@weber.edu) 04 June 2014
#************************************************************************************
#************************************************************************************
# Output files.

hzdat = open('HZs.dat', 'w')
hzcoeff = open('HZ_coefficients.dat', 'w')

#************************************************************************************
# Coeffcients to be used in the analytical expression to calculate habitable zone flux 
# boundaries

seff = [0,0,0,0,0,0]
seffsun  = [1.776,1.107, 0.356, 0.320, 1.188, 0.99] 
a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]


#************************************************************************************
# Writing coefficients into 'HZ_coefficients.dat' file

hzcoeff.write('# The coefficients are as follows. The columns, i, are arranged according to\n')
hzcoeff.write('# the HZ limits given in the paper.\n')
hzcoeff.write('#\n')
hzcoeff.write('# i = 1 --> Recent Venus\n')
hzcoeff.write('# i = 2 --> Runaway Greenhouse\n')
hzcoeff.write('# i = 3 --> Maximum Greenhouse\n')
hzcoeff.write('# i = 4 --> Early Mars\n')
hzcoeff.write('# i = 5 --> Runaway Greenhouse for 5 ME\n')
hzcoeff.write('# i = 6 --> Runaway Greenhouse for 0.1 ME\n')

hzcoeff.write('# First row: S_effSun(i)\n')
hzcoeff.write('# Second row: a(i)\n')
hzcoeff.write('# Third row:  b(i)\n')
hzcoeff.write('# Fourth row: c(i)\n')
hzcoeff.write('# Fifth row:  d(i)\n')
hzcoeff.write('   '+'{:6.6E}'.format(seffsun[0]) + '  ' +
              '{:6.6E}'.format(seffsun[1]) + '  ' +
              '{:6.6E}'.format(seffsun[2]) + '  ' +
              '{:6.6E}'.format(seffsun[3]) + '  ' +
              '{:6.6E}'.format(seffsun[4]) + '  ' +
              '{:6.6E}'.format(seffsun[5]) + '  ' +
              '\n')
hzcoeff.write('   '+'{:6.6E}'.format(a[0]) + '  ' +
              '{:6.6E}'.format(a[1]) + '  ' +
              '{:6.6E}'.format(a[2]) + '  ' +
              '{:6.6E}'.format(a[3]) + '  ' +
              '{:6.6E}'.format(a[4]) + '  ' +
              '{:6.6E}'.format(a[5]) + '  ' +
              '\n')
hzcoeff.write('   '+'{:6.6E}'.format(b[0]) + '  ' +
              '{:6.6E}'.format(b[1]) + '  ' +
              '{:6.6E}'.format(b[2]) + '  ' +
              '{:6.6E}'.format(b[3]) + '  ' +
              '{:6.6E}'.format(b[4]) + '  ' +
              '{:6.6E}'.format(b[5]) + '  ' +
              '\n')
hzcoeff.write('   '+'{:6.5e}'.format(c[0]) + '  ' +
              '{:6.5e}'.format(c[1]) + '  ' +
              '{:6.5e}'.format(c[2]) + '  ' +
              '{:6.5e}'.format(c[3]) + '  ' +
              '{:6.5e}'.format(c[4]) + '  ' +
              '{:6.5e}'.format(c[5]) + '  ' +
              '\n')
hzcoeff.write('   '+'{:6.5e}'.format(d[0]) + '  ' +
              '{:6.5e}'.format(d[1]) + '  ' +
              '{:6.5e}'.format(d[2]) + '  ' +
              '{:6.5e}'.format(d[3]) + '  ' +
              '{:6.5e}'.format(d[4]) + '  ' +
              '{:6.5e}'.format(d[5]) + '  ' +
              '\n')

#************************************************************************************
# Calculating HZ fluxes for stars with 2600 K < T_eff < 7200 K. The output file is
# 'HZ_fluxes.dat'

teff  = 2600.0
hzdat.write('#  Teff(K)        Recent        Runaway        Maximum        Early        5ME Runaway   0.1ME Runaway\n')
hzdat.write('#                 Venus         Greenhouse     Greenhouse     Mars         Greenhouse    Greenhouse\n')

starTemp = []
recentVenus = []
runawayGreenhouse = []
maxGreenhouse = []
earlyMars = []
fivemeRunaway = []
tenthmeRunaway = []

while(teff <= 7201.0): 
  tstar = teff - 5780.0
  for i in range(len(a)):
     seff[i] = seffsun[i] + a[i]*tstar + b[i]*tstar**2 + c[i]*tstar**3 + d[i]*tstar**4

  starTemp.append(teff)
  recentVenus.append(seff[0])
  runawayGreenhouse.append(seff[1])
  maxGreenhouse.append(seff[2])
  earlyMars.append(seff[3])
  fivemeRunaway.append(seff[4])
  tenthmeRunaway.append(seff[5])
     
  hzdat.write('   '+'{:6.0f}'.format(teff) + '      ' +
              '{:6.6E}'.format(seff[0]) + '      ' +
              '{:6.6E}'.format(seff[1]) + '     ' +
              '{:6.6E}'.format(seff[2]) + '   ' +
              '{:6.6E}'.format(seff[3]) + ' ' +
              '{:6.6E}'.format(seff[4]) + ' ' +
              '{:6.6E}'.format(seff[5]) + '  ' +
              '\n')
  teff = teff + 200.0

print '************************************************************'
print ''
print 'The HZ coefficients are printed in HZ_coefficients.dat file.'
print 'HZs for stars with 2600 K <= Teff <=7200 K is in HZs.dat file.'
print ''
print '************************************************************'

hzdat.close()
hzcoeff.close()
