import astropy.io.fits as fits
from astropy.table import Table, vstack
from astropy.io import ascii

"""
Saves SED tables as MUSCLES fits files
"""

def sed_to_ecsv(sed_table):
    """
    quickly saving the whole sed_table to ecsv for testing
    
    """
    star = sed_table.meta['TARGNAME']
    ascii.write(sed_table, 'quicksaves/'+star+'_sed_test.ecsv', format='ecsv', overwrite=True)
    
    
#def make_primary(sed_table, instrument_list)









"""
SIMPLE  =                    T / conforms to FITS standard                      
BITPIX  =                    8 / array data type                                
NAXIS   =                    0 / number of array dimensions                     
EXTEND  =                    T                                                  
TELESCOP= 'MULTI   '                                                            
INSTRUME= 'MULTI   '                                                            
GRATING = 'MULTI   '                                                            
TELESC00= 'MODEL   '                                                            
INSTRU00= 'PHX     '                                                            
GRATIN00= 'NA      '                                                            
TELESC01= 'MODEL   '                                                            
INSTRU01= 'EUV-SCALING'                                                         
GRATIN01= 'NA      '                                                            
TELESC02= 'MODEL   '                                                            
INSTRU02= 'APEC    '                                                            
GRATIN02= 'NA      '                                                            
TELESC03= 'HST     '                                                            
INSTRU03= 'STIS    '                                                            
GRATIN03= 'G230L   '                                                            
TELESC04= 'HST     '                                                            
INSTRU04= 'COS     '                                                            
GRATIN04= 'G160M   '                                                            
TELESC05= 'HST     '                                                            
INSTRU05= 'COS     '                                                            
GRATIN05= 'G130M   '                                                            
TELESC06= 'HST     '                                                            
INSTRU06= 'STIS    '                                                            
GRATIN06= 'G430L   '                                                            
TELESC07= 'HST     '                                                            
INSTRU07= 'COS     '                                                            
GRATIN07= 'G230L   '                                                            
TELESC08= 'MODEL   '                                                            
INSTRU08= 'LYA-RECONSTRUCTION'                                                  
GRATIN08= 'NA      '                                                            
TARGNAME= 'GJ176   '                                                            
RA_TARG =    70.73239599999999                                                  
DEC_TARG=            18.958163                                                  
PROPOSID=                13650                                                  
HLSPNAME= 'Measurements of the Ultraviolet Spectral Characteristics of &'       
CONTINUE  'Low-mass Exoplanet Host Stars'                                       
HLSPACRN= 'MUSCLES '                                                            
HLSPLEAD= 'R. O. Parke Loyd'                                                    
PR_INV_L= 'France  '                                                            
PR_INV_F= 'Kevin   '                                                            
WAVEMIN =    5.000441551208496                                                  
WAVEMAX =             54999.75                                                  
WAVEUNIT= 'ang     '                                                            
AIRORVAC= 'vac     '                                                            
FLUXMIN = -7.1025795946309E-15                                                  
FLUXMAX = 1.10545481750486E-12                                                  
FLUXUNIT= 'erg/s/cm2/ang'                                                       
BOLOFLUX= 1.25681772600649E-08                                                  
LNZ_NORM= 3.96067166323457E-14                                                  
LNZ_GAM =   0.5015959287107679 

"""