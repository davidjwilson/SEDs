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