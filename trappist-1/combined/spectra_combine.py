import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from scipy.io.idl import readsav
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import instruments as inst #pl's istrument code, might just make a function here
from astropy.time import Time
from linsky_euv import euv_estimator
from astropy.units import cds
from craftroom import resample
from scipy.interpolate import interp1d
cds.enable()

"""
v3 20190724 trappist 1 only, depreciated
Script to make Mega-Muscles spectra. Will develop over time.

objectives for this version
output: a file with everything that MUSCLES did i.e.
'WAVELENGTH' 
'WAVELENGTH0' = wavelength bin edges
'WAVELENGTH1'
'FLUX'
'ERROR'
'EXPTIME' exptime of that bit of spectrum
'DQ',
'EXPSTART' exposure dates, mjd
'EXPEND',
'INSTRUMENT', bitwise- could we do this with words or is it what hlsp requires?
'NORMFAC',
'BOLOFLUX' =flux normalised to bolometric flux
'BOLOERR'




"""


#totals = [w_full, f_full, e_full, dq_full, n_full,i_full]
#additions =[w, f, e, dq, n, instrument]

def add_spec(totals, additions):
    """running combinations of wavelength, flux, error, data quality flags, normalisation and instrument"""
    totals =  np.concatenate((totals, additions), axis=1)
    return totals

def sort_totals(totals):
    """sorts totals array by wavelength"""
    arr1inds = totals[0].argsort()
    totals = [t[arr1inds] for t in totals]
    #for t in totals:
     #   t = t[arr1inds]
    return totals
    
def totals_setup():
    """set up arrays to keep running totals of data"""
    names = ['WAVELENGTH','WAVELENGTH0', 'WAVELENGTH1','FLUX','ERROR','EXPTIME',
             'DQ','EXPSTART','EXPEND','INSTRUMENT','NORMFAC']#,'BOLOFLUX' ,'BOLOERR'] #bolos are done when spectrum is compleate
    totals = []
    i = 0
    while i < len(names):
        totals.append([])
        i += 1
    totals = np.array(totals, dtype = float)
    return names, totals

def make_section(totals, data, mask =[], normfac=1.0, clip=[], plot=True):
    """
    runs all the functions needed to add a section to the spectrum
    """
    w, f, additions = build_additions(data, mask=mask, normfac=normfac, clip=clip)
    if plot==True:
        plot_spectrum(w,f)
    totals = add_spec(totals, additions)
    w_start, w_end = w[0],w[-1]
    return w_start, w_end, totals

def build_additions(col, mask=[],  normfac=1.0, clip=[]):
    """
    makes the "addition" array if a mask is applied. "col" is one of the collection dicts generated by the read functions
    """
    if len(mask) > 0:
        if len(clip) > 1: #clip if there are options to clip
            w0, w1, w, f, e, dq, exptime = col['w0'][mask][clip[0]:clip[1]], col['w1'][mask][clip[0]:clip[1]],col['w'][mask][clip[0]:clip[1]],col['f'][mask][clip[0]:clip[1]],col['e'][mask][clip[0]:clip[1]],col['dq'][mask][clip[0]:clip[1]],col['exptime'][mask][clip[0]:clip[1]]
        else:
            w0, w1, w, f, e, dq, exptime = col['w0'][mask], col['w1'][mask],col['w'][mask],col['f'][mask],col['e'][mask],col['dq'][mask],col['exptime'][mask]
    else:
        if len(clip) > 1: #clip if there are options to clip
            w0, w1, w, f, e, dq, exptime = col['w0'][clip[0]:clip[1]], col['w1'][clip[0]:clip[1]],col['w'][clip[0]:clip[1]],col['f'][clip[0]:clip[1]],col['e'][clip[0]:clip[1]],col['dq'][clip[0]:clip[1]],col['exptime'][clip[0]:clip[1]]
        else:
            w0, w1, w, f, e, dq, exptime = col['w0'], col['w1'],col['w'],col['f'],col['e'],col['dq'],col['exptime']        
    f, e = f*normfac, e*normfac #spectrum scaling
    expstart = np.full(len(w),col['expstart'], dtype=float)
    expend = np.full(len(w),col['expend'],dtype=float)
    instrument = np.full(len(w),col['instrument'],dtype=int)
    normfac = np.full(len(w),normfac,dtype=float)
    additions = [w, w0, w1, f, e, exptime,dq, expstart, expend, instrument, normfac]
    return w, f, additions
    
    
def plot_spectrum(w, f):
    """plots a spectrum section"""
    plt.step(w,f)


def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff/2.
    w1 = w + diff/2.
    return w0, w1

def dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code): 
    """
    puts all the bits for each instrument into a dictionary. 
    w0, w1, w, f, e, dq, exptime are arrays of the same length
    expstart, expend are floats
    instrument_code is an integer
    """
    collection = {'w0':w0, 'w1':w1, 'w':w, 'f':f, 'e':e, 'dq':dq, 'exptime':exptime, 'expstart':expstart, 'expend':expend, 'instrument': instrument_code}
    return collection

def read_idl(filepath):
    """
    Extracts data from KF's IDL files
    """
    data = readsav(filepath)
    w, f, e, exptime = data['wave'], data['flux'], data['err'], data['exptime']
    dq = np.zeros(len(w)) #consult with kf on this
    w0, w1 = wavelength_edges(w)
    instrument_list = np.unique(data['grating'])
    instruments = np.array([('hst_cos_'+str(i)[2:-1].lower()) for i in instrument_list])
    instrument_code = np.sum([inst.getinsti(i) for i in instrument_list])
    expstart, expend = 0, 0 #not in idl, placeholder.
    idl_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    return idl_collection

def read_stis_x1d(filepath):
    """
    Extracts data from a hst x1d file
    """
    hdul = fits.open(filepath)
    data = hdul[1].data[0]
    w, f, e, dq = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
    w0, w1 = wavelength_edges(w)
    hdr = hdul[0].header
    exptime = np.full(len(w), hdr['TEXPTIME'])
    expstart, expend = hdr['TEXPSTRT'], hdr['TEXPEND']
    instrument_name = hdr['TELESCOP']+'_sts_'+hdr['OPT_ELEM']  #nb STIS=sts in Pl code 
    instrument_code = inst.getinsti(instrument_name.lower())
    stis_x1d_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    hdul.close()
    return stis_x1d_collection

def read_ecsv(filepath):
    """
    reads an escv with a coadded spectrum. Not much in these yet, need to improve.
    """
    data = Table.read(filepath)#
    w, f, e, dq = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
    w0, w1 = wavelength_edges(w)
    #needs everything else!
    ecsv_collection =  dict_builder(w0, w1, w, f, e, dq, np.zeros(len(w)), 0.,0, 0)
    return ecsv_collection
    
def read_xmm(filepath):
    """
    collects data from CS's xmm outputs and for the apec models
    """
    hdul = fits.open(filepath)
    #spectrum
    data = hdul[1].data
    w, f, e = data['Wave'], data['CFlux'], data['CFLux_Err']
    w0, w1 = w - (data['bin_width']/2), w+(data['bin_width']/2)
    hdr = hdul[0].header
    exptime = np.full(len(w), hdr['HIERARCH pn_DURATION'])
    expstart = Time(hdr['HIERARCH pn_DATE-OBS']).mjd
    expend = Time(hdr['HIERARCH pn_DATE-END']).mjd
    instrument_list = ['xmm_epc_multi','xmm_epc_pn---']
    instrument_code = np.sum([inst.getinsti(i) for i in instrument_list])
    xmm_collection = dict_builder(w0, w1, w, f, e, np.zeros(len(w)),exptime, expstart, expend, instrument_code)
    #model
    data = hdul[2].data
    w, f = data['Wave'], data['Flux']
    w0, w1 = w - (data['bin_width']/2), w+(data['bin_width']/2)
   # hdr = hdul[0].header
    instrument_code = inst.getinsti('mod_apc_-----')
    apec_collection = dict_builder(w0, w1, w, f, np.zeros(len(w)), np.zeros(len(w)),np.zeros(len(w)), 0., 0., instrument_code)
    hdul.close()
    return xmm_collection, apec_collection

def read_lyamod(filepath):
    """
    gather information from a lya model
    """
    lya = Table.read(filepath, format='ascii')
    w, f = lya['WAVELENGTH'], lya['FLUX']
    w0, w1 = wavelength_edges(w)
    instrument_code = inst.getinsti('mod_lya_young')
    lya_collection = dict_builder(w0,w1,w,f, np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w)), 0., 0., instrument_code)
    return lya_collection

def read_scaled_phoenix(filepath):
    """
    gather information from a pre-scaled phoenix model
    """
    data = Table.read(filepath)
    w, f = data['WAVELENGTH'], data['FLUX']
    w0, w1 = wavelength_edges(w)
    instrument_code = inst.getinsti('mod_phx_-----')
    phx_collection = dict_builder(w0,w1,w,f, np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w)), 0., 0., instrument_code)
    return phx_collection

def read_ecsv_phoenix(filepath):
    """
    gather information from a phoenix model saved as ecsv
    """
    data = Table.read(filepath)
    w, f = data['WAVELENGTH'], data['FLUX']
    w0, w1 = wavelength_edges(w)
    instrument_code = inst.getinsti('mod_phx_-----')
    phx_collection = dict_builder(w0,w1,w,f, np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w)), 0., 0., instrument_code)
    return phx_collection
    
def read_phoenix(filepath, wavepath):
    """
    gather information from a phoenix model. Note that flux and wavelegth are in different files
    """
    w = fits.getdata(wavepath, 0)
    f = fits.getdata(filepath)
    w0, w1 = wavelength_edges(w)
    instrument_code = inst.getinsti('mod_phx_-----')
    phx_collection = dict_builder(w0,w1,w,f, np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w)), 0., 0., instrument_code)
    return phx_collection
    
def read_stis_ccd(filepath):
    """
    collects data for a STIS ccd observation
    """
    data = fits.getdata(filepath)[0]
    hdr = fits.getheader(filepath,0)
    w, f, e, dq = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
    w0, w1 = wavelength_edges(w)
    exptime = np.full(len(w), hdr['TEXPTIME'])
    expstart, expend = hdr['TEXPSTRT'], hdr['TEXPEND']
    instrument_name = hdr['TELESCOP']+'_sts_'+hdr['OPT_ELEM']  #nb STIS=sts in Pl code 
    instrument_code = inst.getinsti(instrument_name.lower())
    stis_ccd_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    return stis_ccd_collection

def read_cos_nuv(filepath, fillgap=True):
    """
    collects data from a COS nuv observation. If fill gap =True, fits a polynomial to the gap and trims the zero-flux inner ends off the spectra
    """
    data = fits.getdata(filepath,1)[0:2] #not using reflected order
    hdr0 = fits.getheader(filepath,0)
    hdr1 = fits.getheader(filepath,1)
    w = np.array([], dtype=float)
    f = np.array([], dtype=float)
    e = np.array([], dtype=float)
    dq = np.array([], dtype=int)
    if fillgap == True:
        gap_w, gap_f = nuv_fill(data)
        for dt in data:
            mask = (dt['WAVELENGTH'] < gap_w[0]) | (dt['WAVELENGTH'] > gap_w[-1])
            w= np.concatenate((w, dt['WAVELENGTH'][mask]))
            f = np.concatenate((f, dt['FLUX'][mask]))
            e = np.concatenate((e, dt['ERROR'][mask]))
            dq = np.concatenate((dq, dt['DQ'][mask]))
        gap_w0, gap_w1 = wavelength_edges(gap_w)
        gap_code = inst.getinsti('oth_---_other')
        gap_collection = dict_builder(gap_w0, gap_w1, gap_w, gap_f, np.zeros(len(gap_w)), np.zeros(len(gap_w)),np.zeros(len(gap_w)), 0., 0., gap_code)        
    else:
        for dt in data:
            w= np.concatenate((w, dt['WAVELENGTH']))
            f = np.concatenate((f, dt['FLUX']))
            e = np.concatenate((e, dt['ERROR']))
            dq = np.concatenate((dq, dt['DQ']))
        gap_collection = {}
    w0, w1 = wavelength_edges(w)
    exptime = np.full(len(w), hdr1['EXPTIME'])
    expstart, expend = hdr1['EXPSTART'], hdr1['EXPEND']
    instrument_name = hdr0['TELESCOP']+'_cos_'+hdr0['OPT_ELEM']   
    instrument_code = inst.getinsti(instrument_name.lower())
    cos_nuv_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    return cos_nuv_collection, gap_collection
    
def nuv_fill(data):
    """
    fill in the gap in a cos nuv specturm with a polynomial
    """
    w1, f1 = data[0]['WAVELENGTH'][data[0]['DQ']==0], data[0]['FLUX'][data[0]['DQ']==0]
    w2, f2 = data[1]['WAVELENGTH'][data[1]['DQ']==0], data[1]['FLUX'][data[1]['DQ']==0]
    mgii_mask = (w2 < 2790)|(w2 > 2805) #remove mgii line
    end_w, end_f = np.concatenate((w1,w2[mgii_mask])), np.concatenate((f1,f2[mgii_mask]))
    gap_w = np.arange(w1[-1],w2[0], 1.0)
    gap_f = np.polyval((np.polyfit(end_w,end_f,2)), gap_w)
    return gap_w, gap_f

def make_euv(lya, distance):
    """
    collects the bits for an euv estimation
    """
    w, f = euv_estimator(lya, distance)
    w0, w1 = wavelength_edges(w)
    e, dq, exptime = np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w))
    expstart, expend = 0., 0.
    instrument_code = inst.getinsti('mod_euv_young')
    euv_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    return euv_collection

def dem_to_1A(w,f):
    """
    Converts a DEM model at 5A resolution to 1A resolution
    """
    w1 = np.arange(w[0], w[-1], 1.)
    f1 = interp1d(w, f, fill_value='extrapolate', kind='nearest')(w1)
    return w1, f1


def read_dem(filepath):
    """
    adds a dem model to the EUV
    """
    data = fits.getdata(filepath,1)
    w, bin_f, bin_el, bin_eu = data['Wavelength'], data['Bin-Integrated Flux'], data['Lower_Error_16'], data['Upper_Error_84']
    w0, w1 = wavelength_edges(w)
    f = bin_f/(w1-w0) #convert from bin-intergrated flux to flux 
    w, f = dem_to_1A(w,f)
    w0, w1 = wavelength_edges(w) #remake bin edges
    e, dq, exptime = np.zeros(len(w)), np.zeros(len(w)), np.zeros(len(w))
    expstart, expend = 0., 0.
    instrument_code = inst.getinsti('oth_---_other')
    dem_collection = dict_builder(w0, w1, w, f, e, dq, exptime, expstart, expend, instrument_code)
    return dem_collection
                                                                
def save_to_ecsv(totals, names, star, version, save_path = '', save_1A = False):
    """
    saves the completed spectrum to an ecsv file. No bolflux for now.
    names = ['WAVELENGTH','WAVELENGTH0', 'WAVELENGTH1','FLUX','ERROR','EXPTIME',
             'DQ','EXPSTART','EXPEND','INSTRUMENT','NORMFAC']#,'BOLOFLUX' ,'BOLOERR']
    """
    data = Table([totals[0]*u.AA, totals[1]*u.AA, totals[2]*u.AA, 
                  totals[3]*u.erg/u.s/u.AA/u.cm**2, totals[4]*u.erg/u.s/u.AA/u.cm**2,
                  totals[5]*u.s, totals[6], totals[7]*cds.MJD, totals[8]*cds.MJD,
                  totals[9], totals[10]],
                names=names)
    ascii.write(data, star+'_sed_var_res_'+version+'.ecsv', format = 'ecsv', overwrite=True)
    if save_1A == True:
        save_1A(totals, names, star, version, save_path = '')
        
def save_basic(totals, names, star, version, save_path = '', save_1A = False):
    """
    save just wavelength, flux, error to an ascii file.
    """
    data = Table([totals[0]*u.AA, 
                  totals[3]*u.erg/u.s/u.AA/u.cm**2, totals[4]*u.erg/u.s/u.AA/u.cm**2],
                names=[names[0], names[3], names[4]])
    ascii.write(data, star+'_basic_sed_var_res_'+version+'.ecsv', format = 'ecsv', overwrite=True)
   
    

def save_1A(totals, names, star, version, save_path = ''):
    """
    Saves a spectrum binned to 1A using ZBT's craftroom.resample module
    """
    w, f = resample.bintogrid(totals[0],totals[3],  dx=1.0)
    w, e = resample.bintogrid(totals[0],totals[4],  dx=1.0)
    w0, w1 = wavelength_edges(w)
    #need to dig around in Parke's scripts for the rest of it

"""
Information I need to provide for each spectrum:
-where to get the data
-how to extract it
-what range to plot
-what scaling to use

Information I need to extract:
-wavelength
-wavelength edges (for some)
-flux
-error
-dq
-exptime
-expstart
-expend
-instrument

Information I need to caculate:
-wavelength edges (for some)
-boloflux
-boloerror

"""

