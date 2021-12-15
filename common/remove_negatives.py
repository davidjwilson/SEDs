import numpy as np
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.units import cds
cds.enable()
from scipy.interpolate import interpolate



"""
Removes negative values by iterativly replacing them with the sum of the nearest two points.

20211111
"""
def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff0 = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    diff1 = np.concatenate((diff, np.array([diff[-1]]))) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff0/2.
    w1 = w + diff1/2.
    return w0, w1



def remove_negatives(w, f, e):
    """
    Iteratvly removes negative values by combining them with adjecent bins
    """
    wn, fn, en = w, f, e  
    nz = len(fn[fn <0.0]) 
#     print(nz)
    minfi = np.argmin(fn) # most negative point
    fi = fn[minfi]
    while fi < 0:
        if len(fn) > 2: #can't handle an array less than 2 long
            w0, w1 = wavelength_edges(wn)
            delinds = [] #indices to delete when we're done
            start, end = minfi-1, minfi+2
            if minfi == 0: 
                start = minfi
            elif minfi == len(fn)-1:
                end = minfi +1
            delinds.append(start)
            delinds.append(end-1)
            fi = np.sum(fn[start:end]*(w1[start:end]-w0[start:end])) / (w1[end-1] - w0[start])
            ei = (np.sum(en[start:end]**2 * (w1[start:end]-w0[start:end])**2))**0.5                
            wi = (w0[start]+w1[end-1])/2
            wn[minfi], fn[minfi], en[minfi] = wi, fi, ei
            delinds = np.array(delinds)
            delinds = np.unique(delinds[(delinds != minfi) & (delinds >= 0) & (delinds < len(fn)-1)])
            wn, fn, en  = np.delete(wn, delinds), np.delete(fn, delinds), np.delete(en, delinds)
            minfi = np.argmin(fn) # most negative point
            fi = fn[minfi]
            nz = len(fn[fn <0.0])
#             print('len', len(fn))
#             print(nz)
        else:
            fi = 1e9
    return(wn[fn >0], fn[fn >0], en[fn >0])

def get_line_groups():
    """
    Collections of nearby strong emission lines
    """
    line_groups = np.array([
        [1174.935,1175.265,1175.592,1175.713,1175.713,1175.989,1176.372],
        [1206.499],[1264.737,1265.001],[1238.821], [1242.804],[1294.543],
        [1298.918],[1323.952],[1334.524],[1335.709],[1393.755],[1402.77],
        [1548.201],[1550.772],
        [1640.332,1640.345,1640.375,1640.391,1640.474,1640.49,1640.533],
        [1657.268],
        [1656.267,1656.926,1657.008,1657.379,1657.907,1658.122],
        [1670.787],[2796.35], [2803.53]], dtype='object')
    return line_groups

def make_clean_spectrum(spectrum, dv=0*u.km/u.s):
    """
    Divides a spectrum up into chunks around lines, removes negative values, then sticks them back together.
    Expecting a MM spectrum with 
    'WAVELENGTH',
 'WAVELENGTH0',
 'WAVELENGTH1',
 'FLUX',
 'ERROR',
 'EXPTIME',
 'DQ',
 'EXPSTART',
 'EXPEND'
    """
    w, f, e = spectrum['WAVELENGTH'], spectrum['FLUX'], spectrum['ERROR']
    line_groups = get_line_groups()
    starts = np.array([group[-1]+0.5 for group in line_groups]) #start of each range to remove negatives from
    starts = starts[(starts > w[0]) & (starts < w[-1])]
    starts = dv.to(u.AA, equivalencies=u.doppler_optical(starts*u.AA)).value
    starts = np.insert(starts, 0, w[0])
    ends = np.array([group[0]-0.5 for group in line_groups]) #end of each range to remove negatives from
    end  = dv.to(u.AA, equivalencies=u.doppler_optical(ends*u.AA)).value
    ends = ends[(ends > w[0]) & (ends < w[-1])]
    ends = np.append(ends, w[-1])
    chunks = np.concatenate((starts, ends))
    chunks = chunks[np.argsort(chunks)]
    wnew, fnew, enew = np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
    for i in range(len(chunks)-1):
        start, end = chunks[i], chunks[i+1]
        mask = (w >= start) & (w <= end)
        if len(w[mask]) > 0:
            wi, fi, ei = w[mask], f[mask], e[mask]
            wn, fn, en = remove_negatives(wi, fi, ei)
            wnew= np.concatenate((wnew, wn))
            fnew= np.concatenate((fnew, fn))
            enew= np.concatenate((enew, en))
    args = np.argsort(wnew)
    new_wavelength, new_flux, new_error = wnew[args], fnew[args], enew[args]
    
    new_w0, new_w1 = wavelength_edges(new_wavelength)
    
    new_exptime = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPTIME'])(new_wavelength)
    
    #dq - interploate, then look for unusual values and correct them, summing if the values to either side are different.
   
    new_dq = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['DQ'], kind='previous')(new_wavelength)
    new_dq = new_dq.astype(int)
    
    #expstart - minumum expstart in each bin
    startups = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPSTART'], kind='next')(new_wavelength)
    startdowns = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPSTART'], kind='previous')(new_wavelength)
    new_expstart = np.min([startups, startdowns], axis=0)
    
    #expends - maximum expend in each bin
    endups = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPEND'], kind='next')(new_wavelength)
    enddowns = interpolate.interp1d(spectrum['WAVELENGTH'], spectrum['EXPEND'], kind='previous')(new_wavelength)
    new_expend = np.max([endups, enddowns], axis=0)
    
    names = spectrum.dtype.names
    new_spectrum = Table([new_wavelength*u.AA, new_w0*u.AA, new_w1*u.AA, new_flux*u.erg/u.s/u.cm**2/u.AA, new_error*u.erg/u.s/u.cm**2/u.AA, new_exptime*u.s, 
                           new_dq,new_expstart*cds.MJD, new_expend*cds.MJD], names=names, meta= spectrum.meta)
    return new_spectrum


    