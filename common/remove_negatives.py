import numpy as np

"""
Removes negative values by iterativly replacing them with the sum of the nearest two points.

20211111
"""

def remove_negatives(w0, w1,w, f, e):
    """
    Iteritavly removes negative values by combining them with adjecent bins
    """
    wo, fo, eo = w, f, e
    nz = len(fo[fo <0.0])
    if nz <= 1: 
        wn, fn, en = w, f, e
    while nz > 1: #hangs up on last negative value
        fn = []
        wn = []
        en = []
        inds = []
        last = -2
        for i in range(len(fo)):
            if fo[i] <= 0.0 and i > last+2:
                inds.append(i-1)
                inds.append(i+1)
                start, end = i-1, i+2
                if i == 0: 
                    start = i
                if i == len(fo)-1:
                    end = i+1
                weights = 1 / (eo[start:end]**2)
                fi = np.sum(f[start:end]*(w1[start:end]-w0[start:end])) / (w1[end-1] - w0[start])
                ei = (np.sum(e[start:end]**2 * (w1[start:end]-w0[start:end])**2))**0.5                
                fn.append(fi)
                wn.append((w0[start]+w1[end-1])/2)
                en.append(ei)
                last = i
            else:
                wn.append(wo[i])
                fn.append(fo[i])
                en.append(eo[i])
        inds = np.array(inds)
        inds = np.unique(inds[(inds >= 0) & (inds < len(fo)-1)])
        wn, fn, en = np.array(wn), np.array(fn), np.array(en)
        wn, fn, en  = np.delete(wn, inds), np.delete(fn, inds), np.delete(en, inds)
        nz = len(fn[fn <0.0])
        wo, fo, eo =wn, fn, en
    return(wn[fn >0], fn[fn >0], en[fn >0])

def fill_gaps():
    """
    Fills any gaps larger than 1A
    """
    