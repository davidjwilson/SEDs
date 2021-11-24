import numpy as np

"""
Removes negative values by iterativly replacing them with the sum of the nearest two points.

20211111
"""

def remove_negatives(w, f, e):
    """
    Iteratvly removes negative values by combining them with adjecent bins
    """
    wn, fn, en = w, f, e
    
#     w0, w1 = wavelength_edges(w)
    nz = len(fn[fn <0.0])
#     if nz <= 0: 
#         wn, fn, en = w, f, e
#     while nz > 1: 
    minfi = np.argmin(fn) # most negative point
    fi = fn[minfi]
        
    while fi < 0:
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

        wn[minfi] = wi
        fn[minfi] = fi
        en[minfi] = fi
        delinds = np.array(delinds)
        delinds = np.unique(delinds[(delinds != minfi) & (delinds >= 0) & (delinds < len(fn)-1)])
        wn, fn, en  = np.delete(wn, delinds), np.delete(fn, delinds), np.delete(en, delinds)
        minfi = np.argmin(fn) # most negative point
        fi = fn[minfi]
#         plt.plot(wn, fn)
#         plt.scatter(wn[minfi], fi, marker='x', c='C3')
#         plt.show()
        nz = len(fn[fn <0.0])
#         print(nz)
    return(wn[fn >0], fn[fn >0], en[fn >0])






# def fill_gaps():
#     """
#     Fills any gaps larger than 1A
#     """
    