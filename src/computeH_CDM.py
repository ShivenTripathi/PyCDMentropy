import numpy as np
import scipy.special as special
import scipy.stats as stats

def computeH_CDM(nn, ocnts, ncells, opts):
    # if nargin< 2 
    ocnts = ocnts.reshape(-1,1)
    P = np.dot(ocnts,nn).np.sum(nn)/ncells
    Hmax = -(-P*np.log(P) + (1-P)*np.log(1-P))*ncells
    
    return Hbls, Hvar, CIhandle, internal,Hsamples, opts