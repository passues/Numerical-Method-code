#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 13:21:52 2017

@author: Charles
"""
import numpy as np
from matplotlib.mlab import stineman_interp
import pandas as pd
import scipy.optimize
from sviFitQuarticRoots import sviFitQR
#from optionMetricsTolVols import generateOptionMatricsIvols
from plotIvols import plotIvols

def sviSqrt(sviSqrtParams, k, w0):
    rho = sviSqrtParams[0]
    eta = sviSqrtParams[1]
    
    w = w0/2.*(1+rho*eta/np.sqrt(w0)*k + np.sqrt((eta/np.sqrt(w0)*k+rho)**2+1.-rho**2))
    return w


def sviSqrtFit(ivolData):
    bidVols = ivolData.Bid
    askVols = ivolData.Ask
    expDates = ivolData.Texp.unique()
    nSlices  = len(expDates)
    
    nrows = ivolData.shape[0]
    midV = np.array([float('nan')]*nrows)
    kk   = np.array([float('nan')]*nrows)
    ww0  = np.array([float('nan')]*nrows)
    
    for index in range(nSlices):
        t       = expDates[index]
        texp    = ivolData.Texp
        pick = texp ==t
        pick &= (bidVols>0)&(askVols>0)
        pick = np.array([not x for x in np.isnan(pick)]) & pick
        pick = pick & (bidVols < askVols)
        pick = np.array(pick)
        midVar = (bidVols[pick]**2+askVols[pick]**2)/2.
        f      = (ivolData.Fwd[texp==t])[0]
        k      = np.log(ivolData.Strike[pick]/f)
        w0     = t*stineman_interp(0,k, midVar)
        
        ww0[pick]  = w0
        midV[pick] = midVar
        kk[pick]   = k
          
    tcutoff = min(0.1, max(expDates))
    
    def obj(sviSqrtParams):
        
        sviSqrtVar = sviSqrt(sviSqrtParams,kk,ww0)/texp
        tmp = (midV - sviSqrtVar)**2
        tmp = np.nansum(tmp[texp >= tcutoff])
        
        return tmp
    
    sviSqrtGuess = (-0.7, 1.)
    
    fit = scipy.optimize.minimize(obj, x0=sviSqrtGuess, bounds=((-0.999,0.999), (-np.Inf, np.Inf)), method='L-BFGS-B')
    res = fit.x
    
    
    sel = np.array([not x for x in np.isnan(ww0)])
    w0r = np.unique(ww0[sel])
    rho = np.array([res[0]]*nSlices)
    a   = w0r/2.*(1-rho**2)
    gg  = res[1]/np.sqrt(w0r)
    b   = w0r/2.*gg
    m   = -rho/gg
    sig = np.sqrt(1-rho**2)/gg
                   
    tmp = pd.DataFrame(dict(a=a, b=b, sig=sig, rho=rho, m=m))
    return tmp 

if __name__=='__main__':
#    spxData = pd.read_csv('spxData050915.csv')
#    spxData['strike_price'] = spxData['strike_price']/1000.
#    spxIvols = generateOptionMatricsIvols(spxData)
    fitSqrt050915 = sviSqrtFit(spxIvols)
  #  fit = sviFitQR(spxIvols, fitSqrt050915)
    plotIvols(spxIvols, sviMatrix=fitSqrt050915)