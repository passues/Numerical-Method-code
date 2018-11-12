#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:55:18 2017

@author: Charles
"""
import numpy as np
from svi import svi
import scipy
#import pandas as pd
from optionMetricsTolVols import generateOptionMatricsIvols
from plotIvols import plotIvols
#from sviSqrtFit import sviSqrtFit

def sviFit(ivolData):
    bidVols = ivolData.Bid
    askVols = ivolData.Ask
    vegas   = ivolData.Vega
    expDates = ivolData.Texp.unique()
    nSlices  = len(expDates)
    sviMatrix = []
    
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
        vega   = vegas[pick]
        
        sviGuess = (np.nanmean(midVar), 0.1, 0.1,-0.7, 0)
        def obj(params):
            a = params[0]
            b = params[1]
            sig = params[2]
            rho = params[3]
            m   = params[4]
            sviVar = svi(a,b,sig,rho,m,k)
            minVar = a + b*sig*np.sqrt(np.abs(1-rho**2))
            negVarPenalty = np.min([100, np.exp(-1/minVar)])
            tmp = np.nansum((midVar-sviVar)**2*vega) + negVarPenalty
            
            return tmp*10000
        fit = scipy.optimize.minimize(obj, x0=sviGuess)
        if np.abs(fit.x[3])>0.999:
            fit = scipy.optimize.minimize(obj, x0=sviGuess, bounds=((-10,10),(0,100),(1e-8,100),(-0.999,0.999),(-10,10)), method='L-BFGS-B')
        res = fit.x*np.array([t,t,1,1,1])
        sviMatrix.append(res.tolist())
    
    sviMatrix = pd.DataFrame(sviMatrix, columns=['a','b','sig','rho','m'])
    return sviMatrix

if __name__=='__main__':
#    spxData = pd.read_csv('spxData050915.csv')
#    spxData['strike_price'] = spxData['strike_price']/1000.
#    spxIvols = generateOptionMatricsIvols(spxData)
#    fitSqrt050915 = sviFit(spxIvols)
#    plotIvols(spxIvols, sviMatrix=fitSqrt050915)
    Ivols = pd.read_csv('hsi_20170228.csv')
    is_call = True
    if is_call == False:
        fit_put = ['Strike', 'call_Bid','Bid','call_Ask', 'Ask','Fwd', 'v_c', 'v_p', 'Texp']
        Ivols.columns = fit_put
        Ivols = Ivols.loc[Ivols.Strike<Ivols.Fwd,:]
        Ivols.index = list(range(Ivols.shape[0]))
    else:
        fit_call  = ['Strike', 'Bid','put_Bid','Ask', 'put_Ask','Fwd', 'Vega', 'v_p', 'Texp']
        Ivols.columns = fit_call
        include = np.array(Ivols.Strike>Ivols.Fwd)
        not_include = np.logical_not(include)
        Ivols.Bid[not_include] = Ivols.put_Bid[not_include]
        Ivols.Ask[not_include] = Ivols.put_Ask[not_include]
        
        
      #  Ivols.index = list(range(Ivols.shape[0]))
         
    Ivols.Bid = Ivols.Bid/100.
    Ivols.Ask = Ivols.Ask/100.
    Ivols.Strike = Ivols.Strike
    Ivols.Fwd  = Ivols.Fwd
    fit = sviSqrtFit(Ivols)
    res = plotIvols(Ivols, sviMatrix=fit)
    print('mae', res)
    
