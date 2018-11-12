# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 20:25:20 2017

@author: charles
"""
import numpy as np
import math
import pandas as pd
import scipy.optimize
from svi import svi
from BlackScholes import BSFormula
from sviRoots import sviRoots

def sviFitQR(ivolData, sviGuess, penaltyFactor=100):
    callVals = np.array(ivolData.CallMid)
    expDates = ivolData.Texp.unique()
    nSlices  = len(expDates)
    sviMatrix = sviGuess
    
    for index in range(nSlices-1, -1, -1):
        t = expDates[index]
        texp = np.array(ivolData.Texp)
        
        midVal = np.array(callVals[texp==t])
        pick   = np.logical_not(np.isnan(midVal))
        midVal = midVal[pick]
        f      = (ivolData.Fwd[texp==t])[0]
        k      = np.log(ivolData.Strike[texp==t]/f)[pick]
        
        def sqDist(a, b, sig, rho, m):
            sviVar = svi(a, b, sig, rho, m, k)/t
            outVal = np.array([BSFormula(f, f*math.exp(_k), t, 0, math.sqrt(abs(_sviVar)) ) for _k, _sviVar in zip(k,sviVar)])
            tmp = np.nansum((midVal-outVal)**2)
            return tmp
            
        def sqDistN(a, b, sig, rho, m):
            params = sviGuess.iloc[index,:].to_dict()
            return sqDist(a, b, sig, rho, m)/sqDist(**params)
        
        def crossPenalty(a, b, sig, rho, m):

            cPenalty = 0
            
            if index > 0:
                slicePlusPrevious = pd.concat([pd.DataFrame([dict(a=a, b=b, sig=sig, rho=rho, m=m)]), pd.DataFrame([sviMatrix.iloc[index-1,:].to_dict()])], ignore_index=True)
                cPenalty = sviRoots(slicePlusPrevious)['crossedness']
            else:
                minVar = a + b * sig *np.sqrt(np.abs(1-rho**2))
                negVarPenalty = min(100, math.exp(-1./minVar))
                cPenalty = negVarPenalty
                
            if index < nSlices-1:
                slicePlusNext = pd.concat([pd.DataFrame([dict(a=a, b=b, sig=sig, rho=rho, m=m)]), pd.DataFrame([sviMatrix.iloc[index+1,:].to_dict()])], ignore_index=True)
                cPenalty = cPenalty + sviRoots(slicePlusNext)['crossedness']
                
            return cPenalty*1000
        def obj(params):
            a = params[0]
            b = params[1]
            sig = params[2]
            rho = params[3]
            m   = params[4]
            return sqDistN(a, b, sig, rho, m) + crossPenalty(a, b, sig, rho, m)
            
        _x0 = sviGuess.loc[index,:].tolist()
        fit = scipy.optimize.minimize(obj, x0 =_x0)
        sviMatrix.iloc[index,:] = fit.x
        if abs(fit.x[3])>1:
            fit = scipy.optimize.minimize(obj, x0=sviGuess.iloc[index,:], method='L-BFGS-B', bounds=((-1000,1000),(0,100),(1e-8,100),(-0.999,0.999),(-10,10)))
            sviMatrix.iloc[index,:] = fit.x
    
    sviMatrix = pd.DataFrame(sviMatrix)
    
    return sviMatrix
    

if __name__=='__main__':
    res = sviFitQR(spxIvols, fitSqrt050915)
    