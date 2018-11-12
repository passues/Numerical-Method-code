# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 18:18:03 2017

@author: charles
"""
import math
from scipy.stats import norm
import numpy as np
import pandas as pd 

def BSFormula(S0, K, T, r, sigma):
    x = math.log(S0*1./K) +r*T
    sig = sigma*math.sqrt(T)
    d1 = x/sig + sig/2.
    d2 = d1 - sig
    pv = math.exp(-r*T)
    return S0*norm.cdf(d1) - pv*K*norm.cdf(d2)
    
def BSFormulaPut(S0, K, T, r, sigma):
    x = math.log(S0*1./K) + r*T
    sig = sigma*math.sqrt(T)
    d1 = x/sig + sig/2.
    d2 = d1 -sig
    pv = math.exp(-r*T)
    return S0*norm.cdf(d1) - pv*K*norm.cdf(d2) + pv*K-S0
    
def BSImpliedVolCall(S0, K, T, r, C):
    '''This function works with vectors of strikes and option values'''
    nK     = len(K)
    sigmaL = np.array([1e-10]*nK)
    CL     = np.array([BSFormula(S0, k, T, r, sig) for k, sig in zip(K, sigmaL)])
    sigmaH = np.array([10.]*nK)
    CH     = np.array([BSFormula(S0, k, T, r, sig) for k, sig in zip(K, sigmaH)])
    
    while np.mean(sigmaH-sigmaL) > 1e-10:
        
        sigma  = (sigmaL+sigmaH)/2
        CM     = np.array([BSFormula(S0, k, T, r, sig) for k, sig in zip(K, sigma)])
        CL     = CL + (CM < C) * (CM - CL)
        sigmaL = sigmaL + (CM < C) * (sigma - sigmaL)
        CH     = CH + (CM >= C) * (CM - CH)
        sigmaH = sigmaH + (CM >= C) * (sigma -sigmaH)
        
    return sigma
    
def BSImpliedVolPut(S0, K, T, r, P):
    pv        = math.exp(-r * T)
    intrinsic = K - pv * S0
    nK        = len(K)
    sigmaL    = np.array([1e-10] * nK)
    PL        = np.array([BSFormula(S0, k, T, r, sig) for k, sig in zip(K,sigmaL)])+intrinsic
                        
    sigmaH    = np.array([10.]*nK)    
    PH        = np.array([BSFormula(S0, k, T, r, sig) for k, sig in zip(K, sigmaH)])+ intrinsic
                
    while np.mean(sigmaH-sigmaL) > 1e-10:
        sigma = (sigmaL + sigmaH)/2
        PM    = np.array([BSFormula(S0, k, T, r, sig) for k, sig in zip(K, sigma)]) + intrinsic
        PL    = PL + (PM < P) * (PM - PL)
        sigmaL = sigmaL + (PM < P)*(sigma -sigmaL)
        PH     = PH + (PM >= P)*(PM-PH)
        sigmaH = sigmaH + (PM>=P) *(sigma-sigmaH)
        
    return sigma

def bsOut(xf, T, AK):
    '''function to compute option prices and implied vols given list of final values of underlying'''
    
    nK = len(AK)
    N  = len(xf)
    xfbar = np.mean(xf)
    CAV = [0.] * nK
    BSV = [0.] * nK
    BSVL = [0.] * nK
    BSVH = [0.] * nK
    
    for j in range(nK):
        #bad python, but just copy Jim
        payoff = (xf - AK[j]) * (xf > AK[j])
        CAV[j] = sum(payoff)/N
        err    = math.sqrt(np.var(payoff)/N)
        BSV[j] = BSImpliedVolCall(xfbar, AK[j], T, 0, CAV[j])
        BSVL[j] = BSImpliedVolCall(xfbar, AK[j], T, 0, CAV[j]-err)
        BSVH[j] = BSImpliedVolCall(xfbar, AK[j], T, 0, CAV[j] +err)
    
    return pd.DataFrame(dict(AK=AK, CAV=CAV, BSV=BSV, BSVL=BSVL, BSVH=BSVH))

def analyticOut(callFormula, AK, T):
    nK = len(AK)
    
    callPrice = [0]*nK
    BSV       = [0]*nK
    for j in range(nK):
        callPrice[j] = callFormula(AK[j])
        BSV[j] = BSImpliedVolCall(1, AK[j], T, 0, callPrice[j])
    
    return pd.DataFrame(dict(AK=AK, calllPrice=callPrice, BSV=BSV))

        

    
    