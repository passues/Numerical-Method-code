#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:12:32 2017

@author: Charles
"""

import math
import numpy as np

def sviToJw(sviMatrix, texp):
    '''This is the function to convert the svi params to SVI-JW params
    SIV-JW params is more stable with different inputs
    '''
    
    a   = sviMatrix.get('a')
    b   = sviMatrix.get('b')
    sig = sviMatrix.get('sig')
    rho = sviMatrix.get('rho')
    m   = sviMatrix.get('m')
    
    vt   = 1.0*(a+b*(-rho*m+math.sqrt(m**2+sig**2)))/texp
    bhat = math.sqrt(1./(vt*texp))*b
    psit = bhat/2.*(-m/math.sqrt(m**2+sig**2)+rho)
    pt   = bhat*(1-rho)
    ct   = bhat*(1+rho)
    
    varmint = (a+b*abs(sig)*math.sqrt(1-rho**2))/texp*1.
    tmp     = dict(vt=vt, psit=psit, pt=pt, ct=ct, varmint=varmint, texp=texp)
    
    return tmp

def jwToSvi(jwMatrix):
    '''This is the function to convert SVI-JW params to SVI params'''
    
    vt      = jwMatrix.get('vt')
    psit    = jwMatrix.get('psit')
    pt      = jwMatrix.get('pt')
    ct      = jwMatrix.get('ct')
    varmint = jwMatrix.get('varmint')
    texp    = jwMatrix.get('texp')
    
    sqrtw = math.sqrt(vt*texp)
    bhat  = (pt + ct)/2.
    b     = bhat*sqrtw
    rho   = 1. - pt*1./bhat
    bet   = (rho-2.*psit/bhat)
    alpha = np.sign(bet)*math.sqrt(1./bet**2-1)
    
    m   = (vt -varmint)*texp/b/(-rho+np.sign(alpha)*math.sqrt(1+alpha**2) - alpha*math.sqrt(1-rho**2))
    sig = alpha*m
    a   = varmint*texp - b*sig*math.sqrt(1-rho**2)
    
    tmp = dict(a=a, b=b, sig=sig, rho=rho, m=m)
    return tmp


#if __name__=='__main__':
#    svi_param = dict(a = 0.04, b = 0.4, sig = 0.1, rho = - 0.4, m = 0.1)
#    print(jwToSvi(sviToJw(svi_param, 1)))
    