#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:23:26 2017

@author: Charles
"""

import numpy as np
import matplotlib.pyplot as plt
import functools



def svi(a, b, sig, rho, m, k):
    '''SVI Parametrization'''
    
    return (a+1.*b*(rho*(k-m)+np.sqrt((k-m)*(k-m)+sig*sig)))

    
#if __name__=='__main__':
#    sviparams = dict(a = 0.04, b = 0.4, sig = 0.1, rho = - 0.4, m = 0.1)
#    k = 1
#    callback = functools.partial(svi, sviparams)
#    ks = [x/100. for x in range(-100, 100)]
#    vols = [callback(k) for k in ks]
#    plt.plot(ks, vols)
            
        