#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:37:17 2017

@author: Charles
"""
import numpy as np
import matplotlib.pyplot as plt
from svi import svi

def plotIvols(ivolData,sviMatrix=None, slices=None):
    
    bidVols = np.array(ivolData.Bid)
    askVols = np.array(ivolData.Ask)
    expDates = np.array(ivolData.Texp.unique())
    
    nSlices = len(expDates)
    
    if slices:
        nSlices = len(slices)
    else:
        slices = range(nSlices)
    
    column = np.sqrt(nSlices*2)
    rows = round(column/2.,0)
    columns = round(column,0)
    while rows*columns < nSlices:
        rows = rows + 1 
    
    atmVol = np.array([0.]*nSlices)
    atmSkew = np.array([0.]*nSlices)
    
    #plot all the slices
 #   fig, ax = plt.subplots(rows,columns)
    row_idx = 0
    col_idx = 0
    for index in slices:
        t = expDates[index]
        texp = np.array(ivolData.Texp)
        bidVol = np.array(bidVols[texp==t])
        askVol = np.array(askVols[texp==t])
        midVol = (bidVol + askVol)/2.
        f = ivolData.Fwd[texp==t][0]
        k = np.log(np.array(ivolData.Strike[texp==t])/f)
        include = np.logical_not(np.isnan(bidVol))
        kmin = np.min(k[include])
        kmax = np.max(k[include])
        ybottom = 0.8*np.min(bidVol[include])
        ytop = 1.2*np.nanmax(bidVol[include])
        x_range = [kmin, kmax]  
        y_range = [ybottom, ytop]
        plt.figure(index)
        plt.plot(k, bidVol,'x', color='red')
        plt.plot(k, askVol,'x',color='blue')
        
        def vol(k):
            params = sviMatrix.iloc[index,:].to_dict()
            params['k'] = k
            k.sort()
            return np.sqrt(svi(**params)/t)
        fit_vol = vol(k)
        plt.plot(k, fit_vol)
    
    mid_vol = (bidVol + askVol)/2.
    mae = np.sum(np.abs(mid_vol-fit_vol))/len(mid_vol)
    return mae
            
                
            
            
        
if __name__=='__main__':
#    spxData = pd.read_csv('spxData050915.csv')
#    spxData['strike_price'] = spxData['strike_price']/1000.
#    spxIvols = generateOptionMatricsIvols(spxData)
#    fitSqrt050915 = sviFit(spxIvols)
    plotIvols(spxIvols, sviMatrix=fitSqrt050915)
    
    