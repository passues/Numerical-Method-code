#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:03:03 2017

@author: Charles
"""

import pandas as pd
import datetime
import numpy as np
import scipy.optimize
from BlackScholes import BSImpliedVolCall, BSImpliedVolPut, BSFormula
from sviSqrtFit import sviSqrtFit

def getDays(starts, ends):
    start = [datetime.datetime.strptime(str(s), '%Y%m%d') for s in starts]
    expiry = [datetime.datetime.strptime(str(e), '%Y%m%d') for e in ends]
    return [(e - s).days for s,e in zip(start, expiry)]

def generateOptionMatricsIvols(spxData):
    '''assuming spxData is a pandas dataframe'''
    
    spx_date = spxData.get('date')
    spx_exdate = spxData.get('exdate')
    
    spxData['days'] = getDays(spx_date, spx_exdate)
    out2 = pd.DataFrame()
    
    #get the different maturities
    days = spxData.get('days').copy()
    days.sort(inplace=True)
    days = days.unique()
    
    for numdays in days:
        vdc = spxData[(spxData.days==numdays) & (spxData.cp_flag=='C')]
        vdp = spxData[(spxData.days==numdays) & (spxData.cp_flag=='P')]
        expiration =  spxData[spxData.days==numdays].get('exdate').unique()

        #get the put, call strikes and then the unique of the two
        callStrikes = vdc.get('strike_price').unique()
        callStrikes.sort()
        putStrikes = vdp.get('strike_price').unique()
        putStrikes.sort()
        strikes = callStrikes[np.in1d(callStrikes,putStrikes)]

        #setup for calibration
        nK = len(strikes)
        vols = np.array([0.]*nK)
        imid = np.array([0.]*nK)
        ca   = np.array([0.]*nK)
        cb   = np.array([0.]*nK)
        pa   = np.array([0.]*nK)
        pb   = np.array([0.]*nK)

        if nK >=6 and numdays >0:
            cbb = vdc.best_bid
            pbb = vdp.best_bid
            cba = vdc.best_offer
            pba = vdp.best_offer

            for i in range(nK):
                k = strikes[i]
                cb[i] = cbb[vdc.strike_price==k].mean()
                pb[i] = pbb[vdp.strike_price==k].mean()
                ca[i] = cba[vdc.strike_price==k].mean()
                pa[i] = pba[vdp.strike_price==k].mean()

                ##this is so that we can run put-call parity
                ##C-P = PV*(F-K)
                ibid = cb[i] - pb[i]
                iask = ca[i] - pa[i]
                imid [i] = (ibid+iask)/2.

            pvGuess = 1
            fGuess  = (np.array(imid) + np.array(strikes)).mean()
            nearTheMoneyStrikes = [ strikes[index] for index in np.argsort(abs(np.array(imid)))][:6]

            include = np.in1d(strikes, nearTheMoneyStrikes)
            def obj(params):
                f = params[0]
                pv = params[1]
                ifit = pv*(f-strikes)
                ermid = (ifit-imid)*include
                return sum(np.array(ermid)**2)
            
            fit = scipy.optimize.minimize(obj,x0=(fGuess, pvGuess), bounds=((min(strikes), max(strikes)), (0.5,2)), method='L-BFGS-B')
            ffit = fit.x[0]
            pvfit = fit.x[1]
            
            texp = numdays/365.25
            ivolcbid = BSImpliedVolCall(ffit, strikes, texp, 0, cb/pvfit)
            ivolcask = BSImpliedVolCall(ffit, strikes, texp, 0, ca/pvfit)
            ivolpbid = BSImpliedVolPut(ffit, strikes, texp, 0, pb/pvfit)
            ivolpask = BSImpliedVolPut(ffit, strikes, texp, 0, pa/pvfit)
            
            ivolbid = ivolcbid*(strikes>ffit) + ivolpbid*(strikes<=ffit)
            ivolask = ivolcask*(strikes>ffit) + ivolpask*(strikes<=ffit)
            
            callBid = np.array([BSFormula(ffit, strike, texp, 0, ivol) for strike, ivol in zip(strikes, ivolbid)])
            callAsk = np.array([BSFormula(ffit, strike, texp, 0, ivol) for strike, ivol in zip(strikes, ivolask)])
            exclude = (cb==0)|(pb==0)
            callMid = (callBid + callAsk)/2.0
        
            out = pd.DataFrame(dict(Expiry=expiration.tolist()*nK, Texp=[texp]*nK, Strike=strikes, Bid=ivolbid, Ask=ivolask, Fwd=[ffit]*nK, CallMid=callMid))
            out['CallMid'][exclude] = float('nan')
            out['Bid'][exclude] = float('nan')
            out2 = pd.concat([out2,out])
    
    out2['Bid'][out2['Bid']<(1e-8)] = float('nan')
    out3 = out2.iloc[np.argsort(out2['Texp']),:]
    return out3
                      
    
#    
#if __name__=='__main__':
#    spxData = pd.read_csv('spxData050915.csv')
#    spxData['strike_price'] = spxData['strike_price']/1000.
#    spxIvols = generateOptionMatricsIvols(spxData)
#    #fitSqrt110915 = sviSqrtFit(spxIvols)