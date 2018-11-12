# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 21:16:39 2017

@author: charles
"""
import numpy as np
import cmath
from svi import svi

def sviRoots(sviData):

    a1 = sviData.a[0] 
    b1 = sviData.b[0] 
    s1 = sviData.sig[0] 
    r1 = sviData.rho[0] 
    m1 = sviData.m[0] 
    
    a2 = sviData.a[1] 
    b2 = sviData.b[1] 
    r2 = sviData.rho[1] 
    m2 = sviData.m[1] 
    s2 = sviData.sig[1] 

# The standard form of the quartic is q4 x**4 + q3 x**3 +q2 x**2 +q1 x + q0 == 0
# Note that multiplying all the coefficients qi by a factor should have no effect on the roots.

    q2 = 1000000 * -2 * (-3 * b1 ** 4 * m1 ** 2 + b1 ** 2 * b2 ** 2 * m1 ** 2 + 4 * b1 ** 2 * b2 ** 2 * m1 * m2 + 
                b1 ** 2 * b2 ** 2 * m2 ** 2 - 3 * b2 ** 4 * m2 ** 2 + 6 * b1 ** 4 * m1 ** 2 * r1 ** 2 + 
                b1 ** 2 * b2 ** 2 * m1 ** 2 * r1 ** 2 + 4 * b1 ** 2 * b2 ** 2 * m1 * m2 * r1 ** 2 + 
                b1 ** 2 * b2 ** 2 * m2 ** 2 * r1 ** 2 - 3 * b1 ** 4 * m1 ** 2 * r1 ** 4 - 6 * b1 ** 3 * b2 * m1 ** 2 * r1 * r2 - 
                6 * b1 ** 3 * b2 * m1 * m2 * r1 * r2 - 6 * b1 * b2 ** 3 * m1 * m2 * r1 * r2 - 
                6 * b1 * b2 ** 3 * m2 ** 2 * r1 * r2 + 6 * b1 ** 3 * b2 * m1 ** 2 * r1 ** 3 * r2 + 
                6 * b1 ** 3 * b2 * m1 * m2 * r1 ** 3 * r2 + b1 ** 2 * b2 ** 2 * m1 ** 2 * r2 ** 2 + 
                4 * b1 ** 2 * b2 ** 2 * m1 * m2 * r2 ** 2 + b1 ** 2 * b2 ** 2 * m2 ** 2 * r2 ** 2 + 6 * b2 ** 4 * m2 ** 2 * r2 ** 2 - 
                3 * b1 ** 2 * b2 ** 2 * m1 ** 2 * r1 ** 2 * r2 ** 2 - 12 * b1 ** 2 * b2 ** 2 * m1 * m2 * r1 ** 2 * r2 ** 2 - 
                3 * b1 ** 2 * b2 ** 2 * m2 ** 2 * r1 ** 2 * r2 ** 2 + 6 * b1 * b2 ** 3 * m1 * m2 * r1 * r2 ** 3 + 
                6 * b1 * b2 ** 3 * m2 ** 2 * r1 * r2 ** 3 - 3 * b2 ** 4 * m2 ** 2 * r2 ** 4 - 
                a1 ** 2 * (b1 ** 2 * (-1 + 3 * r1 ** 2) - 6 * b1 * b2 * r1 * r2 + b2 ** 2 * (-1 + 3 * r2 ** 2)) - 
                a2 ** 2 * (b1 ** 2 * (-1 + 3 * r1 ** 2) - 6 * b1 * b2 * r1 * r2 + b2 ** 2 * (-1 + 3 * r2 ** 2)) - 
                2 * a2 * (3 * b1 ** 3 * m1 * r1 * (-1 + r1 ** 2) - b1 ** 2 * b2 * (2 * m1 + m2) * (-1 + 
                    3 * r1 ** 2) * r2 - 3 * b2 ** 3 * m2 * r2 * (-1 + r2 ** 2) + b1 * b2 ** 2 * (m1 + 2 * m2) * 
                   r1 * (-1 + 3 * r2 ** 2)) + 2 * a1 * (3 * b1 ** 3 * m1 * r1 * (-1 + r1 ** 2) - 
                  b1 ** 2 * b2 * (2 * m1 + m2) * (-1 + 3 * r1 ** 2) * r2 - 3 * b2 ** 3 * m2 * r2 * (-1 + 
                    r2 ** 2) + b1 * b2 ** 2 * (m1 + 2 * m2) * r1 * (-1 + 3 * r2 ** 2) + 
                  a2 * (b1 ** 2 * (-1 + 3 * r1 ** 2) - 6 * b1 * b2 * r1 * r2 + b2 ** 2 * (-1 + 3 * r2 ** 2))) - 
                b1 ** 4 * s1 ** 2 + b1 ** 2 * b2 ** 2 * s1 ** 2 + b1 ** 4 * r1 ** 2 * s1 ** 2 - 2 * b1 ** 3 * b2 * r1 * r2 * 
                 s1 ** 2 + b1 ** 2 * b2 ** 2 * r2 ** 2 * s1 ** 2 + b1 ** 2 * b2 ** 2 * s2 ** 2 - b2 ** 4 * s2 ** 2 + 
                b1 ** 2 * b2 ** 2 * r1 ** 2 * s2 ** 2 - 2 * b1 * b2 ** 3 * r1 * r2 * s2 ** 2 + b2 ** 4 * r2 ** 2 * s2 ** 2)
    
    q4 = 1000000 * (b1 ** 4 * (-1 + r1 ** 2) ** 2 - 4 * b1 ** 3 * b2 * r1 * (-1 + r1 ** 2) * r2 - 
                4 * b1 * b2 ** 3 * r1 * r2 * (-1 + r2 ** 2) + b2 ** 4 * (-1 + r2 ** 2) ** 2 + 
                2 * b1 ** 2 * b2 ** 2 * (-1 - r2 ** 2 + r1 ** 2 * (-1 + 3 * r2 ** 2)))
    
    q0 = 1000000 * (a1 ** 4 + a2 ** 4 + b1 ** 4 * m1 ** 4 - 2 * b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 ** 2 + b2 ** 4 * m2 ** 4 - 
                2 * b1 ** 4 * m1 ** 4 * r1 ** 2 - 2 * b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 ** 2 * r1 ** 2 + b1 ** 4 * m1 ** 4 * r1 ** 4 + 
                4 * b1 ** 3 * b2 * m1 ** 3 * m2 * r1 * r2 + 4 * b1 * b2 ** 3 * m1 * m2 ** 3 * r1 * r2 - 
                4 * b1 ** 3 * b2 * m1 ** 3 * m2 * r1 ** 3 * r2 - 2 * b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 ** 2 * r2 ** 2 - 
                2 * b2 ** 4 * m2 ** 4 * r2 ** 2 + 6 * b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 ** 2 * r1 ** 2 * r2 ** 2 - 
                4 * b1 * b2 ** 3 * m1 * m2 ** 3 * r1 * r2 ** 3 + b2 ** 4 * m2 ** 4 * r2 ** 4 + 
                4 * a2 ** 3 * (b1 * m1 * r1 - b2 * m2 * r2) - 4 * a1 ** 3 * (a2 + b1 * m1 * r1 - 
                  b2 * m2 * r2) + 2 * b1 ** 4 * m1 ** 2 * s1 ** 2 - 2 * b1 ** 2 * b2 ** 2 * m2 ** 2 * s1 ** 2 - 
                2 * b1 ** 4 * m1 ** 2 * r1 ** 2 * s1 ** 2 + 4 * b1 ** 3 * b2 * m1 * m2 * r1 * r2 * s1 ** 2 - 
                2 * b1 ** 2 * b2 ** 2 * m2 ** 2 * r2 ** 2 * s1 ** 2 + b1 ** 4 * s1 ** 4 - 2 * b1 ** 2 * b2 ** 2 * m1 ** 2 * s2 ** 2 + 
                2 * b2 ** 4 * m2 ** 2 * s2 ** 2 - 2 * b1 ** 2 * b2 ** 2 * m1 ** 2 * r1 ** 2 * s2 ** 2 + 
                4 * b1 * b2 ** 3 * m1 * m2 * r1 * r2 * s2 ** 2 - 2 * b2 ** 4 * m2 ** 2 * r2 ** 2 * s2 ** 2 - 
                2 * b1 ** 2 * b2 ** 2 * s1 ** 2 * s2 ** 2 + b2 ** 4 * s2 ** 4 + 4 * a2 * (b1 * m1 * r1 - b2 * m2 * r2) * 
                 (-2 * b1 * b2 * m1 * m2 * r1 * r2 + b1 ** 2 * (m1 ** 2 * (-1 + r1 ** 2) - s1 ** 2) + 
                  b2 ** 2 * (m2 ** 2 * (-1 + r2 ** 2) - s2 ** 2)) - 4 * a1 * (a2 + b1 * m1 * r1 - 
                  b2 * m2 * r2) * (a2 ** 2 - 2 * b1 * b2 * m1 * m2 * r1 * r2 + 2 * a2 * (b1 * m1 * r1 - 
                    b2 * m2 * r2) + b1 ** 2 * (m1 ** 2 * (-1 + r1 ** 2) - s1 ** 2) + 
                  b2 ** 2 * (m2 ** 2 * (-1 + r2 ** 2) - s2 ** 2)) + 2 * a2 ** 2 * 
                 (-6 * b1 * b2 * m1 * m2 * r1 * r2 + b1 ** 2 * (m1 ** 2 * (-1 + 3 * r1 ** 2) - s1 ** 2) + 
                  b2 ** 2 * (m2 ** 2 * (-1 + 3 * r2 ** 2) - s2 ** 2)) + 
                2 * a1 ** 2 * (3 * a2 ** 2 - 6 * b1 * b2 * m1 * m2 * r1 * r2 + 6 * a2 * (b1 * m1 * r1 - 
                    b2 * m2 * r2) + b1 ** 2 * (m1 ** 2 * (-1 + 3 * r1 ** 2) - s1 ** 2) + 
                  b2 ** 2 * (m2 ** 2 * (-1 + 3 * r2 ** 2) - s2 ** 2)))
    
    
    q3 = 1000000 * -4 * (b1 ** 4 * m1 * (-1 + r1 ** 2) ** 2 - b1 ** 3 * r1 * (-1 + r1 ** 2) * 
                 (a1 - a2 + b2 * (3 * m1 + m2) * r2) + b2 ** 3 * (-1 + r2 ** 2) * 
                 ((a1 - a2) * r2 + b2 * m2 * (-1 + r2 ** 2)) + b1 * b2 ** 2 * r1 * 
                 (a1 - 3 * a1 * r2 ** 2 - b2 * (m1 + 3 * m2) * r2 * (-1 + r2 ** 2) + 
                  a2 * (-1 + 3 * r2 ** 2)) + b1 ** 2 * b2 * ((a1 - a2) * (-1 + 3 * r1 ** 2) * r2 + 
                  b2 * (m1 + m2) * (-1 - r2 ** 2 + r1 ** 2 * (-1 + 3 * r2 ** 2))))
    
    q1 = 1000000 * 4 * (-(b1 ** 4 * m1 ** 3) + b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 + b1 ** 2 * b2 ** 2 * m1 * m2 ** 2 - 
                b2 ** 4 * m2 ** 3 + 2 * b1 ** 4 * m1 ** 3 * r1 ** 2 + b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 * r1 ** 2 + 
                b1 ** 2 * b2 ** 2 * m1 * m2 ** 2 * r1 ** 2 - b1 ** 4 * m1 ** 3 * r1 ** 4 - b1 ** 3 * b2 * m1 ** 3 * r1 * r2 - 
                3 * b1 ** 3 * b2 * m1 ** 2 * m2 * r1 * r2 - 3 * b1 * b2 ** 3 * m1 * m2 ** 2 * r1 * r2 - b1 * b2 ** 3 * m2 ** 3 * r1 * r2 + b1 ** 3 * b2 * m1 ** 3 * r1 ** 3 * r2 + 3 * b1 ** 3 * b2 * m1 ** 2 * m2 * 
                 r1 ** 3 * r2 + b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 * r2 ** 2 + b1 ** 2 * b2 ** 2 * m1 * m2 ** 2 * r2 ** 2 + 
                2 * b2 ** 4 * m2 ** 3 * r2 ** 2 - 3 * b1 ** 2 * b2 ** 2 * m1 ** 2 * m2 * r1 ** 2 * r2 ** 2 - 
                3 * b1 ** 2 * b2 ** 2 * m1 * m2 ** 2 * r1 ** 2 * r2 ** 2 + 3 * b1 * b2 ** 3 * m1 * m2 ** 2 * r1 * r2 ** 3 + 
                b1 * b2 ** 3 * m2 ** 3 * r1 * r2 ** 3 - b2 ** 4 * m2 ** 3 * r2 ** 4 + a1 ** 3 * (b1 * r1 - b2 * r2) + 
                a2 ** 3 * (-(b1 * r1) + b2 * r2) + a2 ** 2 * (b1 ** 2 * (m1 - 3 * m1 * r1 ** 2) + 3 * b1 * b2 * (m1 + m2) * r1 * r2 + b2 ** 2 * m2 * (1 - 3 * r2 ** 2)) + 
                a1 ** 2 * (b1 ** 2 * (m1 - 3 * m1 * r1 ** 2) + 3 * b1 * r1 * (-a2 + b2 * (m1 + m2) * r2) + 
                  b2 * (3 * a2 * r2 + b2 * (m2 - 3 * m2 * r2 ** 2))) - b1 ** 4 * m1 * s1 ** 2 + 
                b1 ** 2 * b2 ** 2 * m2 * s1 ** 2 + b1 ** 4 * m1 * r1 ** 2 * s1 ** 2 - b1 ** 3 * b2 * m1 * r1 * r2 * s1 ** 2 - 
                b1 ** 3 * b2 * m2 * r1 * r2 * s1 ** 2 + b1 ** 2 * b2 ** 2 * m2 * r2 ** 2 * s1 ** 2 + 
                b1 ** 2 * b2 ** 2 * m1 * s2 ** 2 - b2 ** 4 * m2 * s2 ** 2 + b1 ** 2 * b2 ** 2 * m1 * r1 ** 2 * s2 ** 2 - 
                b1 * b2 ** 3 * m1 * r1 * r2 * s2 ** 2 - b1 * b2 ** 3 * m2 * r1 * r2 * s2 ** 2 + 
                b2 ** 4 * m2 * r2 ** 2 * s2 ** 2 + a2 * (b1 ** 2 * b2 * r2 * (m1 ** 2 * (-1 + 3 * r1 ** 2) + 
                    2 * m1 * m2 * (-1 + 3 * r1 ** 2) - s1 ** 2) + b1 ** 3 * r1 * (-3 * m1 ** 2 * 
                     (-1 + r1 ** 2) + s1 ** 2) + b2 ** 3 * r2 * (3 * m2 ** 2 * (-1 + r2 ** 2) - s2 ** 2) + 
                  b1 * b2 ** 2 * r1 * (m1 * m2 * (2 - 6 * r2 ** 2) + m2 ** 2 * (1 - 3 * r2 ** 2) + s2 ** 2)) + 
                a1 * (3 * a2 ** 2 * (b1 * r1 - b2 * r2) + a2 * (2 * b1 ** 2 * m1 * (-1 + 3 * r1 ** 2) - 
                    6 * b1 * b2 * (m1 + m2) * r1 * r2 + 2 * b2 ** 2 * m2 * (-1 + 3 * r2 ** 2)) + 
                  b1 ** 3 * r1 * (3 * m1 ** 2 * (-1 + r1 ** 2) - s1 ** 2) + b1 ** 2 * b2 * r2 * ( 
                    m1 * m2 * (2 - 6 * r1 ** 2) + m1 ** 2 * (1 - 3 * r1 ** 2) + s1 ** 2) + 
                  b1 * b2 ** 2 * r1 * (2 * m1 * m2 * (-1 + 3 * r2 ** 2) + m2 ** 2 * (-1 + 3 * r2 ** 2) - 
                    s2 ** 2) + b2 ** 3 * r2 * (-3 * m2 ** 2 * (-1 + r2 ** 2) + s2 ** 2)))
    
    term16 = (2 * q2 ** 3 + 27 * q3 ** 2 * q0 - 72 * q4 * q2 * q0 - 9 * q3 * q2 * q1 + 27 * q4 * q1 ** 2)
    
    term21 = (q2 ** 2 / 4 + 3 * q4 * q0 - 3 * q3 * q1 / 4)
    
    term1sq = -256 * term21 ** 3 + term16 ** 2
    
    term1 = cmath.sqrt(term1sq+0*1j)  # Note use of complex arithmetic in R
    
    term23 = (term16 + term1) ** (1/3)
    
    term22 = 3*q4*term23 
    
    temp1 = (4 * 2 ** (1 / 3) * term21)
    temp2 = (3 * 2 ** (1 / 3) * q4)
    temp3 = q3 ** 2 / (4 * q4 ** 2) - (2 * q2) / (3 * q4)
    temp4 = temp1/ term22 + term23/temp2 

    rr = cmath.sqrt(temp3+temp4) 
    
    temp5 = q3 ** 2 / (2 * q4 ** 2) - (4 * q2) / (3 * q4)
    temp6 = (-q3 ** 3 / 4 + q4 * q3 * q2 - 2 * q4 ** 2 * q1) / (q4 ** 3)
    
    ee = q3 ** 2 / (2 * q4 ** 2) - (4 * q2) / (3 * q4) - (4 * 2 ** (1 / 3) * term21) / term22 - term23 / (3 * 2 ** (1 / 3) * q4) -(-q3 ** 3 / 4 + q4 * q3 * q2 - 2 * q4 ** 2 * q1) / (q4 ** 3 * rr) 
    
    dd = q3 ** 2 / (2 * q4 ** 2) - (4 * q2) / (3 * q4) - (4 * 2 ** (1 / 3) * term21) / term22 - term23 / (3 * 2 ** (1 / 3) * q4) + (-q3 ** 3 / 4 + q4 * q3 * q2 - 2 * q4 ** 2 * q1) / (q4 ** 3 * rr) 
    
    temp7 = -q3 / (4 * q4) 
    
    # Potential roots are given by
    roots = np.array([
            -q3/(4*q4)+rr/2+cmath.sqrt(dd)/2,
            -q3/(4*q4)+rr/2-cmath.sqrt(dd)/2,
            -q3/(4*q4)-rr/2+cmath.sqrt(ee)/2,
            -q3/(4*q4)-rr/2-cmath.sqrt(ee)/2])

    # Need to check these are really roots
    kr = roots*(np.abs(np.imag(roots)) < 10**(-10)) 
    def test(k):
        return np.array([(a1 + b1 * (r1 * (_k - m1) + cmath.sqrt((_k - m1) ** 2 + s1 ** 2))) - (a2 + b2 * (r2 * (_k - m2) + cmath.sqrt((_k - m2) ** 2 + s2 ** 2))) for _k in k])
        
            
    num = np.where(np.abs(test(kr))<10**(-10)) # Which potential root is actually a root?

    roots = np.sort(np.real(kr[num]))# Roots in ascending order
    nRoots = len(roots)
    
    crossedness = 0

    if nRoots>1:
        midPoints = (roots[0:(nRoots-1)]+roots[1:nRoots])/2 
    else:
        midPoints=np.array([])
    if(nRoots>0):
        samplePoints = np.array([[roots[0]-1] + midPoints.tolist()+[roots[nRoots-1]+1]])
        short_params = sviData.iloc[0,:].to_dict()
        short_params['k'] = samplePoints
        long_params = sviData.iloc[1,:].to_dict()
        long_params['k'] = samplePoints
        sviShort = svi(**short_params) 
        sviLong = svi(**long_params)
        crossedness = np.max(np.array([0,np.max(sviShort-sviLong)]))  # Maximal amount of crossing bounded below by zero
 
    return dict(roots=roots,crossedness=crossedness)

