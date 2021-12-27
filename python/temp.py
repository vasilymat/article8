# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 08:57:05 2021

@author: Lab
"""

import numpy as np

def sqx(x):
    return (x+7)**4

def gradsolv(start_point,f):
    #f - function for optimization
    n1 = 0
    dx = 0.05
    delta = 0.1
    grd = (f(start_point+dx) - f(start_point))    
    xs = start_point
    y0 = f(xs)
    xs = start_point
    while n1 < 500:        
        xs = xs - delta*grd
        y1 = f(xs)
        grd = f(xs+dx) - y1
        n1 += 1
        y0 = y1
        if np.abs(grd) < 0.001:
            break
            
    print(grd)
    print(n1)
    return y0, xs

a,b = gradsolv(1,sqx)
print(a,b)
