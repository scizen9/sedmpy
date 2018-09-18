# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 12:16:50 2015

@author: nadiablago
"""

from astropy.time import Time
import datetime

def utc2mjd(times):
    t = Time(times, format='isot', scale='utc')
    return t.mjd
    
def mjd2utc(mjd, string=False):
    t = Time(mjd+2400000.5, format='jd', scale="utc")
    
    if (string):
        return t.iso
    else:
        return datetime.datetime.strptime(t.iso, "%Y-%m-%d %H:%M:%S.%f")  
    
def jd2utc(jd, string=False):
    t = Time(jd, format='jd', scale="utc")
    
    if (string):
        return t.iso
    else:
        return datetime.datetime.strptime(t.iso, "%Y-%m-%d %H:%M:%S.%f") 
    
def utc2jd(times):
    t = Time(times, format='iso', scale='utc')
    return t.jd
