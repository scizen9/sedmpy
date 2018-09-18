# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:34:47 2016

@author: nadiablago
"""
import glob
import fitsutils
import rcred
from astropy.io import fits as pf
import numpy as np
import datetime
import time_utils
import coordinates_conversor as cc
from matplotlib import pylab as plt
from scipy.optimize import curve_fit
from astropy.table import Table
        
def get_flats_counts(directory):
    '''
    Reads all the images in a directory marked as "twilight" and uses them to infer the average count rate for a given
    number of seconds after the sunset.
    
    It plots the number of counts and fits a 2D polynomial used to interpolate in the future.
    
    '''

    flatlist = []
    corners = {
    "g" : [1, 910, 1, 900],
    "i" : [1, 910, 1060, 2045],
    "r" : [1040, 2045, 1015, 2045],
    "u" : [1030, 2045, 1, 900]
    }
    
    for f in glob.glob(directory + "/rc*fits"):
        if fitsutils.has_par(f, "IMGTYPE") and fitsutils.get_par(f, "IMGTYPE") == "TWILIGHT":
            flatlist.append(f)
            
    if (len(flatlist)==0):
        print "No suitable twilight flats found in directory: %s"%directory
        return
        
    counts = {"u":[], "g":[], "r":[], "i":[]}
    sun_decs = {"u":[], "g":[], "r":[], "i":[]}
    colors = {"u":"purple", "g":"green", "r":"r", "i":"orange"}
    
    for f in flatlist:
        print f
        print fitsutils.get_par(f, "JD")
        
        bias = rcred.get_overscan_bias_rc(f)
        exptime = fitsutils.get_par(f, "EXPTIME")
        sunsettime = fitsutils.get_par(f, "SUNSET")
        utc = time_utils.jd2utc(fitsutils.get_par(f, "JD"))
        st = datetime.datetime.strptime(sunsettime, "%H:%M")
        elapsed = 3600*(utc.hour - st.hour) + 60*(utc.minute - st.minute) + (utc.second - st.second)
        
        if (elapsed > 5000):
            continue
        print elapsed, utc, st
        
        
        data = pf.open(f)[0].data
        for band in corners.keys():
            c = corners[band]
            sf = data.T[c[0]:c[1], c[2]:c[3]]
            if (np.percentile(sf, 90) < 55000):
                counts[band].append( (np.percentile(sf, 90)-bias)/exptime)
                sun_decs[band].append(elapsed)        
    
    
                 
    t = Table(names=('filter', 'c2', 'c1', 'c0'), dtype=('S1', 'f8', 'f8', 'f8'))


    for band in corners.keys():
        print sun_decs[band], np.log10(counts[band])
        coefs = np.polyfit(sun_decs[band], np.log10(counts[band]), deg=2)#, w=1./np.sqrt(np.log10(counts[band])))
        p = np.poly1d(coefs)
        x = np.linspace(np.min(sun_decs[band]), np.max(sun_decs[band]), 1000)
        
        t.add_row([band, coefs[0], coefs[1], coefs[2]])
        
        plt.plot(sun_decs[band], np.log10(counts[band]), "o", label=band, color=colors[band])
        plt.plot(x, p(x), label="Model "+band, color=colors[band])
        
        
    t.write("/tmp/test_flat", format='csv')
    
    plt.xlabel("Elapsed second since Sunset")
    plt.ylabel("log Counts/s")
    plt.legend()
    plt.show()
    

def get_flat_itime(band, elapsed, counts=20000):
    '''
    Provides the optimum integration time to achieve the desired number of counts.
    Its parameters are the filter in which we are observing, and the time elapsed since sunset.
    
    '''
    
    logcounts = np.log10(counts)
    
    t = Table.read('/tmp/test_flat', format='csv')
    tf = t[t["filter"]==band][0]
    coefs = np.array([tf["c2"], tf["c1"], tf["c0"]])

    p = np.poly1d(coefs)
    
    pint = p.integ()
    
    p0 = pint(elapsed)
    
    print "Initial counts", logcounts
    
    #Line up the different exposure times. Maximum 100s after the elapsed time
    l = np.arange(400) + elapsed
    
    exptime = 0
    myint = 0
    while myint < counts:
        myint = np.trapz(10**p(l[0:exptime]), l[0:exptime])
        exptime = exptime + 1
    
    print "Number of counts",myint, "Exposure time",exptime
    
    
    
    return exptime
    
    
def test_prediction(f):
    
    corners = {
    "g" : [1, 910, 1, 900],
    "i" : [1, 910, 1060, 2045],
    "r" : [1040, 2045, 1015, 2045],
    "u" : [1030, 2045, 1, 900]
    }

    data = pf.open(f)[0].data
    for band in corners.keys():
        c = corners[band]
        sf = data.T[c[0]:c[1], c[2]:c[3]]
        p = np.percentile(sf, 90) < 55000
        # counts[band].append( (np.percentile(sf, 90)-bias)/exptime)

