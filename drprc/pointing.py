# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 17:59:15 2015

@author: nadiablago
"""
from astropy.wcs import WCS
import glob, os
import numpy as np
from astropy.io import fits as pf
import coordinates_conversor as cc
from matplotlib import pylab as plt
import matplotlib

def plot_offset_shift(dirname):
    dras = []
    ddecs = []
    d = []
    r = []
    m = []
    for f in glob.glob(os.path.join(dirname,"*new")):
        image = pf.open(f)
        wcs = WCS(image[0].header)
        rra, rdec = cc.hour2deg(image[0].header['RA'],image[0].header['DEC'] )
        pra, pdec = wcs.wcs_pix2sky(np.array([[1293., 1280.]] , np.float_), 1)[0]
        dra, ddec = cc.get_offset(pra, pdec, rra, rdec)
        
        if np.abs(dra) > 100 or np.abs(ddec)>100:
            continue
        print f, image[0].data.shape , "(",rra, rdec, ")  vs. (",  pra, pdec, ")", dra, ddec
        dras.append(dra)
        ddecs.append(ddec)
        d.append(rdec)
        r.append(rra)
        m.append(image[0].header['JD'])
        
    dras = np.array(dras)
    ddecs = np.array(ddecs)
    
    plt.scatter(dras, ddecs, c=np.array(r), cmap=matplotlib.cm.jet, s=130)
    plt.xlabel('dRA [arcsec]')
    plt.ylabel('dDEC [arcsec]')
    cb = plt.colorbar(label='RA [deg]')

    f = plt.figure()
    plt.scatter(dras, ddecs, c=np.array(d), cmap=matplotlib.cm.jet, s=120)
    plt.xlabel('dRA [arcsec]')
    plt.ylabel('dDEC [arcsec]')
    cb = plt.colorbar(label='DEC [deg]')

    
    f = plt.figure()
    plt.scatter(dras, ddecs, c=(np.array(m)-np.min(m))*24, cmap=matplotlib.cm.jet, s=130)
    plt.xlabel('dRA [arcsec]')
    plt.ylabel('dDEC [arcsec]')
    cb = plt.colorbar(label='JD [hours since first image]')
    
    plt.show()
