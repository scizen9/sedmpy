# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 18:48:13 2015

@author: nadiablago
"""
import numpy as np
from matplotlib import pylab as plt
from scipy import optimize


def doGauss(p, w, data):
    mod = np.ones((w.size,2))
    mod[:,1] = np.exp(-0.5*(p[0]-w)**2/p[1]**2)
    fit,chi = optimize.nnls(mod,data)
    return (fit*mod).sum(1) - data
    
def getGauss(p, w, data):
    mod = np.ones((w.size,2))
    mod[:,1] = np.exp(-0.5*(p[0]-w)**2/p[1]**2)
    fit,chi = optimize.nnls(mod,data)
    return (fit*mod).sum(1) - data

def doModel(p,w,d,getMod=False):
    l1,s1 = p
    mod = np.ones((w.size,3))
    s1sq = s1**2
    mod[:,1] = -1.*np.exp(-0.5*(l1-w)**2/s1sq)
    fit,chi = optimize.nnls((mod.T).T,d)
    if getMod:
        return (fit*mod).sum(1)
    return ((fit*mod).sum(1)-d)
    
def showspec(npyfile):
    s = np.load(npyfile)[0]
    
    hwl =np.array( [3970.07, 4101.76, 4340.47, 4861.33, 6562.80])/10.
    tll = np.array([6875, 7610])/10.
    
    #for si in s['sky_spaxel_ids_A'][0]['spectra']:
    for si in s['spectra']:
        plt.plot(s['nm'], si)
        
    mmax = np.max(s['spectra'])
    mmin = np.min(s['spectra'])
    
    for w in hwl:
        plt.vlines(w, mmin, mmax)
        
    for w in tll:
        plt.vlines(w, mmin, mmax, 'b')
    plt.show()
    
def rms_wl(npyfile, plot=False):
    
    s = np.load(npyfile)[0]
    w = s['nm']*10.
    
    hwl =np.array( [4101.76, 4340.47, 4861.33, 6562.80])
    #hwl =np.array( [3970.07, 4101.76, 4340.47, 4861.33, 6562.80])

    tflux = np.sum(s['spectra'], axis=1)
    
    order = np.array([np.arange(len(tflux)), tflux])
    order = np.sort(order, axis=1).T
    
    allwave = []
        
    for wh in hwl:
        for i in order[0:35, 0]:
            spec = s['spectra'][i]
            #plt.plot(w, spec)
            c = abs(w-wh) < 150
            fit,ier = optimize.leastsq(doModel,[wh,50.],(w[c],spec[c]))
            mod = doModel(fit, w[c], spec[c], True)
            if (plot):
                plt.plot(w[c], spec[c])
                plt.plot(w[c], mod)
            wlInt = fit[0]            
            sigInt = fit[1]
            #print "Wavelength", wlInt, "Sigma", sigInt, "Fit", fit
            if (sigInt) < 60:
                allwave.append(wlInt-wh)
        if (plot):
            plt.show()
        
        print "Wavelength", wh, "Median deviation", np.median(allwave), "Std", np.std(allwave)
        plt.hist(allwave)
        print npyfile.replace('.npy', '%.1d.png'%wh)
        plt.savefig(npyfile.replace('.npy', '%.1d.png'%wh).replace('data/spec', 'plots/wavelengths'))
        plt.clf()

    
    