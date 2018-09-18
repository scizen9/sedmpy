# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 12:16:14 2015

@author: nadiablago
"""
import numpy as np
import datetime, time
import pickle 
import os

from matplotlib import pylab as plt

def get_spectral_window(t, nu=None):
    '''
    Computes the spectral window for the set of observations t, for each frequency in nu.
    Source: equation 2, Eyer & Bartholdi, 1998 (http://arxiv.org/pdf/astro-ph/9808176v1.pdf)
    '''
    
    if(nu==None):
        nu = np.linspace(0. ,20., 10000)

    N = len(t)
    
    #Transform datetime to seconds elapsed since 1970.
    t = np.array([time.mktime(ti.timetuple())/(24*3600.) for ti in t])
    
    T, NU = np.meshgrid(t, nu)
    GNnu = np.abs(np.sum(np.exp(+1j * 2 * np.pi* NU * T), axis=1))**2/N**2
    
    return nu, GNnu 
    
def get_time_lag(tobs):
    '''
    Computes the time lag for a given set of observations.
    '''
    timelag = []
    
    #Compute lag for each field
    for tcurr in tobs:
        latert = tobs[np.array(tobs>tcurr)]
        diffs = latert - tcurr
        diffsd = np.array([np.log10(d.total_seconds()/(24*3600.)) for d in diffs])
        timelag.extend(diffsd)
        
    return timelag
        
def mkdate(text):
    return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S.%f') 
    
def get_field_stats(myfile="/scr2/nblago/kpy/SEDMrph/data/VariableStars/iPTF_exposures.tbl", myfields=[]):
    '''
    Writes the observed dates and coordinates for each PTF field.
    '''
    basedir = os.path.dirname(myfile)
    
    t = np.genfromtxt(myfile, converters={'dateobs':mkdate}, names=True, dtype=None)
    c = np.genfromtxt(os.path.join(basedir, "ptf_corners.tbl"), names=True, 
                      dtype=[('id', '<i8'), ('ra0', '<f5'), ('dec0', '<f5'), ('ra1', '<f5'), ('dec1', '<f5'), ('ra2', '<f5'), ('dec2', '<f5'), ('ra3', '<f5'), ('dec3', '<f5'), ('ra4', '<f5'), ('dec4', '<f5')])


    #print set(t['field'])
    
    if len(myfields) == 0:
        myfields = set(t['field'])
        myfields = list(myfields)
        myfields.sort()
        
        
    fields = np.zeros(len(set(t['field'])), dtype=[('field', '<i8'), 
                                              ('ra0', '<f8'), ('dec0', '<f8'), ('ra1', '<f8'), ('dec1', '<f8'),
                                                ('ra2', '<f8'), ('dec2', '<f8'),('ra3', '<f8'), ('dec3', '<f8'),
                                                ('ra4', '<f8'), ('dec4', '<f8'), 
                                                ('obsdate', np.object)])#, ('obssec', np.object), ('timelags', np.object)])
    for i, f in enumerate(myfields):
        mask = np.array(t['field']==f)
        tf = t[mask]

        
        #Fill in the output vector
        fields[i]['field'] = f
        for j in range(5):
            fields[i]['ra%d'%j] = c[c['id']==f]['ra%d'%j]        
            fields[i]['dec%d'%j] = c[c['id']==f]['dec%d'%j] 
        fields[i]['obsdate'] = tf['dateobs']

        
    pickle.dump(fields, open(os.path.join(basedir, 'field_stats.pickle'), "w"))

    return fields
    
def plot_stats(stats, myfields=[]):
    """
    Plots the time lag and spectral window for all fields indicated into the folder "plots".
    """
    if len(myfields) == 0:
        myfields = set(stats['field'])
        myfields = list(myfields)
        myfields.sort()
        
    timelag = {}
        
    for i, f in enumerate(myfields):
        mask = np.array(stats['field']==f)
        tf = stats[mask][0]
        
        #Compute time lag.
        timelag[f] = get_time_lag(tf['obsdate'])
        
        #Compute spectral window
        nu = np.linspace(0. ,20., 10000)
        nu, w = get_spectral_window(tf['obsdate'], nu)

        if (not os.path.isdir("plots")):
            os.makedirs("plots")
            
        if (len(tf['obsdate']) >1):
            plt.hist(timelag[f], bins=100, range=(-3.5,4))
            plt.xlabel('log(Time$_i$ - Time$_j$ [days])')
            plt.ylabel('Counts')
            plt.savefig('plots/field_lag_%d.png'%f)
            plt.clf()
            
            plt.plot(nu, w)
            plt.xlabel('Frequency [1/day]')
            plt.ylabel('Spectral window')
            plt.savefig('plots/field_swindow_%d.png'%f)
            plt.clf()
    

    
def load_file(myfile):
    """
    Loads the pickled file passed as a parameter.
    """
    f = pickle.load(open(myfile, "r"))
    return f
    