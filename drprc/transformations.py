# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 15:25:46 2015

@author: nadiablago
"""
import numpy as np

def johnson2sdss(U, B, V, R, I):
    '''global transformations between UBVRI and ugriz'''
    
    '''
    #Color Color Term Zeropoint Range
    "gV": (0.630 ± 0.002) (B − V) −(0.124 ± 0.002)
    "ri": (1.007 ± 0.005) (R − I) −(0.236 ± 0.003)
    "rz": (1.584 ± 0.008) (R − I) −(0.386 ± 0.005)
    "rR": (0.267 ± 0.005) (V − R) +(0.088 ± 0.003) V − R ≤ 0.93
    "rR": (0.77 ± 0.04) (V − R) −(0.37 ± 0.04) V − R > 0.93
    "ug": (0.750 ± 0.050) (U − B) + (0.770 ± 0.070) (B − V) +(0.720 ± 0.040)
    "gB": −(0.370 ± 0.002) (B − V) −(0.124 ± 0.002)
    "gr": (1.646 ± 0.008) (V − R) −(0.139 ± 0.004)
    "iI": [0.247, 0.329]'''
    

def sdss2johnson(ref_sdss, savefile=None):
    '''
    Jordi et. al 2006
    
    ugriz -> UBVRcIc
    ================
    
            Transformation
        U-B   =     (0.79 ± 0.02)*(u-g)    - (0.93 ± 0.02)
        U-B   =     (0.52 ± 0.06)*(u-g)    + (0.53 ± 0.09)*(g-r) - (0.82 ± 0.04)
        B-g   =     (0.175 ± 0.002)*(u-g)  + (0.150 ± 0.003)
        B-g   =     (0.313 ± 0.003)*(g-r)  + (0.219 ± 0.002)
        V-g   =     (-0.565 ± 0.001)*(g-r) - (0.016 ± 0.001)
        V-I   =     (0.675 ± 0.002)*(g-i)  + (0.364 ± 0.002) if  g-i <= 2.1
        V-I   =     (1.11 ± 0.02)*(g-i)    - (0.52 ± 0.05)   if  g-i >  2.1
        R-r   =     (-0.153 ± 0.003)*(r-i) - (0.117 ± 0.003)
        R-I   =     (0.930 ± 0.005)*(r-i)  + (0.259 ± 0.002)
        I-i   =     (-0.386 ± 0.004)*(i-z) - (0.397 ± 0.001)
    '''
    
    ref_sdss = np.genfromtxt(ref_sdss, dtype=None, names=True, delimiter=',')

    bands = "BVRI"
    john = np.zeros(len(ref_sdss), dtype=[('id', '<i8'), ('ra', '<f8'), ('dec', '<f8'), \
    ('U', '<f4'), ('B', '<f4'), ('V', '<f4'), ('R', '<f4'), ('I', '<f4'),\
    ('dU', '<f4'), ('dB', '<f4'), ('dV', '<f4'), ('dR', '<f4'), ('dI', '<f4')])

    band_dic = {"B":"g", "V":"g", "R":"r", "I":"i"}
    coldic = {"U":"ug", "B":"gr", "V":"gr", "R":"ri", "I":"iz"}
    coefs = {"U": [np.array([0.79, 0.93]), np.array([0.02, 0.02])],
            "B": [np.array([0.313, 0.219]), np.array([0.003, 0.002])],
             "V": [np.array([-0.565, 0.016]), np.array([0.001, 0.001])],
            "R": [np.array([-0.153, 0.117]), np.array([0.003, 0.003])],
            "I": [np.array([-0.386, 0.397]), np.array([0.004, 0.001])] }
            
    for b in bands:
        col = ref_sdss[coldic[b][0]] - ref_sdss[coldic[b][1]]
        john[b] = np.sum(np.array([col, 1]) * coefs[b][0]) + ref_sdss[band_dic[b]]
        john["d"+b] = np.sum(np.array([col, 1]) * coefs[b][1])

    #U band a bit different
    b = "U"
    col = ref_sdss[coldic[b][0]] - ref_sdss[coldic[b][1]]
    john[b] = np.sum(np.array([col, 1]) * coefs[b][0]) + john["B"]
    john["d"+b] = np.sum( np.array([col, 1]) * coefs[b][1] )
    
    john["ra"] = ref_sdss["ra"]
    john["dec"] = ref_sdss["dec"]
    john["id"] = ref_sdss["objid"]
        
        
    if (not savefile is None):
        np.savetxt(savefile, john, header="id,ra,dec,U,B,V,R,I,dU,dB,dV,dR,dI", fmt="%d,%.5f,%.5f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f")
    return john