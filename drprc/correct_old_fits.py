# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 18:09:41 2016

@author: nadiablago
"""

import glob, os, subprocess
import numpy as np
import fitsutils


def correct_files(dirname):
    
    cmd = "gethead OBJECT *fits | grep -v -i calib | grep PTF | grep -v -E \"\[A|\[B|findi\" | grep ifu | awk '{print $1}' > update_science"
    subprocess.all(cmd, shell=True)
    
    filelist = os.path.join(dirname, "update_science")
    
    for f in np.genfromtxt(filelist, dtype=None):
        fitsutils.update_par(f, "OBJDEC", fitsutils.get_par(f, "DEC"))
        
    for f in np.genfromtxt(filelist, dtype=None):
        fitsutils.update_par(f, "OBJRA", fitsutils.get_par(f, "RA"))
    
    for f in np.genfromtxt(filelist, dtype=None):
        fitsutils.update_par(f, "FILTER", fitsutils.get_par(f, "OBJECT").split()[0].replace("]","").replace("[", ""))
        
    for f in np.genfromtxt(filelist, dtype=None):
        fitsutils.update_par(f, "NAME", fitsutils.get_par(f, "OBJECT").split()[1])
        
    for f in np.genfromtxt(filelist, dtype=None):
        fitsutils.update_par(f, "IMGTYPE", "SCIENCE")