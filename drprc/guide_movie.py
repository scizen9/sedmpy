# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:01:50 2015

@author: nadiablago
"""
import subprocess
from astropy.io import fits as pf

import numpy as np
import aplpy
try:
    import coordinates_conversor
except ImportError:
    import drprc.coordinates_conversor as coordinates_conversor
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
import datetime
import os
import sys
import glob
import argparse
from configparser import ConfigParser
import codecs

from matplotlib import pylab as plt
plt.switch_backend('Agg')

cfg_parser = ConfigParser()

try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                              '../config/sedmconfig.cfg')

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    cfg_parser.read_file(f)

# Default paths
_photpath = cfg_parser.get('paths', 'photpath')

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    
    Creates an animated gif of guider images used for the science image
    in the folder specified as a parameter.
        
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--imfile', type=str, dest="imfile",
                        help='IFU image that requires a finder',
                        default=None)
    
    args = parser.parse_args()

    imfile = args.imfile

    if imfile:
        timestamp = imfile.split('/')[-2]
        rcdir = os.path.join(_photpath, timestamp)
        reduxdir = '/'.join(imfile.split('/')[0:-1])
        objnam = fitsutils.get_par(imfile, "OBJECT").split()[0]
        if 'STD' in objnam:
            objnam = objnam.split('STD-')[-1].split()[0]

        os.chdir(reduxdir)

        print("Making finder for object: %s" % objnam)
        print("Changed to directory where the reduced data is: %s" % reduxdir)
        print("Getting acquisition images from directory: %s" % rcdir)

        # We gather all RC images to locate the Guider ones.
        files = glob.glob(os.path.join(rcdir, "rc*[0-9].fit*"))
        files.sort()
        filesguide = []
        pngdir = os.path.join(rcdir, 'pngraw/guider')

        for f in files:
            try:
                ff = pf.open(f)
            except OSError:
                print("WARNING - corrupt fits file: %s" % f)
                continue
            if "IMGTYPE" in ff[0].header:
                imgtype = ff[0].header["IMGTYPE"]
            else:
                imgtype = ''
            if "OBJECT" in ff[0].header:
                obj = ff[0].header["OBJECT"]
            else:
                obj = ''

            ff.close()

            if 'GUIDE' in imgtype.upper() and objnam in obj:
                pngf = os.path.basename(f).split('.fits')[0] + '_all.png'
                filesguide.append(os.path.join(pngdir, pngf))

        n_guide = len(filesguide)
        print("Found %d files used for quiding %s" % (n_guide, objnam))
        outmov = os.path.join(
            pngdir, os.path.basename(imfile).split('.fits')[0] + '.gif')
        cmd = 'convert -delay 20 ' + ' '.join(filesguide) + ' -loop 1 ' + outmov
        print(cmd)
        subprocess.run(cmd, shell=True)
    else:
        print("Please specify an input IFU image that was guided.")


