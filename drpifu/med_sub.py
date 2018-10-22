# -*- coding: utf-8 -*-
"""

@author: neill
"""
from astropy.io import fits as pf
from scipy import ndimage
import shutil
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    
    Subtract a median filtered background.
        
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--imfile', type=str, dest="imfile",
                        help='IFU image that requires a finder',
                        default=None)
    parser.add_argument('-s', '--size', type=int, dest="size",
                        help='Median filter size (default = 15)', default=15)
    
    args = parser.parse_args()

    imfile = args.imfile
    fsize = args.size

    if os.path.islink(imfile):
        print("Make a local copy...")
        realpath = os.path.realpath(imfile)
        os.unlink(imfile)
        shutil.copy2(realpath, imfile)

    ff = pf.open(imfile, mode='update')
    image = ff[0].data * 1.0    # convert to float
    print("Create median filtered image")
    medfilt = ndimage.median_filter(image, fsize, mode='constant')
    ff[0].data = image - medfilt
    ff[0].scale('float32', 'minmax')
    ff.close()
    print("Done.")
