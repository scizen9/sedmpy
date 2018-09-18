# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 12:22:15 2016

@author: nadiablago
"""
import argparse
import fitsutils
import glob, os, time
import recenter_ifu


def compute_offsets(fitsdir):
    ra_prev = "00:00:00"
    dec_prev = "00:00:00"
    
    lfiles = glob.glob(fitsdir + "/rc*fits")
    lfiles.sort()
    
    
    for f in lfiles:
        f = os.path.basename(f)
        try:

            if ("[r]" in fitsutils.get_par(f, "OBJECT")  or fitsutils.get_par(f, "OBJTYPE")=="ACQUISITION"):
                ra = fitsutils.get_par(f, "OBJRA")
                dec = fitsutils.get_par(f, "OBJDEC")
                obj = fitsutils.get_par(f, "OBJECT")
                isAB = False
                if (fitsutils.has_par(f, "ABPAIR") and fitsutils.get_par(f, "ABPAIR")=="True"):
                    isAB = True
                if ((ra!="" and dec !="" and (ra != ra_prev or dec != dec_prev)) or (fitsutils.get_par(f, "OBJTYPE")=="ACQUISITION")):
                    #print "Found image %s as first acquisition image after the slew. Computing offset for IFU..."%f
                    recenter_ifu.main(os.path.abspath(f), isAB, astro=True, plot=True)
                ra_prev = ra
                dec_prev = dec
        except KeyError:
            print "File %s does not have a good header."%f
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--directory', dest='directory', type=str, help='Directory containing the fits for the night.')

    args = parser.parse_args()
    
    fitsdir = args.directory
    fitsdir = os.path.abspath(fitsdir)

    os.chdir(fitsdir)
    
    print os.getcwd()
    
    while (True):
        compute_offsets(fitsdir)
        print "Brief nap before scanning for more new images..."
        time.sleep(5)
    
