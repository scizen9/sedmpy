#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 10:30:30 2018

@author: nblago
"""
from __future__ import print_function


import glob, os
import fitsutils
import numpy as np
import coordinates_conversor as cc
from pyraf import iraf 
import subprocess
import shutil
import argparse

def solve_astrometry(img, outimage=None, radius=3, with_pix=True, overwrite=True, tweak=3):
    '''
    img: fits image where astrometry should be solved.
    outimage: name for the astrometry solved image. If none is provided, the name will be "a_"img.
    radius: radius of uncertainty on astrometric position in image.
    with_pix: if we want to include the constraint on the pixel size for the RCCam.
    overwrite: wether the astrometrically solved image should go on top of the old one.
    tewak: parameter for astrometry.net
    '''

    from astropy.wcs import InconsistentAxisTypesError
    
    curdir = os.getcwd()
    imgdir = os.path.dirname(img)
    
    os.chdir(imgdir)
    
    img = os.path.abspath(img)
    
    ra = fitsutils.get_par(img, 'RA')
    dec = fitsutils.get_par(img, 'DEC')
    #logger.info( "Solving astrometry on field with (ra,dec)=%s %s"%(ra, dec))
    
    astro = os.path.join( os.path.dirname(img), "a_" + os.path.basename(img))


    #If astrometry exists, we don't run it again.
    if (os.path.isfile(astro) and not overwrite):
        return astro
        

    cmd = "solve-field --ra %s --dec %s --radius %.4f -p --new-fits %s -W none -B none -P none -M none -R none -S none -t %d --overwrite %s "%(ra, dec, radius, astro, tweak, img)
    if (with_pix):
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.4"
    
    print (cmd)

    subprocess.call(cmd, shell=True)
    
    print ("Finished astrometry")
    
    #Cleaning after astrometry.net
    if (os.path.isfile(img.replace(".fits", ".axy"))):
        os.remove(img.replace(".fits", ".axy"))
    if (os.path.isfile(img.replace(".fits", "-indx.xyls"))):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if (os.path.isfile("none")):
        try:
            os.remove("none")
        except:
            print ("Could not remove file none.")
        
    os.chdir(curdir)

    if (not outimage is None and overwrite and os.path.isfile(astro)):
        shutil.move(astro, outimage)
        return outimage
    elif (outimage is None and overwrite and os.path.isfile(astro)):
        shutil.move(astro, img)
        return img
    else:
        return astro
    
def create_masterguide(lfiles, out=None):
    '''
    Receives a list of guider images for the same object.
    It will remove the bias from it, combine them using the median, and comput the astrometry for the image.
    
    '''
    
    curdir = os.getcwd()
    
    if len(lfiles) == 0:
        return
    else:
        os.chdir(os.path.abspath(os.path.dirname(lfiles[0])))
        

    fffile ="/tmp/l_guide"
    np.savetxt(fffile, np.array(lfiles), fmt="%s")

    #If the bias file exists in the directory, we use it, otherwise we pass
    bias_fast = "Bias_rc_fast.fits"
    debias = os.path.isfile(bias_fast)
        
    if debias:
        debiased = ["b_" + os.path.basename(img) for img in lfiles]
        bffile ="/tmp/lb_guider"
        np.savetxt(bffile, np.array(debiased), fmt="%s")


    if (out is None):
        obj = fitsutils.get_par(img, "OBJECT")
        out = os.path.join( os.path.dirname(img), obj.replace(" ", "").replace(":", "")+".fits")


    # Running IRAF
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    #Remove bias from the guider images    
    if debias:
        try:
            iraf.imarith("@"+fffile, "-", bias_fast, "@"+bffile)    
        except:
            print ("Error when debiasing file.")
    else:
        bffile = fffile

        
    #Combine flats
    iraf.imcombine(input = "@"+bffile, \
                    output = out, \
                    combine = "median",\
                    scale = "mode",
                    reject = "sigclip", lsigma = 2., hsigma = 2, gain=1.7, rdnoise=4.)
    #iraf.imstat(out, fields="image,npix,mean,stddev,min,max,mode", Stdout="guide_stats")
    #st = np.genfromtxt("guide_stats", names=True, dtype=None)
    
    #Do some cleaning
    if debias:
        print ( 'Removing from lfiles')
        for f in debiased:
            if os.path.isfile(f): os.remove(f)

    if os.path.isfile(fffile):
        os.remove(fffile)
    if os.path.isfile(bffile):
        os.remove(bffile)
        
    solve_astrometry(out, overwrite=True)
    
    os.chdir(curdir)
    

def __fill_ifu_dic(ifu_img, ifu_dic = {}):
    '''
    Fills some of the parameters for the IFU image that will be used to gather its guider images later.
    
    '''
    imgtype = fitsutils.get_par(ifu_img, "IMGTYPE")
    if not imgtype is None:
        imgtype = imgtype.upper()

    #In Richard's pipeline, the JD is the beginning of the exposure,
    #in Nick's one is the end.
    pipeline_jd_end = fitsutils.get_par(ifu_img, "TELESCOP") == '60'
    
    if imgtype == "SCIENCE" or imgtype =="STANDARD":
        if pipeline_jd_end:
            jd_ini = fitsutils.get_par(ifu_img, "JD") - fitsutils.get_par(ifu_img, "EXPTIME")/ (24*3600.)
            jd_end = fitsutils.get_par(ifu_img, "JD")
        else:
            jd_ini = fitsutils.get_par(ifu_img, "JD") 
            jd_end = fitsutils.get_par(ifu_img, "JD") + fitsutils.get_par(ifu_img, "EXPTIME")/ (24*3600.)
            
        #We only fill the dictionary if the exposure is a valid science or standard image.
        name = fitsutils.get_par(ifu_img, "OBJECT")
        ra = fitsutils.get_par(ifu_img, "RA")
        dec = fitsutils.get_par(ifu_img, "DEC")
        rad, decd = cc.hour2deg(ra, dec)
        exptime = fitsutils.get_par(ifu_img, "EXPTIME")
        ifu_dic[ifu_img] = (name, jd_ini, jd_end, rad, decd, exptime)
    else:
        print ("Image %s is not SCIENCE or STANDARD."%ifu_img)
    
def __combine_guiders(ifu_dic, abspath, outdir):
    '''
    Receives an IFU dictionary with the images that need a solved guider.
    '''        
    rc = np.array(glob.glob(abspath+"/rc*fits"))
        
    rcjd = np.array([fitsutils.get_par(r, "JD") for r in rc])
    imtypes = np.array([fitsutils.get_par(r, "IMGTYPE").upper() for r in rc])
    objnames = np.array([fitsutils.get_par(r, "OBJECT").upper() for r in rc])
    ras = np.array([cc.getDegRaString( fitsutils.get_par(r, "RA")) for r in rc])  
    decs = np.array([ cc.getDegDecString( fitsutils.get_par(r, "DEC")) for r in rc])
    
    for ifu_i in ifu_dic.keys():
        name, jd_ini, jd_end, rad, decd, exptime = ifu_dic[ifu_i]
        #guiders = rc[(imtypes=="GUIDER") * (rcjd >= ifu_dic[ifu_i][1]) * (rcjd <= ifu_dic[ifu_i][2]) ]
        mymask = (rcjd >= jd_ini) * (rcjd <= jd_end) *\
            (np.abs(ras - rad)*np.cos(np.deg2rad(decd))<0.5/60 ) * (np.abs(decs - decd)<0.5/60 )
        guiders = rc[mymask]
        im = imtypes[mymask]
        names = objnames[mymask]
        print ( "For image %s on object %s with exptime %d found guiders:\n %s"%(ifu_i, name, exptime, zip(names, im)))
        out = os.path.join(outdir, "guider_" + os.path.basename(ifu_i))
        create_masterguide(guiders, out=out)
        fitsutils.update_par(out, "IFU_IMG", os.path.basename(ifu_i))
        
def make_guider(ifu_img, outdir):
    '''
    Creates the stacked and astrometry solved RC image for the IFU image passed as a parameter.
    
    ifu_img: string which is the path to the IFU image we want to get the guiders for.
    '''
    ifu_dic = {}
    myabspath = os.path.dirname(os.path.abspath(ifu_img))
    __fill_ifu_dic(ifu_img, ifu_dic)
    __combine_guiders(ifu_dic, myabspath, outdir)
    
    
def make_guiders(ifu_dir, outdir):
    '''
    Looks for all the IFU images in the directory and assembles all the guider images taken with RC for that time interval
    at around the coordinates of the IFU.
    
    ifu_dir: directory where the IFU images are that we want to greated the guiders for.
    '''
    
    abspath = os.path.abspath(ifu_dir)
        
    ifu = np.array(glob.glob(abspath+"/ifu*fits"))
    
    ifu_dic = {}
    for i in ifu:
        __fill_ifu_dic(i, ifu_dic)
        
    __combine_guiders(ifu_dic, abspath, outdir)
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Stacks the photometric rc images from SEDM Rainbow Camera.
        Requires the name of the IFU images to be reduced the guiders for.
        
        %run solve_guider.py -f ifuimage -d OUTDIR 
                
        Reduced images are stored in directory called "OUTDIR".

            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-f', '--ifufile', type=str, help='File containing the ifu image.',default=None)
    parser.add_argument('-d', '--outdir', type=str, help='Directory containing the output stacked guider.', default=None)

    args = parser.parse_args()
    
    ifufile = args.ifufile
    outdir = args.outdir
    
    make_guider(ifufile, outdir)
