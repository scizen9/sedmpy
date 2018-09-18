# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 13:11:52 2015

@author: nadiablago
"""

import matplotlib
matplotlib.use('Agg')

import fitsutils
import os, glob, shutil, sys
import numpy as np
try:
    from pyraf import iraf 
    from pyraf import IrafError
except:
    pass
from astropy.io import fits
from astropy.wcs import WCS

from matplotlib import pylab as plt
import subprocess
import argparse
import time
import datetime
import logging
import sextractor
import zeropoint
from astropy.io import fits
import coordinates_conversor as cc

from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')
_db = parser.get('persistence', 'db')


FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
now = datetime.datetime.utcnow()
timestamp=datetime.datetime.isoformat(now)
creationdate = timestamp
timestamp=timestamp.split("T")[0]

try:
    #Log into a file
    root_dir = _logpath
    logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "rcred_{0}.log".format(timestamp)), level=logging.INFO)
    logger = logging.getLogger('rcred')
except:
    logging.basicConfig(format=FORMAT, filename=os.path.join("/tmp", "rcred_{0}.log".format(timestamp)), level=logging.INFO)
    logger= logging.getLogger("rcred")
    


    
def get_xy_coords(image, ra, dec):
    '''
    Uses the wcs-rd2xy routine to compute the proper pixel number where the target is.
    Sometime the wcs module does not seem to be providing the correct answer, as it does not seem
    to be using the SIP extension.
    
    '''
    import re
    import subprocess
    cmd = "wcs-rd2xy -w %s -r %.5f -d %.5f"%(image, ra, dec)
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output = proc.stdout.read()    
    output = output.split("->")[1]
    
    coords = []
    for s in output.split(","):
        coords.append(float(re.findall("[-+]?\d+[\.]?\d*", s)[0]))
        
    return coords

def get_rd_coords(image, x, y):
    '''
    Uses the wcs-xy2rd routine to compute the ra, dec where the target is for a given pixel.
    Sometime the wcs module does not seem to be providing the correct answer, as it does not seem
    to be using the SIP extension.
    
    '''
    import re
    import subprocess
    cmd = "wcs-xy2rd -w %s -x %.5f -y %.5f"%(image, x, y)
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output = proc.stdout.read()    
    output = output.split("->")[1]
    
    numbers = re.findall("[-+]?\d+[\.]?\d*", output)
    coords = []
    for n in numbers:
        coords.append(float(n))
        
    return coords

    
def create_masterbias(biasdir=None, channel='rc'):
    '''
    Combines slow and fast readout mode biases for the specified channel.
    '''
    
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    if (biasdir is None) or biasdir=="": biasdir = "."
        
    outs = "Bias_%s_slow.fits"%channel
    outf = "Bias_%s_fast.fits"%channel

    doslow = True
    dofast = True
    if (os.path.isfile(os.path.join(biasdir,outs))): 
        logger.warn( "%s master Bias exists!"%outs)
        doslow = False
    if ( os.path.isfile(os.path.join(biasdir,outf))):
        logger.warn("%s master Bias exists!"%outs)
        dofast = False
    if(doslow or dofast):
        logger.info("Starting the Master Bias creation!")
    else:
        return

    os.chdir(biasdir)        
        
    lfastbias = []
    lslowbias = []
    
    #Select all filts that are Bias with same instrument
    for f in glob.glob("rc*fits"):
        try:
            if ( "BIAS" in str.upper(fitsutils.get_par(f, "IMGTYPE").upper()) ):
                if (fitsutils.get_par(f, "ADCSPEED")==2):
                    lfastbias.append(f)
                else:
                    lslowbias.append(f)
        except:
            pass
                
    logger.info("Files for bias SLOW mode: %s"% lslowbias)
    logger.info( "Files for bias FAST mode: %s"% lfastbias)
    
    if len(lfastbias) > 0 and dofast:
        bfile_fast ="lbias_fast_"+channel
        np.savetxt(bfile_fast, np.array(lfastbias), fmt="%s")
        if (os.path.isfile("Bias_stats_fast")): 
            os.remove("Bias_stats_fast")
        iraf.imstat("@"+bfile_fast, Stdout="Bias_stats_fast")
        
        st = np.genfromtxt("Bias_stats_fast", names=True, dtype=None)
        logger.info("%s"%st)
        
        iraf.imcombine(input = "@"+bfile_fast, \
                    output = outf, \
                    combine = "median",\
                    scale = "mode")
        os.remove(bfile_fast)
        
        #copy into the reference folder with current date
        newdir = os.path.join("../../refphot/", os.path.basename(os.path.abspath(biasdir)))
        if (not os.path.isdir(newdir)):
            os.makedirs(newdir)
            shutil.copy(outf, os.path.join(newdir, os.path.basename(outf)) )
        else:
            copy_ref_calib(biasdir, "Bias")


    if len(lslowbias) > 0 and doslow:

        bfile_slow ="lbias_slow_"+channel
        np.savetxt(bfile_slow, np.array(lslowbias), fmt="%s")
        if (os.path.isfile("Bias_stats_slow")): os.remove("Bias_stats_slow")
        iraf.imstat("@"+bfile_slow, Stdout="Bias_stats_slow")
        
        st = np.genfromtxt("Bias_stats_slow", names=True, dtype=None)
        logger.info("%s"%st)
        
        iraf.imcombine(input = "@"+bfile_slow, \
                    output = outs, \
                    combine = "median",\
                    scale = "mode")
        os.remove(bfile_slow)
        
        #copy into the reference folder with current date
        newdir = os.path.join("../../refphot/", os.path.basename(os.path.abspath(biasdir)))
        if (not os.path.isdir(newdir)):
            os.makedirs(newdir)
        shutil.copy(outs, os.path.join(newdir, os.path.basename(outs)) )  
    else:
        copy_ref_calib(biasdir, outs)

def create_masterflat(flatdir=None, biasdir=None, channel='rc', plot=True):
    '''
    Creates a masterflat from both dome flats and sky flats if the number of counts in the given filter
    is not saturated and not too low (between 3000 and 40000). 
    '''
    
    
    if (flatdir == None or flatdir==""): flatdir = "."
        
    if (biasdir == None or biasdir==""): biasdir = flatdir
        
    os.chdir(flatdir)
    
    if (plot and not os.path.isdir("reduced/flats")):
        os.makedirs("reduced/flats")
    
    if (len(glob.glob("Flat_%s*norm.fits"%channel)) == 4):
        logger.info( "Master Flat exists!")
        return 
    if (len(glob.glob("Flat_%s*norm.fits"%channel)) > 0):
        logger.info( "Some Master Flat exist!")
    else:
        logger.info( "Starting the Master Flat creation!")

    bias_slow = "Bias_%s_slow.fits"%channel
    bias_fast = "Bias_%s_fast.fits"%channel
    
    if (not os.path.isfile(bias_slow) and not os.path.isfile(bias_fast) ):
        create_masterbias(biasdir)
     
    lsflat = []
    lfflat = []
    
    obj = ""
    imtype = ""
    
    #Select all filts that are Flats with same instrument
    for f in glob.glob(channel+"*fits"):
        try:
            if fitsutils.has_par(f, "OBJECT"):
                obj = str.upper(fitsutils.get_par(f, "OBJECT"))
            else:
                continue

            if fitsutils.has_par(f, "IMGTYPE"):
                imtype = str.upper(fitsutils.get_par(f, "IMGTYPE"))
            else:
                continue
        
            #if ("RAINBOW CAM" in str.upper(fitsutils.get_par(f, "CAM_NAME")) and  ("DOME" in  obj or "FLAT" in obj or "Twilight" in obj or "TWILIGHT" in imtype or "DOME" in imtype)):
            if ( "twilight" in imtype.lower()):

                if (fitsutils.get_par(f, "ADCSPEED")==2):
                    lfflat.append(f)
                else:
                    lsflat.append(f)
        except:
            logger.error( "Error with retrieving parameters for file %s"% f)
            pass
                
    logger.info( "Files for slow flat %s"% lsflat)
    logger.info( "Files for fast flat %s"% lfflat)
    
    fsfile ="lflat_slow_"+channel
    np.savetxt(fsfile, np.array(lsflat), fmt="%s")
    fffile ="lflat_fast_"+channel
    np.savetxt(fffile, np.array(lfflat), fmt="%s")


    # Running IRAF
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    #Remove bias from the flat
    if len(lsflat) >0:
        iraf.imarith("@"+fsfile, "-", bias_slow, "b_@"+fsfile)
    
    if len(lfflat) >0:
        iraf.imarith("@"+fffile, "-", bias_fast, "b_@"+fffile)    
    
    #Remove the list files
    os.remove(fsfile)
    os.remove(fffile)
    
    #Slices the flats.
    debiased_flats = glob.glob("b_*.fits")
    for f in debiased_flats:
        logger.info( "Slicing file %s"% f)
        try:
            slice_rc(f)
        except:
            logger.error("Error when slicing file... deleting the unsliced one...")
        #Remove the un-sliced file
        os.remove(f)
        
    #Selects the ones that are suitable given the number of counts and combines them.
    bands = ['u', 'g', 'r', 'i']
    for b in bands:
        out = "Flat_%s_%s.fits"%(channel, b)
        out_norm = out.replace(".fits","_norm.fits")

        if (os.path.isfile(out_norm)):
            logger.error( "Master Flat for filter %s exists. Skipping..."%b)
            continue
        
        lfiles = []
        for f in glob.glob('b_*_%s.fits'%b):
            fi = fits.open(f)
	    d = fi[0].data
            status = "rejected" 
            if np.percentile(d, 90)>4000 and np.percentile(d, 90)<45000:
                lfiles.append(f)
		mymode = 1. * np.median(d.flatten())
		d[d>45000] = mymode
		fi[0].data = d
		fi.writeto(f, clobber=True)
		status = "accepted"

            if (plot):
                plt.title("Flat filter %s. %s"%(b, status))
                plt.imshow(d.T, cmap=plt.get_cmap("nipy_spectral"))
		plt.colorbar()
                plt.savefig("reduced/flats/%s"%(f.replace(".fits", ".png")))
                plt.close()
            #Make sure that the optimum number of counts is not too low and not saturated.
        if len(lfiles) == 0:
            logger.error( "WARNING!!! Could not find suitable flats for band %s"%b)
            continue
        if len(lfiles) < 3:
            logger.error( "WARNING!!! Could find less than 3 flats for band %s. Skipping, as it is not reliable..."%b)
            continue
        ffile ="lflat_"+b
        np.savetxt(ffile, np.array(lfiles), fmt="%s")
    
        
        #Cleaning of old files
        if(os.path.isfile(out)): os.remove(out)
        if(os.path.isfile(out_norm)): os.remove(out_norm)
        if(os.path.isfile("Flat_stats")): os.remove("Flat_stats")
        
        
        #Combine flats
        iraf.imcombine(input = "@"+ffile, \
                        output = out, \
                        combine = "median",\
                        scale = "mode",
                        reject = "sigclip", lsigma = 2., hsigma = 2, gain=1.7, rdnoise=4.)
        iraf.imstat(out, fields="image,npix,mean,stddev,min,max,mode", Stdout="Flat_stats")
        st = np.genfromtxt("Flat_stats", names=True, dtype=None)
        #Normalize flats
        iraf.imarith(out, "/", st["MODE"], out_norm)
        
        #Do some cleaning
        logger.info( 'Removing from lfiles')
        for f in glob.glob('b_*_%s.fits'%b):
            os.remove(f)

        os.remove(ffile)
        
        
        if os.path.isfile(fsfile):
            os.remove(fsfile)
        if os.path.isfile(fffile):
            os.remove(fffile)

        #copy into the reference folder with current date
        newdir = os.path.join("../../refphot/", os.path.basename(os.path.abspath(flatdir)))
        if (not os.path.isdir(newdir)):
            os.makedirs(newdir)
        shutil.copy(out_norm, os.path.join(newdir, os.path.basename(out_norm)) )   
    copy_ref_calib(flatdir, "Flat")	

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
        except IrafError:
            iraf.imarith("@"+fffile, "-", bias_fast, "@"+bffile)
    else:
        bffile = fffile

        
    #Combine flats
    iraf.imcombine(input = "@"+bffile, \
                    output = out, \
                    combine = "median",\
                    scale = "mode",
                    reject = "sigclip", lsigma = 2., hsigma = 2, gain=1.7, rdnoise=4.)
    iraf.imstat(out, fields="image,npix,mean,stddev,min,max,mode", Stdout="guide_stats")
    #st = np.genfromtxt("guide_stats", names=True, dtype=None)
    
    #Do some cleaning
    if debias:
        logger.info( 'Removing from lfiles')
        for f in debiased:
            if os.path.isfile(f): os.remove(f)

    if os.path.isfile(fffile):
        os.remove(fffile)
    if os.path.isfile(bffile):
        os.remove(bffile)
        
    solve_astrometry(out, overwrite=False)
    
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
        logger.warn("Image %s is not SCIENCE or STANDARD."%ifu_img)
    
def __combine_guiders(ifu_dic, abspath):
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
            (np.abs(ras - rad)*np.cos(np.deg2rad(decd))<1./60 ) * (np.abs(decs - decd)<1./60 )
        guiders = rc[mymask]
        im = imtypes[mymask]
        names = objnames[mymask]
        logger.info( "For image %s on object %s with exptime %d found guiders:\n %s"%(ifu_i, name, exptime, zip(names, im)))
        out = os.path.join(os.path.dirname(ifu_i), "guider_" + os.path.basename(ifu_i))
        create_masterguide(guiders, out=out)
        fitsutils.update_par(out, "IFU_IMG", os.path.basename(ifu_i))
        
def make_guider(ifu_img):
    '''
    Creates the stacked and astrometry solved RC image for the IFU image passed as a parameter.
    
    ifu_img: string which is the path to the IFU image we want to get the guiders for.
    '''
    ifu_dic = {}
    myabspath = os.path.dirname(os.path.abspath(ifu_img))
    __fill_ifu_dic(ifu_img, ifu_dic)
    __combine_guiders(ifu_dic, myabspath)
    
    
def make_guiders(ifu_dir):
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
        
    __combine_guiders(ifu_dic, abspath)
            
    
    
def mask_stars(image, sexfile, plot=False, overwrite=False):
    ''' 
    Finds the stars in the sextrated file and creates a mask file.
    '''
       
    maskdir = os.path.join(os.path.abspath(os.path.dirname(image)), "masks")

    if not os.path.isdir(maskdir):
        os.makedirs(maskdir)
        
    maskname = os.path.join(maskdir, os.path.basename(image).replace(".fits", ".im.fits"))

    if (os.path.isfile(maskname) and not overwrite):
        return maskname

    print "Creating mask %s"%maskname
    
    hdulist = fits.open(image)
    header = hdulist[0].header
    data = np.ones_like(hdulist[0].data)
    
    stars = np.genfromtxt(sexfile)
    fwhm = stars[:,7]
    mag = stars[:,4]
    flags = np.array(stars[:, 10], dtype=np.int)

    starmask = ( mag < np.percentile(mag, 90)) | (np.bitwise_and(flags, np.repeat(0x004, len(flags)) )>1)
    stars = stars[ starmask]

    fwhm = stars[:,7]

    x = stars[:,0]
    y = stars[:,1]
    
    lenx = data.shape[0]
    leny = data.shape[1]
    X, Y = np.meshgrid( np.arange(leny), np.arange(lenx))
    
    for i in range(len(stars)):
        data[ np.sqrt((X-x[i])**2 + (Y-y[i])**2)< fwhm[i]*5] = 0
            
            
    if (plot):
        print image
        plt.imshow(np.log10(np.abs(hdulist[0].data)), alpha=0.9)
        plt.imshow(data, alpha=0.5)
        plt.show()

    hdu = fits.PrimaryHDU(data)
    hdu.header = header
    newhdulist = fits.HDUList([hdu])
    newhdulist.writeto(maskname)

    return maskname
    
def create_superflat(imdir, filters=["u", "g", "r", "i"]):
    #Locate images for each filter
    imlist = glob.glob("rc*fits")
    
  
    #Run sextractor to locate bright sources
  
    sexfiles = sextractor.run_sex(imlist, overwrite=False)
    maskfiles = []
    
    for i, im in enumerate(imlist): 
        #Create a mask and store it int he mask directory
        maskfile = mask_stars(im, sexfiles[i])        
        maskfiles.append(maskfile)
        fitsutils.update_par(im, "BPM", os.path.relpath(maskfile))
        
    
        
    for filt in filters:
        fimlist = [im for im in imlist if fitsutils.get_par(im, "FILTER") == filt]
        fmasklist = [im for im in maskfiles if fitsutils.get_par(im, "FILTER") == filt]
        
        if len(fimlist) == 0:
            continue
        
        fsfile ="lflat_%s"%filt
        msfile = "lmask_%s"%filt
        np.savetxt(fsfile, np.array(fimlist), fmt="%s")
        np.savetxt(msfile, np.array(fmasklist), fmt="%s")
        
        
        '''masklist = []
        
        for m in fmasklist:
            hdulist = fits.open(m)
            data = hdulist[0].data
            masklist.append(data)
            
            
        masklist = np.array(masklist)
        
        hdu = fits.PrimaryHDU(masklist)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto("mastermask_%s.fits"%filt)'''           
                
        # Running IRAF
        iraf.noao(_doprint=0)
        iraf.imred(_doprint=0)
        iraf.ccdred(_doprint=0)
        
        iraf.imarith("@"+fsfile, "*", "@"+msfile, "m_@"+fsfile)


        #Combine flats
        iraf.imcombine(input = "m_@"+fsfile, \
                        output = "superflat_%s.fits"%filt, \
                        combine = "median",\
                        scale = "mode", \
                        masktype="badvalue",\
                        maskvalue = 0)
                        
        iraf.imstat("superflat_%s.fits"%filt, fields="image,npix,mean,stddev,min,max,mode", Stdout="Flat_stats")
        time.sleep(0.1)
        st = np.genfromtxt("Flat_stats", names=True, dtype=None)
        #Normalize flats
        iraf.imarith("superflat_%s.fits"%filt, "/", st["MODE"], "superflat_%s_norm.fits"%filt)
                        

    
def get_median_bkg(img):
    '''
    Computes the median background.
    '''
    hdu = fits.open(img)
    header = hdu[0].header
    bkg = np.median(hdu[0].data[hdu[0].data > 0])
    return bkg
    
def solve_edges(flat, band):
    '''
    Finds the edges of the flat, and extends the respons of the average to the edges.
    '''
    
    edgedic = {"u":[], "g":[], "i":[], "r":[]}
    
    
    
def copy_ref_calib(curdir, calib="Flat"):
    '''
    Reference master Bias and master Flat are stored in the refphot folder.
    The files are copied if they are not found in the folder where the photometry is being reduced.
    '''

    if calib == "Bias": 
   	calib_dic = { 	"Bias_rc_fast.fits" : False,
			"Bias_rc_slow.fits" : False}
    else:
    	calib_dic = { 	"Flat_rc_u_norm.fits" : False,
			"Flat_rc_g_norm.fits" : False,
			"Flat_rc_r_norm.fits" : False,
			"Flat_rc_i_norm.fits" : False}

    #Get the date of the current directory
    curdir = os.path.abspath(curdir)

    #Check if we have all the calibrations we need.
    for cal in calib_dic.keys():
	calib_dic[cal] = os.path.isfile(os.path.join(curdir, cal))

    print calib_dic

    #If all the calibration files are in place, nothing to do.
    if np.all(calib_dic.values()):
	return

    #Check which dates have the calibration files that we need.
    listcalib = glob.glob(os.path.join(curdir, "../../refphot/*/", calib+"*"))
    #Obtain the folder name
    listdirs = [os.path.dirname(l) for l in listcalib]    
    #Unique name
    calibdates = np.array(list(set(listdirs)))
    

    #compile the dates in datetime format
    dates = []
    for f in calibdates:
        dateobs = os.path.basename(f)
        dates.append(datetime.datetime.strptime(dateobs, "%Y%m%d"))

    dates = np.array(dates)    
    curdate = datetime.datetime.strptime(os.path.basename(curdir), "%Y%m%d")
    
    #Select the folder that is closer to the date of the current directory
    lastdir = np.array(calibdates)[np.argmin(np.abs(curdate-dates))]
   
    #Try to go 100 days before the data was taken.
    i = 0
    while i < 1000 and not np.all(calib_dic.values()):
	i = i+1
	newdate = curdate - datetime.timedelta(i)
	newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
	print "Checking day %s"%newdatedir
	for cal in np.array(calib_dic.keys())[~np.array(calib_dic.values())]:
		c = os.path.join(curdir, "../../refphot/", newdatedir, cal)
		#print "testing %s"%c
		if os.path.isfile(c):
			print "Copying calibration file %s to directory %s"%(c, curdir)
        		shutil.copy(c, os.path.join(curdir, os.path.basename(c)))     
			calib_dic[cal] = True
			

        

def solve_astrometry(img, outimage=None, radius=3, with_pix=True, overwrite=False, tweak=3):
    '''
    img: fits image where astrometry should be solved.
    outimage: name for the astrometry solved image. If none is provided, the name will be "a_"img.
    radius: radius of uncertainty on astrometric position in image.
    with_pix: if we want to include the constraint on the pixel size for the RCCam.
    overwrite: wether the astrometrically solved image should go on top of the old one.
    tewak: parameter for astrometry.net
    '''

    from astropy.wcs import InconsistentAxisTypesError
    
    img = os.path.abspath(img)
    
    ra = fitsutils.get_par(img, 'RA')
    dec = fitsutils.get_par(img, 'DEC')
    #logger.info( "Solving astrometry on field with (ra,dec)=%s %s"%(ra, dec))
    
    astro = os.path.join( os.path.dirname(img), "a_" + os.path.basename(img))


    #If astrometry exists, we don't run it again.
    if (os.path.isfile(astro) and not overwrite):
        return astro
        
    #Store a temporary file with the multiplication of the mask
    mask = "mask.fits"

    cmd = "solve-field --ra %s --dec %s --radius %.4f -p --new-fits %s \
      -W none -B none -P none -M none -R none -S none -t %d --overwrite %s "%(ra, dec, radius, astro, tweak, img)
    if (with_pix):
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.4"
    #logger.info( cmd)

    subprocess.call(cmd, shell=True)
    
    #Cleaning after astrometry.net
    if (os.path.isfile(img.replace(".fits", ".axy"))):
        os.remove(img.replace(".fits", ".axy"))
    if (os.path.isfile(img.replace(".fits", "-indx.xyls"))):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if (os.path.isfile("none")):
        os.remove("none")
        

    if os.path.isfile(astro):
        try:
		is_on_target(img)
        except InconsistentAxisTypesError as e:
		fitsutils.update_par(img, "ONTARGET", 0)
		print "Error detected with WCS when reading file %s. \n %s"%(img, e)
        
    if (not outimage is None and overwrite and os.path.isfile(astro)):
        shutil.move(astro, outimage)
	return outimage
    elif (outimage is None and overwrite and os.path.isfile(astro)):
        shutil.move(astro, img)
	return img
    else:
    	return astro



    

def make_mask_cross(img):  
    
    maskname = os.path.join(os.path.dirname(img), "mask.fits")
    
    if (os.path.isfile(maskname)):
        return maskname
        
    f = fits.open(img)
    data = f[0].data
    
    corners = {
    "g" : [1, 850, 1, 850],
    "i" : [1, 850, 110, 2045],
    "r" : [1150, 2045, 1150, 2045],
    "u" : [1130, 2045, 1, 850]
    }
    
    newdata = np.ones_like(data)*np.nan
    
    for band in corners.keys():
        ix = corners[band]
        newdata[ix[0]:ix[1], ix[2]: ix[3]] = 1
        
    f[0].data = newdata
    f.writeto(maskname, clobber=True)
    
    return maskname
    
def get_masked_image(img):

    # Running IRAF
    iraf.noao(_doprint=0)  
    mask = make_mask_cross(img)    
    masked = img.replace(".fits", "_masked.fits")
    
    if (os.path.isfile(masked)):
        os.remove(masked)
    iraf.imarith(img, "*", mask, masked)
    
    return masked

    
def slice_rc(img):
    '''
    Slices the Rainbow Camera into 4 different images and adds the 'filter' keyword in the fits file.
    '''
    fname = os.path.basename(img)
    fdir = os.path.dirname(img)    
    
    # Running IRAF
    iraf.noao(_doprint=0)    
    
    corners = {
    "g" : [1, 910, 1, 900],
    "i" : [1, 910, 1060, 2045],
    "r" : [1040, 2045, 1015, 2045],
    "u" : [1030, 2045, 1, 900]
    }
    
    
    filenames = []
        
    for i, b in enumerate(corners.keys()):
        logger.info( "Slicing for filter %s"% b)
        name = fname.replace(".fits", "_%s.fits"%b)
        
        #Clean first
        if (os.path.isfile(name)):
            os.remove(name)
            
        iraf.imcopy("%s[%d:%d,%d:%d]"%(img, corners[b][0], corners[b][1], corners[b][2], corners[b][3]), name)
         
        fitsutils.update_par(name, 'filter', b)
        is_on_target(name)
        
        filenames.append(name)
    
    return filenames

def is_on_target(image):
    '''
    Add as a parameter whether the image is on target or not.
    
    '''
    import coordinates_conversor as cc
    
    
    ra, dec = cc.hour2deg(fitsutils.get_par(image, 'OBJRA'), fitsutils.get_par(image, 'OBJDEC'))

    impf = fits.open(image)
    w = WCS(impf[0].header)
    
    filt = fitsutils.get_par(image, "FILTER")
    pra, pdec = get_xy_coords(image, ra,dec)

    shape = impf[0].data.shape
    
    if (pra > 0)  and (pra < shape[0]) and (pdec > 0) and (pdec < shape[1]):
        fitsutils.update_par(image, "ONTARGET", 1)
        return True
    else:
        fitsutils.update_par(image, "ONTARGET", 0)
        return False
        
    
def clean_cosmic(f):
    '''
    From lacosmic.
    '''
    import cosmics
    
    out = f.replace('.fits',  '_clean.fits')
    
    #If it does already exist, just return the name.
    if (os.path.isfile(out)):
        return out
    
    #Otherwise, run the cosmic ray rejection based on LA Cosmic.
    g = fitsutils.get_par(f, "GAIN")
    if (fitsutils.has_par(f, "RDNOISE")):
        rn = fitsutils.get_par(f, "RDNOISE")
    else:
        rn = 20
    array, header = cosmics.fromfits(f)
    
    try:
        c = cosmics.cosmicsimage(array, gain=g, readnoise=rn, sigclip = 8.0, sigfrac = 0.3, satlevel = 64000.0)
        c.run(maxiter = 3)
        out = f.replace('.fits',  '_clean.fits')
    
        cosmics.tofits(out, c.cleanarray, header)
        fitsutils.update_par(out, "CRREJ", 1)

        #os.remove(f)
    except:
        pass
    
    return out
    

    
                
def get_overscan_bias_rc(img):
    '''
    Bias from overscan region.
    '''
    f = fits.open(img)
    bias = np.nanmedian(f[0].data[990-100:990+100,970-100:970+100].flatten())
    
    return bias
    
def get_sequential_name(target_dir, name, i=0):
        '''
        Gets a sequential name if we have a file with the same object name imaged several times.
        '''
        newname = os.path.join(target_dir, name).replace(".fits", "_%d.fits"%i)
        #If the destination file does not exist...
        if (os.path.isfile(name) and not os.path.isfile(newname)):
            return newname
        #If it does exist, but it is another exposure
        elif(os.path.isfile(name) and os.path.isfile(newname)):
            #If it is not the same exposure, add with different name. Otherwise, replace.
            if fitsutils.get_par(name, "JD") != fitsutils.get_par(newname, "JD"):
                newname = get_sequential_name(target_dir, name, i=i+1)
        
        return newname      

def init_header_reduced(image):
    '''
    IQWCS = 1 or 0 / Indicates astrometry has been solved for the field
    IQZEROPT =  1 or 0 /indicates if the zero point was calculated for the image    
    SKYBKG = FLOAT / Average sky background given in counts    
    SEEPIX = FLOAT  / Seeing expressed in pixels    
    ZPCAT = 'String' / Catalog used to calculate zero point    
    ZEROPTU = FLOAT  / Zero point uncertainty    
    ZEROPT  = 'FLOAT'  / Zero point by comparison to catalog
    '''
    
    pardic = {"IQWCS" : 0,\
                "IQZEROPT" : 0,\
                "SKYBKG": 0,\
                "SEEPIX": 0,\
                "ZPCAT" : "none",\
                "ZEROPTU" : 0.,\
                "ZEROPT" : 0.,\
                "CRREJ" : 0}
    fitsutils.update_pars(image, pardic)
    
def plot_image(image):
    '''
    Plots the reduced image into the png folder.

    '''
    logger.info("Plotting image %s"% image)    

    print "Plotting image ", image    

    
    image = os.path.abspath(image)

    
    #Change to image directory
    imdir, imname = os.path.split(image)

    #Create destination directory

    png_dir = os.path.join(imdir, "png")

    if (not os.path.isdir(png_dir)):
        os.makedirs(png_dir)

    f = fits.open(image)[0]
    d = f.data
    h = f.header
    exptime = h.get('EXPTIME', 0)
    name = h.get('OBJECT', 'None')
    filt = h.get('FILTER', 'NA')
    ontarget = h.get('ONTARGET', False)

    plt.imshow(d, origin="lower", vmin=np.percentile(d.flatten(), 5), vmax=np.percentile(d, 95), cmap=matplotlib.cm.cubehelix)
    plt.title("%s %s-band [%ds]. On target=%s"%(name, filt, exptime, ontarget))
    plt.savefig( os.path.join(png_dir, imname.replace(".fits", ".png")))
    plt.close()


def reduce_image(image, flatdir=None, biasdir=None, cosmic=False, astrometry=True, channel='rc', target_dir='reduced', overwrite=False):
    '''
    Applies Flat field and bias calibrations to the image.
    
    Steps:
    
    1. - Solve astrometry on the entire image.
    2. - Computes cosmic ray rejectionon the entire image.
    3. - Compute master bias (if it does not exist) and de-bias the image.
    4. - Separate the image into 4 filters.
    5. - Compute flat field for each filter (if it does not exist) and apply flat fielding on the image.
    6. - Compute the image zeropoint.

    '''
    
    logger.info("Reducing image %s"% image)    

    print "Reducing image ", image    

    
    image = os.path.abspath(image)
    imname = os.path.basename(image).replace(".fits", "")
    try:
        objectname = fitsutils.get_par(image, "NAME").replace(" ","")+"_"+fitsutils.get_par(image, "FILTER")
    except:
        logger.error( "ERROR, image "+ image + " does not have a NAME or a FILTER!!!")
        return

    print "For object", objectname
    logger.info( "For object %s"% objectname)


    
    
    #Change to image directory
    mydir = os.path.dirname(image)
    if mydir=="": mydir = "."
    mydir = os.path.abspath(mydir)
    os.chdir(mydir)
    #Create destination directory
    if (not os.path.isdir(target_dir)):
        os.makedirs(target_dir)

    #If we don't want to overwrite the already extracted images, we check wether they exist.
    if (not overwrite):
        existing = True
        for band in ['u', 'g', 'r', 'i']:
            destfile = os.path.join(target_dir, imname + "_f_b_a_%s_%s_0.fits"%(objectname, band))
            logger.info( "Looking if file %s exists: %s"%( destfile, \
                (os.path.isfile(destfile) ) ) )
            existing = existing and (os.path.isfile( destfile ) )
        if existing:
            return []


        
    #Initialize the basic parameters.
    init_header_reduced(image)

    astro = ""
    if (astrometry):
        logger.info( "Solving astometry for the whole image...")
        img = solve_astrometry(image)
        if (os.path.isfile(img)):
            astro="a_"
            fitsutils.update_par(img, "IQWCS", 1)
        else:
            logger.error( "ASTROMETRY DID NOT SOLVE ON IMAGE %s"% image)
            img = image

    #Update noise parameters needed for cosmic reection
    if (fitsutils.get_par(img, "ADCSPEED")==2):
        fitsutils.update_par(img, "RDNOISE", 20.)
    else:
        fitsutils.update_par(img, "RDNOISE", 4.)

    if (cosmic):
        logger.info( "Correcting for cosmic rays...")
        # Correct for cosmics each filter
        cleanimg = clean_cosmic(os.path.join(os.path.abspath(mydir), img))
        img = cleanimg
    
    
    
    #Compute BIAS
    if (biasdir is None or biasdir==""): biasdir = "."
    create_masterbias(biasdir)
    
    bias_slow = os.path.join(biasdir, "Bias_%s_%s.fits"%(channel, 'slow'))
    bias_fast = os.path.join(biasdir, "Bias_%s_%s.fits"%(channel, 'fast'))
    
    # Running IRAF to DE-BIAS
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    #Compute flat field
    if (flatdir is None or flatdir==""): flatdir = "."
    create_masterflat(flatdir, biasdir)
    
    #New names for the object.
    debiased = os.path.join(os.path.dirname(img), "b_" + os.path.basename(img))
    logger.info( "Creating debiased file, %s"%debiased)
    
    if ( (fitsutils.get_par(img, "ADCSPEED")==0.1 and not os.path.isfile(bias_slow)) \
        or (fitsutils.get_par(img, "ADCSPEED")==2 and not os.path.isfile(bias_fast)) ):
        logger.warn( "Master bias not found! Tryting to copy from reference folder...")
        copy_ref_calib(mydir, "Bias")
        if ( (fitsutils.get_par(img, "ADCSPEED")==0.1 and not os.path.isfile(bias_slow)) \
        or (fitsutils.get_par(img, "ADCSPEED")==2 and not os.path.isfile(bias_fast)) ):
            logger.error( "Bias not found in reference folder")
            return


    #Clean first
    if (os.path.isfile(debiased)):
        os.remove(debiased)
        
    #Debias
    if (fitsutils.get_par(img, "ADCSPEED")==2):
        iraf.imarith(img, "-", bias_fast, debiased)
        fitsutils.update_par(debiased, "BIASFILE", bias_fast)
        fitsutils.update_par(debiased, "RDNOISE", 20.)

    else:
        iraf.imarith(img, "-", bias_slow, debiased)
        fitsutils.update_par(debiased, "BIASFILE", bias_slow)
        fitsutils.update_par(debiased, "RDNOISE", 4.)

    #Set negative counts to zero
    hdu = fits.open(debiased)
    header = hdu[0].header
    hdu[0].data[hdu[0].data<0] = 0
    hdu.writeto(debiased, clobber=True)

    #Slicing the image for flats  
    slice_names = slice_rc(debiased)
    print "Creating sliced files, ", slice_names

    
    #Remove un-sliced image
    os.remove(debiased)

    # DE-flat each filter and store under object name
    for i, debiased_f in enumerate(slice_names):
        b = fitsutils.get_par(debiased_f, 'filter')
        
        deflatted = os.path.join(os.path.dirname(image), target_dir, imname + "_f_b_" + astro + objectname + "_%s.fits"%b)

        #Flat to be used for that filter
        flat = os.path.join(flatdir, "Flat_%s_%s_norm.fits"%(channel, b))

        if (not os.path.isfile(flat)):
            logger.warn( "Master flat not found in %s"% flat)
            copy_ref_calib(mydir, "Flat")
            continue
        else:
            logger.info( "Using flat %s"%flat)
            
        #Cleans the deflatted file if exists
        if (os.path.isfile(deflatted)):
            os.remove(deflatted)

        if (os.path.isfile(debiased_f) and os.path.isfile(flat)):
            logger.info( "Storing de-flatted %s as %s"%(debiased_f, deflatted))
            time.sleep(1)
            iraf.imarith(debiased_f, "/", flat, deflatted)
        else:
            logger.error( "SOMETHING IS WRONG. Error when dividing %s by the flat field %s!"%(debiased_f, flat))
        
        #Removes the de-biased file
        os.remove(debiased_f)
        
        logger.info( "Updating header with original filename and flat field used.")
        fitsutils.update_par(deflatted, "ORIGFILE", os.path.basename(image))
        fitsutils.update_par(deflatted, "FLATFILE", flat)

        slice_names[i] = deflatted
            
                
    #Moving files to the target directory
    for image in slice_names:
        bkg = get_median_bkg(image)
        fitsutils.update_par(image, "SKYBKG", bkg)
        
        #Get basic statistics for the image
        nsrc, fwhm, ellip, bkg = sextractor.get_image_pars(image)
	plot_image(image)
        
        logger.info( "Sextractor statistics: nscr %d, fwhm (arcsec) %.2f, ellipticity %.2f"% (nsrc, fwhm, ellip))
        print "Sextractor statistics: nscr %d, fwhm (arcsec) %.2f, ellipticity %.2f"% (nsrc, fwhm, ellip)
    
        
        dic = {"FWHM": np.round(fwhm, 3) , "FWHMPIX": np.round(fwhm/0.394, 3) , "NSRC":nsrc, "ELLIP": np.round(ellip, 3)}
        #Update the seeing information from sextractor
        fitsutils.update_pars(image, dic)
    

        
    #Compute the zeropoints
    for image in slice_names:
        zeropoint.calibrate_zeropoint(image)
        
    return slice_names


    
def log_db_phot(myfile):
    '''
    Logs the science files.
    '''
    import archive_mod_sedmdb  # this file requires being able to submit observations with pre-determined id values, which isn't allowed in the original


    sdb = archive_mod_sedmdb.SedmDB(dbname='sedmdbtest', host='localhost')

    cmd = "git version %s"%__file__
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    version = p.communicate()[0].replace('\n',"")
    
    pardic = {"id"              : sdb._id_from_time(),
              "phot_calib_id"   : "",
              "observation_id"  : fitsutils.get_par(myfile, "OBS_ID"),
              "astrometry"      : bool(fitsutils.get_par(myfile, "IQWCS")),
              "filter"          : fitsutils.get_par(myfile, "FILTER").strip(),
              "reducedfile"     : os.path.abspath(myfile),
              "sexfile"         : sextractor.run_sex([myfile])[0],
              "maskfile"        : "",
              "pipeline"        : version,
              "marshal_phot_id" : ""}
              
    sdb.add_phot(pardic)
              
             
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Reduces the photometric images from SEDM Rainbow Camera.
        Requires either a list of images to be reduced, or a directory name where the night photometry is.
        As an option, can run lacosmic to remove cosmic rays.
        By default it invokes astrometry.net before the reduction.
        
        %run rcred.py -d PHOTDIR 
                
        Reduced images are stored in directory called "reduced", within the main directory.
        
        Optionally, it can be used to clean the reduction products generated by this pipeline within PHOTDIR diretory using -c option (clean).
        
        %run rcred.py -d PHOTDIR -c
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-l', '--filelist', type=str, help='File containing the list of fits for the night.',default=None)
    parser.add_argument('-d', '--photdir', type=str, help='Directory containing the science fits for the night.', default=None)
    parser.add_argument('-c', '--clean', action="store_true", help='Clean the reduced images?', default=False)
    parser.add_argument('-o', '--overwrite', action="store_true", help='re-reduce and overwrite the reduced images?', default=False)
    parser.add_argument('-p', '--copy', action="store_false", help='disable the copy the reduced folder to transient', default=True)

    parser.add_argument('--cosmic', action="store_true", default=False, help='Whether cosmic rays should be removed.')

    args = parser.parse_args()
    
    filelist = args.filelist
    photdir = args.photdir
    cosmic = args.cosmic
    clean =  args.clean
    overwrite = args.overwrite
    copy = args.copy
    
    myfiles = []


    if (not photdir is None and clean):
        for f in glob.glob("Flat*"):
            os.remove(f)
        for f in glob.glob("Bias*"):
            os.remove(f)
        for f in glob.glob("a_*fits"):
            os.remove(f)
        if (os.path.isdir(os.path.join(photdir, "reduced"))):
            shutil.rmtree(os.path.join(photdir, "reduced"))

    if (not filelist is None ):
        mydir = os.path.dirname(filelist)
        if (mydir==""):
            mydir = "."
        os.chdir(mydir)

        photdir = mydir

        myfiles = np.genfromtxt(filelist, dtype=None)
        myfiles = [os.path.abspath(f) for f in myfiles]
        
    else:

    	if (photdir is None):
        	timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
        	timestamp = timestamp.split("T")[0].replace("-","")
        	photdir = os.path.join(_photpath, timestamp)

		print '''WARNING! You did not specify the directory or the list:\n
        	- A filelist name with the images you want to reduce [-l]
	        OR
        	- The name of the directory which you want to reduce [-d].
	
		A default name for the directory will be assumed on today's date: %s
		'''%photdir

        mydir = os.path.abspath(photdir)
        #Gather all RC fits files in the folder with the keyword IMGTYPE=SCIENCE
        for f in glob.glob(os.path.join(mydir, "rc*fits")):
            try:
            	if (fitsutils.has_par(f, "IMGTYPE") and ( (fitsutils.get_par(f, "IMGTYPE").upper()=="SCIENCE") or ( "ACQ" in fitsutils.get_par(f, "IMGTYPE").upper())) ):
                	myfiles.append(f)
            except:
            	print "problems opening file %s"%f

    create_masterbias(mydir)
    print "Create masterflat",mydir
    create_masterflat(mydir) 
 
    if (len(myfiles)== 0):
        print "Found no files to process"
	sys.exit()
    else:
        print "Found %d files to process"%len(myfiles)

    
    #Reduce them
    reducedfiles = []
    for f in myfiles:
        print f
        #make_mask_cross(f)
        if ( fitsutils.has_par(f, "IMGTYPE") and (fitsutils.get_par(f, "IMGTYPE").upper() == "SCIENCE" or ("ACQUI" in fitsutils.get_par(f, "IMGTYPE").upper() )) ):
		try:
            		reduced = reduce_image(f, cosmic=cosmic, overwrite=overwrite)
            		reducedfiles.extend(reduced)
		except:
            		print "Error when reducing image %s"%f
            		pass

    #If copy is requested, then we copy the whole folder or just the missing files to transient.
    try:
    	zeropoint.lsq_zeropoint(os.path.join(mydir, "reduced/allstars_zp.log"), os.path.join(mydir, "reduced/zeropoint"))
    except ValueError:
	print "Could not calibrate zeropoint for file %s"%(os.path.join(mydir, "reduced/allstars_zp.log"))

    dayname = os.path.basename(os.path.dirname(os.path.abspath(myfiles[0])))
    reducedname = os.path.join(os.path.dirname(os.path.abspath(myfiles[0])), "reduced")
    if (not photdir is None and copy):
        cmd = "rcp -r %s grbuser@transient.caltech.edu:/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s"%(reducedname, dayname)
        subprocess.call(cmd, shell=True)
    elif (not filelist is None and copy):
        for f in reducedfiles:
            cmd = "rcp %s grbuser@transient.caltech.edu:/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s/."%(f, dayname)
            subprocess.call(cmd, shell=True)

    
