# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:18:46 2015

@author: nadiablago
"""

import StringIO
import urllib, base64
from matplotlib import image
from matplotlib import pylab as plt
import glob
from astropy.io import fits as pf
import os, sys
from optparse import OptionParser
import matplotlib
import numpy as np
from astropy.wcs import WCS
import zscale
import time_utils


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-d", "--dir", action="store", type="string", dest="dir", metavar="Directory to be plotted", default=".")
    parser.add_option("-a", "--ra", action="store", type="float", dest="ra", metavar="ra", default=0.0)
    parser.add_option("-b", "--dec", action="store", type="float", dest="dec", metavar="dec", default=0.0)
    parser.add_option("-p", "--pattern", action="store", type="string", dest="pattern", metavar="pattern", default="*.fits")


    (options, args) = parser.parse_args()
    
    if (None in options.__dict__.values()):
    	parser.print_help()
    	sys.exit(0)
    	
    mydir = options.__dict__['dir']
    plot_dir = os.path.join(mydir, "png")
    

    ra = options.__dict__['ra']
    dec = options.__dict__['dec']
    pattern = options.__dict__['pattern']

    print ra, dec
    
    if (not os.path.isdir(plot_dir)):
        os.makedirs(plot_dir)
    for f in glob.glob(mydir + "/" + pattern):
        print "Plotting",os.path.basename(f).replace(".fits", ".png")
        hdulist = pf.open(f)
        if len(hdulist)>1:
            indices = np.arange(len(hdulist)-1)+1
        else:
            indices = np.array([0])
        for i in indices:
        
            prihdr = hdulist[i].header
            img = hdulist[i].data * 1.
            
            nx, ny = img.shape
            if (ra * dec != 0):
    
                # Get pixel coordinates of SN
                wcs = WCS(prihdr)
                try:
                    target_pix = wcs.wcs_sky2pix([(np.array([ra,dec], np.float_))], 1)[0]
                except:
                    print "ERROR when converting sky to wcs. Is astrometry in place? Default coordinates assigned."
                    target_pix = [+nx/2., ny/2.]
    
                print target_pix
            else:
                target_pix = [+nx/2., ny/2.]
               
            img = img - np.nanmin(img)
            av = np.median(img.flatten())
            mi, ma = zscale.zscale(img)
            im = plt.imshow(plt.log10(img), aspect="equal", extent=(0, ny, 0, nx), \
            origin="lower", cmap=matplotlib.cm.gray_r, interpolation="none", vmin=np.log10(av), vmax=np.log10(3*av)) #, interpolation="lanczos")
            plt.scatter(target_pix[0], target_pix[1], marker="x", s=10, c="red")
            plt.colorbar(im)
            filename = os.path.basename(f)
            plt.savefig(os.path.join(plot_dir, filename.replace("."+filename.split(".")[-1], "_{:}.png".format(i))), dpi=200)
            plt.clf()
    
def move_to_discarded(mydir, myfilter, ra, dec):
    import shutil
    
    for f in glob.glob(os.path.join(mydir, myfilter)):
        frames_with_target = get_frames_with_target(f, ra, dec)
        if len(frames_with_target) == 0:
            discarddir = os.path.join(mydir, "discarded")
            if (not os.path.isdir(discarddir)):
                print "Creating directory for discarded files (with no target)"
                os.makedirs(discarddir)
            print "Moving file ",f," to discarded directory",discarddir
            shutil.move(f, os.path.join(discarddir, os.path.basename(f)))
        else:
            print "Object found in frames", frames_with_target, " in file ",f
            print "Extracting this field"
            extract_field_from_moscaic(mydir, os.path.basename(f), frames_with_target, origname=True)
            
def get_frames_with_target(myfile, ra, dec, debug=False):
    
    hdulist = pf.open(myfile)
    if len(hdulist)>1:
        indices = np.arange(len(hdulist)-1)+1
    else:
        indices = np.array([0])
        
    frames = []
    for i in indices:
    
        prihdr = hdulist[i].header
        img = hdulist[i].data * 1.
        
        ny, nx = img.shape
        if (ra * dec != 0):

            # Get pixel coordinates of SN
            wcs = WCS(prihdr)
            try:
                target_pix = wcs.wcs_sky2pix([(np.array([ra,dec], np.float_))], 1)[0]
            except:
                print "ERROR when converting sky to wcs. Is astrometry in place? Default coordinates assigned."
                target_pix = [+nx/2., ny/2.]

            if debug: print i, target_pix
        else:
            target_pix = [+nx/2., ny/2.]
        
        if (target_pix[0] > 0 and target_pix[0]<nx) and (target_pix[1] > 0 and target_pix[1]<ny):
            frames.append(i)
            
    return np.array(frames)
    
def cut_frame_with_target(myfile, ra, dec, h=4000, w=2000, debug=False):
    
    hdulist = pf.open(myfile)[0]

    img = hdulist.data * 1.
    img = img.T
    
    nx, ny = img.shape
    if (ra * dec != 0):

        # Get pixel coordinates of SN
        wcs = WCS(hdulist.header)
        try:
            target_pix = wcs.wcs_sky2pix([(np.array([ra,dec], np.float_))], 1)[0]
        except:
            print "ERROR when converting sky to wcs. Is astrometry in place? Default coordinates assigned."
            target_pix = [+nx/2., ny/2.]

        if debug: print i, target_pix
    else:
        target_pix = [+nx/2., ny/2.]
    
    #If contained in the frame
    if (target_pix[0] > 0 and target_pix[0]<nx) and (target_pix[1] > 0 and target_pix[1]<ny):
        xmin = np.maximum(0, target_pix[0]-w/2.)
        xmax = np.minimum(nx, target_pix[0]+w/2.)
        ymin = np.maximum(0, target_pix[1]-h/2.)
        ymax = np.minimum(ny, target_pix[1]+h/2.)
        
        print "Target", target_pix, xmin, xmax, ymin, ymax
        
        newhdu = pf.PrimaryHDU()
        newhdu.header = hdulist.header
        newhdu.data = hdulist.data[ymin:ymax,xmin:xmax]
        
        newname = os.path.join(os.path.dirname(myfile),"out.fits")
        newhdu.writeto(newname, output_verify="fix", clobber=False)
        print "Extracted region around target and stored to ",newname

    else:
        print "Target not in the frame!"
    
def extract_field_from_moscaic(mydir, myfilter, nfields=None, origname=False):

    if np.isscalar(nfields):
        nfields = np.array([nfields])
    for i, f in enumerate(glob.glob(mydir + "/" + myfilter)):
        hdulist = pf.open(f)
        if nfields==None:
            nfields = np.arange(len(hdulist)-1)+1
        for n in nfields:
            hdu = hdulist[n]
            hdu.header = hdulist[0].header + hdulist[n].header
            hduheader = pf.PrimaryHDU()
            hduheader.header = hdu.header
            hduheader.data = hdu.data
            hdulist1 = pf.HDUList([hduheader, hdu])
            if origname:
                name = os.path.basename(f)
                name = name.replace(".fits", "_%d.fits")
                hdulist1.writeto(name%(n), output_verify="fix", clobber=True)
            else:
                hdulist1.writeto("out%d_%d.fits"%(i,n), output_verify="fix", clobber=True)
        
        
def unify_header(prihdr):
    '''
    Reads the different fields from different telescopes and unifies them under common names.
    
    '''
    
    dic = {"GAIN":0, "RDNOISE":0, "AIRMASS":0, "EXPTIME":0, "AIRMASS":0, "FILTER":0, "MJD-OBS":0}
    filters = {"SDSS-G":"g", "SDSS-R":"r", "SDSS-I":"i", "SDSS-Z":"z"}    
    
    try:
        telescope = prihdr["TELESCOP"]
    except:
        try:
            telescope = prihdr["HIERARCH FPA.TELESCOPE"]
        except:
            print "Could not locate the telescope field"
            return
        
    if telescope == "NOT":
        dic["EXPTIME"] = prihdr["EXPTIME"]
        dic["AIRMASS"] = prihdr["AIRMASS"]
        dic["MJD-OBS"] = time_utils.utc2mjd(prihdr["DATE-OBS"])
        if prihdr["INSTRUME"]=="StanCam":
            dic["FILTER"] = prihdr["STFLTNM"][0]
            dic["GAIN"] = prihdr["GAIN"]
            dic["RDNOISE"] = prihdr["RDNOISE"]
        elif prihdr["INSTRUME"]=="NOTCAM":
            dic["GAIN"] = prihdr["GAIN1"]
            dic["RDNOISE"] = prihdr["RDNOISE1"]
            dic["FILTER"] = prihdr["NCFLTNM2"]
            
    elif telescope == "UKIRT":
        dic["GAIN"] = prihdr["GAIN"]
        dic["RDNOISE"] = prihdr["READNOIS"]
        dic["EXPTIME"] = prihdr["EXP_TIME"]
        dic["AIRMASS"] = prihdr["AMSTART"]
        dic["FILTER"] = prihdr["FILTER"]
        dic["MJD-OBS"] = prihdr["MJD-OBS"]
    elif telescope == "PS1":
        dic["GAIN"] = prihdr["HIERARCH CELL.GAIN"]
        dic["RDNOISE"] = prihdr["HIERARCH CELL.READNOISE"]
        dic["EXPTIME"] = prihdr["EXPTIME"]
        dic["AIRMASS"] = prihdr["AIRMASS"]
        dic["FILTER"] = prihdr["HIERARCH FPA.FILTERID"]
        dic["MJD-OBS"] = prihdr["MJD-OBS"]
    elif telescope == "INT":
            
        dic["GAIN"] = prihdr["GAIN"]
        if ("RDNOISE" in prihdr.keys()):
            dic["RDNOISE"] = prihdr["RDNOISE"]
        else:
            dic["RDNOISE"] = prihdr["READNOIS"]
        dic["EXPTIME"] = prihdr["EXPTIME"]
        dic["AIRMASS"] = prihdr["AIRMASS"]
        dic["FILTER"] = prihdr["WFFBAND"]
        dic["MJD-OBS"] = prihdr["MJD-OBS"]
    elif telescope == "WHT":
        if prihdr["INSTRUME"]=="ACAM":
            dic["GAIN"] = prihdr["GAIN"]
            dic["RDNOISE"] = prihdr["READNOIS"]
            dic["EXPTIME"] = prihdr["EXPTIME"]
            dic["AIRMASS"] = prihdr["AIRMASS"]
            dic["FILTER"] = prihdr["ACAMFILT"]
            dic["MJD-OBS"] = prihdr["MJD-OBS"]
        if "ISIS" in prihdr["INSTRUME"]:
            dic["GAIN"] = prihdr["GAIN"]
            dic["RDNOISE"] = prihdr["READNOIS"]
            dic["EXPTIME"] = prihdr["EXPTIME"]
            dic["AIRMASS"] = prihdr["AIRMASS"]
            dic["FILTER"] = prihdr["ISIFILTA"].replace(".","")
            dic["FILTER2"] = prihdr["ISIFILTB"].replace(".","")
            dic["MJD-OBS"] = prihdr["MJD-OBS"]
            dic["GRISM"] = prihdr["ISIGRAT"]
            dic["SLIT"] = prihdr["ISISLITW"]
    elif telescope == "CFHT 3.6m":
        dic["GAIN"] = prihdr["GAIN"]
        dic["RDNOISE"] = prihdr["RDNOISE"]
        dic["EXPTIME"] = prihdr["EXPTIME"]
        dic["AIRMASS"] = prihdr["AIRMASS"]
        dic["FILTER"] = prihdr["FILTER"]
        dic["MJD-OBS"] = prihdr["MJD-OBS"]
    elif telescope == "Liverpool Telescope":
        dic["GAIN"] = prihdr["GAIN"]
        dic["RDNOISE"] = prihdr["READNOIS"]
        dic["EXPTIME"] = prihdr["EXPTIME"]
        dic["AIRMASS"] = prihdr["AIRMASS"]
        dic["FILTER"] = filters.get(prihdr["FILTER1"], prihdr["FILTER1"])
        dic["MJD-OBS"] = prihdr["MJD"]
        
    elif telescope == "GTC":
        if (prihdr["INSTRUME"]=="OSIRIS"):
            dic["GAIN"] = prihdr["GAIN"]
            dic["EXPTIME"] = prihdr["EXPTIME"]
            dic["AIRMASS"] = prihdr["AIRMASS"]
            dic["FILTER"] = filters.get(prihdr["FILTER1"], prihdr["FILTER1"])
            dic["MJD-OBS"] = prihdr["MJD-OBS"]
            dic["GRISM"] = prihdr["GRISM"]
            dic["SLIT"] = prihdr["SLITW"]
            
            speed = prihdr["RSPEED"]
            #From http://www.gtc.iac.es/instruments/osiris/osiris.php
            if (speed == 200):
                dic["RDNOISE"] = 4.5
            elif (speed == 100):
                dic["RDNOISE"] = 3.5
            else:
                dic["RDNOISE"] = 9
        
    else:
        print "Telescope unknown!"
    
    for i, k in enumerate(dic.keys()):
        if not k in prihdr.keys():
            prihdr[k] = dic[k]
    return prihdr
   
def unify_fits(myfits, overwrite=False, field=-1):
    '''
    Creates a unified standard version of the fits file header, creating standard fields for the most used fields.
    '''
    hdulist = pf.open(myfits)
    
    print hdulist.info()
    
    nfields = len(hdulist)
    datafields = 0
    
    for n in range(nfields):
        if (hdulist[n].data != None):
            datafields +=1
            
    if nfields > 1:
        if (field==-1 and datafields>1):
            print "WARNING, there are several images stored. Plsease, specify which number you want to extract."
            print "Exiting..."
            return 
            
        nfields = np.arange(len(hdulist)-1)+1
        hdu = hdulist[0]
        for n in nfields:
            hdu.header = hdu.header + hdulist[n].header
            
        hdu.data = hdulist[field].data
    else:
        hdu = hdulist[0]
    new_header = unify_header(hdu.header)
    
    primhdu = pf.PrimaryHDU(header = new_header, data = hdu.data)
    hdulist = pf.HDUList([primhdu])
            
    if (overwrite):
        hdulist.writeto(myfits, clobber=True)
    else:
        obj = hdulist[0].header["OBJECT"].replace(" ","")
        filt = hdulist[0].header["FILTER"] 
        inst = hdulist[0].header["INSTRUME"].replace(" ","")
        try:
            slit = "_%.1f" % hdulist[0].header["SLIT"]
            grism = "_" + hdulist[0].header["GRISM"].replace(".", "")

        except:
            slit = ""
            grism = ""
        name = "%s_%s%s%s_%s_1.fits"%(obj, inst,grism, slit, filt)
        
        while os.path.isfile(name):
            seq = int(name.split("_")[-1].replace(".fits", ""))
            name = name.replace("_%d.fits"%seq, "_%d.fits"%(seq+1))
        hdulist.writeto(name, clobber=False)


def get_par(myfits, par):
    '''
    Returns the header parameter from the fits.
    '''
    
    try:
        hdu = pf.open(myfits, ignore_missing_end=True)
        header = hdu[0].header
	if str.upper(par) in header.keys():
	        return header[str.upper(par)]
	else:
		return None
    except IOError:
        return None       
	 
    
def update_par(myfits, par, value):
    '''
    Updates the fits files with the new parameter.
    '''
    hdu = pf.open(myfits, ignore_missing_end=True)
    header = hdu[0].header
    header.set(par, value)
    hdu.writeto(myfits, clobber=True)
    
def update_pars(myfits, pardic):
    '''
    Updates the fits files with the new parameter.
    '''
    hdu = pf.open(myfits, ignore_missing_end=True)
    header = hdu[0].header
    
    for key, value in pardic.iteritems():
        header.set(key, value)
        hdu.writeto(myfits, clobber=True)
    
def has_par(myfits, par):
    '''
    Updates the fits files with the new parameter.
    '''
    hdu = pf.open(myfits, ignore_missing_end=True)
    header = hdu[0].header
    return par in header.keys()
    
def arrange_fits_in_directories(mydir, myfilter, destdir):
    '''
    Automatically sorts the fits files ina directory structure with the filter name as the fierst level, and the integer of the MJD of the observation at
    the second level.
    mydir: where the files to be sorted are.
    myfilter: what filter applies to find them.
    destdir: whhich is the main level directory where fiels should be moved to.
    '''
    import shutil
    for f in glob.glob(os.path.join(mydir, myfilter)):
        header = pf.open(f)[0].header
        filt = header["FILTER"]
        mjd = int(header["MJD-OBS"])
        instrument = header["INSTRUME"]
        destination = os.path.join(destdir, "%s/%s/%d"%(instrument,filt, mjd)) 
        if (not os.path.isdir(destination) ):
            os.makedirs(destination)
            print "Creating directory", destination
        shutil.move(f, os.path.join(destination, os.path.basename(f)))
        
    
def get_gain_ron(fitsfile):

    hdulist = pf.open(fitsfile, ignore_missing_end=True)
    prihdr = hdulist[0].header

    uheader = unify_header(prihdr)
    gain = uheader["GAIN"]
    ron = uheader["RDNOISE"]

    return gain, ron
    
def get_gain_ron_combined_i(gain, ron, n, mode="median"):
    '''
    Returns the combined gain and read out noise for different modes of combining images.
    Avaialable modes: (sum, average, median)
    '''
    
    if (n==1):
        return gain, ron
        
    if (mode=="sum"):
            gain=gain
            ron = np.sqrt(n)*ron
    elif(mode=="average"):
            gain = 1.*n*gain
            ron = np.sqrt(n)*ron
    elif(mode=="median"):
            gain = 2.*n*gain/3.
            ron = np.sqrt(2.*n/3)*ron
    else:
            print "Unknown selected mode. Avaialable values: (sum, average, median)"
    
    print np.round(gain, 5), np.round(ron, 5)
    return np.round(gain, 5), np.round(ron, 5)
    

def get_gain_ron_combined(fitsfile, n, mode="median"):
    gain, ron = get_gain_ron(fitsfile)
    return get_gain_ron_combined_i(gain,ron,n,mode)


def align_combine(fitsdir, myfilter, examine=True):
    from pyraf import iraf 
    
    iraf.noao(_doprint=0)
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    
    os.chdir(fitsdir)
    listfiles = glob.glob(myfilter)
    listfiles.sort()
    
    if (examine):
        print "Opening ",listfiles[0]," to examine."
        iraf.imexamine(input=listfiles[0], \
                    logfile="coords.dat", \
                    keeplog="yes")
        
        with open("align.list",'w') as f:
            for i in listfiles:
                f.write(i+"\n")
    
    print "Aligning with reference:",listfiles[0]
    iraf.imalign( input   =  "@align.list", referenc= listfiles[0], coords  =  "coords.dat", output  = "a_@align.list")  
    
    listfiles = glob.glob("a_"+myfilter)
    listfiles.sort()
    with open("comb.list",'w') as f:
        for i in listfiles:
            f.write(i+"\n")
            
    print "Combining"        
    iraf.imcombine(input = "@comb.list",\
                   output = "out.fits",\
                   combine= "median")
    
    

