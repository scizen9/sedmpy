# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:01:50 2015

@author: nadiablago
"""

from __future__ import print_function

import matplotlib
matplotlib.use("Agg")
from astropy.io import fits as pf
from astropy.wcs import WCS
from astropy import units as u

import numpy as np
import aplpy
import coordinates_conversor
import fitsutils
import datetime
import os, sys
import glob
import argparse

from matplotlib import pylab as plt

import rcred
 
from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')

def finder(myfile, findername, searchrad=0.2/60.):
    
    ra, dec = coordinates_conversor.hour2deg(fitsutils.get_par(myfile, "OBJRA"), fitsutils.get_par(myfile, "OBJDEC"))
    

    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            
    img = img.T

    wcs = WCS(hdulist.header)

    target_pix = wcs.wcs_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_world2pix([(np.array([ra,dec+searchrad], np.float_))], 1)[0]
    dx = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))
    
    imgslice = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx, int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]
    imgslice_target = img[int(target_pix[0])-dx:int(target_pix[0])+dx, int(target_pix[1])-dx:int(target_pix[1])+dx]

    #Maybe the object has moved out of this frame. In this case, make the finder larger.
    x, y = imgslice_target.shape

    print (img.shape, target_pix, corner_pix, dx, int(target_pix[0])-dx, int(target_pix[0])+dx, int(target_pix[1])-dx, int(target_pix[1])+dx)

    if  (x < 2*dx-1) or (y< 2*dx-1):
    	imgslice_target = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx, int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]


    #zmin, zmax = zscale.zscale()
    zmin = np.percentile(imgslice_target.flatten(), 5)
    zmax = np.percentile(imgslice_target.flatten(), 98.5)
   
    print ("Min: %.1f, max: %.1f"%(zmin, zmax) )
    gc = aplpy.FITSFigure(myfile, figsize=(10,9), north=True)
    gc.show_grayscale(vmin=zmin, vmax=zmax, smooth=1, kernel="gauss")
    gc.show_scalebar(0.1/60.)
    gc.scalebar.set_label('10 arcsec')
    gc.scalebar.set_color('white')
    gc.recenter(ra, dec, searchrad)

    #gc.show_markers(ra,dec+searchrad/20.,edgecolor='red',facecolor='none',marker="|",s=250, lw=10)
    #gc.show_markers(ra-(searchrad/20.)/np.cos(np.deg2rad(dec)),dec,edgecolor='red',facecolor='none',marker="_",s=250, lw=10)

    ras = np.array([ra , ra])
    decs = np.array([dec, dec])
    dxs = np.array([0, searchrad/10 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/10, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="red", facecolor="red", head_width=0)

    ras = np.array([ra+searchrad*0.7/ np.cos(np.deg2rad(dec)), ra+searchrad*0.7/ np.cos(np.deg2rad(dec))])
    decs = np.array([dec-searchrad*0.9, dec-searchrad*0.9])
    dxs = np.array([0, searchrad/5 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/5, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="k", facecolor="k")
    gc.add_label(ras[0]+dxs[0]*1.1, decs[0]+dys[0]*1.1, 'N', relative=False, color="k", horizontalalignment="center")
    gc.add_label(ras[1]+dxs[1]*1.1, decs[1]+dys[1]*1.1, 'E', relative=False, color="k", horizontalalignment="center")


    name = fitsutils.get_par(myfile, "NAME").strip()
    filter = fitsutils.get_par(myfile, "FILTER")
    gc.add_label(0.05, 0.95, 'Object: %s'%(name), relative=True, color="white", horizontalalignment="left")                   
    gc.add_label(0.05, 0.9, 'Coordinates: RA=%s DEC=%s'%(coordinates_conversor.deg2hour(ra, dec)), relative=True, color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.84, 'Filter: SDSS %s'%filter, relative=True, color="white", horizontalalignment="left")
    
    gc.save(findername)
    


def simple_finder(myfile, findername, searchrad=0.2/60.):  

    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            
    img = img.T
    img = img[1165:2040, 1137:2040]
    newimg = img #ndimage.filters.gaussian_filter(img, 1, order=0, mode='constant', cval=0.0, truncate=20.0)

    name = fitsutils.get_par(myfile, "NAME")
    filter = fitsutils.get_par(myfile, "FILTER")

    #zmin, zmax = zscale.zscale()
    zmin = np.percentile(newimg.flatten(), 10)
    zmax = np.percentile(newimg.flatten(), 99)
    plt.figure(figsize=(10,9))
    plt.imshow(newimg, origin="lower", cmap=plt.get_cmap('gray'), vmin=zmin, vmax=zmax)


    plt.savefig(findername)

    print ("Created ", findername)

def simple_finder_astro(myfile, findername, searchrad=28./3600):  

    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            

    name = fitsutils.get_par(myfile, "NAME")
    filter = fitsutils.get_par(myfile, "FILTER")

    ra, dec = coordinates_conversor.hour2deg(fitsutils.get_par(myfile, "OBJRA"), fitsutils.get_par(myfile, "OBJDEC"))
    
    wcs = WCS(hdulist.header)

    target_pix = wcs.wcs_world2pix([(np.array([ra,dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_world2pix([(np.array([ra+searchrad,dec+searchrad], np.float_))], 1)[0]
    X = int(target_pix[0])
    Y = int(target_pix[1])
    #Size of the finder in pixels
    
    dx = int(np.abs(np.ceil(corner_pix[0] - target_pix[0])))
    dy = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))
    
    #size = int( (searchrad/0.394)/2)
    
    #zmin, zmax = zscale.zscale()
    newimg = img[X-dx : X+dx, Y-dy : Y+dy]

    zmin = np.percentile(newimg.flatten(), 5)
    zmax = np.percentile(newimg.flatten(), 98.5)
    
    print ("X %d Y %d Size %d, %d zmin=%.2f zmax=%.2f. Size = %s"%(X,Y,dx,dy,zmin,zmax, newimg.shape))

    from astropy.visualization.wcsaxes import SphericalCircle

    plt.figure(figsize=(10,9))
    ax = plt.subplot(projection=wcs)
    ax.imshow(np.flip(newimg, axis=0), \
        origin="lower", cmap=plt.get_cmap('gray'), vmin=zmin, vmax=zmax)
    r = SphericalCircle((ra * u.deg, dec * u.deg), 5./3600 * u.degree,
                     edgecolor='red', facecolor='none',
                     transform=ax.get_transform('fk5'))
    ax.add_patch(r)

    #ax = plt.gca()
    #ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=20,
    #       edgecolor='red', facecolor='none')

    #plt.plot(dy, dx, "+", color="r", ms=20, mfc=None, mew=2)
    #plt.plot(Y, X, "+", color="r", ms=20, mfc=None, mew=2)
    #plt.xlim(Y-dy, Y+dy, X-dx, X+dx)

    plt.savefig(findername)

    print ("Created ", findername)

    
if __name__=="__main__":  
    parser = argparse.ArgumentParser(description=\
    '''
    
    Creates a finder chart for every acquisition image in the folder specified as a parameter.
    As a final step, it copies the acquisition image to the "agn" machine to visualize it.
        
    ''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.', default=None)
    
    args = parser.parse_args()
    
    photdir = args.photdir
    
    
    if (photdir is None):
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-","")
        photdir = os.path.join(_photpath, timestamp)
    else:
        timestamp = os.path.basename(os.path.abspath(photdir))
    
    os.chdir(photdir)
    print ("Changed to directory where the data is. %s"%photdir)
    
    if not (os.path.isdir("finders")):
        os.makedirs("finders")
    
    #We gather all RC images to locate the Acquisition ones.
    files = glob.glob("rc*fits")
    files.sort()
    filesacq = []
    
    for f in files:
        if ( (fitsutils.get_par(f, "IMGTYPE").upper() == "ACQUISITION"  or "ACQ" in fitsutils.get_par(f, "IMGTYPE").upper()) and ("TEST" not in fitsutils.get_par(f, "IMGTYPE").upper())):
            filesacq.append(f)
    
    print ("Found %d files for finders: %s"%(len(filesacq), filesacq))
    
    for f in filesacq:
        try:
            object = fitsutils.get_par(f, "OBJECT").upper()
        except:
            print ('There is no object in this file %s. Skipping the finder and moving to the next file.'%f)
            continue
    
    	#We generate only one finder for each object.
    	name = f.split(".")[0]
    	findername = 'finder_%s_%s_%s.png'%(name, fitsutils.get_par(f, "NAME"), fitsutils.get_par(f, "FILTER"))
	finderpath = os.path.join(photdir, os.path.join("finders/",findername))
    	if not os.path.isfile(finderpath):
        	print ("Generating finder", finderpath)
        	#Solving for astrometry
        	astrof = rcred.solve_astrometry(f)
    
            	try:
                	finder(astrof, finderpath)
            	except AttributeError:
                	print ("Error when generating the finder for file %s"%f)
                	print (sys.exc_info()[0])
    
                	simple_finder_astro(astrof, finderpath)
    
            	except:
                	print ("Error when generating the finder for file %s. Probably montage is broken."%astrof)
                	print (sys.exc_info()[0])
                	simple_finder_astro(astrof, finderpath)

