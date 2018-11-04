# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:01:50 2015

@author: nadiablago
"""
import subprocess
from astropy.io import fits as pf
from astropy.wcs import WCS
from astropy import units as u

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

# import matplotlib
# matplotlib.use("Agg")
from matplotlib import pylab as plt

cfg_parser = ConfigParser()

try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = '/scr2/sedmdrp/sedmpy/drpifu/config/sedmconfig.cfg'

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    cfg_parser.read_file(f)

_rawpath = cfg_parser.get('paths', 'rawpath')
_reduxpath = cfg_parser.get('paths', 'reduxpath')


def finder(myfile, findername, searchrad=0.2/60.):
    
    ra, dec = coordinates_conversor.hour2deg(fitsutils.get_par(myfile, "OBJRA"),
                                             fitsutils.get_par(myfile,
                                                               "OBJDEC"))
    utc = fitsutils.get_par(myfile, "UTC")
    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            
    img = img.T

    wcs = WCS(hdulist.header)

    target_pix = wcs.wcs_world2pix([(np.array([ra, dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_world2pix([(np.array([ra, dec+searchrad], np.float_))],
                                   1)[0]
    dx = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))

    # check bounds
    if (target_pix[0] < 0 or target_pix[0] > img.shape[0] or
       target_pix[1] < 0 or target_pix[1] > img.shape[1]):
        print("ERROR - target outside finder image")
        return
    
    # imgslice = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx,
    #                int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]
    imgslice_target = img[int(target_pix[0])-dx:int(target_pix[0])+dx,
                          int(target_pix[1])-dx:int(target_pix[1])+dx]

    # Maybe the object has moved out of this frame.
    # In this case, make the finder larger.
    x, y = imgslice_target.shape

    print(img.shape, target_pix, corner_pix, dx, int(target_pix[0])-dx,
          int(target_pix[0])+dx, int(target_pix[1])-dx, int(target_pix[1])+dx)

    if (x < 2*dx-1) or (y < 2*dx-1):
        imgslice_target = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx,
                              int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]

    # zmin, zmax = zscale.zscale()
    zmin = np.percentile(imgslice_target.flatten(), 5)
    zmax = np.percentile(imgslice_target.flatten(), 98.5)
   
    print("Min: %.1f, max: %.1f" % (zmin, zmax))
    gc = aplpy.FITSFigure(myfile, figsize=(10, 9), north=True)
    gc.show_grayscale(vmin=zmin, vmax=zmax, smooth=1, kernel="gauss")
    gc.add_scalebar(0.1/60.)
    gc.scalebar.set_label('10 arcsec')
    gc.scalebar.set_color('white')
    gc.recenter(ra, dec, searchrad)

    ras = np.array([ra, ra])
    decs = np.array([dec, dec])
    dxs = np.array([0, searchrad/10 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/10, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="red", facecolor="red",
                   head_width=0)

    ras = np.array([ra+searchrad*0.7 / np.cos(np.deg2rad(dec)),
                    ra+searchrad*0.7 / np.cos(np.deg2rad(dec))])
    decs = np.array([dec-searchrad*0.9, dec-searchrad*0.9])
    dxs = np.array([0, searchrad/5 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/5, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="white", facecolor="white")
    gc.add_label(ras[0]+dxs[0]*1.1, decs[0]+dys[0]*1.1, 'N', relative=False,
                 color="white", horizontalalignment="center")
    gc.add_label(ras[1]+dxs[1]*1.1, decs[1]+dys[1]*1.1, 'E', relative=False,
                 color="white", horizontalalignment="center")

    img_name = fitsutils.get_par(myfile, "NAME").strip()
    img_filter = fitsutils.get_par(myfile, "FILTER")
    gc.add_label(0.05, 0.95, 'Object: %s' % img_name, relative=True,
                 color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.9, 'Filter: SDSS %s' % img_filter, relative=True,
                 color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.85, 'Coordinates: RA=%s DEC=%s' %
                 (coordinates_conversor.deg2hour(ra, dec)), relative=True,
                 color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.8, "UTC: %s" % utc, relative=True,
                 color="white", horizontalalignment="left")
    
    gc.save(findername)
    print("Created %s" % findername)
    

def simple_finder(myfile, findername):

    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            
    img = img.T
    img = img[1165:2040, 1137:2040]
    newimg = img
    # ndimage.filters.gaussian_filter(img, 1, order=0, mode='constant',
    # cval=0.0, truncate=20.0)

    # name = fitsutils.get_par(myfile, "NAME")
    # filter = fitsutils.get_par(myfile, "FILTER")

    # zmin, zmax = zscale.zscale()
    zmin = np.percentile(newimg.flatten(), 10)
    zmax = np.percentile(newimg.flatten(), 99)
    plt.figure(figsize=(10, 9))
    plt.imshow(newimg, origin="lower", cmap=plt.get_cmap('gray'),
               vmin=zmin, vmax=zmax)

    plt.savefig(findername)

    print("Created ", findername)


def simple_finder_astro(myfile, findername, searchrad=28./3600):  

    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            

    # name = fitsutils.get_par(myfile, "NAME")
    # filter = fitsutils.get_par(myfile, "FILTER")

    ra, dec = coordinates_conversor.hour2deg(fitsutils.get_par(myfile, "OBJRA"),
                                             fitsutils.get_par(myfile,
                                                               "OBJDEC"))
    
    wcs = WCS(hdulist.header)

    target_pix = wcs.wcs_world2pix([(np.array([ra, dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_world2pix([(np.array([ra+searchrad, dec+searchrad],
                                              np.float_))], 1)[0]
    targ_x = int(target_pix[0])
    targ_y = int(target_pix[1])
    # Size of the finder in pixels
    
    dx = int(np.abs(np.ceil(corner_pix[0] - target_pix[0])))
    dy = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))
    
    # size = int( (searchrad/0.394)/2)
    
    # zmin, zmax = zscale.zscale()
    newimg = img[targ_x-dx: targ_x+dx, targ_y-dy: targ_y+dy]

    zmin = np.percentile(newimg.flatten(), 5)
    zmax = np.percentile(newimg.flatten(), 98.5)
    
    print("X %d Y %d Size %d, %d zmin=%.2f zmax=%.2f. Size = %s" %
          (targ_x, targ_y, dx, dy, zmin, zmax, newimg.shape))

    from astropy.visualization.wcsaxes import SphericalCircle

    plt.figure(figsize=(10, 9))
    ax = plt.subplot(projection=wcs)
    ax.imshow(np.flip(newimg, axis=0),
              origin="lower", cmap=plt.get_cmap('gray'), vmin=zmin, vmax=zmax)
    r = SphericalCircle((ra * u.deg, dec * u.deg), 5./3600 * u.degree,
                        edgecolor='red', facecolor='none',
                        transform=ax.get_transform('fk5'))
    ax.add_patch(r)

    # ax = plt.gca()
    # ax.scatter(ra, dec, transform=ax.get_transform('fk5'), s=20,
    #       edgecolor='red', facecolor='none')

    # plt.plot(dy, dx, "+", color="r", ms=20, mfc=None, mew=2)
    # plt.plot(Y, X, "+", color="r", ms=20, mfc=None, mew=2)
    # plt.xlim(Y-dy, Y+dy, X-dx, X+dx)

    plt.savefig(findername)

    print("Created ", findername)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    
    Creates a finder chart for every acquisition image in the folder
    specified as a parameter.
    As a final step, it copies the acquisition image to the "agn" machine
    to visualize it.
        
    """, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-d', '--rcdir', type=str, dest="rcdir",
                        help='Directory with rc images from tonight.',
                        default=None)
    parser.add_argument('-r', '--reduxdir', type=str, dest="reduxdir",
                        help='Directory with reduced ifu images from tonight.',
                        default=None)
    parser.add_argument('-i', '--imfile', type=str, dest="imfile",
                        help='IFU image that requires a finder',
                        default=None)
    
    args = parser.parse_args()

    imfile = args.imfile

    if imfile:
        timestamp = imfile.split('/')[-2]
        rcdir = os.path.join(_rawpath, timestamp)
        reduxdir = '/'.join(imfile.split('/')[0:-1])
        objnam = fitsutils.get_par(imfile, "OBJECT").split()[0]
        if 'STD' in objnam:
            objnam = objnam.split('STD-')[-1].split()[0]
    else:
        rcdir = args.rcdir
        reduxdir = args.reduxdir
        objnam = None
        if rcdir is None:
            timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
            timestamp = timestamp.split("T")[0].replace("-", "")
            rcdir = os.path.join(_rawpath, timestamp)
        else:
            timestamp = os.path.basename(os.path.abspath(rcdir))

        if reduxdir is None:
            reduxdir = os.path.join(_reduxpath, timestamp)

    os.chdir(reduxdir)
    print("Making finder for object: %s" % objnam)
    print("Changed to directory where the reduced data is: %s" % reduxdir)
    print("Getting acquisition images from directory: %s" % rcdir)

    if not (os.path.isdir("finders")):
        os.makedirs("finders")

    # We gather all RC images to locate the Acquisition ones.
    files = glob.glob(os.path.join(rcdir, "rc*fits"))
    files.sort()
    filesacq = []

    for f in files:
        try:
            ff = pf.open(f)
        except OSError:
            print("WARNING - corrupt fits file: %s" % f)
            continue
        if "IMGTYPE" in ff[0].header:
            imgtype = ff[0].header["IMGTYPE"]
            if ((imgtype.upper() == "ACQUISITION" or
                 "ACQ" in imgtype.upper()) and ("TEST" not in imgtype.upper())):
                if "OBJECT" in ff[0].header:
                    obj = ff[0].header["OBJECT"]
                    if objnam:
                        if objnam in obj:
                            filesacq.append(f)
                    else:
                        filesacq.append(f)

    n_acq = len(filesacq)
    print("Found %d files for finders:\n%s" % (n_acq, filesacq))

    n_find = 0
    for f in filesacq:
        n_find += 1
        print("Trying finder %d of %d: %s" % (n_find, n_acq, f))
        try:
            objnam = fitsutils.get_par(f, "OBJECT")
        except:
            print('There is no object in this file %s. Skipping the finder'
                  ' and moving to the next file.' % f)
            continue

        # We generate only one finder for each object.
        name = f.split('/')[-1].split(".")[0]
        objnam = objnam.split()[0]
        filt = fitsutils.get_par(f, "FILTER")
        finderplotf = 'finder_%s_%s_%s.png' % (name, objnam, filt)
        finderpath = os.path.join(reduxdir, os.path.join("finders/",
                                                         finderplotf))
        # Check if it was already done
        if not os.path.isfile(finderpath):
            # link rc image into reduxdir
            dest = os.path.join(reduxdir, f.split('/')[-1])
            if os.path.isfile(dest):
                print("RC ACQ image already exists: %s" % dest)
            else:
                os.symlink(f, dest)
            print("Solving astrometry", dest)
            # Solving for astrometry
            astrof = dest.replace(".fits", "_astrom.fits")
            if not os.path.exists(astrof):
                returncode = subprocess.call(['/scr2/sedmdrp/bin/do_astrom',
                                              dest])
                if returncode != 0:
                    print("Astrometry failed, perform median subtraction")
                    returncode = subprocess.call(
                        ['/scr2/sedmdrp/spy',
                         '/scr2/sedmdrp/sedmpy/drpifu/med_sub.py', '-i',
                         f.split('/')[-1]])
                    if returncode != 0:
                        print("Astrometry failed for %s, skipping finder %s" %
                              (dest, finderpath))
                        continue
                    returncode = subprocess.call(['/scr2/sedmdrp/bin/do_astrom',
                                                  dest])
                if returncode != 0:
                    print("Astrometry failed for %s, skipping finder %s" %
                          (dest, finderpath))
                    continue
            else:
                print("Astrometry file already exists: %s" % astrof)
            # Check results
            if not os.path.isfile(astrof):
                print("Astrometry results not found %s" % astrof)
                continue

            try:
                finder(astrof, finderpath)
            except ValueError:
                print("Bad astrometry for this file: %s" % astrof)
                continue
            except AttributeError:
                print("Error when generating the finder for file %s" % f)
                print(sys.exc_info()[0])
                simple_finder_astro(astrof, finderpath)

            except:
                print("Error when generating the finder for file %s. "
                      "Probably montage is broken." % astrof)
                print(sys.exc_info()[0])
                simple_finder_astro(astrof, finderpath)
        else:
            print("Finder already exists: %s" % finderpath)
