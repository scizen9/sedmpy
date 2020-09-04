# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:30:32 2016

@author: nadiablago
"""
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
import subprocess
import os
import sys
from astropy.io import fits as pf
from astropy.wcs import WCS
try:
    import coordinates_conversor as cc
except ImportError:
    import drprc.coordinates_conversor as cc
import numpy as np
import argparse
from matplotlib import pylab as plt
try:
    import zscale
except ImportError:
    import drprc.zscale as zscale
try:
    import fit_utils
except ImportError:
    import drprc.fit_utils as fit_utils
import datetime
import logging


from configparser import ConfigParser
import codecs

parser = ConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.read_file(f)

_logpath = parser.get('paths', 'logpath')
vim = parser.get('paths', 'photpath')

# Log into a file
FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
root_dir = _logpath
timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
timestamp = timestamp.split("T")[0]
logging.basicConfig(format=FORMAT,
                    filename=os.path.join(root_dir,
                                          "listener_{0}.log".format(timestamp)),
                    level=logging.INFO)
logger = logging.getLogger('listener')


def solve_astrometry(img, radius=3.0, with_pix=True, tweak=3):
    """
    img: fits image where astrometry should be solved.
    radius: radius of uncertainty on astrometric position in image.
    """

    try:
        ra = fitsutils.get_par(img, 'OBJRA')
        dec = fitsutils.get_par(img, 'OBJDEC')
    except KeyError:
        ra = fitsutils.get_par(img, 'RA')
        dec = fitsutils.get_par(img, 'DEC')

    mydir = os.path.dirname(img)
    if mydir == "":
        mydir = "."
    os.chdir(mydir)
    
    astromf = os.path.join(os.path.dirname(img), "a_" + os.path.basename(img))
    
    print("Solving astrometry on field with (ra,dec)=", ra, dec,
          "Image", img, "New image", astromf)
    
    cmd = " solve-field --ra %s --dec %s --radius %.4f -p --new-fits " \
          "%s --cpulimit 45 -W none -B none -P none -M none -R none -S none " \
          "-t %d --overwrite %s " % (ra, dec, radius, astromf, tweak, img)
    if with_pix:
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 " \
                    "--scale-high 0.425"
    print(cmd)
    logger.info(cmd)
    
    cmd = cmd + " > /tmp/astrometry_fail  2>/tmp/astrometry_fail"

    subprocess.call(cmd, shell=True)
    
    # Cleaning after astrometry.net
    if os.path.isfile(img.replace(".fits", ".axy")):
        os.remove(img.replace(".fits", ".axy"))
    if os.path.isfile(img.replace(".fits", "-indx.xyls")):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if os.path.isfile("none"):
        os.remove("none")
        
    try:
        os.remove("/tmp/tmp.*")
    except OSError:
        pass

    if not os.path.isfile(astromf):
        print("Astrometry FAILED!")
        logger.error("Astrometry FAILED!")
        
    return astromf


def get_offset_center_failed_astro(af, doplot=False, interactive=True):
    """
    For fields where astrometry is challenging, there is a simple solution.
    Find the brightest peak within the pointing error of the telescope.
    As this fields will usually be centered in Standard stars and very short
    exposure time, the fit intends to

    """

    image = pf.open(af)
    wcs = WCS(image[0].header)
    ra, dec = cc.hour2deg(image[0].header['OBJRA'], image[0].header['OBJDEC'])

    # pra, pdec = wcs.wcs_sky2pix(ra, dec, 0)
    pra, pdec = wcs.wcs_world2pix(ra, dec, 0)
    # Get a local image
    # xl, yl = np.array(wcs.wcs_sky2pix(ra+(60./3600)*np.cos(np.deg2rad(dec)),
    # dec-60./3600, 0), dtype=np.int)
    # xu, yu = np.array(wcs.wcs_sky2pix(ra-(60./3600)*np.cos(np.deg2rad(dec)),
    # dec+60./3600, 0), dtype=np.int)
    imageloc = image[0].data.T[1293-150:1293+150, 1280-150:1280+150]
    
    nx = 300
    ny = 300
            
    def_x = np.argmax(np.sum(imageloc, axis=0))
    def_y = np.argmax(np.sum(imageloc, axis=1))
    
    newx = pra-nx/2.+def_x
    newy = pdec-ny/2.+def_y

    # pra, pdec = wcs.wcs_pix2sky(np.array([[newx[0], newy[0]]],
    #                                     np.float_), 0)[0]
    pra, pdec = wcs.wcs_pix2world(np.array([[newx[0], newy[0]]],
                                           np.float_), 0)[0]
    dra, ddec = cc.get_offset(ra, dec, pra, pdec)
    
    print("Offset", dra, ddec, "Position RA,DEC", pra, pdec)

    x, y, fwhmx, fwhmy, bkg, amp = fit_utils.fit_gauss(imageloc)
    
    if doplot:
        plt.figure(figsize=(8, 8))
        obj = fitsutils.get_par(f, "OBJECT")

        plt.suptitle(obj, fontsize=20)
        zmin, zmax = zscale.zscale(imageloc)
        plt.imshow(imageloc, aspect="auto", interpolation="none",
                   vmin=zmin, vmax=zmax, extent=(0, +300,
                                                                 0, +300))
        plt.plot(x, y, "go", ms=20, label="Centroid using gaussiuan fit.")    
        plt.plot(def_x, def_y, "b*", ms=20, label="Centroid using max/min.")
        plt.plot(150, 150, "wo", ms=20, label="Initial pointing")
        plt.legend()
        
        if interactive:
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(f).replace("raw", "phot"),
                                     os.path.basename(f).replace(".fits",
                                                                 "_std.png")))
        plt.clf()
    return 1, ddec, dra
            
    
def get_offset_center(af, doplot=False, interactive=False):
    """
    Given a fits image, returns the offset in Ra, DEC, that needs to be applied
    for the telescope to go from the current pointing position, to the
    coodinates of the object specified in the fits file.
    """
    
    if not os.path.isfile(af):
        print("File %s does not exist! Returning Zero offsets..." % af)
        return -1, 0, 0
    else:
        image = pf.open(af)
        wcs = WCS(image[0].header)
        rra, rdec = cc.hour2deg(image[0].header['OBJRA'],
                                image[0].header['OBJDEC'])
        # x, y = np.round(wcs.wcs_sky2pix(rra, rdec, 0), 0)
        x, y = np.round(wcs.wcs_world2pix(rra, rdec, 0), 0)
        # pra, pdec = wcs.wcs_pix2sky(np.array([[1293., 1280.]] ,
        # np.float_), 0)[0]
        pra, pdec = wcs.wcs_pix2world(np.array([[1293., 1280.]],
                                               np.float_), 0)[0]
        dra, ddec = cc.get_offset(pra, pdec, rra, rdec)
            
        xl, yu = np.round(wcs.wcs_world2pix(rra+90./3600, rdec-90./3600, 0), 0)
        xu, yl = np.round(wcs.wcs_world2pix(rra-90./3600, rdec+90./3600, 0), 0)

        imageloc = image[0].data.T[xl:xu, yl:yu]

        if imageloc.shape[0] == 0 or imageloc.shape[1] == 0:
            logger.warning("Astrometry has FAILED on this! "
                           "The object is outside the frame! "
                           "Resending to the numb astrometric solution")
            logger.error("Astrometry has FAILED on this! "
                         "The object is outside the frame! "
                         "Resending to the numb astrometric solution")
            print("Pixels are", xl, xu, yl, yu)
            try:
                code, dra, ddec = get_offset_center_failed_astro(
                    af, doplot=doplot, interactive=interactive)
                return 2, dra, ddec
            except:
                return -1, 0, 0
        if doplot:
            plt.figure(figsize=(8, 8))
            
            zmin, zmax = zscale.zscale(imageloc)
    
            obj = fitsutils.get_par(f, "OBJECT")
            plt.suptitle(obj, fontsize=20)
            plt.imshow(imageloc.T, extent=(xl[0], xu[0], yl[0], yu[0]),
                       aspect="equal", interpolation="none",
                       vmin=zmin, vmax=zmax)
            plt.plot(1293., 1280., "ws", ms=7, label="Current pointing")
            plt.plot(x, y, "b*", ms=10, label="Target pointing")
            plt.gca().invert_xaxis()
            plt.legend()
            if interactive:
                plt.show()
            else:
                plt.savefig(os.path.join(
                    os.path.dirname(af).replace("raw", "phot"),
                    os.path.basename(f).replace(".fits", "_a.png")))
            plt.clf()

        return 0, dra, ddec


def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed()  # faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

    
def main(infile, astro=True, doplot=False):
    """
    Computes the coffsets for the acquisition image.
    Returns:
    - Code to indicate if the operation has succeeded.
    - A offset only
    - AB offsets


    Code meaning:
        -  0: Success. All offsets are ok.
        -  1: Offsets computed with Standard Star approach.
        -  2: Offsets were found to be out of the image.
        - -1: Something got very wrong. The astrometry file does not exist.
    """
    offset_file = "/tmp/%s_dra_ddec.txt" % (
        os.path.basename(infile).replace(".fits", ""))

    print("Found image %s as first acquisition image after the slew."
          " Computing offset for IFU..." % infile)
    retcode = 0
    dra = 0
    ddec = 0
    
    # Comment whenever we have the new astrometry file.
    if astro:
        try:
            newfile = solve_astrometry(infile)
            retcode, dra, ddec = get_offset_center(newfile, doplot=doplot,
                                                   interactive=False)
        except Exception as e:
            logger.error(str(sys.exc_info()[0]))
            logger.error(e)
            logger.error('Astrometry failed on file %s. '
                         'Computing the "Failed solution option"' % infile)
            newfile = infile.replace("rc_", "a_rc_")
            
        if not os.path.isfile(newfile):
            retcode, dra, ddec = get_offset_center_failed_astro(
                infile, doplot=True, interactive=False)
    else:
        retcode, dra, ddec = get_offset_center(infile, doplot=doplot,
                                               interactive=False)
        np.savetxt(offset_file,
                   np.array([("CENTER", "%.2f" % dra, "%.2f" % ddec)]),
                   fmt="%s")
        logger.info("Offsets computed for A: \n A %.4f %.4f" % (dra, ddec))

    return retcode, dra, ddec


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""

    Runs astrometry.net on the image specified as a parameter and returns 
    the offset needed to be applied in order to center the object coordinates 
    in the reference pixel.
           
    -i image
    -b It is AB shot.
    -a Astrometry is needed.
    -p plot results.
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--image', type=str, dest="image",
                        help='Fits file with acquisition image.')
    parser.add_argument('-a', '--astro', action='store_true', dest="astro",
                        default=False,
                        help='Solve the astrometry for the image?')
    parser.add_argument('-p', '--plot', action='store_true', dest="plot",
                        default=False,
                        help='Whether we plot the astrometry for the image.')

    args = parser.parse_args()
    
    ret = main(args.image, astro=args.astro, doplot=args.plot)
    print(ret)
    logger.info("Returned from offseets with code %s" % ret[0])
    # retcode, offsets
