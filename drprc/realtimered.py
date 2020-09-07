# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 19:23:29 2016

@author: nadiablago
"""
try:
    import rcred
except ImportError:
    import drprc.rcred as rcred
import subprocess
import glob
import os
import time
import argparse
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
import datetime
import logging
from astropy.io import fits
from matplotlib import pylab as plt
import numpy as np
from scipy.stats import sigmaclip

from configparser import ConfigParser

import codecs

parser = ConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.read_file(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')

plt.switch_backend('Agg')

# Log into a file
log_format = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
root_dir = _logpath
now = datetime.datetime.utcnow()
timestamp = datetime.datetime.isoformat(now)
timestamp = timestamp.split("T")[0]
logging.basicConfig(format=log_format,
                    filename=os.path.join(root_dir,
                                          "rcred_{0}.log".format(timestamp)),
                    level=logging.INFO)
logger = logging.getLogger('realtimered')


def plot_image(image, verbose=False):
    """
    Plots the reduced image into the png folder.

    """
    image = os.path.abspath(image)

    try:
        ff = fits.open(image)[0]
    except OSError:
        logger.error("FATAL! Could not open image %s." % image)
        return
    d = ff.data.astype(np.float64)
    h = ff.header
    imtype = h.get('IMGTYPE', 'None')
    exptime = h.get('EXPTIME', 0)
    name = h.get('OBJECT', 'None')
    filt = h.get('FILTER', 'NA')

    # Sub-dir
    subdir = imtype.lower().strip()
    
    # Change to image directory
    imdir, imname = os.path.split(image)

    # Create destination directory
    png_dir = os.path.join(imdir, "pngraw")
    if not os.path.isdir(png_dir):
        os.makedirs(png_dir)

    png_dir = os.path.join(png_dir, subdir)
    if not os.path.isdir(png_dir):
        os.makedirs(png_dir)

    outfig = os.path.join(png_dir, imname.replace(".fits", "_all.png"))

    if not os.path.isfile(outfig):

        logger.info("Plotting raw %s %s image of %s: %s" %
                    (imtype, filt, name, image))

        corners = {
            "g": [1, 1023, 1, 1023],
            "i": [1, 1023, 1024, 2045],
            "r": [1024, 2045, 1024, 2045],
            "u": [1024, 2045, 1, 1023]
        }

        pltstd = 100.
        for b in corners:
            c, lo, hi = sigmaclip(d[corners[b][2]+150:corners[b][3]-150,
                                    corners[b][0]+150:corners[b][1]-150],
                                  low=2.5, high=2.5)
            std = c.std()
            mid = c.mean()
            d[corners[b][2]:corners[b][3], corners[b][0]:corners[b][1]] -= mid
            if 'bias' in subdir and 'r' in b:
                pltstd = std
            elif 'dome' in subdir or 'twilight' in subdir:
                if std > pltstd:
                    pltstd = std
            else:
                if 'r' in b:
                    if std > pltstd:
                        pltstd = std
            if verbose:
                print("%s %.2f %.2f %.2f" % (b, mid, std, pltstd))

        plt.imshow(d, vmin=-pltstd, vmax=2.*pltstd,
                   cmap=plt.get_cmap('Greys_r'))
        plt.title("{%s} %s %s-band [%ds] " % (imtype, name, filt, exptime))
        plt.colorbar()
        logger.info("As %s", outfig)
        plt.savefig(outfig)
        plt.close()
    else:
        if verbose:
            logger.info("Exists: %s", outfig)


def reduce_on_the_fly(photdir, nocopy=False):
    """
    Waits for new images to appear in the directory to trigger their
    incremental reduction as well.
    """
    # Starting time of this run
    time_ini = datetime.datetime.now()
    # Current time to check against starting time
    time_curr = datetime.datetime.now()
    # Delta time
    deltime = time_curr - time_ini
    # Don't wait longer than this (12hr)
    total_wait = 13.*3600.

    # Do we have files yet?
    whatf = os.path.join(photdir, 'rcwhat.list')
    while not os.path.isfile(whatf) and deltime.total_seconds() < total_wait:
        # Wait 10 minutes
        logger.info("No rcwhat.list file yet, waiting 10 min...")
        time.sleep(600)
        # Check our wait time
        time_curr = datetime.datetime.now()
        deltime = time_curr - time_ini
        if deltime.total_seconds() > total_wait:
            logger.warning("Waited 13hr and no rcwhat file appeared!")
            return

    # Wait for an acquisition
    with open(whatf, 'r') as wtf:
        whatl = wtf.readlines()
    acqs = [wl for wl in whatl if 'ACQ' in wl]
    while len(acqs) <= 0 and deltime.total_seconds() < total_wait:
        for wl in whatl:
            fl = wl.split()[0]
            plot_image(os.path.join(photdir, fl))
        logger.info("No acquisition yet (bright or weather), waiting 10 min...")
        time.sleep(600)
        with open(whatf, 'r') as wtf:
            whatl = wtf.readlines()
        acqs = [wl for wl in whatl if 'ACQ' in wl]
        # Check our wait time
        time_curr = datetime.datetime.now()
        deltime = time_curr - time_ini
        if deltime.total_seconds() > total_wait:
            logger.warning("Waited 13hr and no ACQ appeared!")
            return

    logger.info("We have acquired now, so let's reduce some data!")
    # Get the current the number of files

    nfiles = []
    logger.info("Starting the on-the-fly reduction for directory %s." % photdir)

    dayname = os.path.basename(photdir)

    time_curr = datetime.datetime.now()
    deltime = time_curr - time_ini

    if not nocopy:
        # Make destination directory
        cmd = "rsh transient.caltech.edu -l grbuser mkdir " \
              "/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s" % \
              dayname
        logger.info(cmd)
        subprocess.call(cmd, shell=True)
    
    # Run this loop for 12h after the start.
    while deltime.total_seconds() < total_wait:
        nfilesnew = glob.glob(os.path.join(photdir, "rc*[0-9].fits"))

        if len(nfilesnew) == len(nfiles):
            time.sleep(30)
        else:
            new = [ff for ff in nfilesnew if ff not in nfiles]
            new.sort()
            logger.info("Detected %d new incoming files in the last 30s." %
                        len(new))
            for n in new:
                # Make sure imgtype is available
                if not fitsutils.has_par(n, "IMGTYPE"):
                    print("Image", n, "Does not have an IMGTYPE")
                    time.sleep(0.5)
                    if not fitsutils.has_par(n, "IMGTYPE"):
                        print("Image", n, "STILL Does not have an IMGTYPE")
                        continue
                # Make a plot of image
                imtype = fitsutils.get_par(n, "IMGTYPE")
                plot_image(n)
                if "SCIENCE" in imtype.upper() or "ACQ" in imtype.upper() or \
                        "STANDARD" in imtype.upper():
                    if fitsutils.get_par(n, "EXPTIME") > 30.:
                        do_cosmic = True
                    else:
                        do_cosmic = False
                    reduced = rcred.reduce_image(n, cosmic=do_cosmic)
                    if nocopy:
                        logger.info("Skipping copies to transient")
                    else:
                        # Copy them to transient
                        for r in reduced:
                            toks = os.path.basename(r).split('.')[0].split('_')
                            # Do the filters match?
                            if toks[-2] == toks[-1]:
                                cmd = "rcp %s grbuser@transient.caltech.edu:" \
                                      "/scr3/mansi/ptf/p60phot/" \
                                      "fremling_pipeline/sedm/reduced/%s/" % \
                                      (r, dayname)
                                subprocess.call(cmd, shell=True)
                                logger.info(cmd)
                                logger.info("Successfully copied the image: %s"
                                            % r)
        # Get new delta time
        time_curr = datetime.datetime.now()
        deltime = time_curr - time_ini
        # Update file count
        nfiles = nfilesnew
    logger.info("Concluding RC image processing for %s." % dayname)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""
        Performs on-the-fly RC reduction
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--photdir', type=str, dest="photdir",
                        help='Fits directory file with tonight images.',
                        default=None)
    parser.add_argument('-p', '--nocopy', action="store_true",
                        help='do not copy to transient', default=False)

    args = parser.parse_args()

    if args.photdir is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        pdir = os.path.join(_photpath, timestamp)
    else:
        pdir = args.photdir

    reduce_on_the_fly(os.path.abspath(pdir), nocopy=args.nocopy)
