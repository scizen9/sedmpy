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
import json
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
try:
    from target_mag import get_target_mag
except ImportError:
    from drprc.target_mag import get_target_mag
try:
    from fritz_status_update import update_fritz_status
except ImportError:
    from fritz.fritz_status_update import update_fritz_status

try:
    from pysedmpush import slack
except ImportError:
    slack = None
    print("you need to install pysedmpush to be able to push on slack")
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

SLACK_CHANNEL = "pysedm-report"

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


def plot_raw_image(image, verbose=False, ut_id=None):
    """
    Plots the reduced image into the png folder.

    """
    image = os.path.abspath(image)

    utdate = int(os.path.basename(image).split('rc')[-1].split('_')[0])

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
    focpos = h.get('FOCPOS', 0.)

    # Sub-dir
    subdir = 'test'
    if utdate < 20181210:
        if 'dome' in imtype or 'bias' in imtype:
            subdir = imtype.lower().strip()
        else:
            obtype = h.get('OBJTYPE', 'None')
            if 'Science' in imtype:
                if 'Calibration' in obtype:
                    obname = h.get('OBJNAME')
                    if 'twilight' in obname:
                        subdir = 'twilight'
                    elif 'focus' in obname:
                        subdir = 'focus'
                    else:
                        subdir = 'test'
                else:
                    subdir = obtype.lower().strip()
            else:
                subdir = obtype.lower().strip()
    else:
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

    # Handle Bias and Flat images
    if 'Bias' in imname or 'Flat' in imname:
        out_suffix = '.png'
    else:
        out_suffix = '_all.png'
    # Handle gzipped files
    if imname.endswith("gz"):
        outfile = imname.replace(".fits.gz", out_suffix)
    else:
        outfile = imname.replace(".fits", out_suffix)
    outfig = os.path.join(png_dir, outfile)

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
        if 'FOCUS' in imtype.upper():
            if ut_id is not None:
                plt.title("{%s} %.2f %s %s-band [%ds] " %
                          (imtype, focpos, ut_id, filt, exptime))
            else:
                plt.title("{%s} %.2f %s-band [%ds] " %
                          (imtype, focpos, filt, exptime))
        else:
            if ut_id is not None:
                plt.title("{%s} %s %s %s-band [%ds] " %
                          (imtype, ut_id, name, filt, exptime))
            else:
                plt.title("{%s} %s %s-band [%ds] " %
                          (imtype, name, filt, exptime))
        plt.colorbar()
        logger.info("As %s", outfig)
        plt.savefig(outfig)
        plt.close()
        if verbose:
            print(outfig)
    else:
        if verbose:
            logger.info("Exists: %s", outfig)


def reduce_on_the_fly(photdir, nocopy=False, proc_na=False, do_phot=False):
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
            utid = "_".join(fl.split("_")[1:]).split(".")[0]
            plot_raw_image(os.path.join(photdir, fl), ut_id=utid)
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

    phot_zp = {'u': None, 'g': None, 'r': None, 'i': None}
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
                req_id = fitsutils.get_par(n, "REQ_ID")
                imname = os.path.basename(n).replace(".fits", "")
                utid = "_".join(imname.split("_")[1:])
                plot_raw_image(n, ut_id=utid)
                if "SCIENCE" in imtype.upper() or "ACQ" in imtype.upper() or \
                        "STANDARD" in imtype.upper():
                    if fitsutils.get_par(n, "EXPTIME") > 30.:
                        do_cosmic = True
                    else:
                        do_cosmic = False
                    reduced = rcred.reduce_image(n, cosmic=do_cosmic)
                    # perform quick photometry if requested
                    if do_phot:
                        for rf in reduced:
                            if fitsutils.get_par(rf, "ONTARGET"):
                                target_object = fitsutils.get_par(rf, "OBJECT")
                                target_filter = target_object.split()[-1]
                                target_name = target_object.split()[0]
                                logger.info(
                                    "Getting quick %s-band mag for %s in %s" %
                                    (target_filter, target_name, rf))
                                target_mag, target_magerr, std_zp = \
                                    get_target_mag(rf, zeropoint=phot_zp)
                                if target_mag is None or target_magerr is None:
                                    logger.warning("Quick mag failed!")
                                else:
                                    logger.info("Quick MAG = %.3f +- %.3f" %
                                                (target_mag, target_magerr))
                                if std_zp is not None:
                                    logger.info("Quick MAG_ZP: %.3f" % std_zp)
                                    if phot_zp[target_filter] is None:
                                        phot_zp[target_filter] = std_zp
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
                                # push to slack
                                png_dir = os.path.dirname(r) + '/png/'
                                basename = os.path.basename(r).split('.')[0]
                                imgf = png_dir + basename + '.png'
                                title = "RC image: %s | %s" % (basename, imtype)
                                if slack is not None:
                                    try:
                                        slack.push_image(imgf, caption="",
                                                         title=title,
                                                         channel=SLACK_CHANNEL)
                                    except json.decoder.JSONDecodeError:
                                        print("json error, cannot push %s"
                                              % imgf)
                                else:
                                    print("Cannot push: %s" % imgf)
                    if "SCIENCE" in imtype.upper():
                        t_now = datetime.datetime.now()
                        stat_str = "Complete %4d%02d%02d %02d_%02d_%02d" % (
                            t_now.year, t_now.month, t_now.day,
                            t_now.hour, t_now.minute, t_now.second)
                        update_fritz_status(request_id=req_id, status=stat_str)
                elif "POINTING" in imtype.upper():
                    if fitsutils.get_par(n, "EXPTIME") > 30.:
                        do_cosmic = True
                    else:
                        do_cosmic = False
                    reduced = rcred.reduce_image(n, cosmic=do_cosmic)
                    for r in reduced:
                        # push to slack
                        png_dir = os.path.dirname(r) + '/png/'
                        basename = os.path.basename(r).split('.')[0]
                        imgf = png_dir + basename + '.png'
                        title = "RC image: %s | %s" % (basename, imtype)
                        if slack is not None:
                            # only push r-band where ref pixel is
                            if '_r.png' in imgf:
                                slack.push_image(imgf, caption="",
                                                 title=title,
                                                 channel=SLACK_CHANNEL)
                        else:
                            print("Cannot push: %s" % imgf)
                elif "NA" in imtype.upper() and proc_na:
                    if fitsutils.get_par(n, "EXPTIME") > 30.:
                        do_cosmic = True
                    else:
                        do_cosmic = False
                    reduced = rcred.reduce_image(n, cosmic=do_cosmic)
        # Get new delta time
        time_curr = datetime.datetime.now()
        deltime = time_curr - time_ini
        # Update file count
        nfiles = nfilesnew


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""
        Performs on-the-fly RC reduction
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--photdir', type=str, dest="photdir",
                        help='Fits directory file with tonight images.',
                        default=None)
    parser.add_argument('-n', '--nocopy', action="store_true",
                        help='do not copy to transient', default=False)
    parser.add_argument('-p', '--proc_na', action="store_true",
                        help='process NA image types', default=False)

    args = parser.parse_args()

    if args.photdir is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        pdir = os.path.join(_photpath, timestamp)
    else:
        pdir = args.photdir

    photdir = os.path.abspath(pdir)

    reduce_on_the_fly(photdir, nocopy=args.nocopy, proc_na=args.proc_na)

    dayname = os.path.basename(photdir)
    logger.info("Concluding RC image processing for %s." % dayname)
