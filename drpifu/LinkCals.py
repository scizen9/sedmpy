#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Link calibrations from a previous night into the current directory
"""
import glob
import sys
import os
import logging
import argparse

logging.basicConfig(
    format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
    datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)


def bias_ready(caldir='./'):
    """Check for all required bias calibration files in calibration directory.

    Args:
        caldir (str): directory to check

    Returns:
        bool: True if bias files are present, False if they are not

    """

    ret = False

    # Do we have all the calibration files?
    # Check biases first
    fb = os.path.exists(os.path.join(caldir, 'bias0.1.fits'))
    f2 = os.path.exists(os.path.join(caldir, 'bias2.0.fits'))
    logging.info("Biases ready?: bias0.1: %d, bias2.0: %d" % (fb, f2))
    if fb and f2:
        ret = True

    return ret


def find_recent(redd, fname, destdir, dstr):
    """Find the most recent version of fname and copy it to destdir.

    Look through sorted list of redux directories to find most recent
    version of the input file.  Copy (link) it into the destination directory.

    Args:
        redd (str): reduced directory (something like /scr2/sedm/redux)
        fname (str): what file to look for
        destdir (str): where the file should go
        dstr (str): YYYYMMDD date string of current directory

    Returns:
        bool: True if file found and copied, False otherwise.

    """

    # Default return value
    ret = False
    # Make sure the file doesn't already exist in destdir
    local_file = glob.glob(os.path.join(destdir, dstr + fname))
    if len(local_file) == 1:
        logging.warning("%s already exists in %s" % (fname, destdir))
        ret = True
    # Search in redd for file
    else:
        # Get all but the most recent reduced data directories
        fspec = os.path.join(redd, '20??????')
        redlist = sorted([d for d in glob.glob(fspec)
                          if os.path.isdir(d)])[0:-1]
        # Go back in reduced dir list until we find our file
        for d in reversed(redlist):
            src = glob.glob(os.path.join(d, '20??????' + fname))
            if len(src) == 1:
                os.symlink(src[0], os.path.join(destdir, dstr + fname))
                ret = True
                logging.info("Found %s in directory %s, linking to %s" %
                             (fname, d, destdir))
                break
    if not ret:
        logging.warning(dstr + fname + " not found")

    return ret


def find_recent_bias(redd, fname, destdir):
    """Find the most recent version of fname and copy it to destdir.

    Look through sorted list of redux directories to find most recent
    version of the input file.  Copy it to the destination directory.

    Args:
        redd (str): reduced directory (something like /scr2/sedm/redux)
        fname (str): what file to look for
        destdir (str): where the file should go

    Returns:
        bool: True if file found and copied, False otherwise.

    """

    # Default return value
    ret = False
    # Make sure the file doesn't already exist in destdir
    local_file = glob.glob(os.path.join(destdir, fname))
    if len(local_file) == 1:
        logging.warning("%s already exists in %s" % (fname, destdir))
        ret = True
    # Search in redd for file
    else:
        # Get all but the most recent reduced data directories
        fspec = os.path.join(redd, '20??????')
        redlist = sorted([d for d in glob.glob(fspec)
                          if os.path.isdir(d)])[0:-1]
        # Go back in reduced dir list until we find our file
        for d in reversed(redlist):
            src = glob.glob(os.path.join(d, fname))
            if len(src) == 1:
                os.symlink(src[0], os.path.join(destdir, fname))
                ret = True
                logging.info("Found %s in directory %s, linking to %s" %
                             (fname, d, os.path.join(destdir, fname)))
                break
    if not ret:
        logging.warning("%s not found" % fname)
    return ret


def link_cals(redd='/scr2/sedmdrp/redux', outdir=None):
    # Get current date string
    cur_date_str = str(outdir.split('/')[-1])
    # Check status
    logging.error("Bad calibrations from this night!")
    logging.info("Let's get our calibrations from a previous night")
    nct = find_recent(redd, '_TraceMatch.pkl', outdir,
                      cur_date_str)
    nctm = find_recent(redd, '_TraceMatch_WithMasks.pkl', outdir,
                       cur_date_str)
    ncg = find_recent(redd, '_HexaGrid.pkl', outdir, cur_date_str)
    ncw = find_recent(redd, '_WaveSolution.pkl', outdir, cur_date_str)
    ncf = find_recent(redd, '_Flat.fits', outdir, cur_date_str)
    if not bias_ready(outdir):
        ncb = find_recent_bias(redd, 'bias0.1.fits', outdir)
        nc2 = find_recent_bias(redd, 'bias2.0.fits', outdir)
    else:
        ncb = True
        nc2 = True
    # Check for failure
    if not nct or not nctm or not ncg or not ncw or not ncf or not ncb \
            or not nc2:
        msg = "Calibration stage failed: trace = %s, trace/mask = %s" \
              "grid = %s, wave = %s, flat = %s, " \
              "bias0.1 = %s, bias2.0 = %s, " \
              "stopping" % (nct, nctm, ncg, ncw, ncf, ncb, nc2)
        sys.exit(msg)
    # If we get here, we are done
    logging.info("Using previous night calibration files")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Link cals from previous night

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--reduxdir', type=str, default='/scr2/sedmdrp/redux',
                        help='Output reduced directory (/scr2/sedmdrp/redux)')

    args = parser.parse_args()

    curdir = os.path.curdir()

    link_cals(redd=args.reduxdir, outdir=curdir)
