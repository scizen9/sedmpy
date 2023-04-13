"""Conduct automatic reduction of SEDM data on sedmdrp@minar

Functions
    * :func:`go`           outer loop waits for new data directory
    * :func:`obs_loop`     one night observing loop
    * :func:`cpcal`        copies calibration images into redux directory
    * :func:`cpprecal`     copies calibration images from previous day directory
    * :func:`find_recent`  finds the most recent processed calibration file
    * :func:`cpsci`        copies new science images files into redux directory
    * :func:`proc_stds`    processes standard star observations
    * :func:`proc_bias_crrs`  processes biases and CR rejection
    * :func:`proc_bkg_flex`   processes bkg sub and flex calculation
    * :func:`docp`            low level copy routine
    * :func:`cal_proc_ready`  check if all required raw cal images are present
    * :func:`cube_ready`      check if all required cal files are present
    * :func:`bias_ready`    check if master bias files are present

Note:
    This is used as a python script as follows::

        usage: AutoReduce.py [-h] [--rawdir RAWDIR] [--reduxdir REDUXDIR]

        optional arguments:
          -h, --help           show this help message and exit
          --rawdir RAWDIR      Input raw directory (/data/sedmdrp/raw)
          --reduxdir REDUXDIR  Output reduced directory (/data/sedmdrp/redux)
          --wait               Wait for new directory first (False)
          --piggyback          Don't copy data, copied by another script (False)
          --skip_precal        Skip check of previous day for cal files? (False)
          --date YYYYMMDD      Select date to process (None)
          --update YYYYMMDD    UTDate directory to update (None)
          --local              Process data locally, no push to marshal or slack
                               or db update (False)
          --nodb               Do not update SEDM Db (False)

"""
import time
import glob
import sys
import os
import re
import json
import subprocess
# from subprocess import Popen, PIPE
import astropy.io.fits as pf
import logging
import argparse
import smtplib
from email.message import EmailMessage
import db.SedmDb
from datetime import datetime

from astropy.time import Time
from astropy.coordinates import Angle
import astroplan

import sedmpy_version

try:
    import marshal_commenter as mc
except ImportError:
    import growth.marshal_commenter as mc

try:
    import rcimg
except ImportError:
    import drprc.rcimg as rcimg

drp_ver = sedmpy_version.__version__
logging.basicConfig(
    format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
    datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)

# Get pipeline configuration
# Find config file: default is sedmpy/config/sedmconfig.json
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(sedmpy_version.CONFIG_DIR, "sedmconfig.json")
with open(configfile) as config_file:
    sedm_cfg = json.load(config_file)

# Get paths
_rawpath = sedm_cfg['paths']['rawpath']
_reduxpath = sedm_cfg['paths']['reduxpath']
_srcpath = sedm_cfg['paths']['srcpath']
_nomfs = sedm_cfg['nominal_file_size']


def link_refcube(curdir='./', date_str=None):
    # source (ref) and destination
    srcpre = '/data/sedmdrp/redux/ref/20230407_'
    dstpre = curdir + '/%s_' % date_str
    grid = 'HexaGrid.pkl'
    trace = 'TraceMatch.pkl'
    trace_masks = 'TraceMatch_WithMasks.pkl'
    # link in
    os.symlink(srcpre + grid, dstpre + grid)
    os.symlink(srcpre + trace, dstpre + trace)
    os.symlink(srcpre + trace_masks, dstpre + trace_masks)


def cube_ready(caldir='./', cur_date_str=None, avrmslim=75.0):
    """Check for all required calibration files in calibration directory.

    Args:
        caldir (str): directory to check
        cur_date_str (str): current date in YYYYMMDD format
        avrmslim (float): wavelength solution rms limit for success

    Returns:
        bool: True if calibration files are present, False if any are missing.

    """

    ret = False

    # Files to look for
    if cur_date_str is None:
        tmf = 'TraceMatch.pkl'
        tmmf = 'TraceMatch_WithMasks.pkl'
        hgf = 'HexaGrid.pkl'
        wsf = 'WaveSolution.pkl'
        fff = 'Flat.fits'
        statf = 'wavesolution_stats.txt'
    else:
        tmf = cur_date_str + '_TraceMatch.pkl'
        tmmf = cur_date_str + '_TraceMatch_WithMasks.pkl'
        hgf = cur_date_str + '_HexaGrid.pkl'
        wsf = cur_date_str + '_WaveSolution.pkl'
        fff = cur_date_str + '_Flat.fits'
        statf = cur_date_str + '_wavesolution_stats.txt'

    # check stats for wavesolution
    wave_stats_ok = False
    wstatf = os.path.join(caldir, statf)
    if os.path.exists(wstatf):
        with open(wstatf) as infil:
            for line in infil:
                test = re.findall(r'AvgRMS:', line)
                if test:
                    avrms = float(line.split()[-1])
                    wave_stats_ok = (avrms < avrmslim)
        # Does wavelength solution pass?
        if wave_stats_ok:
            logging.info("Wavelength stats passed")
        else:
            # NO: move bad files away
            os.mkdir(os.path.join(caldir, 'bad'))
            os.system("mv %s_* %s" % (os.path.join(caldir, cur_date_str),
                                      os.path.join(caldir, 'bad')))
            logging.warning("Wavelength stats failed, moved cube to 'bad'")
            # Send email
            msg = EmailMessage()
            msg['To'] = "neill@srl.caltech.edu,"
            msg['Subject'] = "SEDM Error - Wave solution failure for %s" \
                             % cur_date_str
            msg['From'] = 'No_reply_sedm_robot@astro.caltech.edu'
            msg.set_content("Wavelength AvgRMS = %.1f > %.1f nm" % (avrms,
                                                                    avrmslim))
            # Send the message via local SMTP server.
            with smtplib.SMTP('smtp-server.astro.caltech.edu') as s:
                s.send_message(msg)
    # Do we have all the calibration files?
    ft = os.path.exists(os.path.join(caldir, tmf))
    ftm = os.path.exists(os.path.join(caldir, tmmf))
    fg = os.path.exists(os.path.join(caldir, hgf))
    fw = os.path.exists(os.path.join(caldir, wsf))
    ff = os.path.exists(os.path.join(caldir, fff))
    logging.info("Cals ready?: trace: %d, trace/mask: %d, grid: %d, wave: %d, "
                 "flat: %d" % (ft, ftm, fg, fw, ff))
    if ft and ftm and fg and fw and ff:
        ret = True

    return ret


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
    f1 = os.path.exists(os.path.join(caldir, 'bias1.0.fits'))
    logging.info("Biases ready?: bias0.1: %d, bias2.0: %d, bias1.0: %d" %
                 (fb, f2, f1))
    if (fb and f2) or f1:
        ret = True

    return ret


def cal_proc_ready(caldir='./', fsize=_nomfs, mintest=False, ncp=0,
                   test_cal_ims=False):
    """Check counts for all required raw cal file types in caldir directory.

    Args:
        caldir (str): directory where raw cal files reside
        fsize (int): size of completely copied file in bytes
        mintest (bool): test for minimum required number of files
        ncp (int): number of cal images most recently copied
        test_cal_ims (bool): test for presence of input cal images

    Returns:
        bool: True if required raw cal files are present, False otherwise

    """

    nbias = 0
    bias_done = False
    nbias2 = 0
    bias2_done = False
    nbias1 = 0
    bias1_done = False
    nxe = 0
    xe_done = False
    nhg = 0
    hg_done = False
    ncd = 0
    cd_done = False
    ndome = 0
    dome_done = False
    ret = False
    bias_done_str = '10 of 10'
    lamp_done_str = '5 of 5'

    if test_cal_ims:
        dof = glob.glob(os.path.join(caldir, 'dome.fits'))
        hgf = glob.glob(os.path.join(caldir, 'Hg.fits'))
        cdf = glob.glob(os.path.join(caldir, 'Cd.fits'))
        xef = glob.glob(os.path.join(caldir, 'Xe.fits'))
        if len(dof) == 1 and len(hgf) == 1 and len(cdf) == 1 and len(xef) == 1:
            ret = True

    else:

        # Get files in calibration directory
        cflist = sorted(glob.glob(os.path.join(caldir, 'ifu*.fits')))
        # Are there any files yet?
        if len(cflist) > 0:
            # Loop over files
            for cal in cflist:
                # Are we complete?
                if os.stat(cal).st_size >= fsize:
                    # Read FITS header
                    ff = pf.open(cal)
                    hdr = ff[0].header
                    ff.close()
                    # Get OBJECT keyword
                    try:
                        obj = hdr['OBJECT']
                    except KeyError:
                        obj = ''
                    # Get ADCSPEED keyword
                    speed = hdr['ADCSPEED']

                    # Check for calibration files
                    if 'Calib' in obj:
                        if 'bias' in obj:
                            if speed == 2.0:
                                nbias2 += 1
                                if bias_done_str in obj:
                                    bias2_done = True
                            if speed == 0.1:
                                nbias += 1
                                if bias_done_str in obj:
                                    bias_done = True
                            if speed == 1.0:
                                nbias1 += 1
                                if bias_done_str in obj:
                                    bias1_done = True
                        if 'Xe' in obj:
                            nxe += 1
                            if lamp_done_str in obj:
                                xe_done = True
                        if 'dome' in obj:
                            ndome += 1
                            if lamp_done_str in obj:
                                dome_done = True
                        if 'Hg' in obj:
                            nhg += 1
                            if lamp_done_str in obj:
                                hg_done = True
                        if 'Cd' in obj:
                            ncd += 1
                            if lamp_done_str in obj:
                                cd_done = True

            # Do we have the ideal number of calibration files?
            if ((((nbias2 >= 10 or bias2_done) and (nbias >= 10 or bias_done))
                    or (nbias1 >= 10 or bias1_done)) and
                    (nxe >= 5 or xe_done) and (ndome >= 5 or dome_done) and
                    (nhg >= 5 or hg_done) and (ncd >= 5 or cd_done)):
                ret = True
            # Do we have the minimum allowed number of calibration files?
            if mintest:
                if (((nbias2 >= 5 and nbias >= 5) or nbias1 >= 5) and
                        nxe >= 3 and ndome >= 3 and
                        nhg >= 3 and ncd >= 3):
                    ret = True
        logging.info("bias2.0: %d, bias0.1: %d, bias1.0: %d, dome: %d, "
                     "Xe: %d, Hg: %d, Cd: %d" %
                     (nbias2, nbias, nbias1, ndome, nxe, nhg, ncd))
        sys.stdout.flush()
        # Should we process biases?
        if ((nbias2 >= 10 and nbias >= 10) or nbias1 >= 10) and ncp > 0:
            proc_bias_crrs(ncp=ncp)

    return ret
    # END: cal_proc_ready


def docp(src, dest, onsky=True, verbose=False, skip_cals=False, nodb=False):
    """Low level copy from raw directory to redux directory.

    Checks for raw ifu files, while avoiding any test and focus images.
    Uses os.symlink to do the copying (linked to conserve disk space).

    Args:
        src (str): source file
        dest (str): destination file
        onsky (bool): test for dome conditions or not
        verbose (bool): output messages?
        skip_cals (bool): skip copying cal images?
        nodb (bool): don't update SEDM db

    Returns:
        (int, int, int): number of images linked, number of standard
                    star images linked, number of sci objects linked

    """
    # Record copies
    ncp = 0
    # Was a standard star observation copied?
    nstd = 0
    # Was a science object copied
    nobj = 0
    # Read FITS header
    ff = pf.open(src)
    hdr = ff[0].header
    ff.close()
    # Get OBJECT and DOMEST keywords
    try:
        obj = hdr['OBJECT']
    except KeyError:
        logging.warning("Could not find OBJECT keyword, setting to Test")
        obj = 'Test'
    try:
        dome = hdr['DOMEST']
    except KeyError:
        logging.warning("Could not find DOMEST keyword, setting to null")
        dome = ''
    # Check if dome conditions are not right
    if onsky and ('CLOSED' in dome or 'closed' in dome):
        if verbose:
            logging.warning('On sky and dome is closed, skipping %s' % src)
    # All other conditions are OK
    else:
        # Skip test and Focus images
        if skip_cals:
            copy_this = ('test' not in obj and 'Focus:' not in obj and
                         'STOW' not in obj and 'Test' not in obj and
                         'Calib:' not in obj)
        else:
            copy_this = ('test' not in obj and 'Focus:' not in obj and
                         'STOW' not in obj and 'Test' not in obj)

        if copy_this:
            # Symlink to save disk space
            os.symlink(src, dest)
            if 'STD-' in obj:
                nstd = 1
                logging.info("Standard %s linked to %s" % (obj, dest))
            else:
                nobj = 1
                logging.info('Target %s linked to %s' % (obj, dest))
            ncp = 1
            if 'REQ_ID' in hdr and 'OBJ_ID' in hdr:
                if not nodb:
                    # Record in database
                    obs_id = update_observation(src)
                    if obs_id > 0:
                        logging.info("SEDM db accepted observation at id %d"
                                     % obs_id)
                    else:
                        logging.warning("SEDM db rejected observation")
                else:
                    logging.warning("Not updating obs in SEDM db")
            else:
                logging.warning("Missing request and/or object ids,"
                                " no db update")
        # Report skipping and type
        else:
            if verbose and 'test' in obj:
                logging.info('test file %s not linked' % src)
            if verbose and 'Focus:' in obj:
                logging.info('Focus file %s not linked' % src)
            if verbose and 'Calib:' in obj:
                logging.info('calib file %s not linked' % src)

    return ncp, nstd, nobj
    # END: docp


def update_observation(input_fitsfile):
    """ Update the SEDM database observation table on minar
        by adding a new observation"""

    header_dict = {
        'object_id': 'OBJ_ID', 'request_id': 'REQ_ID', 'mjd': 'MJD_OBS',
        'airmass': 'AIRMASS', 'airmass_end': 'ENDAIR', 'exptime': 'EXPTIME',
        'lst': 'LST', 'ra': 'RA', 'dec': 'DEC', 'tel_az': 'TEL_AZ',
        'tel_el': 'TEL_EL', 'tel_pa': 'TEL_PA', 'ra_off': 'RA_OFF',
        'dec_off': 'DEC_OFF', 'imtype': 'IMGTYPE', 'camera': 'CAM_NAME',
        'filter': 'FILTER', 'parang': 'TEL_PA', 'parang_end': 'END_PA',
        'time_elapsed': 'ELAPTIME'
    }
    obs_dict = {
        'object_id': 0, 'request_id': 0, 'mjd': 0.,
        'airmass': 0., 'airmass_end': 0., 'exptime': 0.,
        'lst': ' ', 'ra': 0., 'dec': 0., 'tel_az': 0.,
        'tel_el': 0., 'tel_pa': 0., 'ra_off': 0.,
        'dec_off': 0., 'imtype': ' ', 'camera': ' ',
        'filter': ' ', 'parang': 0., 'parang_end': 0.,
        'time_elapsed': 0.,
        'fitsfile': input_fitsfile.split('/')[-1]
    }
    ff = pf.open(input_fitsfile)

    for key in header_dict.keys():
        hk = header_dict[key]
        if hk in ff[0].header:
            if key == 'dec':
                try:
                    obs_dict[key] = Angle(ff[0].header[hk]+' degrees').degree
                except ValueError:
                    logging.warning("Bad Dec kwd value: %s" % ff[0].header[hk])
                    obs_dict[key] = -99.0
            elif key == 'ra':
                try:
                    obs_dict[key] = Angle(ff[0].header[hk]+' hours').degree
                except ValueError:
                    logging.warning("Bad RA kwd value: %s" % ff[0].header[hk])
                    obs_dict[key] = -99.0
            else:
                obs_dict[key] = ff[0].header[hk]
        else:
            logging.warning("Header keyword not found: %s" % hk)
    ff.close()

    sedmdb = db.SedmDb.SedmDB()
    observation_id, status = sedmdb.add_observation(obs_dict)
    req_id, req_update_status = sedmdb.update_request(
        {'id': obs_dict['request_id'], 'status': 'COMPLETED'})
    logging.info("REQ_ID %d: %s" % (req_id, req_update_status))
    logging.info(status)
    return observation_id
    # END: update_observation


def update_calibration(utdate, src_dir=_reduxpath):
    """ Update the SEDM database spec_calib table on minar
        by adding a new spectral calibration"""

    spec_calib_dict = {}

    src = os.path.join(src_dir, utdate)
    if os.path.exists(src):

        dome_master = os.path.join(src, 'dome.fits')
        if os.path.exists(dome_master):
            spec_calib_dict['dome_master'] = dome_master
        else:
            logging.info("spec cal item not found: %s" % dome_master)

        bias_slow_master = os.path.join(src, 'bias0.1.fits')
        if os.path.exists(bias_slow_master):
            spec_calib_dict['bias_slow_master'] = bias_slow_master
        else:
            logging.info("spec cal item not found: %s" % bias_slow_master)

        bias_fast_master = os.path.join(src, 'bias2.0.fits')
        if os.path.exists(bias_fast_master):
            spec_calib_dict['bias_fast_master'] = bias_fast_master
        else:
            logging.info("spec cal item not found: %s" % bias_fast_master)

        bias_andor_master = os.path.join(src, 'bias1.0.fits')
        if os.path.exists(bias_andor_master):
            spec_calib_dict['bias_slow_master'] = bias_andor_master
        else:
            logging.info("spec cal item not found: %s" % bias_andor_master)

        flat = os.path.join(src, utdate + '_Flat.fits')
        if os.path.exists(flat):
            spec_calib_dict['flat'] = flat
        else:
            logging.info("spec cal item not found: %s" % flat)

        hg_master = os.path.join(src, 'Hg.fits')
        if os.path.exists(hg_master):
            spec_calib_dict['hg_master'] = hg_master
        else:
            logging.info("spec cal item not found: %s" % hg_master)

        xe_master = os.path.join(src, 'Xe.fits')
        if os.path.exists(xe_master):
            spec_calib_dict['xe_master'] = xe_master
        else:
            logging.info("spec cal item not found: %s" % xe_master)

        cd_master = os.path.join(src, 'Cd.fits')
        if os.path.exists(cd_master):
            spec_calib_dict['cd_master'] = cd_master
        else:
            logging.info("spec cal item not found: %s" % cd_master)

        hexagrid = os.path.join(src, utdate + '_HexaGrid.pkl')
        if os.path.exists(hexagrid):
            spec_calib_dict['hexagrid'] = hexagrid
        else:
            logging.info("spec cal item not found: %s" % hexagrid)

        tracematch = os.path.join(src, utdate + '_TraceMatch.pkl')
        if os.path.exists(tracematch):
            spec_calib_dict['tracematch'] = tracematch
        else:
            logging.info("spec cal item not found: %s" % tracematch)

        tracematch_withmasks = os.path.join(
            src, utdate + '_TraceMatch_WithMasks.pkl')
        if os.path.exists(tracematch_withmasks):
            spec_calib_dict['tracematch_withmasks'] = tracematch_withmasks
        else:
            logging.info("spec cal item not found: %s" % tracematch_withmasks)

        wavesolution = os.path.join(src, utdate + '_WaveSolution.pkl')
        if os.path.exists(wavesolution):
            spec_calib_dict['wavesolution'] = wavesolution
        else:
            logging.info("spec cal item not found: %s" % wavesolution)

        dispersionmap = os.path.join(
            src, utdate + '_wavesolution_dispersionmap.png')
        if os.path.exists(dispersionmap):
            spec_calib_dict['dispersionmap'] = dispersionmap
        else:
            logging.info("spec cal item not found: %s" % dispersionmap)

        flatmap = os.path.join(src, utdate + '_flat3d.png')
        if os.path.exists(flatmap):
            spec_calib_dict['flatmap'] = flatmap
        else:
            logging.info("spec cal item not found: %s" % flatmap)

        wstats = os.path.join(src, utdate + '_wavesolution_stats.txt')
        if os.path.exists(wstats):
            with open(wstats) as sf:
                stat_line = sf.readline()
                while stat_line:
                    if 'NSpax' in stat_line:
                        spec_calib_dict['nspaxels'] = int(stat_line.split()[-1])
                        logging.info("Wave NSpax = %d" %
                                     spec_calib_dict['nspaxels'])
                    elif 'MinRMS' in stat_line:
                        spec_calib_dict['wave_rms_min'] = \
                            float(stat_line.split()[-1])
                    elif 'AvgRMS' in stat_line:
                        spec_calib_dict['wave_rms_avg'] =\
                            float(stat_line.split()[-1])
                    elif 'MaxRMS' in stat_line:
                        spec_calib_dict['wave_rms_max'] =\
                            float(stat_line.split()[-1])
                    stat_line = sf.readline()
        else:
            logging.info("spec cal item not found: %s" % wstats)

        dstats = os.path.join(src, utdate + '_dome_stats.txt')
        if os.path.exists(dstats):
            with open(dstats) as sf:
                stat_line = sf.readline()
                while stat_line:
                    if 'NSpax' in stat_line:
                        nspax = int(stat_line.split()[-1])
                        logging.info("Dome NSpax = %d" % nspax)
                        if spec_calib_dict['nspaxels'] == 0:
                            spec_calib_dict['nspaxels'] = nspax
                    elif 'MinWid' in stat_line:
                        spec_calib_dict['width_rms_min'] = \
                            float(stat_line.split()[-1])
                    elif 'AvgWid' in stat_line:
                        spec_calib_dict['width_rms_avg'] = \
                            float(stat_line.split()[-1])
                    elif 'MaxWid' in stat_line:
                        spec_calib_dict['width_rms_max'] = \
                            float(stat_line.split()[-1])
                    stat_line = sf.readline()
        else:
            logging.info("spec cal item not found: %s" % dstats)

    else:
        logging.warning("Source dir does not exist: %s" % src)

    spec_calib_dict['utdate'] = utdate
    spec_calib_dict['drpver'] = drp_ver

    sedmdb = db.SedmDb.SedmDB()
    spec_calib_id, status = sedmdb.add_spec_calib(spec_calib_dict)
    logging.info(status)
    return spec_calib_id
    # END: update_calibration


def proc_bias_crrs(ncp=1, piggyback=False):
    """Process biases and CR rejection steps.

    Args:
        ncp (int): number of images to process
        piggyback (bool): are we using another script to process data?

    Returns:
        bool: True if processing was successful, otherwise False

    """

    # Default return value
    ret = False
    if piggyback:
        ret = True
    else:
        # Get new listing
        retcode = subprocess.call("~/spy what ifu*.fits > what.list",
                                  shell=True)
        if retcode == 0:
            # Generate new Makefile
            retcode = subprocess.call("~/spy plan ifu*.fits", shell=True)
            if retcode == 0:
                # Make bias + bias subtraction
                retcode = subprocess.call(("make", "-j", "16", "bias"))
                if retcode != 0:
                    logging.warning("bias failed, try again")
                    retcode = subprocess.call(("make", "bias"))
                if retcode == 0:
                    # Make CR rejection
                    retcode = subprocess.call(("make", "-j", "8", "crrs"))
                    if retcode != 0:
                        logging.warning("crrs failed, try again")
                        retcode = subprocess.call(("make", "-j", "8", "crrs"))
                    # Success on all fronts!
                    if retcode == 0:
                        logging.info("bias, crrs processed for %d new images"
                                     % ncp)
                        ret = True
                    # Report failures
                    else:
                        logging.error("could not make crrs")
                else:
                    logging.error("could not make bias")
            else:
                logging.error("could not make plan")
        else:
            logging.error("could not make what.list")

    return ret
    # END: proc_bias_crrs


def cpsci(srcdir, destdir='./', fsize=_nomfs, datestr=None, nodb=False):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        srcdir (str): source directory (typically in /data/sedmdrp/raw)
        destdir (str): destination directory (typically in /data/sedmdrp/redux)
        fsize (int): size of completely copied file in bytes
        datestr (str): YYYYMMDD date string
        nodb (bool): skip update of SEDM db

    Returns:
        int: Number of ifu images actually copied

    """

    # Get files in destination directory
    dflist = sorted(glob.glob(os.path.join(destdir, 'ifu*.fits')))
    # Record copies and standard star observations
    ncp = 0
    nstd = 0
    nobj = 0
    copied = []
    stds = []
    sciobj = []
    # Get list of source files
    srcfiles = sorted(glob.glob(os.path.join(srcdir, 'ifu*.fits')))
    # Loop over source files
    for fl in srcfiles:
        # get base filename
        fn = fl.split('/')[-1]
        # Is our source file complete?
        if os.stat(fl).st_size >= fsize:
            # has it been previously copied?
            prev = [s for s in dflist if fn in s]
            # No? then copy the file
            if len(prev) == 0:
                # Call copy
                nc, ns, nob = docp(fl, destdir + '/' + fn, skip_cals=True,
                                   nodb=nodb, verbose=True)
                if nc >= 1:
                    copied.append(fn)
                if ns >= 1:
                    stds.append(fn)
                if nob >= 1:
                    sciobj.append(fn)
                # Record copies
                ncp += nc
                nstd += ns
                nobj += nob
    # We copied files
    logging.info("Linked %d files" % ncp)
    # Do bias subtraction, CR rejection
    if ncp > 0:
        if not proc_bias_crrs(ncp):
            logging.error("Error processing bias/crrs")
        if datestr is None:
            logging.error("Illegal datestr parameter")
            return 0, None

    return ncp, copied
    # END: cpsci


def dosci(destdir='./', datestr=None, local=False, nodb=False,
          nopush_marshal=False, nopush_slack=False, oldext=False):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        destdir (str): destination directory (typically in /data/sedmdrp/redux)
        datestr (str): YYYYMMDD date string
        local (bool): set to skip pushing to marshal and slack
        nodb (bool): if True no update to SEDM db
        nopush_marshal (bool): True if no update to marshal
        nopush_slack (bool): True if no update to slack
        oldext (bool): True to use extract_star.py instead of extractstar.py

    Returns:
        int: Number of ifu images actually copied

    """
    # don't make guider movie if we are local
    guider_movie = not local
    # Record copies and standard star observations
    ncp = 0
    copied = []
    # Get list of source files in destination directory
    srcfiles = sorted(glob.glob(os.path.join(destdir, 'crr_b_ifu*.fits')))
    # Loop over source files
    for fl in srcfiles:
        # get base filename
        fn = fl.split('/')[-1]
        procfn = 'spec*auto*' + fn.split('.')[0] + '*.fits'
        proced = glob.glob(os.path.join(destdir, procfn))
        # Is our source file processed?
        if len(proced) == 0:
            # Read FITS header
            ff = pf.open(fl)
            hdr = ff[0].header
            ff.close()
            # Get OBJECT keyword
            try:
                obj = hdr['OBJECT'].replace(" [A]", "").replace(" ", "-")
            except KeyError:
                logging.warning(
                    "Could not find OBJECT keyword, setting to Test")
                obj = 'Test'
            # Get DOMEST keyword
            try:
                dome = hdr['DOMEST']
            except KeyError:
                logging.warning(
                    "Could not find DOMEST keyword, settting to null")
                dome = ''
            # skip Cal files
            if 'Calib:' in obj:
                continue
            # skip if dome closed
            if 'CLOSED' in dome or 'closed' in dome:
                continue
            # make finder
            make_finder(fl)
            # record action
            copied.append(fn)
            ncp += 1
            # are we a standard star?
            if 'STD-' in obj:
                e3d_good = make_e3d(fnam=fl, destdir=destdir, datestr=datestr,
                                    nodb=nodb, sci=False, hdr=None,
                                    guider_movie=guider_movie)
                if e3d_good:
                    # Get seeing
                    seeing = rcimg.get_seeing(imfile=fn, destdir=destdir,
                                              save_fig=True)
                    # Use auto psf extraction for standard stars
                    if seeing > 0:
                        logging.info("seeing measured as %f" % seeing)
                    else:
                        logging.info("seeing not measured for %s" % fn)
                    if not oldext:
                        cmd = ("extractstar.py", datestr, "--auto", fn,
                               "--std", "--tag", "robot",
                               "--centroid", "brightest", "--seeing", "2.0")
                        # "--seeing", "%.2f" % seeing)
                    else:
                        logging.info("Old extraction method used")
                        cmd = ("extract_star.py", datestr, "--auto", fn,
                               "--std", "--tag", "robot", "--maxpos")
                    logging.info("Extracting std star spectra for " + fn)
                    logging.info(" ".join(cmd))
                    retcode = subprocess.call(cmd)
                    if retcode != 0:
                        logging.error("Error extracting std star spectra for "
                                      + fn)
                        badfn = "spec_auto_notfluxcal_" + fn.split('.')[0] + \
                                "_failed.fits"
                        cmd = ("touch", badfn)
                        subprocess.call(cmd)
                    else:
                        if local or nopush_slack:
                            cmd = ("pysedm_report.py", datestr, "--contains",
                                   fn.split('.')[0])
                        else:
                            cmd = ("pysedm_report.py", datestr, "--contains",
                                   fn.split('.')[0], "--slack")
                        logging.info(" ".join(cmd))
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            logging.error("Error running pysedm_report for " +
                                          fn.split('.')[0])
                        # run Verify.py
                        cmd = "~/sedmpy/drpifu/Verify.py %s --contains %s" % \
                              (datestr, fn.split('.')[0])
                        logging.info(cmd)
                        subprocess.call(cmd, shell=True)
                        # run make report
                        cmd = ("make", "report")
                        logging.info(" ".join(cmd))
                        subprocess.call(cmd)
                        # check if extraction succeeded
                        proced = glob.glob(os.path.join(destdir, procfn))[0]
                        if os.path.exists(proced):
                            if nodb:
                                logging.warning("Not updating spec in SEDM db")
                            else:
                                # Update SedmDb table spec
                                spec_id = update_spec(proced)
                                if spec_id > 0:
                                    logging.info("update of %s with spec_id %d"
                                                 % (proced, spec_id))
                                else:
                                    logging.warning("failed to update spec %s" %
                                                    proced)
                        else:
                            logging.error("Not found: %s" % proced)
                        # Did we generate a flux calibration?
                        flxcal = glob.glob(
                            os.path.join(destdir,
                                         "fluxcal_auto_robot_lstep1__%s_*.fits"
                                         % fn.split('.')[0]))
                        if flxcal:
                            # Generate effective area and efficiency plots
                            cmd = "~/sedmpy/drpifu/Eff.py %s --contains %s" % \
                                  (datestr, fn.split('.')[0])
                            logging.info(cmd)
                            subprocess.call(cmd, shell=True)
                        else:
                            logging.info("No flux calibration generated")
                else:
                    logging.error("Cannot perform extraction for %s" % fn)
            else:
                # Build cube for science observation
                e3d_good = make_e3d(fnam=fl, destdir=destdir, datestr=datestr,
                                    nodb=nodb, sci=True, hdr=hdr,
                                    guider_movie=guider_movie)

                if e3d_good:
                    # Get seeing
                    seeing = rcimg.get_seeing(imfile=fn, destdir=destdir,
                                              save_fig=True)
                    # Use forced psf for science targets
                    if seeing > 0:
                        logging.info("seeing measured as %f" % seeing)

                    else:
                        logging.info("seeing not measured for %s" % fn)
                    if not oldext:
                        cmd = ("extractstar.py", datestr, "--auto", fn,
                               "--autobins", "6", "--tag", "robot",
                               "--centroid", "auto", "--byecr",
                               "--seeing", "2.0")
                        # "--seeing", "%.2f" % seeing)
                    else:
                        logging.info("Old extraction method used")
                        cmd = ("extract_star.py", datestr, "--auto", fn,
                               "--autobins", "6", "--tag", "robot")
                    logging.info("Extracting object spectra for " + fn)
                    logging.info(" ".join(cmd))
                    retcode = subprocess.call(cmd)
                    if retcode != 0:
                        logging.error("Error extracting object spectrum for "
                                      + fn)
                        badfn = "spec_auto_notfluxcal_" + fn.split('.')[0] + \
                                "_failed.fits"
                        cmd = ("touch", badfn)
                        subprocess.call(cmd)
                    else:
                        # Run SNID, SNIascore, and NGSF
                        logging.info("Running SNID, SNIascore, NGSF for " + fn)
                        cmd = ("make", "classify")
                        logging.info(" ".join(cmd))
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            logging.error("Error running SNID, SNIascore, "
                                          "or NGSF")
                        if local or nopush_slack:
                            cmd = ("pysedm_report.py", datestr, "--contains",
                                   fn.split('.')[0])
                        else:
                            cmd = ("pysedm_report.py", datestr, "--contains",
                                   fn.split('.')[0], "--slack")
                        logging.info(" ".join(cmd))
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            logging.error("Error running report for " +
                                          fn.split('.')[0])
                        # Upload spectrum to marshal
                        if local or nopush_marshal:
                            logging.warning("nopush_marshal or local: "
                                            "skipping ztfupload")
                        else:
                            # fritz upload
                            cmd = ("make", "fritzupload")
                            retcode = subprocess.call(cmd)
                            if retcode != 0:
                                logging.error("Error uploading spectra to"
                                              " fritz marshal")
                            # growth upload
                            cmd = ("make", "ztfupload")
                            retcode = subprocess.call(cmd)
                            if retcode != 0:
                                logging.error("Error uploading spectra to"
                                              " growth marshal")
                        # run Verify.py
                        cmd = "~/sedmpy/drpifu/Verify.py %s --contains %s" % \
                              (datestr, fn.split('.')[0])
                        subprocess.call(cmd, shell=True)
                        # notify user that followup successfully completed
                        proced = glob.glob(os.path.join(destdir, procfn))[0]
                        if os.path.exists(proced):
                            if local or nopush_marshal:
                                logging.warning("nopush_marshal or local: "
                                                "skipping email")
                            else:
                                email_user(proced, datestr, obj)
                            if nodb:
                                logging.warning("No update of spec in SEDM db")
                            else:
                                # Update SedmDb table spec
                                spec_id = update_spec(
                                    proced, nopush_marshal=nopush_marshal)
                                logging.info("update of %s with spec_id %d" %
                                             (proced, spec_id))
                        else:
                            logging.error("Not found: %s" % proced)
                        # contsep extraction
                        cmd = ("extractstar.py", datestr, "--auto", fn,
                               "--autobins", "6", "--tag", "contsep",
                               "--centroid", "auto", "--contsep", "--byecr")
                        logging.info("Extracting contsep spectra for " + fn)
                        logging.info(" ".join(cmd))
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            logging.error("Error extracting contsep spectrum "
                                          "for" + fn)
                        else:
                            # Run SNID
                            logging.info("Running SNID for contsep " + fn)
                            cmd = ("make", "classify")
                            logging.info(" ".join(cmd))
                            retcode = subprocess.call(cmd)
                            if retcode != 0:
                                logging.error("Error running SNID")
                            # run Verify.py
                            cmd = "~/sedmpy/drpifu/Verify.py %s " \
                                  "--contains contsep_lstep1__%s" \
                                  % (datestr, fn.split('.')[0])
                            subprocess.call(cmd, shell=True)
                            # run pysedm_report
                            if local or nopush_slack:
                                cmd = ("pysedm_report.py", datestr,
                                       "--contains",
                                       "contsep_lstep1__" + fn.split('.')[0])
                            else:
                                cmd = ("pysedm_report.py", datestr,
                                       "--contains",
                                       "contsep_lstep1__" + fn.split('.')[0],
                                       "--slack")
                            logging.info(" ".join(cmd))
                            retcode = subprocess.call(cmd)
                            if retcode != 0:
                                logging.error("Error running report for " +
                                              "contsep_lstep1__" +
                                              fn.split('.')[0])
                else:
                    logging.error("Cannot perform extraction for %s" % fn)
    return ncp, copied
    # END: dosci


def make_e3d(fnam=None, destdir=None, datestr=None, nodb=False, sci=False,
             hdr=None, guider_movie=False):
    """ Make the e3d cube"""
    cube_good = False
    if fnam is None:
        logging.error("Need input file name")
        fn = None
    else:
        fn = fnam.split('/')[-1]
    if destdir is None:
        logging.error("Need destination dir")
    if datestr is None:
        logging.error("Need date string")
    if sci:
        lab = "Science"
    else:
        lab = "STD"
    cmd = ("ccd_to_cube.py", datestr, "--build", fn, "--solvewcs",
           "--ncore", "8")
    if hdr:
        # Check for moving target: no guider image for those
        if 'RA_RATE' in hdr and 'DEC_RATE' in hdr:
            if hdr['RA_RATE'] != 0. or hdr['DEC_RATE'] != 0.:
                logging.info("Non-sidereal object")
                cmd = ("ccd_to_cube.py", datestr, "--build", fn,
                       "--ncore", "8")

    proccubefn = "e3d_%s_*.fits" % fn.split('.')[0]
    procedcube = glob.glob(os.path.join(destdir, proccubefn))
    if len(procedcube) == 0:
        # Make guider movie?
        if guider_movie:
            subprocess.call("$SEDMPY/bin/guide_movie.py -i %s" % fnam,
                            shell=True)
        # Build cube for observation
        logging.info("Building %s cube for %s" % (lab, fn))
        logging.info(" ".join(cmd))
        retcode = subprocess.call(cmd)
        # Check results
        if retcode != 0:
            logging.error("Error generating cube for " + fn)
            # Do this to prevent constant re-try
            badfn = "spec_auto_notfluxcal_" + fn.split('.')[0] + \
                    "_failed.fits"
            cmd = ("touch", badfn)
            subprocess.call(cmd)
        else:
            logging.info("%s cube generated for %s" % (lab, fn))
            cube_good = True
            if nodb:
                logging.warning("Not updating cube in SEDM db")
            else:
                # Update SedmDb cube table
                cube_id = update_cube(fnam)
                if cube_id > 0:
                    logging.info("SEDM db accepted cube at id %d"
                                 % cube_id)
                else:
                    logging.warning("SEDM db rejected cube")
    else:
        logging.info("%s cube already exists for %s" % (lab, fn))
        cube_good = True
    return cube_good
    # END: make_e3d


def update_spec(input_specfile, update_db=False, nopush_marshal=False):
    """ Update the SEDM database on minar by adding a new spec entry"""

    header_dict = {
        'imgset': 'IMGSET', 'quality': 'QUALITY', 'cubefile': 'SOURCE',
        'reducer': 'REDUCER', 'airmass': 'AIRMASS',
        'atmcorr': 'ATMSCALE', 'pos_ok': 'POSOK', 'srcpos': 'SRCPOS',
        'pos_x_spax': 'XPOS', 'pos_y_spax': 'YPOS', 'psf_model': 'PSFMODEL',
        'psf_fwhm': 'PSFFWHM', 'psf_ell': 'PSFELL', 'psf_adr_pa': 'PSFADRPA',
        'psf_adr_z': 'PSFADRZ', 'psf_adr_c2': 'PSFADRC2',
        'fluxcal': 'FLUXCAL', 'fluxcalfile': 'CALSRC', 'extr_type': 'EXTRTYPE'
    }
    spec_dict = {
        'spec_calib_id': 0, 'observation_id': 0, 'asciifile': '', 'npyfile': '',
        'fitsfile': '', 'imgset': '', 'quality': 0, 'cubefile': '',
        'standardfile': '', 'marshal_spec_id': 0, 'skysub': False,
        'fwhm': 0., 'background': 0., 'line_fwhm': 0., 'extract_x': 0.,
        'extract_y': 0., 'extract_pa': 0., 'extract_a': 0., 'extract_b': 0.,
        'ad_red': 0., 'ad_blue': 0., 'prlltc': 0., 'reducer': '', 'airmass': 0.,
        'atmcorr': 0., 'pos_ok': False, 'srcpos': '', 'pos_x_spax': 0.,
        'pos_y_spax': 0., 'psf_model': '', 'psf_fwhm': 0., 'psf_ell': 0.,
        'psf_adr_pa': 0., 'psf_adr_z': 0., 'psf_adr_c2': 0., 'fluxcal': False,
        'fluxcalfile': '', 'extr_type': '', 'cube_id': 0
    }
    class_header_dict = {
        'classification': 'SNIDTYPE', 'redshift': 'SNIDZMED',
        'redshift_err': 'SNIDZERR', 'score': 'SNIDRLAP', 'phase': 'SNIDAMED',
        'phase_err': 'SNIDAERR', 'class_template': 'SNIDTEMP'
    }
    class_ia_header_dict = {
        'score': 'SNIASCOR', 'score_err': 'SNIASCER',
        'redshift': 'SNIASCZ', 'redshift_err': 'SNIASCZE'
    }
    class_ngsf_header_dict = {
        'classification': 'NGSFTYPE', 'redshift': 'NGSFZ', 'score': 'NGSFCHI2',
        'phase': 'NGSFPHAS', 'class_template': 'NGSFTEMP'
    }
    class_dict = {
        'spec_id': 0, 'object_id': 0, 'classification': '', 'auto': True,
        'redshift': 0., 'redshift_err': 0., 'classifier': 'SNID', 'score': 0.,
        'phase': 0., 'phase_err': 0., 'score_type': 'RLAP',
        'class_source': '', 'class_template': ''
    }
    class_ia_dict = {
        'spec_id': 0, 'object_id': 0, 'classification': 'SNIa', 'auto': True,
        'redshift': 0., 'redshift_err': 0., 'classifier': 'SNIascore',
        'score': 0., 'score_err': 0., 'score_type': 'SNIa',
        'class_source': ''
    }
    class_ngsf_dict = {
        'spec_id': 0, 'object_id': 0, 'classification': '', 'auto': True,
        'redshift': 0., 'clasifier': 'NGSF', 'score': 0.,
        'score_type': 'Chi2/dof', 'class_source': '', 'class_template': ''
    }

    # Get utdate
    indir = '/'.join(input_specfile.split('/')[:-1])
    utdate = indir.split('/')[-1]

    # Read header
    ff = pf.open(input_specfile)

    # Get object name
    objnam = ff[0].header['OBJECT'].split()[0]
    class_dict['class_source'] = objnam
    class_ia_dict['class_source'] = objnam

    # Open database connection
    sedmdb = db.SedmDb.SedmDB()

    # Get marshal number
    marsh = None
    if 'REQ_ID' in ff[0].header:
        request_id = ff[0].header['REQ_ID']
        # Search for target in the database
        try:
            res = sedmdb.get_from_request(["external_id"],
                                          {"id": request_id})[0]
            marsh = res[0]
        except IndexError:
            print("Unable to retrieve marshal number from database")
        except TypeError:
            print("Not a valid REQ_ID")

    # Get header keyword values
    for key in header_dict.keys():
        hk = header_dict[key]
        if hk in ff[0].header:
            spec_dict[key] = ff[0].header[hk]
        else:
            if 'PSFELL' in hk and 'PSFAB' in ff[0].header:
                try:
                    b_a = ff[0].header['PSFAB']
                    smaja = ff[0].header['PSFFWHM']
                    smina = smaja * b_a
                    spec_dict['psf_ell'] = (smaja - smina)/smaja
                except TypeError:
                    logging.warning("PSF values contain NaN")
            else:
                logging.warning("Header keyword not found: %s" % hk)

    # Check for classification info
    # SNID results
    good_class = False
    if 'SNIDTYPE' in ff[0].header:
        if 'NONE' in ff[0].header['SNIDTYPE']:
            logging.info("SNID was unable to type %s" % input_specfile)
        else:
            good_class = True
            for key in class_header_dict.keys():
                hk = class_header_dict[key]
                if hk in ff[0].header:
                    class_dict[key] = ff[0].header[hk]
                else:
                    logging.warning("Header keyword not found: %s" % hk)
            if 'SNIDSUBT' in ff[0].header:
                if '-' not in ff[0].header['SNIDSUBT']:
                    class_dict['classification'] = \
                        class_dict['classification'] + ' ' + \
                        ff[0].header['SNIDSUBT']
    else:
        logging.info("No SNID info in %s" % input_specfile)
    # SNIascore results
    good_ia_class = False
    if 'SNIASCOR' in ff[0].header:
        if ff[0].header['SNIASCOR'] < 0.:
            logging.info("SNIascore was unable to score %s" % input_specfile)
        else:
            good_ia_class = True
            for key in class_ia_header_dict.keys():
                hk = class_ia_header_dict[key]
                if hk in ff[0].header:
                    class_ia_dict[key] = ff[0].header[hk]
                else:
                    logging.warning("Header keyword not found: %s" % hk)
    else:
        logging.info("No SNIascore info in %s" % input_specfile)
    # NGSF results
    good_ngsf_class = False
    if 'NGSFTYPE' in ff[0].header:
        if 'NONE' in ff[0].header['NGSFTYPE']:
            logging.info("NGSF was unable to type %s" % input_specfile)
        else:
            good_ngsf_class = True
            for key in class_ngsf_header_dict.keys():
                hk = class_ngsf_header_dict[key]
                if hk in ff[0].header:
                    class_ngsf_dict[key] = ff[0].header[hk]
                else:
                    logging.warning("Header keyword not found: %s" % hk)
            if 'NGSFSUBT' in ff[0].header:
                class_dict['classification'] = \
                    class_dict['classification'] + ' ' + \
                    ff[0].header['NGSFSUBT']
    else:
        logging.info("No NGSF info in %s" % input_specfile)

    ff.close()

    # Add fitsfile
    spec_dict['fitsfile'] = input_specfile
    # Add asciifile
    spec_dict['asciifile'] = input_specfile.split('.fit')[0] + '.txt'

    # Get marshal spec id
    if 'STD' not in objnam and spec_dict['quality'] <= 2 and marsh != 2:
        srcid = None
        specid = None
        if nopush_marshal:
            logging.warning("nopush_marshal: skipping marshal retrieval")
        else:
            srcid, specid = mc.get_missing_info(objnam, utdate, srcid, specid)
        if srcid is None:
            logging.info("Not an object on the growth marshal: %s" % objnam)
        else:
            if specid is None:
                logging.info("No spectrum found on the growth marshal: %s, %s"
                             % (utdate, objnam))
            else:
                spec_dict['marshal_spec_id'] = specid

    # Check if we've already added this spectrum
    search_fits = input_specfile.replace('+', '\+')
    spec_id = sedmdb.get_from_spec(['id'], {'fitsfile': search_fits},
                                   {'fitsfile': '~'})
    if spec_id:
        logging.info("Spectrum already in db: %s" % input_specfile)
        if update_db:
            logging.info("Updating from %s" % input_specfile)
            spec_dict['id'] = spec_id[0][0]
        else:
            return spec_id[0][0]

    # Get observation and object ids
    ifufile = 'ifu' + '_'.join(
        input_specfile.split('_ifu')[-1].split('_')[:4]) + '.fits'
    observation_id = sedmdb.get_from_observation(['id', 'object_id'],
                                                 {'fitsfile': ifufile},
                                                 {'fitsfile': '~'})
    if observation_id:
        spec_dict['observation_id'] = observation_id[0][0]
        class_dict['object_id'] = observation_id[0][1]
        class_ia_dict['object_id'] = observation_id[0][1]
    else:
        logging.error("No observation_id for %s" % ifufile)
        return -1

    # Get cube id
    cube_id = sedmdb.get_from_cube(['id'], {'observation_id':
                                            spec_dict['observation_id']})
    if cube_id:
        spec_dict['cube_id'] = cube_id[0][0]
    else:
        logging.warning("No cube id for %s" % input_specfile)

    # Get spec_calib id for this utdate
    spec_calib_id = sedmdb.get_from_spec_calib(['id'], {'utdate': utdate})
    if spec_calib_id:
        spec_dict['spec_calib_id'] = spec_calib_id[0][0]

        # Add into database
        spec_id, status = sedmdb.add_spec(spec_dict, update=update_db)
        # update classification
        class_dict['spec_id'] = spec_id
        class_ia_dict['spec_id'] = spec_id
        logging.info(status)

        if good_class:
            class_id, cstatus = sedmdb.add_classification(class_dict)
            if class_id < 0 and update:
                logging.info("SNID Classification already exists for this spec")
            else:
                logging.info("SNID Classification accepted with id %d,"
                             " and status %s" % (class_id, cstatus))
        else:
            logging.info("No SNID classification found in input spectrum")

        if good_ia_class:
            class_ia_id, ciastatus = sedmdb.add_classification(class_ia_dict)
            if class_ia_id < 0 and update:
                logging.info("SNIascore record already exists for this spec")
            else:
                logging.info("SNIascore record accepted with id %d,"
                             " and status %s" % (class_ia_id, ciastatus))
        else:
            logging.info("No SNIascore record found in input spectrum")

        if good_ngsf_class:
            class_ngsf_id, cngsfstatus = sedmdb.add_classification(
                class_ngsf_dict)
            if class_ngsf_id < 0 and update:
                logging.info("NGSF record already exists for this spec")
            else:
                logging.info("NGSF record accepted with id %d,"
                             " and status %s" % (class_ngsf_id, cngsfstatus))
        else:
            logging.info("No NGSF record found in input spectrum")

        return spec_id
    else:
        logging.error("ERROR: no spec_calib_id found for %s" % utdate)
        return -1
    # END: update_spec


def update_cube(input_fitsfile):
    """ Update the SEDM database on minar by adding a new cube entry. """

    header_dict = {
        'ccd_x_flex_corr': 'IFLXCORR', 'ccd_x_flex_px': 'CCDIFLX',
        'ccd_y_flex_corr': 'JFLXCORR', 'ccd_y_flex_px': 'CCDJFLX',
        'atm_corr': 'ATMCORR', 'atm_source': 'ATMSRC',
        'atm_mean_corr': 'ATMSCALE'
    }
    cube_dict = {
        'observation_id': 0,
        'ccd_x_flex_corr': False, 'ccd_x_flex_px': 0.,
        'ccd_y_flex_corr': False, 'ccd_y_flex_px': 0.,
        'atm_corr': False, 'atm_source': '',
        'atm_mean_corr': 1.,
        'spec_calib_id': 0
    }

    # Get cube file
    root = input_fitsfile.split('/')[-1].split('.fits')[0]
    indir = '/'.join(input_fitsfile.split('/')[:-1])
    utdate = indir.split('/')[-1]
    cube_list = glob.glob(os.path.join(indir, 'e3d_'+root+'_*.fits'))

    # Check list
    if len(cube_list) != 1:
        logging.error("ERROR: Ambiguous cube list")
        return -1

    # Read header
    ff = pf.open(cube_list[0])

    # Get header keyword values
    for key in header_dict.keys():
        hk = header_dict[key]
        if hk in ff[0].header:
            cube_dict[key] = ff[0].header[hk]
        else:
            logging.warning("Header keyword not found: %s" % hk)
    ff.close()

    # Open database connection
    sedmdb = db.SedmDb.SedmDB()

    # Get observation id
    fitsfile = 'ifu'+input_fitsfile.split('_ifu')[-1]
    observation_id = sedmdb.get_from_observation(['id'],
                                                 {'fitsfile': fitsfile},
                                                 {'fitsfile': '~'})
    if observation_id:
        cube_dict['observation_id'] = observation_id[0][0]

    # Get spec_calib id for this utdate
    spec_calib_id = sedmdb.get_from_spec_calib(['id'], {'utdate': utdate})
    if spec_calib_id:
        cube_dict['spec_calib_id'] = spec_calib_id[0][0]

        # Add into database
        cube_id, status = sedmdb.add_cube(cube_dict)
        logging.info(status)
        return cube_id
    else:
        logging.error("ERROR: no spec_calib_id found for %s" % utdate)
        return -1
    # END: update_cube


def email_user(spec_file, utdate, object_name):
    """ Send e-mail to requestor indicating followup completed"""
    # get request id
    ff = pf.open(spec_file)
    if 'REQ_ID' in ff[0].header:
        request = int(ff[0].header['REQ_ID'])
    else:
        logging.error("No keyword REQ_ID found in %s\nNo email sent"
                      % spec_file)
        return
    if 'redo' in spec_file:
        quality = 1
    else:
        if 'QUALITY' in ff[0].header:
            quality = int(ff[0].header['QUALITY'])
        else:
            logging.error("No keyword QUALITY found in %s\nNo email sent"
                          % spec_file)
            return
    if quality == 1:
        status = 'IFU auto-extraction has been manually recovered.'
    elif quality == 3:
        status = 'IFU auto-extraction failed: target outside IFU.'
    elif quality == 4:
        status = 'IFU auto-extraction failed: > 20% of flux is negative.'
    elif quality == 5:
        status = 'IFU auto-extraction failed: guider astrometry failure;\n no' \
                 'guarantee that the target is in the IFU or well exposed.\n' \
                 'Manual recovery may be possible.'
    else:
        status = 'IFU extraction succeeded.'
    link = 'http://minar.caltech.edu/data_access/ifu?obsdate=%s' % utdate
    subj = 'SEDM followup status report for %s on %s' % (object_name, utdate)

    sedmdb = db.SedmDb.SedmDB()
    sedmdb.send_email_by_request(requestid=request, template='report_send',
                                 subject=subj,
                                 template_dict={'object_name': object_name,
                                                'link': link, 'status': status})


def make_finder(ffile):
    """ Spawn a process that makes a finder for the ifu file """
    spy = os.path.join(os.getenv("HOME"), "spy")
    prog = os.path.join(os.getenv("HOME"), "sedmpy", "drpifu", "acq_finder.py")
    cmd = (spy, prog, "--imfile", ffile)
    subprocess.Popen(cmd)


def find_recent(redd, fname, destdir, dstr):
    """Find the most recent version of fname and copy it to destdir.

    Look through sorted list of redux directories to find most recent
    version of the input file.  Copy (link) it into the destination directory.

    Args:
        redd (str): reduced directory (something like /data/sedmdrp/redux)
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
        srtlist = sorted([d for d in glob.glob(fspec)
                          if os.path.isdir(d)])
        redlist = srtlist[0:srtlist.index(destdir)]
        logging.info("Looking backwards for %s starting at %s" %
                     (fname, redlist[-1]))
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
        redd (str): reduced directory (something like /data/sedmdrp/redux)
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
        srtlist = sorted([d for d in glob.glob(fspec)
                          if os.path.isdir(d)])
        redlist = srtlist[0:srtlist.index(destdir)]
        logging.info("Looking backwards for %s starting at %s" %
                     (fname, redlist[-1]))
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


def find_recent_fluxcal(redd, fname, destdir):
    """Find the most recent version of fname and copy it to destdir.

    Look through sorted list of redux directories to find most recent
    version of the input file.  Copy it to the destination directory.

    Args:
        redd (str): reduced directory (something like /data/sedmdrp/redux)
        fname (str): what file to look for
        destdir (str): where the file should go

    Returns:
        bool: True if file found and copied, False otherwise.

    """

    # Default return value
    ret = False
    # Make sure the file doesn't already exist in destdir
    local_file = glob.glob(os.path.join(destdir, fname))
    if len(local_file) >= 1:
        logging.warning("%s already exists in %s" % (fname, destdir))
        ret = True
    # Search in redd for file
    else:
        # Get all but the most recent reduced data directories
        fspec = os.path.join(redd, '20??????')
        srtlist = sorted([d for d in glob.glob(fspec)
                          if os.path.isdir(d)])
        redlist = srtlist[0:srtlist.index(destdir)]
        logging.info("Looking backwards for %s starting at %s" %
                     (fname, redlist[-1]))
        # Go back in reduced dir list until we find our file
        for d in reversed(redlist):
            src = sorted(glob.glob(os.path.join(d, fname)))
            for s in src:
                # Skip sym-links
                if os.path.islink(s):
                    continue
                # Read FITS header
                ff = pf.open(s)
                hdr = ff[0].header
                ff.close()
                # Skip if not Telluric corrected
                if 'TELLFLTR' not in hdr:
                    continue
                newfile = os.path.join(destdir, s.split('/')[-1])
                try:
                    os.symlink(s, newfile)
                except OSError:
                    logging.warning("File already exists: %s" % newfile)
                ret = True
                logging.info("Found %s in directory %s, linking to %s" %
                             (fname, d, newfile))
                break
            if ret:
                break
    if not ret:
        logging.warning("%s not found" % fname)
    return ret


def cpprecal(dirlist, destdir='./', fsize=_nomfs, nodb=False):
    """Copy raw cal files from previous date's directory

    Make sure we only look in previous day directory for files created
    within four hours of the day changeover for possible raw calibration
    files required for the current night's calibration.  Copy any such
    files into the destination directory.

    Args:
        dirlist (list): a list of raw dirs (typically in /data/sedmdrp/raw)
        destdir (str): where to put the files
        fsize (int): size of completely copied file in bytes
        nodb (bool): skip update of SEDM db

    Returns:
        int: number of images actually copied

    """

    # Get current and previous dates
    # current source dir
    cdate = dirlist[-1].split('/')[-1]
    # convert to JD
    ctime = Time(cdate[0:4]+'-'+cdate[4:6]+'-'+cdate[6:])
    cjd = ctime.jd
    #
    # previous source dir
    pdate = dirlist[-2].split('/')[-1]
    # convert to JD
    ptime = Time(pdate[0:4]+'-'+pdate[4:6]+'-'+pdate[6:])
    pjd = ptime.jd
    # Record how many images copied
    ncp = 0
    # If there is a previous night, get those files
    if (int(cjd) - int(pjd)) <= 1:
        # Set the previous night as the source directory
        srcdir = dirlist[-2]
        # Get list of previous night's raw cal files
        # (within four hours of day changeover)
        fspec = os.path.join(srcdir, "ifu%s_2*.fits" % pdate)
        flist = sorted(glob.glob(fspec))
        # Loop over file list
        for src in flist:
            if os.stat(src).st_size >= fsize:
                # Read FITS header
                ff = pf.open(src)
                hdr = ff[0].header
                ff.close()
                # Get OBJECT keyword
                obj = hdr['OBJECT']
                # Filter Calibs and avoid test images
                if 'Calib' in obj and 'of' in obj and 'test' not in obj and \
                        'Test' not in obj:
                    # Copy cal images
                    imf = src.split('/')[-1]
                    destfil = os.path.join(destdir, imf)
                    exptime = hdr['EXPTIME']
                    lampcur = hdr['LAMPCUR']
                    # Check for dome exposures
                    if 'dome' in obj:
                        if exptime >= 60. and ('dome' in obj and
                                               'Xe' not in obj and
                                               'Hg' not in obj and 
                                               'Cd' not in obj):
                            if lampcur > 0.0:
                                # Copy dome images
                                if not os.path.exists(destfil):
                                    nc, ns, nob = docp(src, destfil,
                                                       onsky=False,
                                                       verbose=True,
                                                       nodb=nodb)
                                    ncp += nc
                            else:
                                logging.warning("Bad dome - lamp not on: %s"
                                                % src)
                    # Check for arcs
                    elif 'Xe' in obj or 'Cd' in obj or 'Hg' in obj:
                        if exptime > 25.:
                            # Copy arc images
                            if not os.path.exists(destfil):
                                nc, ns, nob = docp(src, destfil, onsky=False,
                                                   verbose=True, nodb=nodb)
                                ncp += nc
                    # Check for biases
                    elif 'bias' in obj:
                        if exptime <= 0.:
                            # Copy bias images
                            if not os.path.exists(destfil):
                                nc, ns, nob = docp(src, destfil, onsky=False,
                                                   verbose=True, nodb=nodb)
                                ncp += nc
            else:
                logging.warning("Truncated file: %s" % src)

    return ncp
    # END: cpprecal


def cpcal(srcdir, destdir='./', fsize=_nomfs, nodb=False):
    """Copy raw cal files from srcdir into destdir.

    Find calibration files taken within 10 hours of the day changeover
    and copy them to the destination directory.

    Args:
        srcdir (str): source for raw cal images
        destdir (str): place to put the cal images
        fsize (int): size of completely copied file in bytes
        nodb (bool): skip update of SEDM db

    Returns:
        int: number of images actually copied

    """

    # Get current date
    sdate = srcdir.split('/')[-1]
    # Get list of current raw calibration files
    # (within 10 hours of day changeover)
    fspec = os.path.join(srcdir, "ifu%s_0*.fits" % sdate)
    flist = sorted(glob.glob(fspec))
    # Record number copied
    ncp = 0
    # Loop over file list
    for src in flist:
        # Get destination filename
        imf = src.split('/')[-1]
        destfil = os.path.join(destdir, imf)
        # Does our local file already exist?
        if glob.glob(destfil):
            # Get the size to compare with source
            loc_size = os.stat(destfil).st_size
        else:
            # Doesn't yet exist locally
            loc_size = 0
        # Get source size to compare
        src_size = os.stat(src).st_size
        # Copy only if source complete or larger than local file
        if src_size >= fsize and src_size > loc_size:
            # Read FITS header
            ff = pf.open(src)
            hdr = ff[0].header
            ff.close()
            # Get OBJECT keyword
            try:
                obj = hdr['OBJECT']
            except KeyError:
                obj = ''
            # Filter Calibs and avoid test images and be sure it is part of
            # a series.
            if 'Calib' in obj and 'of' in obj and 'test' not in obj and \
                    'Test' not in obj:
                exptime = hdr['EXPTIME']
                lampcur = hdr['LAMPCUR']
                # Check for dome exposures
                if 'dome' in obj:
                    if exptime > 30. and ('dome' in obj and
                                          'Xe' not in obj and
                                          'Hg' not in obj and
                                          'Cd' not in obj):
                        if lampcur > 0.0:
                            # Copy dome images
                            nc, ns, nob = docp(src, destfil, onsky=False,
                                               verbose=True, nodb=nodb)
                            ncp += nc
                        else:
                            logging.warning("Bad dome - lamp not on: %s" % src)
                # Check for arcs
                elif 'Xe' in obj or 'Cd' in obj or 'Hg' in obj:
                    if exptime > 15.:
                        # Copy arc images
                        nc, ns, nob = docp(src, destfil, onsky=False,
                                           verbose=True, nodb=nodb)
                        ncp += nc
                # Check for biases
                elif 'bias' in obj:
                    if exptime <= 0.:
                        # Copy bias images
                        nc, ns, nob = docp(src, destfil, onsky=False,
                                           verbose=True, nodb=nodb)
                        ncp += nc

    return ncp
    # END: cpcal


def obs_loop(rawlist=None, redd=None, check_precal=True, indir=None,
             piggyback=False, local=False, nodb=False, use_refcube=False,
             nopush_marshal=False, nopush_slack=False, oldext=False):
    """One night observing loop: processes calibrations and science data

    Copy raw cal files until we are ready to process the night's
    calibrations.  When ready, process them.  If we get to UT = 3:00,
    the sun has set and we then retrieve the processed cal files from
    the most recent night.  Next, enter a loop that waits for new ifu
    images, copies them and performs the basic bias removal and CR
    rejection processing.  If there are standard star observations,
    re-generate the standard star calibration.  If we get no new ifu
    files and we get to UT = 15:00, then the sun is up and we exit
    the loop.

    Args:
        rawlist (list): list of raw data dirs (usually in /data/sedmdrp/raw)
        redd (str): reduced directory (something like /data/sedmdrp/redux)
        check_precal (bool): should we check for images from previous night?
        indir (str): input directory for single night processing
        piggyback (bool): basic processing done by StartObs.py
        local (bool): True if no marshal/slack update required
        nodb (bool): True if no update to SEDM db
        use_refcube (bool): True to use ref traces from 20230407
        nopush_marshal (bool): True if no update to marshal
        nopush_slack (bool): True if no update to slack
        oldext (bool): True to use extract_star.py instead of extracstar.py

    Returns:
        bool: True if night completed normally, False otherwise

    Note:
        KeyboardInterrupt handler exits gracefully with a ctrl-C.

    """
    # Set up Observatory params and night limiting epochs
    p60 = astroplan.Observer.at_site(sedm_cfg['observatory']['name'])
    obstime = Time(datetime.utcnow())
    evening_civil_twilight = p60.twilight_evening_civil(obstime,
                                                        which='nearest')
    morning_civil_twilight = p60.twilight_morning_civil(obstime, which='next')
    print("evening civil twilight: ", evening_civil_twilight.iso)
    print("morning civil twilight: ", morning_civil_twilight.iso)

    # Source directory is most recent raw dir
    if indir is None:
        srcdir = rawlist[-1]
    else:
        srcdir = indir
    # Default return value
    ret = False
    # Output directory is based on source dir
    outdir = os.path.join(redd, srcdir.split('/')[-1])
    # Current date string
    cur_date_str = str(outdir.split('/')[-1])
    # Do we have a new directory?  This tells us we are observing tonight
    if not os.path.exists(outdir):
        # Make it
        os.mkdir(outdir)
    # Go there
    os.chdir(outdir)
    if indir is not None:
        fl = glob.glob("ifu*.fits")
        if len(fl) > 0:
            # Generate new Makefile
            retcode = subprocess.call("~/spy plan ifu*.fits", shell=True)
            if retcode != 0:
                logging.warning("Error making plan in %s" % indir)
        else:
            logging.warning("No fits files in %s yet, so no plan made" % indir)
    # report
    logging.info("Raw files from  : %s" % srcdir)
    logging.info("Reduced files to: %s" % outdir)

    # Check if processed cal files are ready
    if not cube_ready(outdir, cur_date_str):
        # Wait for cal files until sunset
        if piggyback:
            logging.info("Skipping check for raw cal files")
            ncp = 0
        else:
            if check_precal:
                # Copy raw cal files from previous date directory
                npre = cpprecal(rawlist, outdir, nodb=nodb)
                logging.info("Linked %d raw cal files from %s" % (npre,
                                                                  rawlist[-2]))
            # Now check the current source dir for raw cal files
            ncp = cpcal(srcdir, outdir, nodb=nodb)
            logging.info("Linked %d raw cal files from %s" % (ncp, srcdir))
        # Now loop until we have the raw cal files we need or sun is down
        while not cal_proc_ready(outdir, ncp=ncp, test_cal_ims=piggyback):
            # Wait a minute
            logging.info("waiting 60s for more raw cal files...")
            now = Time(datetime.utcnow())
            time.sleep(60)
            if piggyback:
                logging.info("checking for processed cal files")
                ncp = 0
            else:
                if check_precal and now.to_datetime().hour >= 20:
                    logging.info("checking %s for new raw cal files..."
                                 % rawlist[-2])
                    ncp = cpprecal(rawlist, outdir, nodb=nodb)
                    logging.info("Linked %d raw cal files from %s"
                                 % (ncp, rawlist[-2]))
                else:
                    logging.info("checking %s for new raw cal files..."
                                 % srcdir)
                    ncp = cpcal(srcdir, outdir, nodb=nodb)
                logging.info("Linked %d raw cal files from %s" % (ncp, srcdir))
            if ncp <= 0:
                # Check to see if we are still before an hour after sunset
                now = Time(datetime.utcnow())
                if now < evening_civil_twilight:
                    logging.info("UT  = %s < civil twilight (%s),"
                                 " so keep  waiting" %
                                 (now.iso.split()[-1],
                                  evening_civil_twilight.iso.split()[-1]))
                else:
                    logging.info("UT = %s >= civil twilight (%s), "
                                 "time to get a cal set" %
                                 (now.iso.split()[-1],
                                  evening_civil_twilight.iso.split()[-1]))
                    break
            else:
                # Get new listing
                retcode = subprocess.call("~/spy what ifu*.fits > what.list",
                                          shell=True)
                # Link what.txt
                if not os.path.islink(os.path.join('what.txt')):
                    os.symlink('what.list', 'what.txt')
                if retcode != 0:
                    logging.error("what oops!")

        # Process calibrations if we are using them
        if cal_proc_ready(outdir, mintest=True, test_cal_ims=piggyback):
            # bias subtract and CR reject
            start_time = time.time()
            if proc_bias_crrs(20, piggyback=piggyback):
                procb_time = int(time.time() - start_time)
                if not piggyback:
                    # Make cal images
                    subprocess.call(("make", "calimgs"))
                # Process calibration
                start_time = time.time()
                if use_refcube:
                    link_refcube(curdir=outdir, date_str=cur_date_str)
                    logging.info("linked Traces from ref dir into %s" %
                                 cur_date_str)
                else:
                    cmd = ("ccd_to_cube.py", cur_date_str, "--tracematch",
                           "--hexagrid")
                    logging.info(" ".join(cmd))
                    subprocess.call(cmd)
                procg_time = int(time.time() - start_time)
                if os.path.exists(
                   os.path.join(outdir, cur_date_str + '_HexaGrid.pkl')):
                    # Process wavelengths
                    start_time = time.time()
                    # Spawn nsub sub-processes to solve wavelengths faster
                    nsub = 8
                    cmd = ("derive_wavesolution.py", cur_date_str,
                           "--nsub", "%d" % nsub)
                    logging.info(" ".join(cmd))
                    subprocess.Popen(cmd)
                    time.sleep(60)
                    # Get a list of solved spaxels
                    wslist = glob.glob(os.path.join(outdir, cur_date_str +
                                                    '_WaveSolution_range*.pkl'))
                    # Wait until they are all finished
                    nfin = len(wslist)
                    while nfin < nsub:
                        time.sleep(60)
                        wslist = glob.glob(
                            os.path.join(outdir, cur_date_str +
                                         '_WaveSolution_range*.pkl'))
                        if len(wslist) != nfin:
                            print("\nFinished %d out of %d parts"
                                  % (len(wslist), nsub))
                            nfin = len(wslist)
                        else:
                            print(".", end="", flush=True)
                    logging.info("Finished all %d parts, merging..." % nsub)
                    # Merge the solutions
                    subprocess.call(("derive_wavesolution.py", cur_date_str,
                                     "--merge"))
                procw_time = int(time.time() - start_time)
                if os.path.exists(
                   os.path.join(outdir, cur_date_str + '_WaveSolution.pkl')):
                    # Process flat
                    start_time = time.time()
                    cmd = ("ccd_to_cube.py", cur_date_str, "--flat")
                    logging.info(" ".join(cmd))
                    subprocess.call(cmd)
                    if not (os.path.exists(
                            os.path.join(outdir, cur_date_str + '_Flat.fits'))):
                        logging.info("Making of %s_Flat.fits failed!"
                                     % cur_date_str)
                else:
                    logging.error("Making of %s cube failed!" % cur_date_str)
                procf_time = int(time.time() - start_time)
                # Report times
                logging.info("Calibration processing took "
                             "%d s (bias,crrs), %d s (grid),"
                             "%d s (waves),  and %d s (flat)" %
                             (procb_time, procg_time, procw_time, procf_time))
                # Make cube report
                cmd = "~/miniconda3/bin/python " \
                      "~/sedmpy/drpifu/CubeReport.py %s" % cur_date_str
                if not local:
                    cmd += " --slack"   # send to slack pysedm_report channel
                logging.info(cmd)
                subprocess.call(cmd, shell=True)
        # Check status
        if cube_ready(outdir, cur_date_str):
            if nodb:
                logging.warning("Not updating SEDM db")
            else:
                # Update spec_calib table in sedmdb
                spec_calib_id = update_calibration(cur_date_str)
                logging.info("SEDM db accepted spec_calib at id %d" %
                             spec_calib_id)
        else:
            logging.error("These calibrations failed!")
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
                nc1 = find_recent_bias(redd, 'bias1.0.fits', outdir)
            else:
                ncb = True
                nc2 = True
                nc1 = True
            # Check for bias failure
            if not ncb or not nc2:
                if not nc1:
                    msg = "Calibration stage biases failed: bias0.1 = %s, " \
                          "bias2.0 = %s, bias1.0 = %s, stopping" % (ncb,
                                                                    nc2, nc1)
                    sys.exit(msg)
                else:
                    logging.info("Using Andor single speed biases")
            # Check for geom failure
            if not nct or not nctm or not ncg or not ncw or not ncf:
                msg = "Calibration stage geom failed: trace = %s, " \
                      "trace/mask = %s, grid = %s, wave = %s, flat = %s, "\
                      "stopping" % (nct, nctm, ncg, ncw, ncf)
                sys.exit(msg)
            # If we get here, we are done
            logging.info("Using older calibration files")
    else:
        logging.info("Calibrations already present in %s" % outdir)

    logging.info("Calibration stage complete, ready for science!")
    # Link recent flux cal file
    find_recent_fluxcal(redd, 'fluxcal*.fits', outdir)
    # Keep track of no copy
    nnc = 0
    # loop and copy new files
    doit = True
    try:
        while doit:
            # Wait for morning civil twilight

            # Wait a minute
            logging.info("waiting 60s for new ifu images...")
            sys.stdout.flush()
            time.sleep(60)
            # Check for new ifu images
            logging.info("checking %s for new ifu images..." % srcdir)
            sys.stdout.flush()
            # Record starting time for new file processing
            start_time = time.time()
            if piggyback:
                nsci, science = dosci(outdir, datestr=cur_date_str,
                                      local=local, nodb=nodb,
                                      nopush_marshal=nopush_marshal,
                                      nopush_slack=nopush_slack, oldext=oldext)
                ncp = nsci
            else:
                ncp, copied = cpsci(srcdir, outdir, datestr=cur_date_str,
                                    nodb=nodb)
                nsci, science = dosci(outdir, datestr=cur_date_str,
                                      local=local, nodb=nodb,
                                      nopush_marshal=nopush_marshal,
                                      nopush_slack=nopush_slack, oldext=oldext)
            # We copied some new ones so report processing time
            if ncp > 0:
                proc_time = int(time.time() - start_time)
                logging.info("%d new ifu images copied and %d processed "
                             "in %d s" % (ncp, nsci, proc_time))
                sys.stdout.flush()
                nnc = 0
            else:
                nnc += 1
            # Have we been waiting for a while?
            if nnc > 3:
                # Check time
                now = Time(datetime.utcnow())
                if now >= morning_civil_twilight:
                    # No new observations but civil twilight has begun
                    logging.info("No new images for %d minutes and UT = "
                                 "%s > %s so twilight has begun!" %
                                 (nnc, now.iso.split()[-1],
                                  morning_civil_twilight.iso.split()[-1]))
                    logging.info(
                        "Time to wait until we have a new raw directory")
                    doit = False
                    # Normal termination
                    subprocess.call(("make", "report"))
                    ret = True
                else:
                    logging.info("No new image for %d minutes but UT = "
                                 "%s <= %s, so civil twilight has not started, "
                                 "keep waiting" %
                                 (nnc, now.iso.split()[-1],
                                  morning_civil_twilight.iso.split()[-1]))
                if indir is not None:
                    logging.info("Done processing images from %s", indir)
                    doit = False
                    subprocess.call(("make", "report"))
                    ret = True

    # Handle a ctrl-C
    except KeyboardInterrupt:
        sys.exit("Exiting")

    return ret
    # END: obs_loop


def update(red_dir=_reduxpath, ut_dir=None):
    """Update re-extracted spectra"""

    # Check inputs
    if not ut_dir:
        return
    # First update ztf marshal
    cmd = ("make", "ztfupload")
    retcode = subprocess.call(cmd)
    if retcode:
        logging.warning("Not all spectra uploaded to marshal")
    # Now update sedmdb
    flist = glob.glob(os.path.join(red_dir, ut_dir, "spec*.fits"))
    for file in flist:
        if 'failed' in file:
            logging.info("Skipping failed extraction: %s" % file)
            continue
        update_spec(file)


def clean_post_redux(outdir, utdstr):
    """Remove/compress unused files after reduction"""
    ndel = 0
    # Remove raw file links
    flist = glob.glob(os.path.join(outdir, "ifu%s_*.fits" % utdstr))
    flist.extend(glob.glob(os.path.join(outdir, 'rc%s_*.fits' % utdstr)))
    for fl in flist:
        if os.path.islink(fl):
            os.remove(fl)
            ndel += 1
    # Remove intermediate processing files
    flist = glob.glob(os.path.join(outdir, "b_ifu*.fits"))
    for fl in flist:
        os.remove(fl)
        ndel += 1
    # Remove intermediate calib files, compress others
    n_gzip = 0
    flist = glob.glob(os.path.join(outdir, "*crr_b_ifu*.fits"))
    for fl in flist:
        if 'failed' in fl:
            continue
        ff = pf.open(fl)
        hdr = ff[0].header
        ff.close()
        try:
            obj = hdr['OBJECT']
        except KeyError:
            obj = 'Test'
        if 'Calib' in obj:
            os.remove(fl)
            ndel += 1
        else:
            rute = fl.split('/')[-1]
            if rute.startswith('maskcrr_') or rute.startswith('crr_b_ifu') or \
                    rute.startswith('bkgd_crr_b_ifu') or \
                    rute.startswith('forcepsf') or \
                    rute.startswith('guider_crr_b'):
                subprocess.call(["gzip", fl])
                n_gzip += 1
    # Compress calib files
    flist = glob.glob(os.path.join(outdir, "*dome.fits"))
    flist.extend(glob.glob(os.path.join(outdir, '??.fits')))
    flist.extend(glob.glob(os.path.join(outdir, 'bias*.fits')))
    # Compress rainbow cam images
    flist.extend(glob.glob(os.path.join(outdir, 'rc*.fits')))
    for fl in flist:
        if os.path.islink(fl):
            continue
        subprocess.call(["gzip", fl])
        n_gzip += 1
    # append dir to backup file
    back_file = sedm_cfg['backup']['redux_backup_file']
    try:
        with open(back_file, 'w') as bf:
            bf.writelines(utdstr + "\n")
        print("%s written to %s, ready for rsync" % (utdstr, back_file))
    except OSError:
        print("Cannot open backup file for update: %s" % back_file)

    return ndel, n_gzip


def go(rawd=_rawpath, redd=_reduxpath, wait=False, use_refcube=False,
       check_precal=True, indate=None, piggyback=False, local=False,
       nopush_marshal=False, nopush_slack=False, nodb=False, oldext=False):
    """Outermost infinite loop that watches for a new raw directory.

    Keep a list of raw directories in `redd` and fire off
    the obs_loop procedure when a new directory appears.  Check for
    a new raw directory every 10 minutes.

    Args:
        rawd (str): raw directory, should be like /data/sedmdrp/raw
        redd (str): reduced directory, should be like /data/sedmdrp/redux
        wait (bool): wait for new directory, else start right away
        use_refcube (bool): Use reference traces from 20230407
        check_precal (bool): should we check previous night for cals?
        indate (str): input date to process: YYYYMMDD (e.g. 20180626)
        piggyback (bool): True if using other script to copy data
        local (bool): True if no marshal/slack update required
        nopush_marshal (bool): True if no marshal update required
        nopush_slack (bool): True if no slack update required
        nodb (bool): True if no update of SEDM Db
        oldext (bool): True to use extract_star.py instead of extractstar.py

    Returns:
        None

    Note:
        KeyboardInterrupt handler exits gracefully with a ctrl-C.

    """

    # Infinite loop
    dobs = True
    stat = True
    # Keep track of iterations
    its = 0
    # Get all raw directories
    fspec = os.path.join(rawd, '20??????')
    rawlist = sorted([d for d in glob.glob(fspec) if os.path.isdir(d)])
    nraw = len(rawlist)
    logging.info("Found %d raw directories in %s: putting reduced data in %s" %
                 (nraw, rawd, redd))
    if indate is None:
        logging.info("Latest raw directory is %s" % rawlist[-1])

        if not wait:
            stat = obs_loop(rawlist, redd, check_precal=check_precal,
                            piggyback=piggyback, local=local, nodb=nodb,
                            nopush_marshal=nopush_marshal,
                            use_refcube=use_refcube, nopush_slack=nopush_slack,
                            oldext=oldext)
            its += 1
            logging.info("Finished SEDM observing iteration %d in raw dir %s" %
                         (its, rawlist[-1]))
        try:
            while dobs:
                if stat:
                    logging.info("Now we wait until we get a new raw directory")
                    waiting = True
                    while waiting:
                        logging.info("waiting 10min for new raw directory...")
                        sys.stdout.flush()
                        time.sleep(600)
                        # Get all raw directories
                        new_rawlist = sorted([d for d in glob.glob(fspec)
                                              if os.path.isdir(d)])
                        new_nraw = len(new_rawlist)
                        if new_nraw > nraw:
                            waiting = False
                            sys.stdout.flush()
                            rawlist = new_rawlist
                            nraw = new_nraw
                            logging.info("Starting next SEDM observing "
                                         "iteration with raw dir %s" %
                                         rawlist[-1])
                        else:
                            now = datetime.utcnow()
                            logging.info("UT = %02d:%02d No new directories "
                                         "yet, so keep waiting" %
                                         (now.hour, now.minute))
                            sys.stdout.flush()
                logging.info("Found %d raw directories in %s: "
                             "putting reduced data in %s" % (nraw, rawd, redd))
                logging.info("Latest raw directory is %s" % rawlist[-1])
                stat = obs_loop(rawlist, redd, check_precal=check_precal,
                                piggyback=piggyback, local=local, nodb=nodb,
                                nopush_marshal=nopush_marshal,
                                nopush_slack=nopush_slack,
                                use_refcube=use_refcube, oldext=oldext)
                its += 1
                logging.info("Finished SEDM observing iteration %d in "
                             "raw dir %s" % (its, rawlist[-1]))
        # Handle a ctrl-C
        except KeyboardInterrupt:
            sys.exit("Exiting")
    else:
        indir = os.path.join(rawd, indate)
        logging.info("Processing raw data from %s" % indir)
        stat = obs_loop(rawlist, redd, check_precal=check_precal, indir=indir,
                        piggyback=piggyback, local=local, nodb=nodb,
                        nopush_marshal=nopush_marshal,
                        nopush_slack=nopush_slack,
                        use_refcube=use_refcube, oldext=oldext)
        its += 1
        logging.info("Finished SEDM processing in raw dir %s with status %d" %
                     (indir, stat))
    # END: go


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Start SEDM pipeline

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--rawdir', type=str, default=_rawpath,
                        help='Input raw directory (%s)' % _rawpath)
    parser.add_argument('--reduxdir', type=str, default=_reduxpath,
                        help='Output reduced directory (%s)' % _reduxpath)
    parser.add_argument('--wait', action="store_true", default=False,
                        help='Wait for new directory first')
    parser.add_argument('--piggyback', action="store_true", default=False,
                        help='Do not copy data, copied by another script')
    parser.add_argument('--skip_precal', action="store_true", default=False,
                        help='Skip check of previous day for cal files?')
    parser.add_argument('--date', type=str, default=None,
                        help='Select date to process')
    parser.add_argument('--update', type=str, default=None,
                        help='UTDate directory to update')
    parser.add_argument('--local', action="store_true", default=False,
                        help='Process data locally only (no push to marshal or '
                             'slack)')
    parser.add_argument('--nopush_marshal', action="store_true", default=False,
                        help='Do not push to marshal')
    parser.add_argument('--nopush_slack', action="store_true", default=False,
                        help='Do not push to slack')
    parser.add_argument('--nodb', action="store_true", default=False,
                        help='Do not update SEDM Db')
    parser.add_argument('--oldext', action="store_true", default=False,
                        help='Use extract_star.py instead of extractstar.py')
    parser.add_argument('--clean', action="store_true", default=False,
                        help='Clean UTDate directory')
    parser.add_argument('--use_refcube', action="store_true", default=False,
                        help="Use reference traces from 20230407")

    args = parser.parse_args()

    if args.update:
        update(red_dir=args.reduxdir, ut_dir=args.update)
    elif args.clean:
        if args.date is not None:
            odir = os.path.join(args.reduxdir, args.date)
            print("Cleaning %s" % odir)
            nrm, ngzip = clean_post_redux(odir, args.date)
            print("%d intermediate files removed, %d files gzipped" %
                  (nrm, ngzip))
        else:
            print("Error: Must provide a UTDate to clean with --date")
    else:
        if args.local:
            arg_nodb = True
            arg_nopush_slack = True
            arg_nopush_marshal = True
        else:
            arg_nodb = args.nodb
            arg_nopush_slack = args.nopush_slack
            arg_nopush_marshal = args.nopush_marshal
        go(rawd=args.rawdir, redd=args.reduxdir, wait=args.wait,
           check_precal=(not args.skip_precal), indate=args.date,
           piggyback=args.piggyback, local=args.local, nodb=arg_nodb,
           nopush_marshal=arg_nopush_marshal, nopush_slack=arg_nopush_slack,
           oldext=args.oldext, use_refcube=args.use_refcube)
