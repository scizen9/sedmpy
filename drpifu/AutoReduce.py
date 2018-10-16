"""Conduct automatic reduction of SEDM data in sedmdrp@pharos

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

        usage: StartObs.py [-h] [--rawdir RAWDIR] [--reduxdir REDUXDIR]

        optional arguments:
          -h, --help           show this help message and exit
          --rawdir RAWDIR      Input raw directory (/scr2/sedm/raw)
          --reduxdir REDUXDIR  Output reduced directory (/scr2/sedm/redux)

"""
import time
import glob
import sys
import os
import subprocess
import astropy.io.fits as pf
import argparse
import ephem
import db.SedmDb

from astropy.time import Time
from astropy.coordinates import Angle

try:
    import Version
except ImportError:
    import drpifu.Version as Version

drp_ver = Version.ifu_drp_version()


def cube_ready(caldir='./', cur_date_str=None):
    """Check for all required calibration files in calibration directory.

    Args:
        caldir (str): directory to check
        cur_date_str (str): current date in YYYYMMDD format

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
    else:
        tmf = cur_date_str + '_TraceMatch.pkl'
        tmmf = cur_date_str + '_TraceMatch_WithMasks.pkl'
        hgf = cur_date_str + '_HexaGrid.pkl'
        wsf = cur_date_str + '_WaveSolution.pkl'
        fff = cur_date_str + '_Flat.fits'

    # Do we have all the calibration files?
    ft = os.path.exists(os.path.join(caldir, tmf))
    ftm = os.path.exists(os.path.join(caldir, tmmf))
    fg = os.path.exists(os.path.join(caldir, hgf))
    fw = os.path.exists(os.path.join(caldir, wsf))
    ff = os.path.exists(os.path.join(caldir, fff))
    print("Cals ready?: trace: %d, trace/mask: %d, grid: %d, wave: %d, "
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
    print("Biases ready?: bias0.1: %d, bias2.0: %d" % (fb, f2))
    if fb and f2:
        ret = True

    return ret


def cal_proc_ready(caldir='./', fsize=8400960, mintest=False, ncp=0,
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
                    f = pf.open(cal)
                    hdr = f[0].header
                    f.close()
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
            if ((nbias2 >= 10 or bias2_done) and (nbias >= 10 or bias_done) and
                    (nxe >= 5 or xe_done) and (ndome >= 5 or dome_done) and
                    (nhg >= 5 or hg_done) and (ncd >= 5 or cd_done)):
                ret = True
            # Do we have the minimum allowed number of calibration files?
            if mintest:
                if (nbias2 >= 5 and nbias >= 5 and nxe >= 3 and ndome >= 3 and
                        nhg >= 3 and ncd >= 3):
                    ret = True
        print("bias2.0: %d, bias0.1: %d, dome: %d, Xe: %d, Hg: %d, Cd: %d" %
              (nbias2, nbias, ndome, nxe, nhg, ncd))
        sys.stdout.flush()
        # Should we process biases?
        if nbias2 >= 10 and nbias >= 10 and ncp > 0:
            proc_bias_crrs(ncp=ncp)

    return ret
    # END: cal_proc_ready


def docp(src, dest, onsky=True, verbose=False):
    """Low level copy from raw directory to redux directory.

    Checks for raw ifu files, while avoiding any test and focus images.
    Uses os.symlink to do the copying (linked to conserve disk space).

    Args:
        src (str): source file
        dest (str): destination file
        onsky (bool): test for dome conditions or not
        verbose (bool): print messages?

    Returns:
        (int, int, int): number of images linked, number of standard
                    star images linked, number of sci objects linked

    """

    # Read FITS header
    f = pf.open(src)
    hdr = f[0].header
    f.close()
    # Get OBJECT and DOMEST keywords
    obj = hdr['OBJECT']
    dome = hdr['DOMEST']
    # Record copies
    ncp = 0
    # Was a standard star observation copied?
    nstd = 0
    # Was a science object copied
    nobj = 0
    # Check if dome conditions are not right
    if onsky and ('CLOSED' in dome or 'closed' in dome):
        if verbose:
            print('On sky and dome is closed, skipping %s' % src)
    # All other conditions are OK
    else:
        # Skip test and Focus images
        if 'test' not in obj and 'Focus:' not in obj and 'STOW' not in obj and \
                'Test' not in obj:
            # Symlink to save disk space
            os.symlink(src, dest)
            if 'STD-' in obj:
                nstd = 1
                print("Standard %s linked to %s" % (obj, dest))
            else:
                nobj = 1
                print('Target %s linked to %s' % (obj, dest))
            ncp = 1
            # Record in database
            obs_id = update_observation(src)
            if obs_id > 0:
                print("SEDM db accepted observation at id %d" % obs_id)
            else:
                print("SEDM db rejected observation")
        # Report skipping and type
        else:
            if verbose and 'test' in hdr['OBJECT']:
                print('test file %s not linked' % src)
            if verbose and 'Focus:' in hdr['OBJECT']:
                print('Focus file %s not linked' % src)

    return ncp, nstd, nobj
    # END: docp


def update_observation(input_fitsfile):
    """ Update the SEDM database observation table on pharos
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
            if key is 'dec':
                obs_dict[key] = Angle(ff[0].header[hk]+' degrees').degree
            elif key is 'ra':
                obs_dict[key] = Angle(ff[0].header[hk]+' hours').degree
            else:
                obs_dict[key] = ff[0].header[hk]
        else:
            print("Header keyword not found: %s" % hk)
    ff.close()

    sedmdb = db.SedmDb.SedmDB()
    observation_id, status = sedmdb.add_observation(obs_dict)
    print(status)
    return observation_id
    # END: update_observation


def update_calibration(utdate, src_dir='/scr2/sedmdrp/redux'):
    """ Update the SEDM database spec_calib table on pharos
        by adding a new spectral calibration"""

    spec_calib_dict = {}

    src = os.path.join(src_dir, utdate)
    if os.path.exists(src):

        dome_master = os.path.join(src, 'dome.fits')
        if os.path.exists(dome_master):
            spec_calib_dict['dome_master'] = dome_master
        else:
            print("spec cal item not found: %s" % dome_master)

        bias_slow_master = os.path.join(src, 'bias0.1.fits')
        if os.path.exists(bias_slow_master):
            spec_calib_dict['bias_slow_master'] = bias_slow_master
        else:
            print("spec cal item not found: %s" % bias_slow_master)

        bias_fast_master = os.path.join(src, 'bias2.0.fits')
        if os.path.exists(bias_fast_master):
            spec_calib_dict['bias_fast_master'] = bias_fast_master
        else:
            print("spec cal item not found: %s" % bias_fast_master)

        flat = os.path.join(src, utdate + '_Flat.fits')
        if os.path.exists(flat):
            spec_calib_dict['flat'] = flat
        else:
            print("spec cal item not found: %s" % flat)

        hg_master = os.path.join(src, 'Hg.fits')
        if os.path.exists(hg_master):
            spec_calib_dict['hg_master'] = hg_master
        else:
            print("spec cal item not found: %s" % hg_master)

        xe_master = os.path.join(src, 'Xe.fits')
        if os.path.exists(xe_master):
            spec_calib_dict['xe_master'] = xe_master
        else:
            print("spec cal item not found: %s" % xe_master)

        cd_master = os.path.join(src, 'Cd.fits')
        if os.path.exists(cd_master):
            spec_calib_dict['cd_master'] = cd_master
        else:
            print("spec cal item not found: %s" % cd_master)

        hexagrid = os.path.join(src, utdate + '_HexaGrid.pkl')
        if os.path.exists(hexagrid):
            spec_calib_dict['hexagrid'] = hexagrid
        else:
            print("spec cal item not found: %s" % hexagrid)

        tracematch = os.path.join(src, utdate + '_TraceMatch.pkl')
        if os.path.exists(tracematch):
            spec_calib_dict['tracematch'] = tracematch
        else:
            print("spec cal item not found: %s" % tracematch)

        tracematch_withmasks = os.path.join(
            src, utdate + '_TraceMatch_WithMasks.pkl')
        if os.path.exists(tracematch_withmasks):
            spec_calib_dict['tracematch_withmasks'] = tracematch_withmasks
        else:
            print("spec cal item not found: %s" % tracematch_withmasks)

        wavesolution = os.path.join(src, utdate + '_WaveSolution.pkl')
        if os.path.exists(wavesolution):
            spec_calib_dict['wavesolution'] = wavesolution
        else:
            print("spec cal item not found: %s" % wavesolution)

        dispersionmap = os.path.join(
            src, utdate + '_wavesolution_dispersionmap.png')
        if os.path.exists(dispersionmap):
            spec_calib_dict['dispersionmap'] = dispersionmap
        else:
            print("spec cal item not found: %s" % dispersionmap)

        flatmap = os.path.join(src, utdate + '_flat3d.png')
        if os.path.exists(flatmap):
            spec_calib_dict['flatmap'] = flatmap
        else:
            print("spec cal item not found: %s" % flatmap)

    else:
        print("Source dir does not exist: %s" % src)

    spec_calib_dict['utdate'] = utdate
    spec_calib_dict['drpver'] = drp_ver

    sedmdb = db.SedmDb.SedmDB()
    spec_calib_id, status = sedmdb.add_spec_calib(spec_calib_dict)
    print(status)
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
        retcode = subprocess.call("~/spy what ifu*.fits > what.list", shell=True)
        if retcode == 0:
            # Generate new Makefile
            retcode = subprocess.call("~/spy plan ifu*.fits", shell=True)
            if retcode == 0:
                # Make bias + bias subtraction
                retcode = subprocess.call(("make", "-j", "16", "bias"))
                if retcode != 0:
                    print("bias failed, try again")
                    retcode = subprocess.call(("make", "bias"))
                if retcode == 0:
                    # Make CR rejection
                    retcode = subprocess.call(("make", "-j", "8", "crrs"))
                    if retcode != 0:
                        print("crrs failed, try again")
                        retcode = subprocess.call(("make", "-j", "8", "crrs"))
                    # Success on all fronts!
                    if retcode == 0:
                        print("bias, crrs processed for %d new images" % ncp)
                        ret = True
                    # Report failures
                    else:
                        print("could not make crrs")
                else:
                    print("could not make bias")
            else:
                print("could not make plan")
        else:
            print("could not make what.list")

    return ret
    # END: proc_bias_crrs


def cpsci(srcdir, destdir='./', fsize=8400960, datestr=None):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        srcdir (str): source directory (typically in /scr2/sedm/raw)
        destdir (str): destination directory (typically in /scr2/sedm/redux)
        fsize (int): size of completely copied file in bytes
        datestr (str): YYYYMMDD date string

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
    for f in srcfiles:
        # get base filename
        fn = f.split('/')[-1]
        # Is our source file complete?
        if os.stat(f).st_size >= fsize:
            # has it been previously copied?
            prev = [s for s in dflist if fn in s]
            # No? then copy the file
            if len(prev) == 0:
                # Call copy
                nc, ns, nob = docp(f, destdir + '/' + fn)
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
    print("Linked %d files" % ncp)
    # Do bias subtraction, CR rejection
    if ncp > 0:
        if not proc_bias_crrs(ncp):
            print("Error processing bias/crrs")
        if datestr is None:
            print("Illegal datestr parameter")
            return 0, None

    return ncp, copied
    # END: cpsci


def dosci(destdir='./', datestr=None):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        destdir (str): destination directory (typically in /scr2/sedm/redux)
        datestr (str): YYYYMMDD date string

    Returns:
        int: Number of ifu images actually copied

    """

    # Record copies and standard star observations
    ncp = 0
    copied = []
    # Get list of source files in destination directory
    srcfiles = sorted(glob.glob(os.path.join(destdir, 'crr_b_ifu*.fits')))
    # Loop over source files
    for f in srcfiles:
        # get base filename
        fn = f.split('/')[-1]
        procfn = 'spec*auto*' + fn.split('.')[0] + '*.fits'
        proced = glob.glob(os.path.join(destdir, procfn))
        # Is our source file processed?
        if len(proced) == 0:
            # Read FITS header
            ff = pf.open(f)
            hdr = ff[0].header
            ff.close()
            # Get OBJECT keyword
            obj = hdr['OBJECT'].split()[0]
            # Get DOMEST keyword
            dome = hdr['DOMEST']
            # skip Cal files
            if 'Calib:' in obj:
                continue
            # skip if dome closed
            if 'CLOSED' in dome or 'closed' in dome:
                continue
            # make finder
            make_finder(f)
            # record action
            copied.append(fn)
            ncp += 1
            # are we a standard star?
            if 'STD-' in obj:
                # Build cube for STD observation
                print("Building STD cube for " + fn)
                # Don't solve WCS for standards (always brightest in IFU)
                cmd = ("ccd_to_cube.py", datestr, "--build", fn, "--noguider")
                print(" ".join(cmd), flush=True)
                retcode = subprocess.call(cmd)
                # Check results
                if retcode != 0:
                    print("Error generating cube for " + fn)
                else:
                    # TODO: update SedmDb cube table
                    # Use auto psf aperture for standard stars
                    print("Extracting std star spectra for " + fn)
                    cmd = ("extract_star.py", datestr, "--auto", fn, "--std")
                    print(" ".join(cmd), flush=True)
                    retcode = subprocess.call(cmd)
                    if retcode != 0:
                        print("Error extracting std star spectra for " + fn)
                        badfn = "spec_auto_notfluxcal_" + fn.split('.')[0] + \
                                "_failed.fits"
                        cmd = ("touch", badfn)
                        subprocess.call(cmd)
                    else:
                        cmd = ("pysedm_report.py", datestr, "--contains",
                               fn.split('.')[0], "--slack")
                        print(" ".join(cmd), flush=True)
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            print("Error running report for " +
                                  fn.split('.')[0])
                        # run Verify.py
                        cmd = "~/sedmpy/drpifu/Verify.py %s --contains %s" % \
                              (datestr, fn.split('.')[0])
                        subprocess.call(cmd, shell=True)
                        # TODO: update SedmDb spec table
            else:
                # Build cube for science observation
                print("Building science cube for " + fn)
                # Solve WCS for science targets
                cmd = ("ccd_to_cube.py", datestr, "--build", fn, "--solvewcs")
                print(" ".join(cmd), flush=True)
                retcode = subprocess.call(cmd)
                # Check results
                if retcode != 0:
                    print("Error generating cube for " + fn)
                else:
                    # TODO: update SedmDb cube table
                    # Use forced psf for science targets
                    print("Extracting object spectra for " + fn)
                    cmd = ("extract_star.py", datestr, "--auto", fn,
                           "--autobins", "6")
                    print(" ".join(cmd), flush=True)
                    retcode = subprocess.call(cmd)
                    if retcode != 0:
                        print("Error extracting object spectrum for " + fn)
                        badfn = "spec_auto_notfluxcal_" + fn.split('.')[0] + \
                                "_failed.fits"
                        cmd = ("touch", badfn)
                        subprocess.call(cmd)
                    else:
                        print("Running SNID for " + fn)
                        cmd = ("make", "classify")
                        print(" ".join(cmd), flush=True)
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            print("Error running SNID")
                        cmd = ("pysedm_report.py", datestr, "--contains",
                               fn.split('.')[0], "--slack")
                        print(" ".join(cmd), flush=True)
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            print("Error running report for " +
                                  fn.split('.')[0])
                        # Upload spectrum to marshal
                        cmd = ("make", "ztfupload")
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            print("Error uploading spectra to marshal")
                        # run Verify.py
                        cmd = "~/sedmpy/drpifu/Verify.py %s --contains %s" % \
                              (datestr, fn.split('.')[0])
                        subprocess.call(cmd, shell=True)
                        # notify user that followup successfully completed
                        proced = glob.glob(os.path.join(destdir, procfn))[0]
                        if os.path.exists(proced):
                            email_user(proced, datestr, obj)
                            # TODO: update spec table in SedmDb
                        else:
                            print("Not found: %s" % proced)
    return ncp, copied
    # END: dosci


def update_cube(input_fitsfile):
    """ Update the SEDM database on pharos by adding a new cube entry. """

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
        print("ERROR: Ambiguous cube list")
        return -1
    # Read header
    ff = pf.open(cube_list[0])
    # Get header keyword values
    for key in header_dict.keys():
        hk = header_dict[key]
        if hk in ff[0].header:
            cube_dict[key] = ff[0].header[hk]
        else:
            print("Header keyword not found: %s" % hk)
    ff.close()
    # Open database connection
    sedmdb = db.SedmDb.SedmDB()
    # Get spec_calib id for this utdate
    spec_calib_id = sedmdb.get_from_spec_calib(['id'], {'utdate': utdate})[0]
    cube_dict['spec_calib_id'] = spec_calib_id[0]
    # Add into database
    cube_id, status = sedmdb.add_cube(cube_dict)
    print(status)
    return cube_id
    # END: update_cube


def email_user(spec_file, utdate, object_name):
    """ Send e-mail to requestor indicating followup completed"""
    # get request id
    ff = pf.open(spec_file)
    request = int(ff[0].header['REQ_ID'])
    quality = int(ff[0].header['QUALITY'])
    if quality == 3:
        status = 'IFU auto-extraction failed: target outside IFU.'
    elif quality == 4:
        status = 'IFU auto-extraction failed: > 20% of flux is negative.'
    elif quality == 5:
        status = 'IFU auto-extraction failed: guider astrometry failure - ' \
                 'may be fixed by manual extraction.'
    else:
        status = 'IFU auto-extraction succeeded.'
    link = 'http://pharos.caltech.edu/data_access/ifu?obsdate=%s' % utdate
    subj = 'SEDM followup status report for %s on %s' % (object_name, utdate)

    sedmdb = db.SedmDb.SedmDB()
    sedmdb.send_email_by_request(requestid=request, template='report_send',
                                 subject=subj,
                                 template_dict={'object_name': object_name,
                                                'link': link, 'status': status})


def make_finder(ffile):
    """ Spawn a process that makes a finder for the ifu file """
    cmd = ("/scr2/sedmdrp/spy", "/scr2/sedmdrp/sedmpy/drpifu/acq_finder.py",
           "--imfile", ffile)
    subprocess.Popen(cmd)


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
        print("%s already exists in %s" % (fname, destdir))
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
                print("Found %s in directory %s, linking to %s" %
                      (fname, d, destdir))
                break
    if not ret:
        print(dstr + fname + " not found")

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
        print("%s already exists in %s" % (fname, destdir))
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
                print("Found %s in directory %s, linking to %s" %
                      (fname, d, os.path.join(destdir, fname)))
                break
    if not ret:
        print("%s not found" % fname)
    return ret


def find_recent_fluxcal(redd, fname, destdir):
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
    if len(local_file) >= 1:
        print("%s already exists in %s" % (fname, destdir))
        ret = True
    # Search in redd for file
    else:
        # Get all but the most recent reduced data directories
        dspec = os.path.join(redd, '20??????')
        redlist = sorted([d for d in glob.glob(dspec)
                          if os.path.isdir(d)])[0:-1]
        # Go back in reduced dir list until we find our file
        for d in reversed(redlist):
            src = sorted(glob.glob(os.path.join(d, fname)))
            for s in src:
                # Skip sym-links
                if os.path.islink(s):
                    continue
                # Read FITS header
                f = pf.open(s)
                hdr = f[0].header
                f.close()
                # Skip if not Telluric corrected
                if 'TELLFLTR' not in hdr:
                    continue
                try:
                    newfile = os.path.join(destdir, s.split('/')[-1])
                    os.symlink(s, newfile)
                except OSError:
                    print("File already exists: %s" % newfile)
                ret = True
                print("Found %s in directory %s, linking to %s" %
                      (fname, d, newfile))
                break
            if ret:
                break
    if not ret:
        print("%s not found" % fname)
    return ret


def cpprecal(dirlist, destdir='./', fsize=8400960):
    """Copy raw cal files from previous date's directory

    Make sure we only look in previous day directory for files created
    within four hours of the day changeover for possible raw calibration
    files required for the current night's calibration.  Copy any such
    files into the destination directory.

    Args:
        dirlist (list): a list of raw directories (typically in /scr2/sedm/raw)
        destdir (str): where to put the files
        fsize (int): size of completely copied file in bytes

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
                f = pf.open(src)
                hdr = f[0].header
                f.close()
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
                                                       verbose=True)
                                    ncp += nc
                            else:
                                print("Bad dome - lamp not on: %s" % src)
                    # Check for arcs
                    elif 'Xe' in obj or 'Cd' in obj or 'Hg' in obj:
                        if exptime > 25.:
                            # Copy arc images
                            if not os.path.exists(destfil):
                                nc, ns, nob = docp(src, destfil, onsky=False,
                                                   verbose=True)
                                ncp += nc
                    # Check for biases
                    elif 'bias' in obj:
                        if exptime <= 0.:
                            # Copy bias images
                            if not os.path.exists(destfil):
                                nc, ns, nob = docp(src, destfil, onsky=False,
                                                   verbose=True)
                                ncp += nc
            else:
                print("Truncated file: %s" % src)

    return ncp
    # END: cpprecal


def cpcal(srcdir, destdir='./', fsize=8400960):
    """Copy raw cal files from srcdir into destdir.

    Find calibration files taken within 10 hours of the day changeover
    and copy them to the destination directory.

    Args:
        srcdir (str): source for raw cal images
        destdir (str): place to put the cal images
        fsize (int): size of completely copied file in bytes

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
            f = pf.open(src)
            hdr = f[0].header
            f.close()
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
                    if exptime > 100. and ('dome' in obj and
                                           'Xe' not in obj and
                                           'Hg' not in obj and 
                                           'Cd' not in obj):
                        if lampcur > 0.0:
                            # Copy dome images
                            nc, ns, nob = docp(src, destfil, onsky=False,
                                               verbose=True)
                            ncp += nc
                        else:
                            print("Bad dome - lamp not on: %s" % src)
                # Check for arcs
                elif 'Xe' in obj or 'Cd' in obj or 'Hg' in obj:
                    if exptime > 25.:
                        # Copy arc images
                        nc, ns, nob = docp(src, destfil, onsky=False,
                                           verbose=True)
                        ncp += nc
                # Check for biases
                elif 'bias' in obj:
                    if exptime <= 0.:
                        # Copy bias images
                        nc, ns, nob = docp(src, destfil, onsky=False,
                                           verbose=True)
                        ncp += nc

    return ncp
    # END: cpcal


def obs_loop(rawlist=None, redd=None, check_precal=True, indir=None,
             piggyback=False):
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
        rawlist (list): list of raw data directories (usually in /scr2/sedm/raw)
        redd (str): reduced directory (something like /scr2/sedm/redux)
        check_precal (bool): should we check for images from previous night?
        indir (str): input directory for single night processing
        piggyback (bool): basic processing done by StartObs.py

    Returns:
        bool: True if night completed normally, False otherwise

    Note:
        KeyboardInterrupt handler exits gracefully with a ctrl-C.

    """
    # Set up Observatory params
    p60 = ephem.Observer()
    p60.lat = '33:21:00'
    p60.lon = '-116:52:00'
    p60.elevation = 1706
    sun = ephem.Sun()

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
    # report
    print("Raw files from  : %s\nReduced files to: %s" % (srcdir, outdir))
    # Check if processed cal files are ready
    if not cube_ready(outdir, cur_date_str):
        # Wait for cal files until sunset
        sunset = p60.next_setting(sun)
        if piggyback:
            print("Skipping check for raw cal files")
            ncp = 0
        else:
            # Copy raw cal files from previous date directory
            npre = cpprecal(rawlist, outdir)
            print("Linked %d raw cal files from %s" % (npre, rawlist[-2]))
            # Now check the current source dir for raw cal files
            ncp = cpcal(srcdir, outdir)
            print("Linked %d raw cal files from %s" % (ncp, srcdir))
        # Now loop until we have the raw cal files we need or sun is down
        while not cal_proc_ready(outdir, ncp=ncp, test_cal_ims=piggyback):
            # Wait a minute
            print("waiting 60s...", flush=True)
            now = ephem.now()
            time.sleep(60)
            if piggyback:
                print("checking for processed cal files")
                ncp = 0
            else:
                if check_precal and now.tuple()[3] >= 20:
                    print("checking %s for new raw cal files..." % rawlist[-2])
                    ncp = cpprecal(rawlist, outdir)
                    print("Linked %d raw cal files from %s" % (ncp,
                                                               rawlist[-2]))
                else:
                    print("checking %s for new raw cal files..." % srcdir)
                    ncp = cpcal(srcdir, outdir)
                print("Linked %d raw cal files from %s" % (ncp, srcdir),
                      flush=True)
            if ncp <= 0:
                # Check to see if we are still before an hour after sunset
                now = ephem.now()
                if now < sunset + ephem.hour:
                    print("UT  = %02d/%02d %02d:%02d < sunset "
                          "(%02d/%02d %02d:%02d) + 1hr, "
                          "so keep waiting" % (now.tuple()[1], now.tuple()[2],
                                               now.tuple()[3], now.tuple()[4],
                                               sunset.tuple()[1],
                                               sunset.tuple()[2],
                                               sunset.tuple()[3],
                                               sunset.tuple()[4]),
                          flush=True)
                else:
                    print("UT = %02d/%02d %02d:%02d >= sunset "
                          "(%02d/%02d %02d:%02d) + 1hr, "
                          "time to get a cal set" % (now.tuple()[1],
                                                     now.tuple()[2],
                                                     now.tuple()[3],
                                                     now.tuple()[4],
                                                     sunset.tuple()[1],
                                                     sunset.tuple()[2],
                                                     sunset.tuple()[3],
                                                     sunset.tuple()[4]),
                          flush=True)
                    break
            else:
                # Get new listing
                retcode = subprocess.call("~/spy what ifu*.fits > what.list",
                                       shell=True)
                if retcode != 0:
                    print("what oops!")

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
                cmd = ("ccd_to_cube.py", cur_date_str, "--tracematch",
                       "--hexagrid")
                print(" ".join(cmd), flush=True)
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
                    print(" ".join(cmd), flush=True)
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
                    print("Finished all %d parts, merging..." % nsub)
                    # Merge the solutions
                    subprocess.call(("derive_wavesolution.py", cur_date_str,
                                     "--merge"))
                procw_time = int(time.time() - start_time)
                if os.path.exists(
                   os.path.join(outdir, cur_date_str + '_WaveSolution.pkl')):
                    # Process flat
                    start_time = time.time()
                    cmd = ("ccd_to_cube.py", cur_date_str, "--flat")
                    print(" ".join(cmd), flush=True)
                    subprocess.call(cmd)
                    if not (os.path.exists(
                            os.path.join(outdir, cur_date_str + '_Flat.fits'))):
                        print("Making of %s_Flat.fits failed!" % cur_date_str)
                else:
                    print("Making of %s cube failed!" % cur_date_str)
                procf_time = int(time.time() - start_time)
                # Report times
                print("Calibration processing took "
                      "%d s (bias,crrs), %d s (grid),"
                      " %d s (waves), and %d s (flat)" %
                      (procb_time, procg_time, procw_time, procf_time))

        # Check status
        if not cube_ready(outdir, cur_date_str):
            print("These calibrations failed!")
            print("Let's get our calibrations from a previous night")
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
            print("Using older calibration files")
    else:
        print("Calibrations already present in %s" % outdir)

    # Update spec_calib table in sedmdb
    spec_calib_id = update_calibration(cur_date_str)
    print("SedmDb table spec_calib updated with id %d" % spec_calib_id)

    print("Calibration stage complete, ready for science!")
    # Link recent flux cal file
    find_recent_fluxcal(redd, 'fluxcal*.fits', outdir)
    # Keep track of no copy
    nnc = 0
    # loop and copy new files
    doit = True
    try:
        while doit:
            # Wait for next sunrise
            sunrise = p60.next_rising(sun)
            # Wait a minute
            print("waiting 60s for new ifu images...")
            sys.stdout.flush()
            time.sleep(60)
            # Check for new ifu images
            print("checking %s for new ifu images..." % srcdir)
            sys.stdout.flush()
            # Record starting time for new file processing
            start_time = time.time()
            if piggyback:
                nsci, science = dosci(outdir, datestr=cur_date_str)
                ncp = nsci
            else:
                ncp, copied = cpsci(srcdir, outdir, datestr=cur_date_str)
                nsci, science = dosci(outdir, datestr=cur_date_str)
            # We copied some new ones so report processing time
            if ncp > 0:
                proc_time = int(time.time() - start_time)
                print("%d new ifu images copied and %d processed in %d s" %
                      (ncp, nsci, proc_time))
                sys.stdout.flush()
                nnc = 0
            else:
                nnc += 1
            # Have we been waiting for a while?
            if nnc > 3:
                # Check time
                now = ephem.now()
                if now >= sunrise:
                    # No new observations and sun is probably up!
                    print("No new images for %d minutes and UT = "
                          "%02d/%02d %02d:%02d > "
                          "%02d/%02d %02d:%02d so sun is up!" %
                          (nnc, now.tuple()[1], now.tuple()[2],
                           now.tuple()[3], now.tuple()[4],
                           sunrise.tuple()[1], sunrise.tuple()[2],
                           sunrise.tuple()[3], sunrise.tuple()[4]))
                    print("Time to wait until we have a new raw directory")
                    doit = False
                    # Normal termination
                    subprocess.call(("make", "report"))
                    ret = True
                else:
                    print("No new image for %d minutes but UT = %02d/%02d "
                          "%02d:%02d <= "
                          "%02d/%02d %02d:%02d, so sun is still down, "
                          "keep waiting" %
                          (nnc, now.tuple()[1], now.tuple()[2],
                           now.tuple()[3], now.tuple()[4],
                           sunrise.tuple()[1], sunrise.tuple()[2],
                           sunrise.tuple()[3], sunrise.tuple()[4]))

    # Handle a ctrl-C
    except KeyboardInterrupt:
        sys.exit("Exiting")

    return ret
    # END: obs_loop


def go(rawd='/scr2/sedm/raw', redd='/scr2/sedmdrp/redux', wait=False,
       check_precal=True, indate=None, piggyback=False):
    """Outermost infinite loop that watches for a new raw directory.

    Keep a list of raw directories in `redd` and fire off
    the obs_loop procedure when a new directory appears.  Check for
    a new raw directory every 10 minutes.

    Args:
        rawd (str): raw directory, should be /scr2/sedm/raw
        redd (str): reduced directory, should be like /scr2/sedm/redux
        wait (bool): wait for new directory, else start right away
        check_precal (bool): should we check previous night for cals?
        indate (str): input date to process: YYYYMMDD (e.g. 20180626)
        piggyback (bool): True if using other script to copy data

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
    print("Found %d raw directories in %s: putting reduced data in %s" %
          (nraw, rawd, redd))
    if indate is None:
        print("Latest raw directory is %s" % rawlist[-1])

        if not wait:
            stat = obs_loop(rawlist, redd, check_precal=check_precal,
                            piggyback=piggyback)
            its += 1
            print("Finished SEDM observing iteration %d in raw dir %s" %
                  (its, rawlist[-1]))
        try:
            while dobs:
                if stat:
                    print("Now we wait until we get a new raw directory")
                    waiting = True
                    while waiting:
                        print("waiting 10min for new raw directory...")
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
                            print("Starting next SEDM observing iteration with "
                                  "raw dir %s" % rawlist[-1])
                        else:
                            gm = time.gmtime()
                            print("UT = %02d:%02d No new directories yet, "
                                  "so keep waiting" % (gm.tm_hour, gm.tm_min))
                            sys.stdout.flush()
                print("Found %d raw directories in %s: "
                      "putting reduced data in %s" % (nraw, rawd, redd))
                print("Latest raw directory is %s" % rawlist[-1])
                stat = obs_loop(rawlist, redd, check_precal=check_precal,
                                piggyback=piggyback)
                its += 1
                print("Finished SEDM observing iteration %d in raw dir %s" %
                      (its, rawlist[-1]))
        # Handle a ctrl-C
        except KeyboardInterrupt:
            sys.exit("Exiting")
    else:
        indir = os.path.join(rawd, indate)
        print("Processing raw data from %s" % indir)
        stat = obs_loop(rawlist, redd, check_precal=check_precal, indir=indir,
                        piggyback=piggyback)
        its += 1
        print("Finished SEDM processing in raw dir %s with status %d" %
              (indir, stat))
    # END: go


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Start SEDM pipeline

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--rawdir', type=str, default='/scr2/sedm/raw',
                        help='Input raw directory (/scr2/sedm/raw)')
    parser.add_argument('--reduxdir', type=str, default='/scr2/sedmdrp/redux',
                        help='Output reduced directory (/scr2/sedmdrp/redux)')
    parser.add_argument('--wait', action="store_true", default=False,
                        help='Wait for new directory first')
    parser.add_argument('--piggyback', action="store_true", default=False,
                        help='Do not copy data, copied by another script')
    parser.add_argument('--skip_precal', action="store_true", default=False,
                        help='Skip check of previous day for cal files?')
    parser.add_argument('--date', type=str, default=None,
                        help='Select date to process')

    args = parser.parse_args()

    go(rawd=args.rawdir, redd=args.reduxdir, wait=args.wait,
       check_precal=(not args.skip_precal), indate=args.date,
       piggyback=args.piggyback)
