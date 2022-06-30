"""Conduct automatic reduction of SEDM data on sedmdrp@minar

Functions
    * :func:`red_loop`     one night observing loop
    * :func:`dosci`        processes new science images files
    * :func:`proc_bias_crrs`  processes biases and CR rejection
    * :func:`cal_proc_ready`  check if all required raw cal images are present
    * :func:`cube_ready`      check if all required cal files are present
    * :func:`bias_ready`    check if master bias files are present

Note:
    This is used as a python script as follows::

        usage: Reduce.py

        optional arguments:
          None

"""
import time
import glob
import sys
import os
import argparse
import subprocess
import astropy.io.fits as pf
try:
    import drpifu.AutoReduce as drp
except ImportError:
    import AutoReduce as drp


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
    print("Cals ready?: trace: %d, trace/mask: %d, grid: %d, "
          "wave: %d, flat: %d" % (ft, ftm, fg, fw, ff))
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
    nbias2 = 0
    nxe = 0
    nhg = 0
    ncd = 0
    ndome = 0
    ret = False

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
                            if speed == 0.1:
                                nbias += 1
                        if 'Xe' in obj:
                            nxe += 1
                        if 'dome' in obj:
                            ndome += 1
                        if 'Hg' in obj:
                            nhg += 1
                        if 'Cd' in obj:
                            ncd += 1

            # Do we have the ideal number of calibration files?
            if (nbias2 >= 10 and nbias >= 10 and nxe >= 5 and ndome >= 5 and
                    nhg >= 5 and ncd >= 5):
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
        retcode = os.system("~/spy what ifu*.fits > what.list")
        if retcode == 0:
            # Generate new Makefile
            # Are we using a previous calibration set?
            retcode2 = os.system("~/spy plan ifu*.fits")
            if retcode2 == 0:
                # Make bias + bias subtraction
                retcode3 = os.system("make -j 16 bias")
                if retcode3 != 0:
                    print("bias failed, try again")
                    retcode3 = os.system("make bias")
                if retcode3 == 0:
                    # Make CR rejection
                    retcode4 = os.system("make -j 8 crrs")
                    if retcode4 != 0:
                        print("crrs failed, try again")
                        retcode4 = os.system("make crrs")
                    # Success on all fronts!
                    if retcode4 == 0:
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


def dosci(destdir='./', datestr=None, ztfupld=False, slack=False,
          dbupdate=False):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        destdir (str): destination directory (typically in /data/sedmdrp/redux)
        datestr (str): YYYYMMDD date string
        ztfupld (bool): upload to ZTF marshal?
        slack (bool): upload to pysedm-report Slack channel?
        dbupdate (bool): update SEDm database?

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
                if retcode > 0:
                    print("Error generating cube for " + fn)
                else:
                    if dbupdate:
                        # Update SedmDb cube table
                        cube_id = drp.update_cube(f)
                        if cube_id > 0:
                            print("SEDM db accepted cube at id %d" % cube_id)
                        else:
                            print("SEDM db rejected cube")
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
                        if slack:
                            cmd = ("pysedm_report.py", datestr, "--contains",
                                   fn.split('.')[0], "--slack")
                        else:
                            cmd = ("pysedm_report.py", datestr, "--contains",
                                   fn.split('.')[0])
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
                    if dbupdate:
                        # Update SedmDb cube table
                        cube_id = drp.update_cube(f)
                        if cube_id > 0:
                            print("SEDM db accepted cube at id %d" % cube_id)
                        else:
                            print("SEDM db rejected cube")
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
                        if retcode > 0:
                            print("Error running SNID")
                        if slack:
                            cmd = ("pysedm_report.py", "--contains",
                                   fn.split('.')[0], "--slack")
                        else:
                            cmd = ("pysedm_report.py", "--contains",
                                   fn.split('.')[0])
                        print(" ".join(cmd), flush=True)
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            print("Error running report for " +
                                  fn.split('.')[0])
                        # Upload spectrum to marshal
                        if ztfupld:
                            cmd = ("make", "ztfupload")
                            retcode = subprocess.call(cmd)
                            if retcode != 0:
                                print("Error uploading spectra to marshal")
                        # run Verify.py
                        cmd = "~/sedmpy/drpifu/Verify.py %s --contains %s" % \
                              (datestr, fn.split('.')[0])
                        subprocess.call(cmd, shell=True)
                        # TODO: update SedmDb spec table
    return ncp, copied
    # END: dosci


def red_loop(outdir=None, upld=False, slack=False, dbup=False):
    """One night observing loop: processes calibrations and science data

    Args:
        outdir (str): directory for single night processing
        upld (bool): upload to ZTF marshal?
        slack (bool): upload to pysedm-report Slack channel?
        dbup (bool): update the SEDm database

    Returns:
        bool: True if night completed normally, False otherwise

    Note:
        KeyboardInterrupt handler exits gracefully with a ctrl-C.

    """

    # Default return value
    ret = False
    # Current date string
    cur_date_str = str(outdir.split('/')[-1])
    # report
    print("Raw/Reduced files from/to: %s" % outdir)
    # Check if processed cal files are ready
    if not cube_ready(outdir, cur_date_str):
        # Get new listing
        retcode = subprocess.call("~/spy what ifu*.fits > what.list",
                                  shell=True)
        if retcode > 0:
            print("what oops!")

        # Process calibrations if we are using them
        if cal_proc_ready(outdir, mintest=True):
            # bias subtract and CR reject
            start_time = time.time()
            if proc_bias_crrs(20):
                procb_time = int(time.time() - start_time)
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
    else:
        print("Calibrations already present in %s" % outdir)

    if dbup:
        # Update spec_calib table in sedmdb
        spec_calib_id = drp.update_calibration(cur_date_str)
        print("SEDM db accepted spec_calib at id %d" % spec_calib_id)

    print("Calibration stage complete, ready for science!")
    # process files
    start_time = time.time()
    nsci, science = dosci(outdir, datestr=cur_date_str, ztfupld=upld,
                          slack=slack, dbupdate=dbup)
    # We copied some new ones so report processing time

    proc_time = int(time.time() - start_time)
    print("%d new ifu images processed in %d s" % (nsci, proc_time))
    sys.stdout.flush()

    return ret
    # END: red_loop


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""Reduce a night's SEDM data

                """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--upload', action="store_true", default=False,
                        help='Upload to ZTF marshal')
    parser.add_argument('--toslack', action="store_true", default=False,
                        help='Upload to pysedm-report Slack channel')
    parser.add_argument('--dbupdate', action="store_true", default=False,
                        help='Update SEDm database')

    args = parser.parse_args()

    # Get current directory
    curdir = os.getcwd()
    reduxdir = '/'.join(curdir.split('/')[:-1])
    reduxenv = os.getenv("SEDMREDUXPATH")
    # Test if environment set correctly
    if reduxenv is None:
        print("ERROR: please set SEDMREDUXPATH to something")
    else:
        if reduxenv not in reduxdir:
            print("ERROR: setenv SEDMREDUXPATH correctly and try again")
        else:
            red_loop(outdir=curdir, upld=args.upload, slack=args.toslack,
                     dbup=args.dbupdate)
