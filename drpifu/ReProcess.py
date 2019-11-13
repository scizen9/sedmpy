"""Conduct automatic reduction of SEDM data in sedmdrp@pharos

Functions
    * :func:`cal_loop`     one night calibration loop
    * :func:`cal_proc_ready`  check if all required raw cal images are present
    * :func:`cube_ready`      check if all required cal files are present
    * :func:`update_calibration`    update cal cube in SEDM db

Note:
    This is used as a python script as follows::

        usage: ReProcess.py [-h] [--rawdir RAWDIR] [--reduxdir REDUXDIR]

        optional arguments:
          -h, --help           show this help message and exit
          --reduxdir REDUXDIR  Output reduced directory (/scr2/sedm/redux)
          --date YYYYMMDD      Select date to process (None)
          --nodb               Do not update SEDM Db (False)

"""
import time
import glob
import os
import shutil
import re
import subprocess
import logging
import argparse
from astropy.io import fits as pf

from configparser import ConfigParser
import codecs
from AutoReduce import update_spec, make_e3d, update_calibration

try:
    import rcimg
except ImportError:
    import drprc.rcimg as rcimg

try:
    import Version
except ImportError:
    import drpifu.Version as Version

drp_ver = Version.ifu_drp_version()
logging.basicConfig(
    format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
    datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)

# Get pipeline configuration
cfg_parser = ConfigParser()
# Find config file: default is sedmpy/drpifu/config/sedmconfig.cfg
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                              'config/sedmconfig.cfg')
# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    cfg_parser.read_file(f)
# Get paths
_rawpath = cfg_parser.get('paths', 'rawpath')
_reduxpath = cfg_parser.get('paths', 'reduxpath')
_srcpath = cfg_parser.get('paths', 'srcpath')


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
                    wave_stats_ok = (float(line.split()[-1]) < 35.0)
        # Does wavelength solution pass?
        if wave_stats_ok:
            logging.info("Wavelength stats passed")
        else:
            # NO: move bad files away
            os.mkdir(os.path.join(caldir, 'bad'))
            os.system("mv %s_* %s" % (os.path.join(caldir, cur_date_str),
                                      os.path.join(caldir, 'bad')))
            logging.warning("Wavelength stats failed, moved cube to 'bad'")
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


def cal_proc_ready(caldir='./'):
    """Check counts for all required raw cal file types in caldir directory.

    Args:
        caldir (str): directory where raw cal files reside

    Returns:
        bool: True if required raw cal files are present, False otherwise

    """

    ret = False

    # Check for gzipped versions first
    dof = glob.glob(os.path.join(caldir, 'dome.fits.gz'))
    if len(dof) == 1:
        subprocess.run(["gunzip", dof[0]])
    hgf = glob.glob(os.path.join(caldir, 'Hg.fits.gz'))
    if len(hgf) == 1:
        subprocess.run(["gunzip", hgf[0]])
    cdf = glob.glob(os.path.join(caldir, 'Cd.fits.gz'))
    if len(cdf) == 1:
        subprocess.run(["gunzip", cdf[0]])
    xef = glob.glob(os.path.join(caldir, 'Xe.fits.gz'))
    if len(xef) == 1:
        subprocess.run(["gunzip", xef[0]])

    dof = glob.glob(os.path.join(caldir, 'dome.fits'))
    hgf = glob.glob(os.path.join(caldir, 'Hg.fits'))
    cdf = glob.glob(os.path.join(caldir, 'Cd.fits'))
    xef = glob.glob(os.path.join(caldir, 'Xe.fits'))
    if len(dof) == 1 and len(hgf) == 1 and len(cdf) == 1 and len(xef) == 1:
        ret = True

    return ret
    # END: cal_proc_ready


def delete_old_pysedm_files(odir, ut_date, keep_spec=False):
    """Remove all the old pysedm output files"""
    if keep_spec:
        # Make archive directory
        archdir = os.path.join(odir, 'pysedm')
        if not os.path.exists(archdir):
            os.mkdir(archdir)
    else:
        archdir = None
    # Get list of pysedm files to remove
    flist = glob.glob(os.path.join(odir, '%s_*' % ut_date))
    flist.extend(glob.glob(os.path.join(odir, '*_crr_b_ifu%s*' % ut_date)))
    flist.extend(glob.glob(os.path.join(odir, 'bkgd_dome.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'e3d_dome.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'pysedm_run.log')))
    flist.extend(glob.glob(os.path.join(odir, 'report.txt')))
    # Remove them
    if len(flist) > 0:
        ndelfile = 0
        nkeepspec = 0
        for fl in flist:
            # Keep spectra?
            if keep_spec:
                if 'spec_' in fl:
                    shutil.move(fl, archdir)
                    nkeepspec += 1
                else:
                    os.remove(fl)
                    ndelfile += 1
            else:
                os.remove(fl)
    else:
        logging.warning("No pysedm files found in %s" % odir)
    # If we keep the spectra, gzip them
    if keep_spec:
        # Now gzip the files in the archive
        flist = os.listdir(archdir)
        for fl in flist:
            if 'gz' not in fl:
                subprocess.run(["gzip", os.path.join(archdir, fl)])
    # END: delete_old_pysedm_files


def archive_old_pysedm_files(odir, ut_date):
    """Move all the old pysedm output files to ./pysedm/"""
    # Make archive directory
    archdir = os.path.join(odir, 'pysedm')
    if not os.path.exists(archdir):
        os.mkdir(archdir)
    # Get list of pysedm files to move
    flist = glob.glob(os.path.join(odir, '%s_*' % ut_date))
    flist.extend(glob.glob(os.path.join(odir, '*_crr_b_ifu%s*' % ut_date)))
    flist.extend(glob.glob(os.path.join(odir, 'bkgd_dome.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'e3d_dome.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'pysedm_run.log')))
    flist.extend(glob.glob(os.path.join(odir, 'report.txt')))
    flist.extend(glob.glob(os.path.join(odir, 'pysedm_run.log')))
    # Move them into the archive
    if len(flist) > 0:
        nfilemv = 0
        for fl in flist:
            rute = fl.split('/')[-1]
            if not os.path.exists(os.path.join(archdir, rute)):
                shutil.move(fl, archdir)
                nfilemv += 1
        logging.info("Moved %d old pysedm files into %s" % (nfilemv, archdir))
    else:
        logging.warning("No pysedm files found in %s" % odir)
    # Now gzip the files in the archive
    flist = os.listdir(archdir)
    for fl in flist:
        if 'gz' not in fl:
            subprocess.run(["gzip", os.path.join(archdir, fl)])
    # END: archive_old_pysedm_files


def archive_old_kpy_files(odir):
    """Move all the old kpy output files to ./kpy/"""
    # Make archive directory
    archdir = os.path.join(odir, 'kpy')
    if not os.path.exists(archdir):
        os.mkdir(archdir)
    # Get list of kpy files to move
    flist = glob.glob(os.path.join(odir, '*.npy'))
    flist.extend(glob.glob(os.path.join(odir, '*SEDM*')))
    flist.extend(glob.glob(os.path.join(odir, 'allspec_*')))
    flist.extend(glob.glob(os.path.join(odir, 'bs_*')))
    flist.extend(glob.glob(os.path.join(odir, 'bgd_*')))
    flist.extend(glob.glob(os.path.join(odir, 'back_*')))
    flist.extend(glob.glob(os.path.join(odir, 'cog_*')))
    flist.extend(glob.glob(os.path.join(odir, 'cube_*')))
    flist.extend(glob.glob(os.path.join(odir, 'flex_bs_*')))
    flist.extend(glob.glob(os.path.join(odir, 'image_*')))
    flist.extend(glob.glob(os.path.join(odir, 's_*')))
    flist.extend(glob.glob(os.path.join(odir, 'seg_*')))
    flist.extend(glob.glob(os.path.join(odir, 'var_*')))
    flist.extend(glob.glob(os.path.join(odir, 'XYs_*')))
    flist.extend(glob.glob(os.path.join(odir, 'cat_Hg.fits.txt')))
    flist.extend(glob.glob(os.path.join(odir, 'default.conv')))
    flist.extend(glob.glob(os.path.join(odir, 'deleteme')))
    flist.extend(glob.glob(os.path.join(odir, 'ds9_dome.fits_segments.reg')))
    flist.extend(glob.glob(os.path.join(odir, 'filtered_Hg.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'flat-field-values.pdf')))
    flist.extend(glob.glob(os.path.join(odir, 'rough.reg')))
    flist.extend(glob.glob(os.path.join(odir, 'Standard_Correction.pdf')))
    flist.extend(glob.glob(os.path.join(odir, 'run.log')))
    # Move them into the archive
    if len(flist) > 0:
        nfilemv = 0
        for fl in flist:
            rute = fl.split('/')[-1]
            if not os.path.exists(os.path.join(archdir, rute)):
                shutil.move(fl, archdir)
                nfilemv += 1
        logging.info("Moved %d old kpy files into %s" % (nfilemv, archdir))
    else:
        logging.warning("No kpy files found in %s" % odir)
    # Now check spec files
    flist = glob.glob(os.path.join(odir, 'spec_*'))
    if len(flist) > 0:
        nspecmv = 0
        for fl in flist:
            if '_crr_b_ifu' not in fl:
                shutil.move(fl, archdir)
                nspecmv += 1
        logging.info("Moved %d kpy spec files into %s" % (nspecmv, archdir))
    else:
        logging.warning("No spec files found in %s" % odir)
    # Now gzip the files in the archive
    flist = os.listdir(archdir)
    for fl in flist:
        if 'gz' not in fl:
            subprocess.run(["gzip", os.path.join(archdir, fl)])
    # END: archive_old_kpy_files


def dosci(destdir='./', datestr=None, nodb=False):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        destdir (str): destination directory (typically in /scr2/sedm/redux)
        datestr (str): YYYYMMDD date string
        nodb (bool): if True no update to SEDM db

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
            # record action
            copied.append(fn)
            ncp += 1
            # are we a standard star?
            if 'STD-' in obj:
                e3d_good = make_e3d(fnam=f, destdir=destdir, datestr=datestr,
                                    nodb=nodb, sci=False, hdr=None)
                if e3d_good:
                    # Get seeing
                    seeing = rcimg.get_seeing(imfile=fn, destdir=destdir,
                                              save_fig=True)
                    # Use auto psf extraction for standard stars
                    if seeing > 0:
                        logging.info("seeing measured as %f" % seeing)
                        cmd = ("extractstar.py", datestr, "--auto", fn,
                               "--std", "--tag", "robot",
                               "--centroid", "brightest",
                               "--seeing", "%.2f" % seeing)
                    else:
                        logging.info("seeing not measured for %s" % fn)
                    cmd = ("extract_star.py", datestr, "--auto", fn,
                           "--std", "--tag", "robot", "--maxpos")  # ,
                    # "--centroid", "brightest")
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
                        cmd = ("pysedm_report.py", datestr, "--contains",
                               fn.split('.')[0])
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
                e3d_good = make_e3d(fnam=f, destdir=destdir, datestr=datestr,
                                    nodb=nodb, sci=True, hdr=hdr)

                if e3d_good:
                    logging.info("e3d science cube generated")
                else:
                    logging.error("Cannot perform extraction for %s" % fn)
    return ncp, copied
    # END: dosci


def cal_loop(redd=None, indir=None, nodb=False,
             arch_kpy=False, arch_pysedm=False):
    """Create calibration files for one night.

    Args:
        redd (str): reduced directory (something like /scr2/sedm/redux)
        indir (str): input directory for single night processing
        nodb (bool): True if no update to SEDM db
        arch_kpy (bool): True to archive old kpy files into ./kpy/
        arch_pysedm (bool): True to archive old pysedm files into ./pysedm/

    Returns:
        bool: True if cals completed normally, False otherwise

    """
    # Default return value
    ret = False
    # Output directory is based on redd and indir
    outdir = os.path.join(redd, indir)
    # Current date string
    cur_date_str = str(outdir.split('/')[-1])
    # change to directory
    os.chdir(outdir)
    # report
    logging.info("Reduced files to: %s" % outdir)
    # Do we archive first?
    if arch_kpy:
        archive_old_kpy_files(outdir)
    if arch_pysedm:
        archive_old_pysedm_files(outdir, cur_date_str)
    else:
        delete_old_pysedm_files(outdir, cur_date_str, keep_spec=True)

    # Check calibration status
    cal_good = False
    if not cube_ready(outdir, cur_date_str):
        # Process calibrations
        if cal_proc_ready(outdir):
            # Process calibration
            start_time = time.time()
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
                        procf_time = int(time.time() - start_time)
                        # Report times
                        logging.info("Calibration processing took %d s (grid), "
                                     "%d s (waves),  and %d s (flat)" %
                                     (procg_time, procw_time, procf_time))
                        if nodb:
                            logging.warning("Not updating SEDM db")
                        else:
                            # Update spec_calib table in sedmdb
                            spec_calib_id = update_calibration(cur_date_str)
                            logging.info(
                                "SEDM db accepted spec_calib at id %d"
                                % spec_calib_id)
                        cal_good = True
                        logging.info(
                            "Calibration stage complete, ready for science!")
                        # Gzip cal images
                        subprocess.run(["gzip", "dome.fits", "Hg.fits",
                                        "Cd.fits", "Xe.fits"])
                else:
                    logging.error("Making of wavesolution failed!")
            else:
                logging.error("Making of tracematch and hexagrid failed!")
        # Check status
        if not cube_ready(outdir, cur_date_str):
            logging.error("These calibrations failed!")
    else:
        logging.info("Calibrations already present in %s" % outdir)
        cal_good = True
    # Proceed to build e3d cubes
    if cal_good:
        # Gunzip input files
        logging.info("Gunzipping crr_b_ifu%s*.fits.gz in dir %s"
                     % (cur_date_str, cur_date_str))
        subprocess.run(["gunzip", "crr_b_ifu%s*.fits.gz" % cur_date_str])
        dosci(outdir, datestr=cur_date_str, nodb=nodb)
        # Re-gzip input files
        # Gunzip input files
        logging.info("Gzipping crr_b_ifu%s*.fits in dir %s"
                     % (cur_date_str, cur_date_str))
        subprocess.run(["gzip", "crr_b_ifu%s*.fits" % cur_date_str])

    return ret
    # END: cal_loop


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Make cals for PYSEDM pipeline

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--reduxdir', type=str, default=_reduxpath,
                        help='Output reduced directory (%s)' % _reduxpath)
    parser.add_argument('--date', type=str, default=None,
                        help='Select date to process')
    parser.add_argument('--nodb', action="store_true", default=False,
                        help='Do not update SEDM Db')
    parser.add_argument('--archive_kpy', action="store_true", default=False,
                        help='Archive any existing kpy files')
    parser.add_argument('--archive_pysedm', action="store_true", default=False,
                        help='Archive any existing pysedm files')

    args = parser.parse_args()

    if not args.date:
        logging.error("Must provide a YYYYMMDD date with --date")
    else:
        cal_loop(redd=args.reduxdir, indir=args.date, nodb=args.nodb,
                 arch_kpy=args.archive_kpy, arch_pysedm=args.archive_pysedm)
