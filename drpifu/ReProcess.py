"""Conduct automatic reduction of SEDM data in sedmdrp@pharos

Functions
    * :func:`reproc`     one night calibration loop
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
import datetime
from astropy.io import fits as pf
from astropy.time import Time, TimeDelta

from configparser import ConfigParser
import codecs
try:
    from AutoReduce import update_spec, make_e3d, update_calibration, \
        proc_bias_crrs
except ImportError:
    from drpifu.AutoReduce import update_spec, make_e3d, update_calibration, \
        proc_bias_crrs

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
# Find config file: default is sedmpy/config/sedmconfig.cfg
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                              '../config/sedmconfig.cfg')
# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    cfg_parser.read_file(f)
# Get paths
_rawpath = cfg_parser.get('paths', 'rawpath')
_reduxpath = cfg_parser.get('paths', 'reduxpath')
_srcpath = cfg_parser.get('paths', 'srcpath')


def docp(src, dest, onsky=True, verbose=False, skip_cals=False, header=None):
    """Low level copy from raw directory to redux directory.

    Checks for raw ifu files, while avoiding any test and focus images.
    Uses os.symlink to do the copying (linked to conserve disk space).

    Args:
        src (str): source file
        dest (str): destination file
        onsky (bool): test for dome conditions or not
        verbose (bool): output messages?
        skip_cals (bool): skip copying cal images?
        header (dict): header from source image

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
    # Read FITS header if needed
    if header is None:
        ff = pf.open(src)
        hdr = ff[0].header
        ff.close()
    else:
        hdr = header
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
            if '.gz' in src:
                subprocess.call("gunzip -c %s > %s" % (src, dest), shell=True)
            else:
                # Symlink to save disk space
                os.symlink(src, dest)
            if 'STD-' in obj:
                nstd = 1
                logging.info("Standard %s linked to %s" % (obj, dest))
            else:
                nobj = 1
                logging.info('Target %s linked to %s' % (obj, dest))
            ncp = 1
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
    # Record copies and standard star observations
    ncp = 0
    nstd = 0
    nobj = 0
    copied = []
    stds = []
    sciobj = []
    # Get list of source files
    srcfiles = sorted(glob.glob(os.path.join(srcdir,
                                             'ifu%s_*.fits*' % datestr)))
    # Loop over source files
    for src in srcfiles:
        # get base filename
        fn = src.split('/')[-1].split('.gz')[0]
        destfl = os.path.join(destdir, fn)
        # Is our source file complete?
        if os.stat(src).st_size >= fsize or '.gz' in src:
            # Call copy
            nc, ns, nob = docp(src, destfl, skip_cals=True)
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
    logging.info("Linked/copied %d files" % ncp)

    return ncp, copied
    # END: cpsci


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
    fspec = os.path.join(srcdir, "ifu%s_0*.fits*" % sdate)
    flist = sorted(glob.glob(fspec))
    # Record number copied
    ncp = 0
    # Loop over file list
    for src in flist:
        # Get destination filename
        imf = src.split('/')[-1].split('.gz')[0]
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
        if (src_size >= fsize and src_size > loc_size) or '.gz' in src:
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
            if 'Calib' in obj and 'test' not in obj and 'Test' not in obj:
                exptime = hdr['EXPTIME']
                if exptime < 0:
                    exptime = 0.0
                    logging.warning("illegal EXPTIME value, setting to 0.0")
                lampcur = hdr['LAMPCUR']
                if lampcur < 0:
                    lampcur = 1.0
                    logging.warning("illegal LAMPCUR value, setting to 1.0")
                # Check for dome exposures
                if 'dome' in obj:
                    if exptime > 100. and ('dome' in obj and
                                           'Xe' not in obj and
                                           'Hg' not in obj and
                                           'Cd' not in obj):
                        if lampcur > 0.0:
                            # Copy dome images
                            nc, ns, nob = docp(src, destfil, onsky=False,
                                               verbose=True, header=hdr)
                            ncp += nc
                        else:
                            logging.warning("Bad dome - lamp not on: %s" % src)
                # Check for arcs
                elif 'Xe' in obj or 'Cd' in obj or 'Hg' in obj:
                    if exptime > 10.:
                        # Copy arc images
                        nc, ns, nob = docp(src, destfil, onsky=False,
                                           verbose=True, header=hdr)
                        ncp += nc
                # Check for biases
                elif 'bias' in obj:
                    if exptime <= 0.:
                        # Copy bias images
                        nc, ns, nob = docp(src, destfil, onsky=False,
                                           verbose=True, header=hdr)
                        ncp += nc

    return ncp
    # END: cpcal


def cpprecal(rawd=None, destdir=None, cdate=None, fsize=8400960):
    """Copy raw cal files from previous date's directory

    Make sure we only look in previous day directory for files created
    within four hours of the day changeover for possible raw calibration
    files required for the current night's calibration.  Copy any such
    files into the destination directory.

    Args:
        rawd (str): top level raw directory (like '/scr2/sedm/raw')
        destdir (str): where to put the files
        cdate (str): date of interest in format 'YYYYMMDD'
        fsize (int): size of completely copied file in bytes

    Returns:
        int: number of images actually copied

    """

    # Get current and previous dates
    # Do Time arithmetic
    ctime = Time(cdate[0:4]+'-'+cdate[4:6]+'-'+cdate[6:])
    one_day = TimeDelta(1.0, format='jd')
    # Previous night
    ptime = ctime - one_day
    #
    # previous source dir
    pdate = ''.join(ptime.iso.split()[0].split('-'))
    # Set the previous night as the source directory
    srcdir = os.path.join(rawd, pdate)

    # Record how many images copied
    ncp = 0
    # If there is a previous night, get those files
    if os.path.exists(srcdir):
        # Get list of previous night's raw cal files
        # (within four hours of day changeover)
        flist = sorted(glob.glob(os.path.join(srcdir,
                                              "ifu%s_2*.fits*" % pdate)))
        # Loop over file list
        for src in flist:
            if os.stat(src).st_size >= fsize or '.gz' in src:
                # Read FITS header
                ff = pf.open(src)
                hdr = ff[0].header
                ff.close()
                # Get OBJECT keyword
                obj = hdr['OBJECT']
                # Filter Calibs and avoid test images
                if 'Calib' in obj and 'test' not in obj and 'Test' not in obj:
                    # Copy cal images
                    imf = src.split('/')[-1].split('.gz')[0]
                    destfil = os.path.join(destdir, imf)
                    exptime = hdr['EXPTIME']
                    if exptime < 0:
                        exptime = 0.0
                        logging.warning("illegal EXPTIME value, setting to 0.0")
                    lampcur = hdr['LAMPCUR']
                    if lampcur < 0:
                        lampcur = 1.0
                        logging.warning("illegal LAMPCUR value, setting to 1.0")
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
                                                       header=hdr)
                                    ncp += nc
                            else:
                                logging.warning("Bad dome - lamp not on: %s"
                                                % src)
                    # Check for arcs
                    elif 'Xe' in obj or 'Cd' in obj or 'Hg' in obj:
                        if exptime > 10.:
                            # Copy arc images
                            if not os.path.exists(destfil):
                                nc, ns, nob = docp(src, destfil, onsky=False,
                                                   verbose=True, header=hdr)
                                ncp += nc
                    # Check for biases
                    elif 'bias' in obj:
                        if exptime <= 0.:
                            # Copy bias images
                            if not os.path.exists(destfil):
                                nc, ns, nob = docp(src, destfil, onsky=False,
                                                   verbose=True, header=hdr)
                                ncp += nc
            else:
                logging.warning("Truncated file: %s" % src)
    else:
        logging.warning("No previous directory: %s" % srcdir)

    return ncp
    # END: cpprecal


def cube_ready(caldir='./', cur_date_str=None, ignore_bad=False):
    """Check for all required calibration files in calibration directory.

    Args:
        caldir (str): directory to check
        cur_date_str (str): current date in YYYYMMDD format
        ignore_bad (bool): should we ignore a bad wave stat?

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
            if ignore_bad:
                logging.warning("Wavelength stats failed, proceeding anyway")
            else:
                if not os.path.exists(os.path.join(caldir, 'bad')):
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
    # Test
    dof = glob.glob(os.path.join(caldir, 'dome.fits'))
    hgf = glob.glob(os.path.join(caldir, 'Hg.fits'))
    cdf = glob.glob(os.path.join(caldir, 'Cd.fits'))
    xef = glob.glob(os.path.join(caldir, 'Xe.fits'))
    if len(dof) == 1 and len(hgf) == 1 and len(cdf) == 1 and len(xef) == 1:
        ret = True
    return ret
    # END: cal_proc_ready


def archive_old_data(odir, ut_date):
    """Remove old files and take appropriate action for archiving"""
    # Extract old positions, if available
    pos_dict = get_extract_pos(odir, ut_date)
    # First archive old kpy files (if any)
    archive_old_kpy_files(odir)
    # Now remove old pysedm cals and cubes
    delete_old_pysedm_files(odir, ut_date, keep_spec=True, keep_cubes=False)
    # Now delete old raw data
    delete_old_raw_files(odir)
    # Output positions, if available
    if len(pos_dict) > 0:
        write_extract_pos(pos_dict=pos_dict, odir=odir)
    else:
        logging.warning("No old pysedm positions found")


def delete_old_raw_files(odir):
    """Remove all the old raw input files"""

    # Get list of raw files to remove
    flist = glob.glob(os.path.join(odir, "*.fits*"))
    flist.extend(glob.glob(os.path.join(odir, '*.lst')))
    ndel = 0
    for fl in flist:
        if os.path.islink(fl):
            continue
        os.remove(fl)
        ndel += 1
    logging.info("Removed %d old raw files" % ndel)


def clean_raw_data(odir):
    """Remove unused raw files after reduction"""
    # Remove intermediate processing files
    flist = glob.glob(os.path.join(odir, "b_ifu*.fits"))
    flist.extend(glob.glob(os.path.join(odir, "ifu*.fits")))
    for fl in flist:
        os.remove(fl)
    ndel = len(flist)
    # Remove calib files
    flist = glob.glob(os.path.join(odir, "*crr_b_ifu*.fits"))
    for fl in flist:
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
            if rute.startswith('maskcrr_'):
                subprocess.call(["gzip", fl])

    return ndel


def delete_old_pysedm_files(odir, ut_date, keep_spec=False, keep_cubes=False):
    """Remove all the old pysedm output files"""
    # Potential archive directory
    archdir = os.path.join(odir, 'pysedm')
    # Keep track of how many files
    ndelfile = 0
    nkeepspec = 0
    nkeepcube = 0
    # Get list of pysedm files to remove
    flist = glob.glob(os.path.join(odir, '%s_*' % ut_date))
    flist.extend(glob.glob(os.path.join(odir, '*_crr_b_ifu%s*' % ut_date)))
    flist.extend(glob.glob(os.path.join(odir, 'bkgd_dome.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'e3d_dome.fit*')))
    flist.extend(glob.glob(os.path.join(odir, 'pysedm_run.log')))
    flist.extend(glob.glob(os.path.join(odir, 'report.txt')))
    # Remove or move them
    if len(flist) > 0:
        if keep_spec:
            # Make archive directory
            if not os.path.exists(archdir):
                os.mkdir(archdir)
        for fl in flist:
            # Skip deleting linked fluxcal file
            if 'fluxcal_' in fl and os.path.islink(fl):
                logging.info("Preserve linked fluxcal file: %s" % fl)
                continue
            # Do we keep cals and cubes?
            if keep_cubes:
                rute = fl.split('/')[-1]
                if rute.startswith(ut_date) or \
                   rute.startswith('e3d_') or rute.startswith('flex'):
                    nkeepcube += 1
                    continue
            # Keep spectra?
            if keep_spec:
                if 'spec_' in fl and 'snid' not in fl and \
                        '_ea.fit' not in fl and '_failed' not in fl and \
                        '.png' not in fl and '.pdf' not in fl or \
                        'report.txt' in fl:
                    shutil.move(fl, archdir)
                    nkeepspec += 1
                else:
                    os.remove(fl)
                    ndelfile += 1
            else:
                os.remove(fl)
                ndelfile += 1
        # If we keep the spectra, gzip them
        if keep_spec:
            # Now gzip the files in the archive
            flist = os.listdir(archdir)
            for fl in flist:
                if 'gz' not in fl:
                    subprocess.run(["gzip", os.path.join(archdir, fl)])
            logging.info("Moved %d old pysedm spectra to %s" %
                         (nkeepspec, archdir))
        if keep_cubes:
            logging.info("Kept %d old pysedm cal and cube files in %s" %
                         (nkeepcube, ut_date))
        logging.info("Deleted %d old pysedm files in %s" %
                     (ndelfile, ut_date))
    else:
        logging.warning("No pysedm files found in %s" % odir)
    # END: delete_old_pysedm_files


def archive_old_pysedm_extractions(redd=None, ut_date=None):
    """Move all the old pysedm extraction output files to ./pysedm/"""
    # Output dir
    odir = os.path.join(redd, ut_date)
    # Check for old positions
    if not os.path.exists(os.path.join(odir, 'old_positions.txt')):
        pos_dic = get_extract_pos(odir, ut_date)
        write_extract_pos(pos_dict=pos_dic, odir=odir)
    # Get list of pysedm files to move
    flist = glob.glob(os.path.join(odir, '*auto*_crr_b_ifu%s*' % ut_date))
    flist.extend(glob.glob(os.path.join(odir,
                                        "fluxcal_crr_b_ifu%s*" % ut_date)))
    flist.extend(glob.glob(os.path.join(odir,
                                        '*notfluxcal*_crr_b_ifu%s*' % ut_date)))
    flist.extend(glob.glob(os.path.join(odir, 'report*.txt')))
    flist.extend(glob.glob(os.path.join(odir, 'rp_extract.log')))
    # Move them into the archive
    if len(flist) > 0:
        # get time stamp
        mdate = datetime.datetime.fromtimestamp(os.path.getmtime(flist[0]))
        stamp = "%4d%02d%02d" % (mdate.year, mdate.month, mdate.day)
        # Archive directory
        archdir = os.path.join(odir, "pysedm_%s" % stamp)
        # Make archive directory
        if not os.path.exists(archdir):
            os.mkdir(archdir)
        nfilemv = 0
        for fl in flist:
            rute = fl.split('/')[-1]
            if not os.path.exists(os.path.join(archdir, rute)):
                shutil.move(fl, archdir)
                nfilemv += 1
        # Now gzip the files in the archive
        flist = os.listdir(archdir)
        for fl in flist:
            if 'gz' not in fl:
                subprocess.run(["gzip", os.path.join(archdir, fl)])
        logging.info("Archived %d old pysedm extraction files into %s" %
                     (nfilemv, archdir))
    else:
        logging.warning("No pysedm extraction files found in %s" % odir)
    # END: archive_old_pysedm_files


def archive_old_kpy_files(odir):
    """Move all the old kpy output files to ./kpy/"""
    # Potential archive directory
    archdir = os.path.join(odir, 'kpy')
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
    nfilemv = 0
    if len(flist) > 0:
        # Add report
        flist.extend(glob.glob(os.path.join(odir, 'report.txt')))
        # Make archive directory
        if not os.path.exists(archdir):
            os.mkdir(archdir)
        for fl in flist:
            rute = fl.split('/')[-1]
            if not os.path.exists(os.path.join(archdir, rute)):
                shutil.move(fl, archdir)
                nfilemv += 1
        logging.info("Moved %d old kpy files into %s" % (nfilemv, archdir))
    else:
        logging.warning("No old kpy files found in %s" % odir)
    # Now check spec files
    flist = glob.glob(os.path.join(odir, 'spec_*'))
    nspecmv = 0
    if len(flist) > 0:
        # Make archive directory
        if not os.path.exists(archdir):
            os.mkdir(archdir)
        for fl in flist:
            if '_crr_b_ifu' not in fl:
                shutil.move(fl, archdir)
                nspecmv += 1
        logging.info("Moved %d old kpy spec files into %s" % (nspecmv, archdir))
    else:
        logging.warning("No spec files found in %s" % odir)
    # Now gzip the files in the archive
    if os.path.exists(archdir) and nfilemv+nspecmv > 0:
        flist = os.listdir(archdir)
        logging.info("gzipping old kpy files...")
        for fl in flist:
            if '.gz' not in fl:
                subprocess.run(["gzip", os.path.join(archdir, fl)])
    # END: archive_old_kpy_files


def dosci(destdir='./', datestr=None, nodb=False, posdic=None, oldext=False):
    """Copies new science ifu image files from srcdir to destdir.

    Searches for most recent ifu image in destdir and looks for and
    copies any ifu images in srcdir that are newer and complete.
    Then bias subtracts and CR rejects the copied images.  If any are standard
    star observations, process them as well.

    Args:
        destdir (str): destination directory (typically in /scr2/sedm/redux)
        datestr (str): YYYYMMDD date string
        nodb (bool): if True no update to SEDM db
        posdic (dict): position dictionary from get_extract_pos function
        oldext (bool): True to use old extract_star instead of new extractstar

    Returns:
        int: Number of ifu images actually copied

    """

    # Record copies and standard star observations
    nextr = 0
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
                obj = hdr['OBJECT'].split("[")[0].strip().replace(" ", "-")
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
            # are we a standard star?
            if 'STD-' in obj:
                e3d_good = make_e3d(fnam=fl, destdir=destdir, datestr=datestr,
                                    nodb=nodb, sci=False, hdr=None)
                if e3d_good:
                    # Get seeing
                    seeing = rcimg.get_seeing(imfile=fn, destdir=destdir,
                                              save_fig=True)
                    # Use auto psf extraction for standard stars
                    if seeing > 0:
                        logging.info("seeing measured as %f" % seeing)
                        if oldext:
                            cmd = ("extract_star.py", datestr, "--auto", fn,
                                   "--std", "--tag", "robot", "--maxpos")
                        else:
                            cmd = ("extractstar.py", datestr, "--auto", fn,
                                   "--std", "--tag", "robot",
                                   "--centroid", "brightest",
                                   "--seeing", "%.2f" % seeing)
                    else:
                        logging.info("seeing not measured for %s" % fn)
                        if oldext:
                            cmd = ("extract_star.py", datestr, "--auto", fn,
                                   "--std", "--tag", "robot", "--maxpos")
                        else:
                            cmd = ("extractstar.py", datestr, "--auto", fn,
                                   "--std", "--tag", "robot",
                                   "--centroid", "brightest")
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
                        nextr += 1
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
                        if oldext:
                            fcspec = "fluxcal_auto_robot_lstep1__%s_*.fits" % \
                                     fn.split('.')[0]
                        else:
                            fcspec = "fluxcal_%s_*.fits" % fn.split('.')[0]
                        flxcal = glob.glob(os.path.join(destdir, fcspec))
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
                                    nodb=nodb, sci=True, hdr=hdr)
                if e3d_good:
                    logging.info("e3d science cube generated")
                    # Do we have a position
                    pos_key = fn.split('_', 2)[-1].split('.')[0]
                    if pos_key not in posdic:
                        logging.info("No position for %s" % fn)
                        continue
                    # Position
                    xpos = posdic[pos_key][0]
                    ypos = posdic[pos_key][1]
                    # Get seeing
                    seeing = rcimg.get_seeing(imfile=fn, destdir=destdir,
                                              save_fig=True)
                    # Use forced psf for science targets
                    if seeing > 0:
                        logging.info("seeing measured as %f" % seeing)
                        if oldext:
                            cmd = ("extract_star.py", datestr, "--auto", fn,
                                   "--autobins", "6", "--tag", "robot",
                                   "--centroid", "%.2f" % xpos, "%.2f" % ypos)
                        else:
                            cmd = ("extractstar.py", datestr, "--auto", fn,
                                   "--autobins", "6", "--tag", "robot",
                                   "--centroid", "%.2f" % xpos, "%.2f" % ypos,
                                   "--seeing", "%.2f" % seeing)
                    else:
                        logging.info("seeing not measured for %s" % fn)
                        if oldext:
                            cmd = ("extract_star.py", datestr, "--auto", fn,
                                   "--autobins", "6", "--tag", "robot",
                                   "--centroid", "%.2f" % xpos, "%.2f" % ypos)
                        else:
                            cmd = ("extractstar.py", datestr, "--auto", fn,
                                   "--autobins", "6", "--tag", "robot",
                                   "--centroid", "%.2f" % xpos, "%.2f" % ypos)
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
                        nextr += 1
                        logging.info("Running SNID for " + fn)
                        cmd = ("make", "classify")
                        logging.info(" ".join(cmd))
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            logging.error("Error running SNID")
                        cmd = ("pysedm_report.py", datestr, "--contains",
                               fn.split('.')[0])
                        logging.info(" ".join(cmd))
                        retcode = subprocess.call(cmd)
                        if retcode != 0:
                            logging.error("Error running report for " + fn)
                        # run Verify.py
                        cmd = "~/sedmpy/drpifu/Verify.py %s --contains %s" % \
                              (datestr, fn.split('.')[0])
                        subprocess.call(cmd, shell=True)
                        # update database
                        proced = glob.glob(os.path.join(destdir, procfn))[0]
                        if os.path.exists(proced):
                            if nodb:
                                logging.warning("No update of spec in SEDM db")
                            else:
                                # Update SedmDb table spec
                                spec_id = update_spec(proced)
                                logging.info("update of %s with spec_id %d" %
                                             (proced, spec_id))
                        else:
                            logging.error("Not found: %s" % proced)
                else:
                    logging.error("Cannot perform extraction for %s" % fn)
    return nextr
    # END: dosci


def write_extract_pos(pos_dict=None, odir=None):
    """Write old positions to old_positions.txt"""
    ofil = os.path.join(odir, 'old_positions.txt')
    with open(ofil, 'w') as fh:
        logging.info("Writing positions to old_positions.txt")
        for k, v in pos_dict.items():
            rec = "%s %.2f %.2f" % (k, v[0], v[1])
            logging.info(rec)
            fh.write(rec+"\n")
    logging.info("%d old positions written to old_positions.txt" %
                 len(pos_dict))


def read_extract_pos(indir):
    """Read old_positions.txt into positions dictionary"""
    pos_dict = {}
    ndic = 0
    pfile = os.path.join(indir, 'old_positions.txt')
    if os.path.exists(pfile):
        with open(pfile, 'r') as pfh:
            for cnt, line in enumerate(pfh):
                parts = line.split()
                try:
                    pos_dict[parts[0]] = (float(parts[1]), float(parts[2]))
                    ndic += 1
                except:
                    continue
    else:
        logging.warning("old_positions.txt not found in %s" % indir)
    return pos_dict


def get_extract_pos(indir, indate):
    """Create a dictionary of A positions from pysedm spec_*.fits files"""
    out_dict = {}
    redo_dict = {}
    # Get input crr_b_ifu files
    flist = glob.glob(os.path.join(indir, 'spec_*_crr_b_ifu%s_*.fit*' % indate))
    # loop over input spec files
    for specfs in flist:
        # Skip unwanted spectra
        if '_aperture_' in specfs:
            continue
        if '_ea' in specfs:
            continue
        if '_failed' in specfs:
            continue
        # Get fits header
        ff = pf.open(specfs)
        header = ff[0].header
        ff.close()
        # Get OBJECT keyword
        try:
            objname = header['OBJECT']
        except KeyError:
            objname = 'Test'
        # Skip cal, STD files
        if 'Calib' in objname:
            continue
        if 'STD-' in objname:
            continue
        # Get image ID
        rute = "ifu" + "_".join(specfs.split('_ifu')[-1].split('_', 4)[:4])

        # Find good positions
        if '_redo' in specfs:
            try:
                redo_xpos = header['XPOS']
                redo_ypos = header['YPOS']
                redo_dict[rute] = (redo_xpos, redo_ypos)
            except KeyError:
                continue
        else:
            try:
                xpos = header['XPOS']
                ypos = header['YPOS']
                out_dict[rute] = (xpos, ypos)
            except KeyError:
                continue
    # END for specfs in flist:
    # Check redo dictionary
    if len(redo_dict) > 0:
        for k, v in redo_dict.items():
            out_dict[k] = v

    logging.info("Found %d positions in %s" % (len(out_dict), indate))
    return out_dict
    # END: get_extract_pos


def re_extract(redd=None, ut_date=None, nodb=False, oldext=False,
               ignore_bad=False):
    """Re-extract spectra from cube files for one night.

    Args:
        redd (str): reduced directory (something like /scr2/sedmdrp/redux)
        ut_date (str): input directory for single night processing
        nodb (bool): True if no update to SEDM db
        oldext (bool): True to use old extract_star instead of new extractstar
        ignore_bad (bool): extract in spite of bad cals?

    Returns:
        bool: True if cals completed normally, False otherwise

    """
    # Default return value
    ret = False
    # Report extraction version
    if oldext:
        logging.info("Using old extract_star.py routine")
    else:
        logging.info("Using new extractstar.py routine")
    # Output directory is based on redd and ut_date
    outdir = os.path.join(redd, ut_date)
    # change to directory
    os.chdir(outdir)
    # report
    logging.info("Extracted spectra to: %s" % outdir)
    # Get kpy position dictionary
    pos_dic = read_extract_pos(outdir)
    if len(pos_dic) <= 0:
        logging.warning("No pysedm positions found")
    else:
        logging.info("%d pysedm positions found" % len(pos_dic))

    # Check calibration status
    if not cube_ready(outdir, ut_date, ignore_bad=ignore_bad):
        # Report status
        logging.error("No calibrations!  Please run ReProcess.py --calibrate.")
    else:
        logging.info("Calibrations present in %s" % outdir)
        # Process observations
        nextr = dosci(outdir, datestr=ut_date, nodb=nodb, posdic=pos_dic,
                      oldext=oldext)
        logging.info("%d spectra extracted." % nextr)
        # Re-gzip input files
        cmd = ["gzip crr_b_ifu%s*.fits" % ut_date]
        logging.info(cmd)
        subprocess.call(cmd, shell=True)

    return ret
    # END: re_extract


def re_cube(redd=None, ut_date=None, nodb=False, ignore_bad=False):
    """Create e3d cube files for one night.

    Args:
        redd (str): reduced directory (something like /scr2/sedmdrp/redux)
        ut_date (str): input directory for single night processing
        nodb (bool): True if no update to SEDM db
        ignore_bad (boo): should we ignore bad wave stats?

    Returns:
        int: Number of cubes made

    """
    # Output directory is based on redd and ut_date
    destdir = os.path.join(redd, ut_date)
    # Report
    logging.info("Cube files to: %s" % destdir)
    # Count cubes attempted and successfully made
    ncube = 0
    ntry = 0
    # Check calibration status
    if not cube_ready(destdir, ut_date, ignore_bad=ignore_bad):
        logging.error("No calibrations!  Please run ReProcess.py --calibrate.")
    else:
        logging.info("Calibrations present in %s" % destdir)
        # Proceed to build e3d cubes
        # Get list of source files in destination directory
        srcfiles = sorted(glob.glob(os.path.join(destdir, 'crr_b_ifu*.fits')))
        # Loop over source files
        for fl in srcfiles:
            # Read FITS header
            ff = pf.open(fl)
            hdr = ff[0].header
            ff.close()
            # Get OBJECT keyword
            try:
                obj = hdr['OBJECT']
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
            # Now we try for a cube
            ntry += 1
            # are we a standard star?
            if 'STD-' in obj:
                if make_e3d(fnam=fl, destdir=destdir, datestr=ut_date,
                            nodb=nodb):
                    ncube += 1
            else:
                # Build cube for science observation
                if make_e3d(fnam=fl, destdir=destdir, datestr=ut_date,
                            nodb=nodb, sci=True, hdr=hdr):
                    ncube += 1
        # END: for fl in srcfiles:
        logging.info("%d cubes made out of %d attempted" % (ncube, ntry))

    return ncube
    # END: re_cube


def re_calib(redd=None, ut_date=None, nodb=False):
    """Create calibration files for one night.

    Args:
        redd (str): reduced directory (something like /scr2/sedmdrp/redux)
        ut_date (str): input directory for single night processing
        nodb (bool): True if no update to SEDM db

    Returns:
        bool: True if cals completed normally, False otherwise

    """
    # Output directory is based on redd and ut_date
    outdir = os.path.join(redd, ut_date)
    # Go there
    os.chdir(outdir)
    # report
    logging.info("Calibration files to: %s" % outdir)

    # Check calibration status
    if not cube_ready(outdir, ut_date):
        # Process calibrations
        if cal_proc_ready(outdir):
            # Process calibration
            start_time = time.time()
            cmd = ("ccd_to_cube.py", ut_date, "--tracematch",
                   "--hexagrid")
            logging.info(" ".join(cmd))
            subprocess.call(cmd)
            procg_time = int(time.time() - start_time)
            if os.path.exists(
               os.path.join(outdir, ut_date + '_HexaGrid.pkl')):
                # Process wavelengths
                start_time = time.time()
                # Spawn nsub sub-processes to solve wavelengths faster
                nsub = 8
                cmd = ("derive_wavesolution.py", ut_date,
                       "--nsub", "%d" % nsub)
                logging.info(" ".join(cmd))
                subprocess.Popen(cmd)
                time.sleep(60)
                # Get a list of solved spaxels
                wslist = glob.glob(os.path.join(outdir, ut_date +
                                                '_WaveSolution_range*.pkl'))
                # Wait until they are all finished
                nfin = len(wslist)
                while nfin < nsub:
                    time.sleep(60)
                    wslist = glob.glob(
                        os.path.join(outdir, ut_date +
                                     '_WaveSolution_range*.pkl'))
                    if len(wslist) != nfin:
                        print("\nFinished %d out of %d parts"
                              % (len(wslist), nsub))
                        nfin = len(wslist)
                    else:
                        print(".", end="", flush=True)
                logging.info("Finished all %d parts, merging..." % nsub)
                # Merge the solutions
                subprocess.call(("derive_wavesolution.py", ut_date,
                                 "--merge"))
                procw_time = int(time.time() - start_time)
                if os.path.exists(
                   os.path.join(outdir, ut_date + '_WaveSolution.pkl')):
                    # Process flat
                    start_time = time.time()
                    cmd = ("ccd_to_cube.py", ut_date, "--flat")
                    logging.info(" ".join(cmd))
                    subprocess.call(cmd)
                    if not (os.path.exists(
                            os.path.join(outdir, ut_date + '_Flat.fits'))):
                        logging.info("Making of %s_Flat.fits failed!"
                                     % ut_date)
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
                            spec_calib_id = update_calibration(ut_date)
                            logging.info(
                                "SEDM db accepted spec_calib at id %d"
                                % spec_calib_id)
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
        if not cube_ready(outdir, ut_date):
            logging.error("These calibrations failed!")
    else:
        logging.info("Calibrations already present in %s" % outdir)
    # END: re_calib


def re_reduce(rawd=None, redd=None, ut_date=None, noprecal=False):
    """Relink/copy raw files and produce reduced cal files"""
    # Input/output directories
    idir = os.path.join(rawd, ut_date)
    odir = os.path.join(redd, ut_date)

    # Make sure we are in our output directory
    os.chdir(odir)

    # Archive/delete any old data
    archive_old_data(odir, ut_date)

    if noprecal:
        logging.info("Skipping pre-cal check")
        npre = 0
    else:
        # Check for precal files
        npre = cpprecal(rawd=rawd, destdir=odir, cdate=ut_date)
        logging.info("%d pre-cal images copied" % npre)

    # Now get any cals from this date
    ncal = cpcal(srcdir=idir, destdir=odir)
    logging.info("%d cal images copied" % ncal)

    # Now get science images from this date
    nsci, sci = cpsci(srcdir=idir, destdir=odir, datestr=ut_date)
    logging.info("%d sci images copied" % nsci)

    # Now do basic reduction
    proc_bias_crrs(npre+ncal+nsci)

    # Make cal images
    logging.info("Calling make calimgs")
    subprocess.call(("make", "calimgs"))

    # Now clean up raw files
    ncln = clean_raw_data(odir)
    logging.info("%d raw files cleaned" % ncln)

    logging.info("Reduction complete.")
    # END: re_reduce


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Make cals for PYSEDM pipeline

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--rawdir', type=str, default=_rawpath,
                        help='Input raw directory (%s)' % _rawpath)
    parser.add_argument('--reduxdir', type=str, default=_reduxpath,
                        help='Output reduced directory (%s)' % _reduxpath)
    parser.add_argument('--date', type=str, default=None,
                        help='Select date to process')
    parser.add_argument('--nodb', action="store_true", default=False,
                        help='Do not update SEDM Db')
    parser.add_argument('--archive_pysedm', action="store_true", default=False,
                        help='Archive any existing pysedm files')
    parser.add_argument('--oldext', action="store_true", default=False,
                        help='Use old extract_star.py for extraction')
    parser.add_argument('--reduce', action="store_true", default=False,
                        help='Re-reduce raw images')
    parser.add_argument('--skip_precal', action="store_true", default=False,
                        help='Skip check for pre-cal data from previous night')
    parser.add_argument('--calibrate', action="store_true", default=False,
                        help='Re-make calibrations')
    parser.add_argument('--cube', action="store_true", default=False,
                        help='Re-make e3d*.fits cubes')
    parser.add_argument('--ignore_bad', action="store_true", default=False,
                        help='Proceed in spite of bad calibration')
    parser.add_argument('--extract', action="store_true", default=False,
                        help='Re-extract targets')

    args = parser.parse_args()

    if not args.date:
        logging.error("Must provide a YYYYMMDD date with --date")
    else:
        if args.archive_pysedm:
            archive_old_pysedm_extractions(redd=args.reduxdir,
                                           ut_date=args.date)
        elif args.reduce:
            re_reduce(rawd=args.rawdir, redd=args.reduxdir, ut_date=args.date,
                      noprecal=args.skip_precal)
        elif args.calibrate:
            re_calib(redd=args.reduxdir, ut_date=args.date, nodb=args.nodb)
        elif args.cube:
            re_cube(redd=args.reduxdir, ut_date=args.date, nodb=args.nodb,
                    ignore_bad=args.ignore_bad)
        elif args.extract:
            re_extract(redd=args.reduxdir, ut_date=args.date, nodb=args.nodb,
                       oldext=args.oldext, ignore_bad=args.ignore_bad)
        else:
            logging.error("Must specify one of "
                          "--[reduce|calibrate|cube|extract|archive_pysedm]")
