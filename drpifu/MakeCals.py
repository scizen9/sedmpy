"""Conduct automatic reduction of SEDM data in sedmdrp@pharos

Functions
    * :func:`cal_loop`     one night calibration loop
    * :func:`cal_proc_ready`  check if all required raw cal images are present
    * :func:`cube_ready`      check if all required cal files are present
    * :func:`update_calibration`    update cal cube in SEDM db

Note:
    This is used as a python script as follows::

        usage: MakeCals.py [-h] [--rawdir RAWDIR] [--reduxdir REDUXDIR]

        optional arguments:
          -h, --help           show this help message and exit
          --reduxdir REDUXDIR  Output reduced directory (/scr2/sedm/redux)
          --date YYYYMMDD      Select date to process (None)
          --nodb               Do not update SEDM Db (False)

"""
import time
import glob
import os
import re
import subprocess
import logging
import argparse
import db.SedmDb

from configparser import ConfigParser
import codecs

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


def cal_loop(redd=None, indir=None, nodb=False):
    """Create calibration files for one night.

    Args:
        redd (str): reduced directory (something like /scr2/sedm/redux)
        indir (str): input directory for single night processing
        nodb (bool): True if no update to SEDM db

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

    # Check if processed cal files are ready
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
                    procf_time = int(time.time() - start_time)
                    # Report times
                    logging.info("Calibration processing took "
                                 "%d s (grid), %d s (waves),  and %d s (flat)" %
                                 (procg_time, procw_time, procf_time))
                    if nodb:
                        logging.warning("Not updating SEDM db")
                    else:
                        # Update spec_calib table in sedmdb
                        spec_calib_id = update_calibration(cur_date_str)
                        logging.info(
                            "SEDM db accepted spec_calib at id %d"
                            % spec_calib_id)

                    logging.info(
                        "Calibration stage complete, ready for science!")
                    # Re-gzip cal images
                    subprocess.run(["gzip", "dome.fits", "Hg.fits", "Cd.fits",
                                    "Xe.fits"])
                else:
                    logging.error("Making of wavesolution failed!")
            else:
                logging.error("Making of tracematch and hexagrid failed!")
        # Check status
        if not cube_ready(outdir, cur_date_str):
            logging.error("These calibrations failed!")
    else:
        logging.info("Calibrations already present in %s" % outdir)

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

    args = parser.parse_args()

    if not args.date:
        logging.error("Must provide a YYYYMMDD date with --date")
    else:
        cal_loop(redd=args.reduxdir, indir=args.date, nodb=args.nodb)
