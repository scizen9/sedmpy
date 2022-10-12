import glob
import os
import sys
import time
import subprocess
import numpy as np
import json
import datetime

import sedmpy_version

# Get pipeline configuration
# Find config file: default is sedmpy/config/sedmconfig.json
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(sedmpy_version.CONFIG_DIR, 'sedmconfig.json')
with open(configfile) as config_file:
    sedm_cfg = json.load(config_file)

_logpath = sedm_cfg['paths']['logpath']
_photpath = sedm_cfg['paths']['photpath']


def seeing_stats(statfile=None):
    """return min, mean, max seeing from guider images"""

    if not statfile:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        statfile = os.path.join(_photpath, timestamp, "stats/stats.log")

    if os.path.isfile(statfile):
        count = 0
        while os.path.getsize(statfile) == 0 and count < 5:
            count += 1
            print("stats file empty, retrying ", count)
            time.sleep(1)
        s = np.genfromtxt(statfile, delimiter=",", dtype=None, encoding=None)
    else:
        return None, None, None, None

    try:
        s.sort(order="f2")
    except ValueError:
        return None, None, None, None

    s = s[s["f3"] > 1]
    s = s[s["f9"] == 'GUIDER']

    see_min = s['f4'].min()
    see_max = s['f4'].max()
    see_med = np.median(s['f4'])
    see_std = np.nanstd(s['f4'])

    return see_min, see_med, see_max, see_std


def report():
    """Generate DRP report using output spec_*.txt files"""
    # comment?
    if len(sys.argv) > 1:
        comment = input("Enter comment: ")
    else:
        comment = None

    obsdate = os.getcwd().split('/')[-1]
    flist = glob.glob("spec_*.txt")
    flist.sort()
    print("\nReport generated on %s" % time.strftime("%c"))
    if comment:
        print("\n%s" % comment)
    print("\nSEDM DRP run in %s\nFound %d spec_*.txt files" %
          (os.getcwd(), len(flist)))
    print("See http://minar.caltech.edu/data_access/ifu?obsdate=%s\n" % obsdate)

    statsfile = os.path.join(_photpath, obsdate, 'stats/stats.log')
    see_min, see_mean, see_max, see_std = seeing_stats(statsfile)
    if see_mean is not None:
        print("RCam Seeing (FWHM) median, StdDev, min, max: %.2f, %.2f, "
              "%.2f, %.2f\n" % (see_mean, see_std, see_min, see_max))
    else:
        print("Seeing data unavailable\n")

    print("UTStart  Object               Exptime Air    Qual"
          "                   method  Allocation           "
          "Type            z        Rlap   SNIaS SNIaSe   SNIaZ SNIaZe")
    recs = []
    for f in flist:
        # Get object name
        tname = f.split('_ifu')[-1].split('_')[4:]
        if len(tname) > 1:
            objname = '_'.join(tname).split('.txt')[0]
        else:
            objname = tname[0].split('.txt')[0]
        # Get time string
        tstr = ':'.join(f.split('_ifu')[-1].split('_')[1:4])
        # check the ascii spectrum file for SNID data
        with open(f, "r") as sfl:
            lines = sfl.readlines()
            # check for SNID classification
            clas = [li for li in lines if "SNIDMATCHTYPE" in li]
            ctype = ""
            if len(clas) > 0:
                for cl in clas:
                    ctype += (" %s" % cl.split()[-1])
            # check for SNID subtype
            clas = [li for li in lines if "SNIDMATCHSUBTYPE" in li]
            stype = ""
            if len(clas) > 0:
                for cl in clas:
                    st = cl.split()[-1]
                    if st != '-':
                        stype += ("%s" % st)
                    else:
                        stype = ""
                if len(stype) > 0:
                    ctype += "-"
                    ctype += stype
            # get SNID redshift
            zmch = [li for li in lines if "REDSHIFT" in li]
            if len(zmch) > 0:
                zmch = ("%.4f" % float(zmch[0].split()[-1]))
            else:
                zmch = ""
            # get SNID rlap
            rlap = [li for li in lines if "RLAP" in li]
            if len(rlap) > 0:
                rlap = ("%.2f" % float(rlap[0].split()[-1]))
            else:
                rlap = ""
            # get SNIascore
            snia_score = [li for li in lines if "SNIASCORE" in li]
            if len(snia_score) > 0:
                snia_score = ("%.3f" % float(snia_score[0].split()[-1]))
            else:
                snia_score = ""
            # get SNIascore err
            snia_score_err = [li for li in lines if "SNIASCORE_ERR" in li]
            if len(snia_score_err) > 0:
                snia_score_err = ("%.3f" % float(snia_score_err[0].split()[-1]))
            else:
                snia_score_err = ""
            # get SNIascore redshift
            snia_z = [li for li in lines if "SNIASCORE_Z" in li]
            if len(snia_z) > 0:
                snia_z = ("%.3f" % float(snia_z[0].split()[-1]))
            else:
                snia_z = ""
            # get SNIascore redshift err
            snia_z_err = [li for li in lines if "SNIASCORE_ZERR" in li]
            if len(snia_z_err) > 0:
                snia_z_err = ("%.3f" % float(snia_z_err[0].split()[-1]))
            else:
                snia_z_err = ""
            # get flux calibration
            flxcal = [li for li in lines if "FLUXCAL" in li]
            if len(flxcal) > 0:
                flxcal = flxcal[0].split()[-1]
            else:
                flxcal = "False"
            # get method
            meth = f.split('__crr')[0].split('spec_')[-1]
            meth = meth.split('auto_')[-1]
            # get exposure time
            expt = [li for li in lines if "EXPTIME" in li]
            if len(expt) > 0:
                expt = ("%.1f" % float(expt[0].split()[-1]))
            else:
                expt = "-"
            # get airmass
            air = [li for li in lines if "AIRMASS" in li]
            if len(air) > 0:
                air = ("%.3f" % float(air[0].split()[-1]))
            else:
                air = "-"
            # get program ID
            prid = [li for li in lines if "P60PRID" in li]
            if len(prid) > 0:
                prid = prid[0].split()[-1]
            else:
                prid = "-"
            # get Quality
            quality = [li for li in lines if "QUALITY" in li]
            if len(quality) > 0:
                quality = int(quality[0].split()[-1])
            else:
                quality = 9
            sfl.close()
        if ctype == "" or quality == 5:
            if "STD" in f:
                ctype = " STD"
            else:
                ctype = " QUALITY_%d" % quality
                stype = ""
                zmch = ""
                rlap = ""
                snia_score = ""
                snia_score_err = ""
                snia_z = ""
                snia_z_err = ""

        recs.append("%8s %-20s %7s %5s     %d %24s  %-18s  %-14s  "
                    "%-8s %5s  %6s %6s  %6s %6s" %
                    (tstr, objname, expt, air, quality, meth, prid, ctype,
                     zmch, rlap, snia_score, snia_score_err,
                     snia_z, snia_z_err))
    recs.sort()
    for r in recs:
        print(r)
    # Check for contsep uploads
    flist = glob.glob("spec_auto_contsep_*.upl")
    if len(flist) > 0:
        print("\nThere were/was %d contsep spectra uploaded to the marshal"
              % len(flist))
        for f in flist:
            tstr = ':'.join(f.split('ifu')[-1].split('_')[1:4])
            tok = f.split("_")[-1].split(".")[0]
            print("%8s %-25s contsep" % (tstr, tok))
    else:
        print("\nThere were no contsep spectra uploaded")
    # Check for failed extractions
    flist = glob.glob("spec_*failed.fits")
    if len(flist) > 0:
        flist.sort()
        print("\nThere were/was %d failed extraction(s)" % len(flist))
        for f in flist:
            tstr = ':'.join(f.split('ifu')[-1].split('_')[1:4])
            tok = f.split("_failed")[0].split("_ifu")[-1]
            try:
                out = subprocess.check_output(('grep', tok, 'what.list'),
                                              universal_newlines=True)
                out = out.split()[3]
            except subprocess.CalledProcessError:
                out = tok
            print("%8s %-25s FAILED" % (tstr, out))
    else:
        print("\nThere were no failed extractions")


if __name__ == '__main__':
    report()
