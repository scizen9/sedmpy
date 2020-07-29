#! /usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import os
import argparse
import logging
import subprocess
import astropy.io.fits as pf
from astropy.time import Time

from configparser import ConfigParser
import codecs

try:
    import Version
except ImportError:
    import drpifu.Version as Version

drp_ver = Version.ifu_drp_version()

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
# Get path
_rawpath = cfg_parser.get('paths', 'rawpath')


header_types = {"TELESCOP": str, "LST": str,
                "MJD_OBS": float, "JD": float,
                "APPEQX": float, "EQUINOX": int, "TEL_HA": str,
                "RA": str, "TEL_RA": str, "DEC": str, "TEL_DEC": str,
                "TEL_AZ": float, "TEL_EL": float, "AIRMASS": float,
                "TEL_PA": float, "RA_RATE": float, "DEC_RATE": float,
                "RA_OFF": float, "DEC_OFF": float,
                "TELHASP": float, "TELDECSP": float,
                "HA_REFR": float, "DEC_REFR": float, "FOCPOS": float,
                "IFUFOCUS": float, "IFUFOC2": float, "DOMEST": str,
                "DOMEMO": str, "DOME_GAP": float, "DOMEAZ": float,
                "WSCRMO": str, "TELCONT": str,
                "LAMPSTAT": str, "LAMPCUR": float,
                "HG_LAMP": str, "XE_LAMP": str, "CD_LAMP": str,
                "TELPOWST": str, "OILSTAT": str, "WEASTAT": str, "SUNSTAT": str,
                "REMOTST": str, "TELRDST": str, "HAAX_ST": str, "FOCSTAT": str,
                "DEC_AX": str, "OBJECT": str, "OBJTYPE": str, "IMGTYPE": str,
                "OBJNAME": str, "OBJEQX": str, "OBJRA": str, "OBJDEC": str,
                "ORA_RAT": float, "ODEC_RAT": float,
                "SUNRISE": str, "SUNSET": str,
                "TEL_MO": str, "WSCR_EL": float, "SOL_RA": str, "SOL_DEC": str,
                "WIND_DIR": str, "WSP_CUR": float, "WSP_AVG": float,
                "OUT_AIR": float, "OUT_HUM": float, "OUT_DEW": float,
                "IN_AIR": float, "IN_HUM": float, "IN_DEW": float,
                "MIR_TEMP": float, "TOP_AIR": float, "PRI_TEMP": float,
                "SEC_TEMP": float, "FLO_TEMP": float, "BOT_TEMP": float,
                "MID_TEMP": float, "TOP_TEMP": float, "WETNESS": float,
                "FILTER": str, "ABPAIR": bool, "IMGSET": str, "NAME": str,
                "UTC": str, "OBSDATE": str, "OBSTIME": str, "P60PRID": str,
                "P60PRNM": str, "P60PRPI": str, "EMAIL": str, "INSTRUME": str,
                "REQ_ID": int, "OBJ_ID": int, "CRVAL1": float, "CRVAL2": float,
                "BARPRESS": float, "SECPRESS": float,
                "END_SHUT": str, "ENDAIR": float, "ENDDOME": str,
                "END_RA": str, "END_DEC": str, "END_PA": float,
                "ENDBARPR": float, "ENDSECPR": float, "ELAPTIME": float}


def date_time_from_filename(fname):
    # get date and time from filename
    if 'ifu' in fname:
        date_str = fname.split('ifu')[-1].split('_')[0]
    elif 'rc' in fname:
        date_str = fname.split('rc')[-1].split('_')[0]
    else:
        date_str = '11111111'
    date = "-".join([date_str[:4], date_str[4:6], date_str[6:]])
    try:
        time_str = fname.split('_', 1)[-1].split('.')[0]
        time = ":".join([time_str.split('_')[0],
                         time_str.split('_')[1],
                         time_str.split('_')[2]])
    except IndexError:
        time = "00:00:00"
    return date, time


def sedm_fix_header(fname):
    """Make sure required keywords are present and correct"""
    try:
        ff = pf.open(fname, 'update')
    except OSError:
        logging.error("Unable to read %s" % fname)
        return False
    # Get root filename
    rute = fname.split('/')[-1]
    # Make sure header is correct
    try:
        if len(ff[0].header) < 30:
            logging.warning("Short header: %s" % rute)
            return False
    except IndexError:
        logging.error("Problem with file: %s" % rute)
        return False
    # Count updated keywords
    nupkey = 1
    # Equinox
    ff[0].header['EQUINOX'] = 2000
    # lamp status
    lamps_dic = {'LAMPSTAT': 'off', 'HG_LAMP': 'off',
                 'CD_LAMP': 'off', 'XE_LAMP': 'off'}
    # image type
    if 'OBJECT' in ff[0].header:
        obj = ff[0].header['OBJECT'].split('[')[0].strip()
        if 'STD-' in obj:
            ff[0].header['IMGTYPE'] = 'Standard'
            ff[0].header['ABPAIR'] = False
            ff[0].header['NAME'] = obj
            nupkey += 3
        elif 'Calib' in obj:
            if 'dome' in obj:
                ff[0].header['IMGTYPE'] = 'dome'
                nupkey += 1
                lamps_dic['LAMPSTAT'] = 'on'
            elif 'bias' in obj:
                ff[0].header['IMGTYPE'] = 'bias'
                nupkey += 1
            else:
                ff[0].header['IMGTYPE'] = 'lamp'
                nupkey += 1
                if 'Hg' in obj:
                    lamps_dic['HG_LAMP'] = 'on'
                elif 'Cd' in obj:
                    lamps_dic['CD-LAMP'] = 'on'
                elif 'Xe' in obj:
                    lamps_dic['XE_LAMP'] = 'on'
                else:
                    logging.warning("Unk. lamp: %s" % rute)
            ff[0].header['NAME'] = obj
            ff[0].header['ABPAIR'] = False
            nupkey += 2
        else:
            ff[0].header['IMGTYPE'] = 'Science'
            ff[0].header['NAME'] = obj.replace(" ", "-")
            nupkey += 2
            # ABPAIR status (assume true unless otherwise)
            if 'ABPAIR' not in ff[0].header:
                ff[0].header['ABPAIR'] = True
                nupkey += 1
            else:
                if type(ff[0].header['ABPAIR']) != bool:
                    if type(ff[0].header['ABPAIR']) == str:
                        ff[0].header['ABPAIR'] = ('T' in ff[0].header['ABPAIR'])
                    else:
                        ff[0].header['ABPAIR'] = True
                    nupkey += 1
    else:
        obj = ''
    # Set lamps
    for k, v in lamps_dic.items():
        # Check current value
        if k in ff[0].header:
            # If null, replace with dic value
            if len(ff[0].header[k]) <= 0:
                ff[0].header[k] = v
                nupkey += 1
        # If not in header, use dic value
        else:
            ff[0].header[k] = v
            nupkey += 1
    # RA rate
    if 'RA_RATE' not in ff[0].header:
        if 'RARATE' in ff[0].header:
            ra_rate = float(ff[0].header['RARATE'])
            ff[0].header['RA_RATE'] = ra_rate
        else:
            logging.warning("No RARATE: %s" % rute)
            ff[0].header['RA_RATE'] = 0.
        nupkey += 1
    # DEC rate
    if 'DEC_RATE' not in ff[0].header:
        if 'DECRATE' in ff[0].header:
            dec_rate = float(ff[0].header['DECRATE'])
            ff[0].header['DEC_RATE'] = dec_rate
        else:
            logging.warning("No DECRATE: %s" % rute)
            ff[0].header['DEC_RATE'] = 0.
        nupkey += 1
    # Humidity
    if 'IN_HUM' not in ff[0].header:
        if 'Inside_Rel_Hum' in ff[0].header:
            rel_hum = ff[0].header['Inside_Rel_Hum']
            ff[0].header['IN_HUM'] = rel_hum
        else:
            ff[0].header['IN_HUM'] = -999.
            logging.warning("No humidity: %s" % rute)
        nupkey += 1
    # Temperature
    if 'IN_AIR' not in ff[0].header:
        if 'Inside_Air_Temp' in ff[0].header:
            in_temp = ff[0].header['Inside_Air_Temp']
            ff[0].header['IN_AIR'] = in_temp
        else:
            ff[0].header['IN_AIR'] = -999.
            logging.warning("No temperature: %s" % rute)
        nupkey += 1
    # Parallactic Angle
    if 'TEL_PA' not in ff[0].header:
        if 'PRLLTC' in ff[0].header:
            tel_pa = ff[0].header['PRLLTC']
            ff[0].header['TEL_PA'] = tel_pa
        else:
            logging.warning("No PA: %s" % rute)
            ff[0].header['TEL_PA'] = -999.
        nupkey += 1
    # MJD Obs
    if 'MJD_OBS' not in ff[0].header:
        if 'JD' in ff[0].header:
            mjd = ff[0].header['JD'] - 2400000.5
            ff[0].header['MJD_OBS'] = mjd
        else:
            logging.warning("No JD: %s" % rute)
            ff[0].header['MJD_OBS'] = 0.
        nupkey += 1
    # Dome status
    if 'DOMEST' not in ff[0].header:
        ff[0].header['DOMEST'] = 'Open'
        nupkey += 1
    # FOCUS keywords
    if 'FOCPOS' not in ff[0].header:
        if 'SECFOCUS' in ff[0].header:
            focus = ff[0].header['SECFOCUS']
            ff[0].header['FOCPOS'] = (focus, "Position of secondary focus")
            nupkey += 1
    # Refraction keywords
    if 'HA_REFR' not in ff[0].header:
        if 'HAREFR' in ff[0].header:
            refr = ff[0].header['HAREFR']
        elif 'RA_REFR' in ff[0].header:
            refr = ff[0].header['RA_REFR']
        else:
            refr = -999.
        ff[0].header['HA_REFR'] = (refr, "Telescope HA Refraction")
        nupkey += 1
    # Dome gap
    if 'DOMEGAP' in ff[0].header:
        domegap = ff[0].header['DOMEGAP']
        ff[0].header['DOME_GAP'] = domegap
        nupkey += 1
    # OBSDATE and OBSTIME
    if 'OBSDATE' not in ff[0].header or 'OBSTIME' not in ff[0].header:
        dt, tm = date_time_from_filename(fname)
        ff[0].header['OBSDATE'] = (dt, "UT Start Date")
        ff[0].header['OBSTIME'] = (tm, "UT Start Time")
        nupkey += 2
    # AIRMASS
    if type(ff[0].header['AIRMASS']) == str:
        try:
            ff[0].header['AIRMASS'] = float(ff[0].header['AIRMASS'])
        except ValueError:
            ff[0].header['AIRMASS'] = (1.0, 'bogus airmass')
    # END values
    if 'ENDAIR' not in ff[0].header:
        if 'AIRMASS' in ff[0].header:
            ff[0].header['ENDAIR'] = ff[0].header['AIRMASS'] * 1.5
            nupkey += 1
    if 'END_PA' not in ff[0].header:
        if 'TEL_PA' in ff[0].header:
            ff[0].header['END_PA'] = ff[0].header['TEL_PA']
            nupkey += 1
    # INSTRUME
    ff[0].header['INSTRUME'] = 'SEDM-P60'
    # TELESCOP
    ff[0].header['TELESCOP'] = '60'
    nupkey += 2
    # OBJNAME
    if 'OBJNAME' in ff[0].header:
        if 'simulated' in ff[0].header['OBJNAME'] or \
                'telinit' in ff[0].header['OBJNAME']:
            ff[0].header['OBJNAME'] = obj
            nupkey += 1
    else:
        ff[0].header['OBJNAME'] = obj
        nupkey += 1
    # SER_NO
    if 'SER_NO' in ff[0].header:
        if 'Demo' in ff[0].header['SER_NO']:
            if 'rc' in fname:
                ff[0].header['SER_NO'] = '04001312'
            elif 'ifu' in fname:
                ff[0].header['SER_NO'] = '05313416'
            nupkey += 1
    # OBJRA, OBJDEC
    if 'OBJRA' not in ff[0].header:
        if 'OBRA' in ff[0].header:
            ff[0].header['OBJRA'] = ff[0].header['OBRA']
    if 'OBJDEC' not in ff[0].header:
        if 'OBDEC' in ff[0].header:
            ff[0].header['OBJDEC'] = ff[0].header['OBDEC']

    # Now verify header types
    for k, v in header_types.items():
        if v == str:
            if k not in ff[0].header:
                ff[0].header[k] = ''
                nupkey += 1
            else:
                if type(ff[0].header[k]) != v:
                    try:
                        newval = str(ff[0].header[k])
                        ff[0].header[k] = newval
                    except ValueError:
                        ff[0].header[k] = ''
                    nupkey += 1
        elif v == float:
            if k not in ff[0].header:
                ff[0].header[k] = -999.
                nupkey += 1
            else:
                if type(ff[0].header[k]) != v:
                    try:
                        newval = float(ff[0].header[k])
                        ff[0].header[k] = newval
                    except ValueError:
                        ff[0].header[k] = -999.
                    nupkey += 1
        elif v == int:
            if k not in ff[0].header:
                ff[0].header[k] = -999
                nupkey += 1
            else:
                if type(ff[0].header[k]) != v:
                    try:
                        newval = int(ff[0].header[k])
                        ff[0].header[k] = newval
                    except ValueError:
                        ff[0].header[k] = -999
                    nupkey += 1
        elif v == bool:
            if k not in ff[0].header:
                ff[0].header[k] = False
                nupkey += 1
            else:
                if type(ff[0].header[k]) != v:
                    try:
                        newval = bool(ff[0].header[k])
                        ff[0].header[k] = newval
                    except ValueError:
                        ff[0].header[k] = False
                    nupkey += 1
        else:
            logging.warning("Illegal type for %s: %s" % (k, rute))
    # Put version in header
    ff[0].header['HFIXVERS'] = (drp_ver, "HdrFix version")
    ff[0].header['HFIXDATE'] = (Time.now().fits, "HdrFix fix date")
    ff[0].header['HFIXNKEY'] = (nupkey, "HdrFix number of updated keywords")
    # Close
    ff.close()

    return True


if __name__ == "__main__":

    logging.basicConfig(
        format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
        datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)

    # setup arguments parser
    parser = argparse.ArgumentParser(
        description="""Fix FITS headers for SEDM observations.""",
        formatter_class=argparse.RawTextHelpFormatter)
    # setup arguments
    parser.add_argument('--date', type=str, default=None,
                        help='file spec for fits files')
    parser.add_argument('--rawdir', type=str, default=_rawpath,
                        help='raw file directory root (%s)' % _rawpath)
    args = parser.parse_args()

    if not args.date:
        logging.error("Must provide a YYYYMMDD date with --date")
    elif not args.rawdir:
        logging.error("Must provide a raw file directory with --rawdir")
    else:
        # Get input directory
        indir = os.path.join(args.rawdir, args.date)
        # Get a place to put bad files, if needed
        badir = os.path.join(indir, 'bad')
        # get a list of files
        flist = glob.glob(os.path.join(indir, "*.fits*"))
        # report
        logging.info("Fixing %d fits headers..." % len(flist))
        # Count
        nbad = 0
        ngood = 0
        # loop over files
        for file_name in flist:
            if sedm_fix_header(file_name):
                ngood += 1
            else:
                nbad += 1
                if not os.path.exists(badir):
                    os.mkdir(badir)
                os.system("mv %s %s" % (file_name, badir))
        if nbad > 0:
            logging.warning("%d bad FITS files in %s" % (nbad, args.date))
            logging.info("gzipping bad files")
            blist = os.listdir(badir)
            for fl in blist:
                subprocess.run(["gzip", os.path.join(badir, fl)])
        if ngood > 0:
            logging.info("%d good FITS files fixed in %s" % (ngood, args.date))
        # Update listings
        os.chdir(indir)
        os.system("rm *what.list")
        os.system("/scr2/sedmdrp/spy what ifu*.fits > what.list")
        os.system("/scr2/sedmdrp/spy what rc*.fits > rcwhat.list")
