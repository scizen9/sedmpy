# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 13:11:52 2015

@author: nadiablago
"""

try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
import os
import glob
import shutil
import sys
import numpy as np

import ccdproc

from matplotlib import pylab as plt
import subprocess
import argparse
import datetime
import logging

try:
    import sextractor
except ImportError:
    import drprc.sextractor as sextractor

from astropy.io import fits

try:
    import coordinates_conversor as cc
except ImportError:
    import drprc.coordinates_conversor as cc

try:
    import cosmics
except ImportError:
    import drprc.cosmics as cosmics

from configparser import ConfigParser
import codecs

parser = ConfigParser()

configfile = os.environ['SEDMCONFIG']

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.read_file(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')
_db = parser.get('persistence', 'db')

FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
now = datetime.datetime.utcnow()
timestamp = datetime.datetime.isoformat(now)
creationdate = timestamp
timestamp = timestamp.split("T")[0]

try:
    # Log into a file
    root_dir = _logpath
    logging.basicConfig(
        format=FORMAT,
        filename=os.path.join(root_dir,
                              "rcred_{0}.log".format(timestamp)),
        level=logging.INFO)
    logger = logging.getLogger('rcred')
except OSError:
    logging.basicConfig(
        format=FORMAT,
        filename=os.path.join("/tmp", "rcred_{0}.log".format(timestamp)),
        level=logging.INFO)
    logger = logging.getLogger("rcred")


def get_xy_coords(image, ra, dec):
    """
    Uses the wcs-rd2xy routine to compute the proper pixel number where the
    target is.  Sometime the wcs module does not seem to be providing the
    correct answer, as it does not seem to be using the SIP extension.

    """
    import re
    import subprocess
    cmd = "wcs-rd2xy -w %s -r %.5f -d %.5f" % (image, ra, dec)
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output = str(proc.stdout.read())
    output = output.split("->")[1]

    coords = []
    for ss in output.split(","):
        coords.append(float(re.findall("[-+]?\d+[\.]?\d*", ss)[0]))

    return coords


def create_masterbias(biasdir=None, channel='rc'):
    """
    Combines slow and fast readout mode biases for the specified channel.
    """

    if (biasdir is None) or biasdir == "":
        biasdir = "."

    outs = "Bias_%s_slow.fits" % channel
    outf = "Bias_%s_fast.fits" % channel

    doslow = True
    dofast = True
    if os.path.isfile(os.path.join(biasdir, outs)):
        logger.warning("%s master Bias exists!" % outs)
        doslow = False
    if os.path.isfile(os.path.join(biasdir, outf)):
        logger.warning("%s master Bias exists!" % outs)
        dofast = False

    if doslow or dofast:
        logger.info("Starting the Master Bias creation!")
    else:
        return

    os.chdir(biasdir)

    lfastbias = []
    lslowbias = []

    # Select all filts that are Bias with same instrument
    for ff in glob.glob("rc*fits"):
        try:
            if "BIAS" in str.upper(fitsutils.get_par(ff, "IMGTYPE").upper()):
                if fitsutils.get_par(ff, "ADCSPEED") == 2:
                    lfastbias.append(ff)
                else:
                    lslowbias.append(ff)
        except (KeyError, OSError, AttributeError):
            pass

    logger.info("Files for bias SLOW mode: %s" % lslowbias)
    logger.info("Files for bias FAST mode: %s" % lfastbias)

    if len(lfastbias) > 0 and dofast:

        fstacked = ccdproc.combine(lfastbias, method="median",
                                   sigma_clip=True,
                                   sigma_clip_low_thresh=None,
                                   sigma_clip_high_thresh=2.0, unit='adu')
        fstacked.header['HISTORY'] = 'Master Bias stacked'
        fstacked.header['NSTACK'] = (len(lfastbias), 'number of images stacked')
        fstacked.header['STCKMETH'] = ("median", 'method used for stacking')
        for ii, fname in enumerate(lfastbias):
            fstacked.header['STACKF%d' % (ii+1)] = (fname, 'stack input file')
        hdulist = fstacked.to_hdu()
        hdulist.writeto(outf)

        # copy into the reference folder with current date
        newdir = os.path.join("../../refphot/",
                              os.path.basename(os.path.abspath(biasdir)))
        if not os.path.isdir(newdir):
            os.makedirs(newdir)
            shutil.copy(outf, os.path.join(newdir, os.path.basename(outf)))
        else:
            copy_ref_calib(biasdir, "Bias")

    if len(lslowbias) > 0 and doslow:

        sstacked = ccdproc.combine(lslowbias, method="median",
                                   sigma_clip=True,
                                   sigma_clip_low_thresh=None,
                                   sigma_clip_high_thresh=2.0, unit='adu')

        sstacked.header['HISTORY'] = 'Master Bias stacked'
        sstacked.header['NSTACK'] = (len(lslowbias), 'number of images stacked')
        sstacked.header['STCKMETH'] = ("median", 'method used for stacking')
        for ii, fname in enumerate(lslowbias):
            sstacked.header['STACKF%d' % (ii + 1)] = (fname, 'stack input file')
        hdulist = sstacked.to_hdu()
        hdulist.writeto(outs)

        # copy into the reference folder with current date
        newdir = os.path.join("../../refphot/",
                              os.path.basename(os.path.abspath(biasdir)))
        if not os.path.isdir(newdir):
            os.makedirs(newdir)
        shutil.copy(outs, os.path.join(newdir, os.path.basename(outs)))
    else:
        copy_ref_calib(biasdir, outs)


def create_masterflat(flatdir=None, biasdir=None, channel='rc', plot=True):
    """
    Creates a masterflat from both dome flats and sky flats if the number of
    counts in the given filter is not saturated and not too low
    (between 3000 and 40000).
    """
    if flatdir is None or flatdir == "":
        flatdir = "."

    if biasdir is None or biasdir == "":
        biasdir = flatdir

    os.chdir(flatdir)

    if plot and not os.path.isdir("reduced/flats"):
        os.makedirs("reduced/flats")

    if len(glob.glob("Flat_%s*norm.fits" % channel)) == 4:
        logger.info("Master Flat exists!")
        return
    if len(glob.glob("Flat_%s*norm.fits" % channel)) > 0:
        logger.info("Some Master Flat exist!")
    else:
        logger.info("Starting the Master Flat creation!")

    bias_slow = "Bias_%s_slow.fits" % channel
    bias_fast = "Bias_%s_fast.fits" % channel

    if not os.path.isfile(bias_slow) and not os.path.isfile(bias_fast):
        create_masterbias(biasdir)

    lsflat = []
    lfflat = []

    # Select all filts that are Flats with same instrument
    for ff in glob.glob(channel + "*fits"):
        try:
            if fitsutils.has_par(ff, "IMGTYPE"):
                imtype = str.upper(fitsutils.get_par(ff, "IMGTYPE"))
            else:
                continue

            # if ("RAINBOW CAM" in str.upper(fitsutils.get_par(f, "CAM_NAME"))
            # and  ("DOME" in  obj or "FLAT" in obj or "Twilight" in obj or
            # "TWILIGHT" in imtype or "DOME" in imtype)):
            if "twilight" in imtype.lower():

                if fitsutils.get_par(ff, "ADCSPEED") == 2:
                    lfflat.append(ff)
                else:
                    lsflat.append(ff)
        except (OSError, KeyError):
            logger.error("Error with retrieving parameters for file %s" % ff)
            pass

    logger.info("Files for slow flat %s" % lsflat)
    logger.info("Files for fast flat %s" % lfflat)

    # Remove bias from the flat
    debiased_flats = []
    if len(lsflat) > 0:
        sbias = ccdproc.fits_ccddata_reader(filename=bias_slow)
        for sfl in lsflat:
            rawf = ccdproc.fits_ccddata_reader(filename=sfl, unit='adu')
            rawf.data = rawf.data.astype(np.float64)
            rawf.data -= sbias.data
            rawf.header['HISTORY'] = 'Bias subtracted'
            rawf.header['MBIAS'] = (bias_slow, 'master bias file')
            redhdul = rawf.to_hdu()
            redhdul.writeto('b_'+sfl)
            debiased_flats.append('b_'+sfl)

    if len(lfflat) > 0:
        fbias = ccdproc.fits_ccddata_reader(filename=bias_fast)
        for ffl in lfflat:
            rawf = ccdproc.fits_ccddata_reader(filename=ffl, unit='adu')
            rawf.data = rawf.data.astype(np.float64)
            rawf.data -= fbias.data
            rawf.header['HISTORY'] = 'Bias subtracted'
            rawf.header['MBIAS'] = (bias_fast, 'master bias file')
            redhdul = rawf.to_hdu()
            redhdul.writeto('b_' + ffl)
            debiased_flats.append('b_'+ffl)

    # Slices the flats.
    for ff in debiased_flats:
        logger.info("Slicing file %s" % ff)
        try:
            slice_rc(ff, calib=True)
        except OSError:
            logger.error("Error when slicing file, "
                         "deleting the unsliced one...")
        # Remove the un-sliced file
        os.remove(ff)

    # Selects the ones that are suitable given
    # the number of counts and combines them.
    bands = ['u', 'g', 'r', 'i']
    for b in bands:
        out = "Flat_%s_%s.fits" % (channel, b)
        out_norm = out.replace(".fits", "_norm.fits")

        if os.path.isfile(out_norm):
            logger.error("Master Flat for filter %s exists. Skipping..." % b)
            continue

        lfiles = []
        status = "none"
        d = None
        for ff in glob.glob('b_*_%s.fits' % b):
            fi = fits.open(ff)
            d = fi[0].data
            status = "rejected"
            if 4000 < np.percentile(d, 90) < 45000:
                lfiles.append(ff)
                mymode = 1. * np.median(d.flatten())
                d[d > 45000] = mymode
                fi[0].header['FLMODE'] = (mymode, 'median of flat level')
                fi[0].data = d
                fi.writeto(ff, overwrite=True)
                status = "accepted"

            if plot:
                plt.title("Flat filter %s. %s" % (b, status))
                plt.imshow(d.T, cmap=plt.get_cmap("nipy_spectral"))
                plt.colorbar()
                plt.savefig("reduced/flats/%s" % (ff.replace(".fits",
                                                             ".png")))
                plt.close()
        # Make sure that the optimum number of counts
        # is not too low and not saturated.
        if len(lfiles) == 0:
            logger.error("WARNING!!! Could not find suitable flats for band %s"
                         % b)
            continue
        if len(lfiles) < 3:
            logger.error("WARNING!!! Could find less than 3 flats for band %s."
                         "Skipping, as it is not reliable..." % b)
            continue

        # Cleaning of old files
        if os.path.isfile(out):
            os.remove(out)
        if os.path.isfile(out_norm):
            os.remove(out_norm)
        if os.path.isfile("Flat_stats"):
            os.remove("Flat_stats")

        # read in flats
        scales = []
        stack = []
        ref_mode = 1.
        for ffl in lfiles:
            ccdd = ccdproc.fits_ccddata_reader(ffl)
            stack.append(ccdd)
            if len(stack) == 1:
                ref_mode = ccdd.header['FLMODE']
            scales.append(ccdd.header['FLMODE'] / ref_mode)

        stacked = ccdproc.combine(stack, method="median", sigma_clip=True,
                                  sigma_clip_low_thresh=2.,
                                  sigma_clip_high_thresh=2., scale=scales)
        stacked.header['HISTORY'] = 'Master Flat combined'
        stacked.header['NSTACK'] = (len(lfiles), 'number of images stacked')
        stacked.header['STCKMETH'] = ("median", 'method used for stacking')
        for ii, fname in enumerate(lfiles):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, 'stack input file')
        hdulist = stacked.to_hdu()
        hdulist.writeto(out)

        # Normalize flat
        mymode = np.nanmedian(stacked.data[100:-100, 100:-100])
        stacked.data /= mymode
        stacked.header['HISTORY'] = 'Master Flat normalized'
        stacked.header['FLSCALE'] = (mymode, 'Divided by this to normalize')
        hdulist = stacked.to_hdu()
        hdulist.writeto(out_norm)

        # Do some cleaning
        logger.info('Removing from lfiles')
        for ff in glob.glob('b_*_%s.fits' % b):
            os.remove(ff)

        # copy into the reference folder with current date
        newdir = os.path.join("../../refphot/",
                              os.path.basename(os.path.abspath(flatdir)))
        if not os.path.isdir(newdir):
            os.makedirs(newdir)
        shutil.copy(out_norm, os.path.join(newdir, os.path.basename(out_norm)))
    copy_ref_calib(flatdir, "Flat")


def get_median_bkg(img):
    """
    Computes the median background.
    """
    hdu = fits.open(img)
    # header = hdu[0].header
    bkg = np.median(hdu[0].data[hdu[0].data > 0])
    return bkg


def copy_ref_calib(curdir, calib="Flat"):
    """
    Reference master Bias and master Flat are stored in the refphot folder.
    The files are copied if they are not found in the folder where the
    photometry is being reduced.
    """

    if calib == "Bias":
        calib_dic = {"Bias_rc_fast.fits": False, "Bias_rc_slow.fits": False}
    else:
        calib_dic = {"Flat_rc_u_norm.fits": False,
                     "Flat_rc_g_norm.fits": False,
                     "Flat_rc_r_norm.fits": False,
                     "Flat_rc_i_norm.fits": False}

    # Get the date of the current directory
    curdir = os.path.abspath(curdir)

    # Check if we have all the calibrations we need.
    for cal in calib_dic:
        calib_dic[cal] = os.path.isfile(os.path.join(curdir, cal))

    # If all the calibration files are in place, nothing to do.
    if np.all(calib_dic.values()):
        return

    # Check which dates have the calibration files that we need.
    listcalib = glob.glob(os.path.join(curdir, "../../refphot/*/", calib + "*"))
    # Obtain the folder name
    listdirs = [os.path.dirname(fl) for fl in listcalib]
    # Unique name
    calibdates = np.array(list(set(listdirs)))

    # compile the dates in datetime format
    dates = []
    for ff in calibdates:
        dateobs = os.path.basename(ff)
        dates.append(datetime.datetime.strptime(dateobs, "%Y%m%d"))

    # dates = np.array(dates)
    curdate = datetime.datetime.strptime(os.path.basename(curdir), "%Y%m%d")

    # Select the folder that is closer to the date of the current directory
    # lastdir = np.array(calibdates)[np.argmin(np.abs(curdate - dates))]

    # Try to go 100 days before the data was taken.
    i = 0
    while i < 1000 and not np.all(calib_dic.values()):
        i = i + 1
        newdate = curdate - datetime.timedelta(i)
        newdatedir = "%d%02d%02d" % (newdate.year, newdate.month, newdate.day)
        print("Checking day %s" % newdatedir)
        for cal in np.array(calib_dic.keys())[~np.array(calib_dic.values())]:
            c = os.path.join(curdir, "../../refphot/", newdatedir, cal)

            if os.path.isfile(c):
                print("Copying calibration file %s to directory %s" % (c,
                                                                       curdir))
                shutil.copy(c, os.path.join(curdir, os.path.basename(c)))
                calib_dic[cal] = True


def solve_astrometry(img, outimage=None, radius=3, with_pix=True,
                     overwrite=False, tweak=3):
    """
    img: fits image where astrometry should be solved.
    outimage: name for the astrometry solved image. Defaults to "a_"img.
    radius: radius of uncertainty on astrometric position in image.
    with_pix: if we want to include the constraint on the pixel size
                for the RCCam.
    overwrite: if the astrometrically solved image should overwrite the old one.
    tweak: parameter for astrometry.net
    """

    from astropy.wcs import InconsistentAxisTypesError

    img = os.path.abspath(img)

    ra = fitsutils.get_par(img, 'RA')
    dec = fitsutils.get_par(img, 'DEC')
    # logger.info( "Solving astrometry on field with (ra,dec)=%s %s"%(ra, dec))

    astro = os.path.join(os.path.dirname(img), "a_" + os.path.basename(img))

    # If astrometry exists, we don't run it again.
    if os.path.isfile(astro) and not overwrite:
        return astro

    # Store a temporary file with the multiplication of the mask
    # mask = "mask.fits"

    cmd = "solve-field --ra %s --dec %s --radius %.4f -p --new-fits %s \
      -W none -B none -P none -M none -R none -S none -t %d --overwrite %s " %\
          (ra, dec, radius, astro, tweak, img)
    if with_pix:
        cmd = cmd + \
              " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.4"
    # logger.info(cmd)

    subprocess.call(cmd, shell=True)

    # Cleaning after astrometry.net
    if os.path.isfile(img.replace(".fits", ".axy")):
        os.remove(img.replace(".fits", ".axy"))
    if os.path.isfile(img.replace(".fits", "-indx.xyls")):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if os.path.isfile("none"):
        os.remove("none")

    if os.path.isfile(astro):
        try:
            is_on_target(img)
        except InconsistentAxisTypesError as e:
            fitsutils.update_par(img, "ONTARGET", 0)
            print("Error detected with WCS when reading file %s. \n %s" % (img,
                                                                           e))
            # If no name is provided, and overwrite is true, we move the
            # astrometry solved file into the image with name outimage.
    if outimage is not None and overwrite and os.path.isfile(astro):
        shutil.move(astro, outimage)
        return outimage
    # If ther is no name, we overwrite the old file.
    elif outimage is None and overwrite and os.path.isfile(astro):
        shutil.move(astro, img)
        return img
    # Otherwise, we just return the name of the astrometry sovled image
    else:
        return astro


def make_mask_cross(img):
    maskname = os.path.join(os.path.dirname(img), "mask.fits")

    if os.path.isfile(maskname):
        return maskname

    ff = fits.open(img)
    data = ff[0].data

    corners = {
        "g": [1, 850, 1, 850],
        "i": [1, 850, 110, 2045],
        "r": [1150, 2045, 1150, 2045],
        "u": [1130, 2045, 1, 850]
    }

    newdata = np.ones_like(data) * np.nan

    for band in corners.keys():
        ix = corners[band]
        newdata[ix[0]:ix[1], ix[2]: ix[3]] = 1

    ff[0].data = newdata
    ff.writeto(maskname, overwrite=True)

    return maskname


def slice_rc(img, calib=False):
    """
    Slices the Rainbow Camera into 4 different images and adds the 'filter'
    keyword in the fits file.
    """
    fname = os.path.basename(img)
    # fdir = os.path.dirname(img)

    corners = {
        "g": [1, 910, 1, 900],
        "i": [1, 910, 1060, 2045],
        "r": [1040, 2045, 1015, 2045],
        "u": [1030, 2045, 1, 900]
    }

    frame = ccdproc.fits_ccddata_reader(img, unit='adu')

    filenames = []

    for i, b in enumerate(corners.keys()):
        logger.info("Slicing for filter %s" % b)
        name = fname.replace(".fits", "_%s.fits" % b)

        # Clean first
        if os.path.isfile(name):
            os.remove(name)

        f_frame = frame.copy()

        f_frame.data = f_frame.data[corners[b][2]:corners[b][3],
                                    corners[b][0]:corners[b][1]]
        f_frame.header['FILTER'] = (b, 'RC filter')
        hdul = f_frame.to_hdu()
        hdul.writeto(name)

        fitsutils.update_par(name, 'filter', b)
        if not calib:
            is_on_target(name)

        filenames.append(name)

    return filenames


def is_on_target(image):
    """
    Add as a parameter whether the image is on target or not.

    """
    ra, dec = cc.hour2deg(fitsutils.get_par(image, 'OBJRA'),
                          fitsutils.get_par(image, 'OBJDEC'))

    impf = fits.open(image)
    # w = WCS(impf[0].header)

    # filt = fitsutils.get_par(image, "FILTER")
    pra, pdec = get_xy_coords(image, ra, dec)

    shape = impf[0].data.shape

    if (pra > 0) and (pra < shape[0]) and (pdec > 0) and (pdec < shape[1]):
        fitsutils.update_par(image, "ONTARGET", 1)
        return True
    else:
        fitsutils.update_par(image, "ONTARGET", 0)
        return False


def clean_cosmic(fl):
    """
    From lacosmic.
    """
    out = fl.replace('.fits', '_clean.fits')

    # If it does already exist, just return the name.
    if os.path.isfile(out):
        return out

    # Otherwise, run the cosmic ray rejection based on LA Cosmic.
    g = fitsutils.get_par(fl, "GAIN")
    if fitsutils.has_par(fl, "RDNOISE"):
        rn = fitsutils.get_par(fl, "RDNOISE")
    else:
        rn = 20
    array, header = cosmics.fromfits(fl)

    try:
        c = cosmics.cosmicsimage(array, gain=g, readnoise=rn, sigclip=8.0,
                                 sigfrac=0.3, satlevel=64000.0)
        c.run(maxiter=3)
        out = fl.replace('.fits', '_clean.fits')

        cosmics.tofits(out, c.cleanarray, header)
        fitsutils.update_par(out, "CRREJ", 1)

        # os.remove(fl)
    except:
        pass

    return out


def get_overscan_bias_rc(img):
    """
    Bias from overscan region.
    """
    ff = fits.open(img)
    bias = np.nanmedian(ff[0].data[990 - 100:990 + 100,
                        970 - 100:970 + 100].flatten())

    return bias


def get_sequential_name(target_dir, name, i=0):
    """
    Gets a sequential name if we have a file with the same object name
    imaged several times.
    """
    newname = os.path.join(target_dir, name).replace(".fits", "_%d.fits" % i)
    # If the destination file does not exist...
    if os.path.isfile(name) and not os.path.isfile(newname):
        return newname
    # If it does exist, but it is another exposure
    elif os.path.isfile(name) and os.path.isfile(newname):
        # If it is not the same exposure, add with different name.
        # Otherwise, replace.
        if fitsutils.get_par(name, "JD") != fitsutils.get_par(newname, "JD"):
            newname = get_sequential_name(target_dir, name, i=i + 1)

    return newname


def init_header_reduced(image):
    """
    IQWCS = 1 or 0 / Indicates astrometry has been solved for the field
    IQZEROPT =  1 or 0 /indicates if the zero point was calculated for the image
    SKYBKG = FLOAT / Average sky background given in counts
    SEEPIX = FLOAT  / Seeing expressed in pixels
    ZPCAT = 'String' / Catalog used to calculate zero point
    ZEROPTU = FLOAT  / Zero point uncertainty
    ZEROPT  = 'FLOAT'  / Zero point by comparison to catalog
    """

    pardic = {"IQWCS": 0,
              "IQZEROPT": 0,
              "SKYBKG": 0,
              "SEEPIX": 0,
              "ZPCAT": "none",
              "ZEROPTU": 0.,
              "ZEROPT": 0.,
              "CRREJ": 0}
    fitsutils.update_pars(image, pardic)


def plot_image(image):
    """
    Plots the reduced image into the png folder.

    """
    logger.info("Plotting image %s" % image)

    print("Plotting image ", image)

    image = os.path.abspath(image)

    # Change to image directory
    imdir, imname = os.path.split(image)

    # Create destination directory

    png_dir = os.path.join(imdir, "png")

    if not os.path.isdir(png_dir):
        os.makedirs(png_dir)

    ff = fits.open(image)[0]
    d = ff.data
    h = ff.header
    exptime = h.get('EXPTIME', 0)
    name = h.get('OBJECT', 'None')
    filt = h.get('FILTER', 'NA')
    ontarget = h.get('ONTARGET', False)

    plt.imshow(d, origin="lower", vmin=np.percentile(d.flatten(), 5),
               vmax=np.percentile(d, 95), cmap=plt.get_cmap('cubehelix'))
    plt.title("%s %s-band [%ds]. On target=%s" % (name, filt,
                                                  exptime, ontarget))
    plt.savefig(os.path.join(png_dir, imname.replace(".fits", ".png")))
    plt.close()


def reduce_image(image, flatdir=None, biasdir=None, cosmic=False,
                 astrometry=True, channel='rc', target_dir='reduced',
                 overwrite=False):
    """
    Applies Flat field and bias calibrations to the image.

    Steps:

    1. - Solve astrometry on the entire image.
    2. - Computes cosmic ray rejectionon the entire image.
    3. - Compute master bias (if it does not exist) and de-bias the image.
    4. - Separate the image into 4 filters.
    5. - Compute flat field for each filter (if it does not exist) and apply
            flat fielding on the image.
    6. - Compute the image zeropoint.

    """

    logger.info("Reducing image %s" % image)

    print("Reducing image ", image)

    image = os.path.abspath(image)
    imname = os.path.basename(image).replace(".fits", "")
    try:
        objectname = fitsutils.get_par(image, "NAME").replace(" ", "") + "_" + \
                     fitsutils.get_par(image, "FILTER")
    except KeyError:
        logger.error("ERROR, image " + image +
                     " does not have a NAME or a FILTER!!!")
        return

    print("For object", objectname)
    logger.info("For object %s" % objectname)

    # Change to image directory
    curdir = os.path.dirname(image)
    if curdir == "":
        curdir = "."
    curdir = os.path.abspath(curdir)
    os.chdir(curdir)
    # Create destination directory
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)

    # If we don't want to overwrite the already extracted images,
    # we check wether they exist.
    if not overwrite:
        existing = True
        for band in ['u', 'g', 'r', 'i']:
            destfile = os.path.join(target_dir, imname + "_f_b_a_%s_%s_0.fits" %
                                    (objectname, band))
            logger.info("Looking if file %s exists: %s" %
                        (destfile, (os.path.isfile(destfile))))
            existing = existing and (os.path.isfile(destfile))
        if existing:
            return []

    # Initialize the basic parameters.
    init_header_reduced(image)

    astro = ""
    if astrometry:
        logger.info("Solving astometry for the whole image...")
        img = solve_astrometry(image)
        if os.path.isfile(img):
            astro = "a_"
            fitsutils.update_par(img, "IQWCS", 1)
        else:
            logger.error("ASTROMETRY DID NOT SOLVE ON IMAGE %s" % image)
            img = image
    else:
        img = image

    # Update noise parameters needed for cosmic reection
    if fitsutils.get_par(img, "ADCSPEED") == 2:
        fitsutils.update_par(img, "RDNOISE", 20.)
    else:
        fitsutils.update_par(img, "RDNOISE", 4.)

    if cosmic:
        logger.info("Correcting for cosmic rays...")
        # Correct for cosmics each filter
        cleanimg = clean_cosmic(os.path.join(os.path.abspath(curdir), img))
        img = cleanimg

    # Compute BIAS
    if biasdir is None or biasdir == "":
        biasdir = "."
    create_masterbias(biasdir)

    bias_slow = os.path.join(biasdir, "Bias_%s_%s.fits" % (channel, 'slow'))
    bias_fast = os.path.join(biasdir, "Bias_%s_%s.fits" % (channel, 'fast'))

    # Compute flat field
    if flatdir is None or flatdir == "":
        flatdir = "."
    create_masterflat(flatdir, biasdir)

    # New names for the object.
    debiased = os.path.join(os.path.dirname(img), "b_" + os.path.basename(img))
    logger.info("Creating debiased file, %s" % debiased)

    if ((fitsutils.get_par(img, "ADCSPEED") == 0.1 and
         not os.path.isfile(bias_slow)) or
            (fitsutils.get_par(img, "ADCSPEED") == 2 and
             not os.path.isfile(bias_fast))):
        logger.warning(
            "Master bias not found! Tryting to copy from reference folder...")
        copy_ref_calib(curdir, "Bias")
        if ((fitsutils.get_par(img, "ADCSPEED") == 0.1 and
             not os.path.isfile(bias_slow)) or
                (fitsutils.get_par(img, "ADCSPEED") == 2 and
                 not os.path.isfile(bias_fast))):
            logger.error("Bias not found in reference folder")
            return

    # Clean first
    if os.path.isfile(debiased):
        os.remove(debiased)

    # Debias
    rawf = ccdproc.fits_ccddata_reader(filename=img, unit='adu')
    rawf.data = rawf.data.astype(np.float64)
    if rawf.header["ADCSPEED"] == 2:
        fbias = ccdproc.fits_ccddata_reader(filename=bias_fast)
        rawf.data -= fbias.data
        rawf.header['BIASFILE'] = (bias_fast, 'Master bias')
        rawf.header['RDNOISE'] = (20., 'Read noise in electrons')
    else:
        sbias = ccdproc.fits_ccddata_reader(filename=bias_slow)
        rawf.data -= sbias.data
        rawf.header['BIASFILE'] = (bias_slow, 'Master bias')
        rawf.header['RDNOISE'] = (20., 'Read noise in electrons')
    # Set negative counts to zero
    rawf.data[rawf.data < 0] = 0
    hdul = rawf.to_hdu()
    hdul.writeto(debiased)

    # Slicing the image for flats
    slice_names = slice_rc(debiased)
    print("Creating sliced files, ", slice_names)

    # Remove un-sliced image
    os.remove(debiased)

    # DE-flat each filter and store under object name
    for i, debiased_f in enumerate(slice_names):
        b = fitsutils.get_par(debiased_f, 'filter')

        deflatted = os.path.join(
            os.path.dirname(image), target_dir,
            imname + "_f_b_" + astro + objectname + "_%s.fits" % b)

        # Flat to be used for that filter
        flat = os.path.join(flatdir, "Flat_%s_%s_norm.fits" % (channel, b))

        if not os.path.isfile(flat):
            logger.warning("Master flat not found in %s" % flat)
            copy_ref_calib(curdir, "Flat")
            continue
        else:
            logger.info("Using flat %s" % flat)

        # Cleans the deflatted file if exists
        if os.path.isfile(deflatted):
            os.remove(deflatted)

        if os.path.isfile(debiased_f) and os.path.isfile(flat):
            logger.info("Storing de-flatted %s as %s" % (debiased_f, deflatted))
            flatf = ccdproc.fits_ccddata_reader(flat)
            debif = ccdproc.fits_ccddata_reader(debiased_f)
            debif.data /= flatf.data
            debif.header['HISTORY'] = 'Flat-field corrected'
            debif.header['FLATFILE'] = (flat, 'Master flat file')
            debif.header['ORIGFILE'] = (os.path.basename(image),
                                        'Original filename')
            hdul = debif.to_hdu()
            hdul.writeto(deflatted)
        else:
            logger.error("SOMETHING IS WRONG. Error when dividing %s by "
                         "the flat field %s!" % (debiased_f, flat))

        # Removes the de-biased file
        os.remove(debiased_f)

        slice_names[i] = deflatted

    # Moving files to the target directory
    for image in slice_names:
        bkg = get_median_bkg(image)
        fitsutils.update_par(image, "SKYBKG", bkg)

        # Get basic statistics for the image
        nsrc, fwhm, ellip, bkg = sextractor.get_image_pars(image)
        plot_image(image)

        logger.info("Sextractor statistics: nscr %d, fwhm (arcsec) "
                    "%.2f, ellipticity %.2f" % (nsrc, fwhm, ellip))
        print("Sextractor statistics: nscr %d, fwhm (arcsec) %.2f, "
              "ellipticity %.2f" % (nsrc, fwhm, ellip))

        dic = {"FWHM": np.round(fwhm, 3),
               "FWHMPIX": np.round(fwhm / 0.394, 3),
               "NSRC": nsrc,
               "ELLIP": np.round(ellip, 3)}
        # Update the seeing information from sextractor
        fitsutils.update_pars(image, dic)

    return slice_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""

        Reduces the photometric images from SEDM Rainbow Camera.
        Requires either a list of images to be reduced, or a directory name
        where the night photometry is. As an option, can run lacosmic to remove
        cosmic rays. By default it invokes astrometry.net before the reduction.
        
        %run rcred.py -d PHOTDIR 
                
        Reduced images are stored in directory called "reduced", within the 
        main directory.
        
        Optionally, it can be used to clean the reduction products generated 
        by this pipeline within PHOTDIR diretory using -c option (clean).
        
        %run rcred.py -d PHOTDIR -c
            
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-l', '--filelist', type=str,
                        help='File containing the list of fits for the night.',
                        default=None)
    parser.add_argument('-d', '--photdir', type=str,
                        help='Directory containing the science fits '
                             'for the night.',
                        default=None)
    parser.add_argument('-c', '--clean', action="store_true",
                        help='Clean the reduced images?', default=False)
    parser.add_argument('-o', '--overwrite', action="store_true",
                        help='re-reduce and overwrite the reduced images?',
                        default=False)
    parser.add_argument('-p', '--nocopy', action="store_true",
                        help='disable the copy the reduced folder to transient',
                        default=False)
    parser.add_argument('--cosmic', action="store_true", default=False,
                        help='Whether cosmic rays should be removed.')

    args = parser.parse_args()

    myfiles = []

    if args.photdir is not None and args.clean:
        for f in glob.glob("Flat*"):
            os.remove(f)
        for f in glob.glob("Bias*"):
            os.remove(f)
        for f in glob.glob("a_*fits"):
            os.remove(f)
        if os.path.isdir(os.path.join(args.photdir, "reduced")):
            shutil.rmtree(os.path.join(args.photdir, "reduced"))

    if args.filelist is not None:
        mydir = os.path.dirname(args.filelist)
        if mydir == "":
            mydir = "."
        os.chdir(mydir)

        photdir = mydir

        myfiles = np.genfromtxt(args.filelist, dtype=None)
        myfiles = [os.path.abspath(f) for f in myfiles]

    else:

        if args.photdir is None:
            timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
            timestamp = timestamp.split("T")[0].replace("-", "")
            photdir = os.path.join(_photpath, timestamp)

            s = ("WARNING! You did not specify the directory or the list:",
                 "- A filelist name with the images you want to reduce [-l] OR",
                 "- The name of the directory which you want to reduce [-d].",
                 "",
                 "A default name for the directory will be assumed on today's "
                 "date: %s" % photdir)

    mydir = os.path.abspath(args.photdir)
    # Gather all RC fits files in the folder with the keyword IMGTYPE=SCIENCE
    for f in glob.glob(os.path.join(mydir, "rc*fits")):
        try:
            if (fitsutils.has_par(f, "IMGTYPE") and
                    ((fitsutils.get_par(f, "IMGTYPE").upper() == "SCIENCE") or (
                        "ACQ" in fitsutils.get_par(f, "IMGTYPE").upper()))):
                myfiles.append(f)
        except (OSError, KeyError):
            print("problems opening file %s" % f)

    create_masterbias(mydir)
    print("Create masterflat", mydir)
    create_masterflat(mydir)

    if len(myfiles) == 0:
        print("Found no files to process")
        sys.exit()
    else:
        print("Found %d files to process" % len(myfiles))

    # Reduce them
    reducedfiles = []
    for f in myfiles:
        print(f)
        # make_mask_cross(f)
        if (fitsutils.has_par(f, "IMGTYPE") and
                (fitsutils.get_par(f, "IMGTYPE").upper() == "SCIENCE" or (
                    "ACQUI" in fitsutils.get_par(f, "IMGTYPE").upper()))):
            try:
                reduced = reduce_image(f, cosmic=args.cosmic,
                                       overwrite=args.overwrite)
                reducedfiles.extend(reduced)
            except OSError:
                print("Error when reducing image %s" % f)
                pass

    # If copy is requested, then we copy the whole folder or just the
    # missing files to transient.

    dayname = os.path.basename(os.path.dirname(os.path.abspath(myfiles[0])))
    reducedname = os.path.join(os.path.dirname(os.path.abspath(myfiles[0])),
                               "reduced")
    if args.photdir is not None and not args.nocopy:
        com = "rcp -r %s grbuser@transient.caltech.edu:" \
              "/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s" % \
              (reducedname, dayname)
        subprocess.call(com, shell=True)
    elif args.filelist is not None and not args.nocopy:
        for f in reducedfiles:
            com = "rcp %s grbuser@transient.caltech.edu:" \
                  "/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s/."\
                  % (f, dayname)
            subprocess.call(com, shell=True)
