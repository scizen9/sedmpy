#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import pil
except ImportError:
    import drpifu.pil as pil
import datetime
from astropy.io import fits
import glob
import os
from configparser import ConfigParser
import codecs

cfg_parser = ConfigParser()

try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                              '../config/sedmconfig.cfg')

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    cfg_parser.read_file(f)

_rawpath = cfg_parser.get('paths', 'rawpath')
_reduxpath = cfg_parser.get('paths', 'reduxpath')


def time_from_fspec(filespec=None, imtype="ifu"):
    """Return time in seconds based on filename"""
    fsec = -1
    tstr = ''
    if filespec is not None:
        try:
            tstr = filespec.split(imtype)[-1].split('_')[1:4]
        except IndexError:
            pass
        if len(tstr) == 3:
            fsec = int(tstr[0]) * 3600 + int(tstr[1]) * 60 + int(tstr[2])

    return fsec


def build_image_report(indir=None, fspec=None):
    """ """
    now = datetime.datetime.now()
    timestring = "made on %d-%02d-%02d at %02d:%02d:%02d" % (now.year,
                                                             now.month,
                                                             now.day,
                                                             now.hour,
                                                             now.minute,
                                                             now.second)
        
    footer = pil.get_buffer([20, 1], timestring, fontsize=15,
                            textprop=dict(color="0.5"), barcolor="k")

    try:
        specfile = glob.glob("spec_*"+fspec+"*.fits")[0]
    except IndexError:
        print("No spec_*%s*.fits file found" % fspec)
        return None, None

    print("Reading %s header" % specfile)
    try:
        header = fits.getheader(specfile)
    except OSError:
        print("This observation failed, fix and remove spec_*_failed.fits")
        return None, None

    filesourcename = specfile.split("spec_")[-1].split(".fits")[0]

    # remove the [A] in 'TARGET [A]' and remove spaces
    object_name = header['OBJECT'].split()[0].replace(" ", "-")
    
    if "STD" in filesourcename:
        is_std = True
    else:
        is_std = False

    # check for quality
    if 'QUALITY' not in header:
        header['QUALITY'] = -1

    # Spectrum ID
    spec_id = filesourcename.split("ifu"+indir+"_")[-1].split("_" +
                                                              object_name)[0]
    # Missing plot format
    prop_missing = dict(fontsize=20, textprop=dict(color="C1"))

    # Spaxels used:
    try:
        spaxel_file = glob.glob("ifu_spaxels*"+fspec+"*.png")[0]
        img_spax = pil.Image.open(spaxel_file)
    except IndexError:
        print("Cannot find ifu spaxel file")
        img_spax = pil.get_buffer([7, 7], "Spaxel IFU image missing",
                                  **prop_missing)

    # Output Spectra
    all_spectra_files = glob.glob("spec*"+fspec+"*.png")
    extention = "%s.png" % object_name
    pysedm_spec_file = glob.glob("spec*"+fspec+"*"+extention)
    if pysedm_spec_file:
        pysedm_spec_file = pysedm_spec_file[0]
        if not is_std:
            typed_spectra = [f for f in all_spectra_files
                             if not f.endswith(extention)]
            used_spec_file = pysedm_spec_file if len(typed_spectra) == 0 \
                else typed_spectra[0]
        else:
            calib_spectra = glob.glob("calibcheck_spec_"+filesourcename + "*.png")
            used_spec_file = pysedm_spec_file if len(calib_spectra) == 0 \
                else calib_spectra[0]
        try:
            img_spec = pil.Image.open(used_spec_file)
        except FileNotFoundError:
            print("Cannot find spectra plot")
            img_spec = pil.get_buffer([13, 7], "Spectra image missing",
                                      **prop_missing)
    else:
        print("No spectrum plot")
        img_spec = pil.get_buffer([13, 7], "Spectra image missing",
                                  **prop_missing)

    # Acquisition finder
    t_obs = time_from_fspec(specfile)
    try:
        if is_std:
            fspec = os.path.join(_reduxpath, "%s/finders/finder_*ACQ-%s_*.png"
                                 % (indir, object_name.split("STD-")[-1]))
        else:
            fspec = os.path.join(_reduxpath, "%s/finders/finder_*ACQ-%s_*.png"
                                 % (indir, object_name))
        finder_file = glob.glob(fspec)
        if not finder_file:
            if is_std:
                fspec = os.path.join(_reduxpath, "%s/finders/finder_*_%s_*.png"
                                     % (indir, object_name.split("STD-")[-1]))
            else:
                fspec = os.path.join(_reduxpath, "%s/finders/finder_*_%s_*.png"
                                     % (indir, object_name))
            finder_file = glob.glob(fspec)
        if finder_file:
            # Get delta time = (t_obs - t_acq)
            off_time_list = []
            for f in finder_file:
                off_time_list.append(abs(t_obs - time_from_fspec(filespec=f,
                                                                 imtype="rc")))
            off_min = min(off_time_list)
            # Use the one that is closest to the observation time
            for f in finder_file:
                t_off = abs(t_obs - time_from_fspec(filespec=f, imtype="rc"))
                if t_off == off_min:
                    finder_file = f
                    break
            # Check time offset
            if off_min > 120:
                print("Warning: time offset = %d > 120s" % off_min)

            img_find = pil.Image.open(finder_file)
        else:
            print("Cannot find %s" % fspec)
            img_find = pil.get_buffer([8, 7], "Finder image missing",
                                      **prop_missing)
    except IndexError:
        print("Cannot find %s" % fspec)
        img_find = pil.get_buffer([8, 7], "Finder image missing",
                                  **prop_missing)

    # ---------
    # PSF Extraction Specials
    # ---------
    # PSF
    try:
        psf_file = glob.glob("psfprofile_"+filesourcename + "*.png")[0]
        img_psf = pil.Image.open(psf_file).crop((50, 0, 995, 500))
    except IndexError:
        print("Cannot find psfprofile image")
        img_psf = pil.get_buffer([8, 4], "PSF Profile image missing",
                                 **prop_missing)

    # ============== #
    #  Combination   #
    # ============== #    
    # title
    title = "%s" % object_name
    title += " | exptime: %.1f s" % header["EXPTIME"]
    title += " | qual: %d" % header["QUALITY"]
    title += " | file ID: %s " % spec_id

    # Auto
    title_img = pil.get_buffer([8, 1], title, fontsize=16, hline=[0.3, 0.7],
                               barprop=dict(lw=1))

    img_upperright = pil.get_image_column([title_img,  img_find])
    img_upperleft = pil.get_image_row([img_spax])
    img_upper = pil.get_image_row([img_upperleft, img_upperright])

    img_lower = pil.get_image_row([img_psf, img_spec])

    img_combined = pil.get_image_column([img_upper, img_lower, footer])

    filesourcename = filesourcename.replace("_robotA_", "_robot_")
    filesourcename = filesourcename.replace("_robotB_", "_robot_")
            
    return img_combined, "verify_"+filesourcename+".png"
        

#################################
#
#   MAIN 
#
#################################
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=""" verify spectral extraction
            """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('infile', type=str, default=None,
                        help='cube filepath (UT date as YYYYMMDD)')

    parser.add_argument('--contains',  type=str, default="*",
                        help='Provide here part of the filename.')

    # ================ #
    # END of Option    #
    # ================ #
    args = parser.parse_args()

    print("Verifying observation %s in directory %s" % (args.contains,
                                                        args.infile))

    img_report, report_filename = build_image_report(indir=args.infile,
                                                     fspec=args.contains)
    if img_report is not None and report_filename is not None:
        img_report.save(report_filename, dpi=(1000, 1000))
