# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 20:19:04 2016
Updated on Thu Feb  6 2020 by neill

@author: nadiablago
"""

import os
import shutil
import subprocess
import numpy as np
from astropy.io import fits as pf
from astropy.io import ascii

import datetime
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils

from matplotlib import pylab as plt

import sedmpy_version
import json

# Get pipeline configuration
# Find config file: default is sedmpy/config/sedmconfig.json
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(sedmpy_version.CONFIG_DIR, 'sedmconfig.json')
with open(configfile) as config_file:
    sedm_cfg = json.load(config_file)

_focuspath = sedm_cfg['paths']['focuspath']

rootdir = _focuspath

if not (os.path.isdir(rootdir)):
    rootdir = "/tmp"


def run_sex(flist, mask=False, cosmics=False, overwrite=False):
    
    d = os.path.dirname(flist[0])
    if d == "":
        d = "."
    os.chdir(d)

    # Create the directory where the sextracted images are going to go.
    sexdir = os.path.join(d, "sextractor")    
    if not os.path.isdir(sexdir):
        os.makedirs(sexdir)
        
    newlist = []
    for ff in flist:
        newimage = os.path.join(sexdir, os.path.basename(ff).replace(".fits",
                                                                     ".sex"))

        if os.path.isfile(newimage) and not overwrite:
            newlist.append(newimage)
            print("Sextracted image %s already exists." % newimage)
        else:
            try:
                ff = os.path.abspath(ff)
                if mask:
                    print("Mask no longer implemented.")
                    out = ff
                else:
                    out = ff
                    
                if cosmics and (not fitsutils.has_par(out, "CRREJ") or
                                fitsutils.get_par(out, "CRREJ") == 0):
                    print("sextractor cosmic ray cleaning not implemented")
                    out = ff
    
                cmd = "sex -c %s/config/daofind.sex %s" % (os.environ["SEDMPY"],
                                                           out)
                subprocess.call(cmd, shell=True)
                print(cmd)
                shutil.move("image.sex", newimage)
                newlist.append(newimage)
            except IOError:
                print("IOError detected reading file", ff)
                pass
        
    return newlist


def analyse_sex(sexfileslist, plot=True, interactive=False):
    """
    Analyses the sextractor filelist to determine
     the best focus for the RC camera.

    #   1 X_IMAGE         Object position along x                 [pixel]
    #   2 Y_IMAGE         Object position along y                 [pixel]
    #   3 ALPHA_J2000     Right ascension of barycenter (J2000)   [deg]
    #   4 DELTA_J2000     Declination of barycenter (J2000)       [deg]
    #   5 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR         [mag]
    #   6 MAGERR_BEST     RMS error for MAG_BEST                  [mag]
    #   7 FWHM_WORLD      FWHM assuming a gaussian core           [deg]
    #   8 FWHM_IMAGE      FWHM assuming a gaussian core           [pixel]
    #   9 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
    #  10 BACKGROUND      Background at centroid position         [count]
    #  11 FLAGS           Extraction flags
    #  12 A_IMAGE         Isophotal image mejor axis
    #  13 B_IMAGE         Isophotal image minor axis
    #  14 THETA_IMAGE     Isophotal image position angle
    #  15 PETRO_RADIUS    Petrosian radius

    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    """
    sexfileslist = list(sexfileslist)    
    sexfileslist.sort()
 
    focpos = []
    fwhms = []
    std_fwhm = []
    for i, ff in enumerate(sexfileslist):
        fits = ff.replace("sextractor/", "").replace(".sex", ".fits")
        hdul = pf.open(fits)
        pos = float(hdul[0].header['focpos'])

        s = np.genfromtxt(ff, comments="#")
        print(ff, "Initial number of sources", len(s))

        s = s[s[:, 1] < 2000]
        s = s[s[:, 10] != 4]

        # Only objects with FWHM less than 40 pixels... but larger than 2
        s = s[s[:, 7] < 60]
        s = s[s[:, 7] > 1]
        
        # Select bright magnitudes
        s = s[s[:, 4] < np.percentile(s[:, 4], 30)]
        # Select round sources (ellipticity is 1-axis_ratio)
        s = s[s[:, 8] < np.percentile(s[:, 8], 30)]
        print(ff, "number of sources", len(s))
 
        focpos.append(pos)
        fwhms.append(np.nanmean(s[:, 7]*0.394))
        mad = np.median(np.abs(s[:, 7] -
                               np.nanmean(s[:,
                                          7]))) / 0.67448975019608171 * 0.394
        std_fwhm.append(mad)
    
    focpos = np.array(focpos)
    fwhms = np.array(fwhms)
    std_fwhm = np.array(std_fwhm)

    n = len(fwhms)
    
    best_seeing_id = int(np.nanargmin(fwhms))
    # We will take 4 datapoints on the left and right of the best value.
    selected_ids = np.arange(-4, 5, 1)
    selected_ids = selected_ids + best_seeing_id
    selected_ids = np.minimum(selected_ids, n-1)
    selected_ids = np.maximum(selected_ids, 0)
    print("FWHMS: %s, focpos: %s, Best seeing id: %d. Selected ids %s"
          % (fwhms, focpos, best_seeing_id, selected_ids))
    selected_ids = np.array(list(set(selected_ids)))

    focpos = focpos[selected_ids]
    fwhms = fwhms[selected_ids]
    std_fwhm = std_fwhm[selected_ids]
    
    std_fwhm = np.maximum(1e-5, np.array(std_fwhm))
    
    coefs = np.polyfit(focpos, fwhms, w=1/std_fwhm, deg=2)
    
    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    print("Best focus:%.2f" % x[np.argmin(p(x))], coefs, std_fwhm)

    if plot:
        plt.title("Best focus:%.2f" % x[np.argmin(p(x))])
        with open(os.path.join(rootdir, "focus"), "w") as ff:
            ff.write(str(focpos))
            ff.write(str(fwhms))
        plt.errorbar(focpos, fwhms, yerr=std_fwhm, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("FWHM (arcsec)")
        if interactive:
            plt.show()
        else:
            plt.savefig(
                os.path.join(os.path.dirname(sexfileslist[0]),
                             "focus_%s.png" %
                             (datetime.datetime.utcnow()).strftime(
                                 "%Y%m%d-%H:%M:%S")))
            plt.clf()
    return x[np.argmin(p(x))], coefs[0]


def analyse_sex_ifu_spectrograph(sexfileslist, plot=True, interactive=False,
                                 debug=False):
    """
    Analyses the sextractor filelist to determine the best focus
     for the IFU spectrograph instrument.

    #   1 X_IMAGE         Object position along x                 [pixel]
    #   2 Y_IMAGE         Object position along y                 [pixel]
    #   3 ALPHA_J2000     Right ascension of barycenter (J2000)   [deg]
    #   4 DELTA_J2000     Declination of barycenter (J2000)       [deg]
    #   5 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR         [mag]
    #   6 MAGERR_BEST     RMS error for MAG_BEST                  [mag]
    #   7 FWHM_WORLD      FWHM assuming a gaussian core           [deg]
    #   8 FWHM_IMAGE      FWHM assuming a gaussian core           [pixel]
    #   9 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
    #  10 BACKGROUND      Background at centroid position         [count]
    #  11 FLAGS           Extraction flags
    #  12 A_IMAGE         Isophotal image mejor axis
    #  13 B_IMAGE         Isophotal image minor axis
    #  14 THETA_IMAGE     Isophotal image position angle
    #  15 PETRO_RADIUS    Petrosian radius


    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    """
    sexfileslist = list(sexfileslist)    
    sexfileslist.sort()
 
    focpos = []
    fwhms = []
    std_fwhm = []
    mags = []
    mags_std = []
    trace_width = []
    trace_width_std = []
    
    for i, ff in enumerate(sexfileslist):
        print(ff)
        fits = ff.replace("sextractor/", "").replace(".sex", ".fits")
        hdul = pf.open(fits)
        pos = float(hdul[0].header['ifufocus'])

        s = np.genfromtxt(ff, comments="#",
                          dtype=[("x", np.float),
                                 ("y", np.float),
                                 ("ra", np.float),
                                 ("dec", np.float),
                                 ("mag", np.float),
                                 ("magerr", np.float),
                                 ("fwhm_world", np.float),
                                 ("fwhm_image", np.float),
                                 ("ellipticity", np.float),
                                 ("background", np.float),
                                 ("flags", np.float),
                                 ("a_image", np.float),
                                 ("b_image", np.float),
                                 ("theta_image", np.float),
                                 ("petro_radius", np.float)])

        sb = s[(s["x"] > 200) * (s["x"] < 1848) *
               (s["y"] > 200) * (s["y"] < 1848)]
        
        sb = sb[np.abs(sb["theta_image"]) < 10]
        sb = sb[np.abs(sb["a_image"]) > 20]

        if debug:
            plt.figure(figsize=(10, 7))
            plt.suptitle('Focus %.2f' % pos, fontsize=10,
                         horizontalalignment="right")
            ax = plt.subplot2grid((2, 3), (0, 0))
            magmap = ax.scatter(sb["x"], sb["y"], c=10**(-0.4*sb["mag"]), s=5,
                                edgecolors='face')
            plt.title("flux / spaxel")
            plt.colorbar(magmap, label="flux", format="%.2e")
            
            ax = plt.subplot2grid((2, 3), (0, 1))
            ax.hist(sb["a_image"], bins=20, range=(20, 80))
            plt.title("A image")

            ax = plt.subplot2grid((2, 3), (0, 2))
            ax.hist(sb["ellipticity"], bins=50, range=(0.9, 1))
            plt.title("Ellipticity")

            ax = plt.subplot2grid((2, 3), (1, 0))
            ax.hist(sb["b_image"], bins=20, range=(0.8, 2))
            plt.title("B image")
            
            ax = plt.subplot2grid((2, 3), (1, 1))
            ax.hist(sb["theta_image"], bins=20, range=(-1, 1))
            plt.title("Theta (A/B)")
            
            ax = plt.subplot2grid((2, 3), (1, 2))
            ax.hist(sb["mag"], bins=20)
            plt.title("mag")
            
            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     os.path.basename(ff).replace(".sex",
                                                                  ".png")))
            plt.clf()
        
        focpos.append(pos)
        fwhms.append(np.median(sb["a_image"]))
        mags.append(np.average(sb["mag"]))
        mags_std.append(np.std(sb["mag"]))
        trace_width.append(np.median(sb["b_image"]))
        trace_width_std.append(np.minimum(0.3, np.std(sb["b_image"])))

        std_fwhm.append(np.std(sb["a_image"]))
    
    focpos = np.array(focpos)
    fwhms = np.array(fwhms)
    std_fwhm = np.minimum(6.5, np.array(std_fwhm))

    mags_std = np.array(mags_std)

    coefs = np.polyfit(focpos, fwhms, w=1/std_fwhm, deg=2)
    coefs_mag = np.polyfit(focpos, mags, w=1/mags_std, deg=2)

    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    pmag = np.poly1d(coefs_mag)
    
    print("Best focus:%.2f" % x[np.argmax(p(x))], coefs, std_fwhm)

    if plot:
        plt.title("Best focus:%.2f" % x[np.argmax(p(x))])
        with open(os.path.join(rootdir, "focus"), "w") as ff:
            ff.write(str(focpos))
            ff.write(str(fwhms))
        plt.errorbar(focpos, fwhms, yerr=std_fwhm, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("Median trace longitude [pixels]")
        if interactive:
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_ifu_tracelength_%s.png" %
                                     (datetime.datetime.utcnow()).strftime(
                                         "%Y%m%d-%H:%M:%S")))
            plt.clf()

    if plot:
        plt.title("Best focus:%.2f" % x[np.argmin(pmag(x))])
        plt.errorbar(focpos, mags, yerr=mags_std, fmt="o")
        plt.plot(x, pmag(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("Mag")
        if interactive:
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_ifu_mag_%s.png" %
                                     (datetime.datetime.utcnow()).strftime(
                                         "%Y%m%d-%H:%M:%S")))
            plt.clf()
            
    return x[np.argmin(pmag(x))], np.abs(coefs_mag[0])


def analyse_sex_ifu(sexfileslist, plot=True, interactive=False, debug=False):
    """
    Analyses the sextractor filelist to determine the best focus
    for the IFU camera stage controller.

    #   1 X_IMAGE         Object position along x                 [pixel]
    #   2 Y_IMAGE         Object position along y                 [pixel]
    #   3 ALPHA_J2000     Right ascension of barycenter (J2000)   [deg]
    #   4 DELTA_J2000     Declination of barycenter (J2000)       [deg]
    #   5 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR         [mag]
    #   6 MAGERR_BEST     RMS error for MAG_BEST                  [mag]
    #   7 FWHM_WORLD      FWHM assuming a gaussian core           [deg]
    #   8 FWHM_IMAGE      FWHM assuming a gaussian core           [pixel]
    #   9 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
    #  10 BACKGROUND      Background at centroid position         [count]
    #  11 FLAGS           Extraction flags
    #  12 A_IMAGE         Isophotal image mejor axis
    #  13 B_IMAGE         Isophotal image minor axis
    #  14 THETA_IMAGE     Isophotal image position angle
    #  15 PETRO_RADIUS    Petrosian radius

    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    """
    sexfileslist = list(sexfileslist)    
    sexfileslist.sort()
 
    focpos = []
    mags = []
    std_mags = []

    for i, ff in enumerate(sexfileslist):
        print(ff)
        fits = ff.replace("sextractor/", "").replace(".sex", ".fits")
        hdul = pf.open(fits)
        pos = float(hdul[0].header['focpos'])

        s = np.genfromtxt(ff, comments="#",
                          dtype=[("x", np.float),
                                 ("y", np.float),
                                 ("ra", np.float),
                                 ("dec", np.float),
                                 ("mag", np.float),
                                 ("magerr", np.float),
                                 ("fwhm_world", np.float),
                                 ("fwhm_image", np.float),
                                 ("ellipticity", np.float),
                                 ("background", np.float),
                                 ("flags", np.float),
                                 ("a_image", np.float),
                                 ("b_image", np.float),
                                 ("theta_image", np.float),
                                 ("petro_radius", np.float)])
        
        # for the focus purpose, we only want traces,
        # meaning that their ellipticity needs to be large.
        sb = s[(s["x"] > 200) * (s["x"] < 1848) *
               (s["y"] > 200) * (s["y"] < 1848)]
        # sb = sb[sb["a_image"]> 30]

        if debug:
            plt.figure(figsize=(10, 7))
            plt.suptitle('Focus %.2f' % pos, fontsize=10,
                         horizontalalignment="right")
            ax = plt.subplot2grid((2, 2), (0, 0))
            magmap = ax.scatter(sb["x"], sb["y"], c=10**(-0.4*sb["mag"]), s=10,
                                edgecolors='face')
            plt.title("flux / spaxel")
            plt.colorbar(magmap, label="flux", format="%.1e")

            ax = plt.subplot2grid((2, 2), (0, 1))
            ax.hist(sb["a_image"], bins=30, range=(0, 60))
            plt.title("A image")

            ax = plt.subplot2grid((2, 2), (1, 0))
            ax.hist(sb["b_image"], bins=20, range=(0, 5))
            plt.title("B image")
            
            ax = plt.subplot2grid((2, 2), (1, 1))
            ax.hist(sb["mag"], bins=25, range=(-15, -4))
            plt.title("mags image")
            
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     os.path.basename(ff).replace(".sex",
                                                                  ".png")))
            plt.clf()
        
        focpos.append(pos)
        avex = np.average(sb["x"])
        avey = np.average(sb["y"])
        mags.append(np.std(np.sqrt((sb["x"]-avex)**2 + (sb["y"]-avey)**2)))
        # mags.append(np.percentile(sb["mag"], 25))
        std_mags.append(np.std(sb["mag"] < np.percentile(sb["mag"], 25)))
    
    focpos = np.array(focpos)
    mags = np.array(mags)
    std_mags = np.maximum(1e-5, np.array(std_mags))

    coefs = np.polyfit(focpos, mags, w=1/std_mags, deg=2)
    
    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    print("Best focus:%.2f" % x[np.argmax(p(x))], coefs, std_mags)

    if plot:
        plt.figure()
        plt.title("Best focus IFU:%.2f" % x[np.argmax(p(x))])
        with open(os.path.join(rootdir, "focus"), "w") as ff:
            ff.write(str(focpos))
            ff.write(str(mags))
        plt.errorbar(focpos, mags, yerr=std_mags, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("Mag brightest spaxel")
        if interactive:
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_ifu_%s.png" %
                                     (datetime.datetime.utcnow()).strftime(
                                         "%Y%m%d-%H:%M:%S")))
            plt.clf()
        plt.close("all")

    return x[np.argmax(p(x))], coefs[0]


def analyze_img(sexfile, arcsecpix=0.394, is_rcam=True, mag_quantile=0.8,
                ellp_quantile=0.25, radius=10, y_min=50, y_max=2000):

    if not os.path.exists(sexfile):
        print("File not found: %s" % sexfile)
        return 0, 0, 0, 0

    data = ascii.read(sexfile)
    df = data.to_pandas()

    mag = df['MAG_BEST'].quantile(mag_quantile)
    ellip = df['ELLIPTICITY'].quantile(ellp_quantile)
    df = df[(df['MAG_BEST'] < mag) & (df['ELLIPTICITY'] < ellip)]
    df = df[(df['Y_IMAGE'] < y_max) & (df['Y_IMAGE'] > y_min)]
    df = df[(df['FLAGS'] <= 0)]
    # cut on size
    size_cut = df['FWHM_IMAGE'].median() + 2.5 * df['FWHM_IMAGE'].std()
    df = df[(df['FWHM_IMAGE'] < size_cut)]
    nextract = len(df)
    avgfwhm = df['FWHM_IMAGE'].median() * arcsecpix
    avgellip = df['ELLIPTICITY'].median()
    avgbkg = df['BACKGROUND'].median()
    return nextract, avgfwhm, avgellip, avgbkg


def analyse_image(sexfile, arcsecpix=0.394, is_rccam=True):
    """
    Analyses the sextractor filelist to determine the best focus.
    If FWHM in pixes is required, set arcsecpix=1

    #   1 X_IMAGE         Object position along x                 [pixel]
    #   2 Y_IMAGE         Object position along y                 [pixel]
    #   3 ALPHA_J2000     Right ascension of barycenter (J2000)   [deg]
    #   4 DELTA_J2000     Declination of barycenter (J2000)       [deg]
    #   5 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR         [mag]
    #   6 MAGERR_BEST     RMS error for MAG_BEST                  [mag]
    #   7 FWHM_WORLD      FWHM assuming a gaussian core           [deg]
    #   8 FWHM_IMAGE      FWHM assuming a gaussian core           [pixel]
    #   9 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
    #  10 BACKGROUND      Background at centroid position         [count]
    #  11 FLAGS           Extraction flags
    #  12 A_IMAGE         Isophotal image mejor axis
    #  13 B_IMAGE         Isophotal image minor axis
    #  14 THETA_IMAGE     Isophotal image position angle
    #  15 PETRO_RADIUS    Petrosian radius

    returns: A tuple containing:
        1. - Number of extracted sources.
        2. - FWHM in arcsecs.
        3. - Ellipticity.
        4. - Background


    """

    s = np.genfromtxt(sexfile, comments="#", dtype=[("x", np.float),
                                                    ("y", np.float),
                                                    ("ra", np.float),
                                                    ("dec", np.float),
                                                    ("mag", np.float),
                                                    ("magerr", np.float),
                                                    ("fwhm_world", np.float),
                                                    ("fwhm_image", np.float),
                                                    ("ellipticity", np.float),
                                                    ("background", np.float),
                                                    ("flags", np.float),
                                                    ("a_image", np.float),
                                                    ("b_image", np.float),
                                                    ("theta_image", np.float),
                                                    ("petro_radius", np.float)])

    if s is None or s.ndim == 0 or len(s) == 0:
        print("Empty content of the file for file %s. " % sexfile)
        return 0, 0, 0, 0

    x = s["x"]
    y = s["y"]
    
    if is_rccam:
        # Avoid cross region
        s = s[((y < 850) | (y > 1125))*((x < 885) | (x > 1540))]

    # Select with good flags only.
    s = s[s["flags"] == 0]

    nsources = len(s) 
    if nsources == 0:
        return 0, 0, 0, 0

    # Select round sources (ellipticity is 1-axis_ratio)
    s = s[s["ellipticity"] < 0.3]
    ellipticity = np.nanmedian(s["ellipticity"])

    # Select FWHM at least 3 pixels and lower than 15 arcsec
    s = s[(s["fwhm_image"] > 3)*(s["fwhm_image"]*arcsecpix < 15)]

    nsources = len(s) 
    if nsources == 0:
        return 0, 0, 0, 0

    # Select bright magnitudes
    s = s[s["mag"] < np.percentile(s["mag"], 20)]

    nsources = len(s)
    if nsources == 0:
        return 0, 0, 0, 0
       
    fwhm = np.nanmedian(s["fwhm_image"]*arcsecpix)
    bkg = np.nanmedian(s["background"])
    
    return nsources, fwhm, ellipticity, bkg


def get_focus(lfiles, plot=True, interactive=False):
    """
    Receives a list of focus files and returns the best focus value.

    """
    sexfiles = run_sex(lfiles)
    focus, sigma = analyse_sex(sexfiles, plot=plot, interactive=interactive)
    return focus, sigma


def get_focus_ifu(lfiles, plot=True, debug=False, interactive=False):
    """
    Receives a list of focus ifu files and returns the best focus.
    """
    sexfiles = run_sex(lfiles)
    res = analyse_sex_ifu(sexfiles, plot=plot, debug=debug,
                          interactive=interactive)
    
    print(res)
    focus, sigma = res
    return focus, sigma


def get_image_pars(image, arcsecpix=0.394):
    """
    Returns a set of statistics for a given image.
    """
    sexfiles = run_sex([image])
    pars = analyse_image(sexfiles[0], arcsecpix=arcsecpix)
    
    return pars[0], pars[1], pars[2], pars[3]
