from __future__ import print_function
# -*- coding: utf-8 -*-
"""
Created on Sat May 23 18:23:02 2015

@author: nadiablago
"""

from matplotlib import pylab as plt

import numpy as np
import datetime

import os
import math
from astropy.io import fits as pf

try:
    import zscale
except ImportError:
    import drprc.zscale as zscale
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
import glob
import argparse
try:
    import coordinates_conversor as cc
except ImportError:
    import drprc.coordinates_conversor as cc
from astropy.wcs import WCS
import matplotlib.lines as mlines
import subprocess
import warnings

from photutils import SkyCircularAperture

from configparser import ConfigParser
import codecs

parser = ConfigParser()

configfile = os.environ['SEDMCONFIG']

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.read_file(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')


def fxn():
    warnings.warn("deprecated", DeprecationWarning)


def get_app_phot(coords, image, plot_only=False, wcsin="world",
                 fwhm=2.5, plotdir=None, box=15, arcsecpix=0.394):
    """
    coords: files:
    wcsin: can be "world", "logic"
    """
    # Load packages; splot is in the onedspec package, which is in noao. 
    # The special keyword _doprint=0 turns off displaying the tasks 
    # when loading a package.

    imdir = os.path.dirname(image)
    imname = os.path.basename(image)
    
    if plotdir is None:
        plotdir = os.path.join(imdir, "photometry")
    
    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)
        
    out_name = os.path.join(plotdir, imname + ".seq.mag")
    clean_name = os.path.join(plotdir, imname + ".app.mag")

    print("Will create output files", out_name, clean_name)

    # Read values from .ec file
    
    fwhm_value = fwhm/arcsecpix
    
    if fitsutils.has_par(image, 'FWHM'):
        fwhm_value = fitsutils.get_par(image, 'FWHM')/arcsecpix
    elif fwhm is None:
        fwhm_value = 3.5/arcsecpix
    # if fitsutils.has_par(image, 'AIRMASS'):
        # airmass_value = fitsutils.get_par(image, 'AIRMASS')
    # else:
        # airmass_value = 1.3
 
    # exptime = fitsutils.get_par(image, 'EXPTIME')
    gain = fitsutils.get_par(image, 'GAIN')
    # noise = fitsutils.get_par(image, 'RDNOISE')

    print("FWHM: %.1f pixels, %.1f arcsec" % (fwhm_value, fwhm_value*arcsecpix))
    # Set aperture radius to three times the PSF radius
    aperture_rad = math.ceil(float(fwhm_value)*1.5)
    sky_rad = math.ceil(aperture_rad)*4

    if not plot_only:

        iraf.noao.digiphot.apphot.qphot(image=image,
                                        cbox=box,
                                        annulus=sky_rad,
                                        dannulus=20.,
                                        aperture=str(aperture_rad),
                                        coords=coords,
                                        output=out_name,
                                        plotfile="",
                                        zmag=0.,
                                        exposure="exptime",
                                        airmass="airmass",
                                        filter="filter",
                                        obstime="DATE",
                                        fwhm=fwhm_value,
                                        epadu=gain,
                                        interactive="no",
                                        radplots="yes",
                                        verbose="no",
                                        graphics="stdgraph",
                                        display="stdimage",
                                        icommands="",
                                        wcsin=wcsin,
                                        wcsout="logical",
                                        gcommands="")

        # iraf.noao.digiphot.apphot.phot(image=image, cbox=5., annulus=12.4,
        # dannulus=10., salgori = "centroid", aperture=9.3,wcsin="world",
        # wcsout="tv", interac = "no", coords=coords, output=out_name)
        iraf.txdump(out_name,
                    "id,image,xcenter,ycenter,xshift,yshift,fwhm,msky,stdev,"
                    "cier,rapert,sum,area,nsky,flux,itime,mag,merr", "yes",
                    Stdout=clean_name)

    ma = np.genfromtxt(clean_name, comments="#", dtype=[("id", "<f4"),
                                                        ("image", "|S20"),
                                                        ("X", "<f4"),
                                                        ("Y", "<f4"),
                                                        ("Xshift", "<f4"),
                                                        ("Yshift", "<f4"),
                                                        ("fwhm", "<f4"),
                                                        ("msky", "<f4"),
                                                        ("stdev", "<f4"),
                                                        ("flags", np.int),
                                                        ("rapert", "<f4"),
                                                        ("sum", "<f4"),
                                                        ("area", "<f4"),
                                                        ("nsky", "<f4"),
                                                        ("flux", "<f4"),
                                                        ("itime", "<f4"),
                                                        ("fit_mag", "<f4"),
                                                        ("fiterr", "<f4")])
    if ma.size > 0:
        m = ma[~np.isnan(ma["fit_mag"])]
    else:
        print("Only one object found!")
        m = np.array([ma])
        
    hdulist = pf.open(image)
    # prihdr = hdulist[0].header
    img = hdulist[0].data * 1.
    # nx, ny = img.shape

    dim_x = int(np.floor(np.sqrt(len(m))))
    dim_y = int(np.ceil(len(m)*1./dim_x))
    outerrad = sky_rad+10
    cutrad = outerrad + 15
    
    print("Cutrad %.1f" % cutrad)

    plt.suptitle("FWHM=%.2f arcsec. %d stars" % (fwhm_value*arcsecpix, len(m)))
    for i in np.arange(dim_x):
        for j in np.arange(dim_y):
            if i*dim_y + j < len(m):
                k = i*dim_y + j
                # print (dimX, dimY, i, j, k)
                ax = plt.subplot2grid((dim_x, dim_y), (i, j))
                y1, y2, x1, x2 = m[k]["X"]-cutrad, m[k]["X"]+cutrad, \
                    m[k]["Y"]-cutrad, m[k]["Y"]+cutrad
                y1, y2, x1, x2 = int(y1), int(y2), int(x1), int(x2)
                y1 = np.maximum(y1, 0)
                y2 = np.maximum(y2, 0)
                x1 = np.maximum(x1, 0)
                x2 = np.maximum(x2, 0)
                zmin = None
                zmax = None
                try:
                    zmin, zmax = zscale.zscale(img[x1:x2, y1:y2], nsamples=1000,
                                               contrast=0.25)
                except ValueError:
                    print(y1, y2, x1, x2)
                    print(img[x1:x2, y1:y2])
                    sh = img[x1:x2, y1:y2].shape
                    if sh[0] > 0 and sh[1] > 0:
                        zmin = np.nanmin(img[x1:x2, y1:y2])
                        zmax = np.nanmax(img[x1:x2, y1:y2])
                ax.imshow(img[x1:x2, y1:y2], aspect="equal",
                          extent=(-cutrad, cutrad, -cutrad, cutrad),
                          cmap=plt.get_cmap('Grays_r'),
                          interpolation="none", vmin=zmin, vmax=zmax)
                c1 = plt.Circle((0, 0), edgecolor="r", facecolor="none",
                                radius=aperture_rad)
                c2 = plt.Circle((0, 0), edgecolor="orange", facecolor="none",
                                radius=sky_rad)
                c3 = plt.Circle((0, 0), edgecolor="yellow", facecolor="none",
                                radius=sky_rad+20)
                plt.gca().add_artist(c1)
                plt.gca().add_artist(c2)
                plt.gca().add_artist(c3)
                ax.set_xticks([])
                ax.set_yticks([])
        
                plt.text(+5, +5, "%d" % m[k]["id"])
                plt.text(-cutrad, -cutrad, "%.2f$\pm$%.2f" % (m[k]["fit_mag"],
                                                              m[k]["fiterr"]),
                         color="b")

    plt.tight_layout()
    plt.savefig(os.path.join(plotdir, imname + "plot.png"))
    plt.clf()


def get_xy_coords(image, ra, dec):
    """
    Uses the wcs-rd2xy routine to compute the proper pixel number where the
    target is. Sometime the wcs module does not seem to be providing the correct
    answer, as it does not seem to be using the SIP extension.

    """
    w = WCS(image)
    pix_all = w.all_world2pix(ra, dec, 0)
    coords = [float(pix_all[0]), float(pix_all[1])]
        
    return coords


def get_app_phot_target(image, ra=None, dec=None, plot=True, wcsin="logical",
                        fwhm=None, box=15, arcsecpix=0.394, app=2):
    """
    coords: files:
    wcsin: can be "world", "logic"
    fwhm: in arcsec

    """
    # Load packages; splot is in the onedspec package, which is in noao. 
    # The special keyword _doprint=0 turns off displaying the tasks 
    # when loading a package.
    
    impf = pf.open(image)
    wcs = WCS(impf[0].header)
    # Check that actually the object is within this frame.
    if ra is None or dec is None:
        if fitsutils.has_par(image, "OBJRA") and fitsutils.has_par(image,
                                                                   "OBJRA"):
            ra, dec = cc.hour2deg(fitsutils.get_par(image, 'OBJRA'),
                                  fitsutils.get_par(image, 'OBJDEC'))
        else:
            ra, dec = cc.hour2deg(fitsutils.get_par(image, 'RA'),
                                  fitsutils.get_par(image, 'DEC'))
        print("Assuming ra=%.5f, dec=%.5f" % (ra, dec))
        pra, pdec = get_xy_coords(image, ra, dec)

    else:
        if "logic" in wcsin:
            pra, pdec = ra, dec
        else:
            # Using new method to derive the X, Y pixel coordinates,
            # as wcs module does not seem to be working well.
            try:
                # pra, pdec = wcs.wcs_sky2pix(ra, dec, 1)
                # print ("Retrieved the pixel number")
                pra, pdec = get_xy_coords(image, ra, dec)
            except IndexError:
                print("Error with astrometry.net."
                      " trying the rudimentary method.")
                pra, pdec = wcs.wcs_world2pix(ra, dec, 1)
            # pra, pdec = wcs.wcs_sky2pix(np.array([ra, dec], ndmin=2), 1)[0]

    shape = impf[0].data.shape
    
    if (pra > 0) and (pra < shape[0]) and (pdec > 0) and (pdec < shape[1]):
        pass
    else:
        print(image, "ERROR! Object coordinates are outside this frame."
                     " Skipping any aperture photometry!!")
        print(pra, pdec, shape)
        return

    imdir = os.path.dirname(image)
    imname = os.path.basename(image)
    plotdir = os.path.join(imdir, "photometry")

    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)
        
    out_name = os.path.join(plotdir, imname + ".seq.mag")
    clean_name = os.path.join(plotdir, imname + ".objapp.mag")
    
    if fwhm is not None:
        fwhm_value = fwhm
    elif fitsutils.has_par(image, 'FWHM'):
        fwhm_value = fitsutils.get_par(image, 'FWHM')
    else:
        # Put some default value for Palomar
        fwhm_value = 1.5
        
    if wcsin == 'logical':
        fwhm_value = fwhm_value / arcsecpix 
        
    # if fitsutils.has_par(image, 'AIRMASS'):
    #    airmass_value = fitsutils.get_par(image, 'AIRMASS')
    # else:
    #    airmass_value = 1.3
        
    if not fitsutils.has_par(image, "EXPTIME"):
        if fitsutils.has_par(image, "ITIME") and fitsutils.has_par(image,
                                                                   "COADDS"):
            exptime = fitsutils.get_par(image, "ITIME") * \
                      fitsutils.get_par(image, "COADDS")
            fitsutils.update_par(image, "EXPTIME", exptime)
    # exptime = fitsutils.get_par(image, 'EXPTIME')
    gain = fitsutils.get_par(image, 'GAIN')
    
    # print ("FWHM", fwhm_value)
    # Set aperture radius to two times the PSF radius
    aperture_rad = math.ceil(float(fwhm_value)*app)
    sky_rad = math.ceil(aperture_rad*app*2)
    
    # print (aperture_rad, sky_rad)

    print("Saving coodinates for the object in pixels", pra, pdec)
    coords = "/tmp/coords.dat"    
    np.savetxt("/tmp/coords.dat", np.array([[pra, pdec]]), fmt="%.4f %.4f")

    if os.path.isfile(out_name):
        os.remove(out_name)
    if os.path.isfile(clean_name):
        os.remove(clean_name)

    iraf.noao.digiphot.apphot.qphot(image=image,
                                    cbox=box,
                                    annulus=sky_rad,
                                    dannulus=20.,
                                    aperture=str(aperture_rad),
                                    coords=coords,
                                    output=out_name,
                                    plotfile="",
                                    zmag=0.,
                                    exposure="exptime",
                                    airmass="airmass",
                                    filter="filter",
                                    obstime="DATE",
                                    epadu=gain,
                                    interactive="no",
                                    radplots="yes",
                                    verbose="no",
                                    graphics="stdgraph",
                                    display="stdimage",
                                    icommands="",
                                    wcsin="logical",
                                    wcsout="logical",
                                    gcommands="")

    # iraf.noao.digiphot.apphot.phot(image=image, cbox=5., annulus=12.4,
    # dannulus=10., salgori = "centroid", aperture=9.3,wcsin="world",
    # wcsout="tv", interac = "no", coords=coords, output=out_name)
    iraf.txdump(out_name, "id,image,xcenter,ycenter,xshift,yshift,fwhm,msky,"
                          "stdev,cier,rapert,sum,area,nsky,flux,itime,mag,merr",
                "yes", Stdout=clean_name)

    ma = np.genfromtxt(clean_name, comments="#", dtype=[("id", "<f4"),
                                                        ("image", "|S20"),
                                                        ("X", "<f4"),
                                                        ("Y", "<f4"),
                                                        ("Xshift", "<f4"),
                                                        ("Yshift", "<f4"),
                                                        ("fwhm", "<f4"),
                                                        ("msky", "<f4"),
                                                        ("stdev", "<f4"),
                                                        ("flags", np.int),
                                                        ("rapert", "<f4"),
                                                        ("sum", "<f4"),
                                                        ("area", "<f4"),
                                                        ("nsky", "<f4"),
                                                        ("flux", "<f4"),
                                                        ("itime", "<f4"),
                                                        ("fit_mag", "<f4"),
                                                        ("fiterr", "<f4")])
    if ma.size > 0:
        ma = np.array([ma])
        # m = ma[~np.isnan(ma["fit_mag"])]
    else:
        m = np.array([ma])
        print("Only one object found!", m)

    insmag = np.round(ma['fit_mag'][0], 3)
    insmagerr = np.round(ma['fiterr'][0], 3)  
    if fitsutils.has_par(image, "ZEROPT") and fitsutils.has_par(image,
                                                                "ZEROPTU"):
        mag = insmag + float(fitsutils.get_par(image, "ZEROPT"))
        magerr = np.sqrt(insmagerr**2 + float(fitsutils.get_par(image,
                                                                "ZEROPTU"))**2)
    else:
        mag = 0
        magerr = 0

    if np.isnan(mag):
        mag, magerr = 0, 0
        insmag, insmagerr = 0, 0

    fitsutils.update_par(image, "INSMAG", "%.3f" % insmag)
    fitsutils.update_par(image, "INSMAGER", "%.3f" % insmagerr)
    fitsutils.update_par(image, "APPMAG", np.round(mag, 3))
    fitsutils.update_par(image, "APPMAGER", np.round(magerr, 3))
         
    if plot:
        
        # zmin, zmax = zscale.zscale(impf[0].data.T[
        # pra-50:pra+50,pdec-50:pdec+50])
        
        # zmin, zmax = zscale.zscale(impf[0].data)
           
        # im = plt.imshow(impf[0].data, vmin=zmin, vmax=zmax)
        print(np.percentile(impf[0].data, 5), np.percentile(impf[0].data, 95))
        impf[0].data[np.isnan(impf[0].data)] = np.nanmedian(impf[0].data)
        print(np.percentile(impf[0].data, 5), np.percentile(impf[0].data, 95))

        im = plt.imshow(impf[0].data, vmin=np.percentile(impf[0].data, 5),
                        vmax=np.percentile(impf[0].data, 95))
       
        xx = int(ma["X"][0])
        yy = int(ma["Y"][0])
        pra = int(pra)
        pdec = int(pdec)
        
        plt.scatter(xx, yy, marker="o", s=100, facecolor="none",
                    edgecolor="red")
        plt.colorbar(im)
        plt.savefig(os.path.join(plotdir, imname+".png"), dpi=200)
        plt.clf()
        
        zmin, zmax = zscale.zscale(impf[0].data.T[xx-50:xx+50, yy-50:yy+50].T)
        im = plt.imshow(impf[0].data.T[pra-50:pra+50, pdec-50:pdec+50].T,
                        vmin=zmin, vmax=zmax, interpolation="none",
                        extent=(-50, 50, -50, 50))
        c1 = plt.Circle((pra-xx, pdec-yy), edgecolor="k", facecolor="none",
                        radius=aperture_rad, label="Initial position")
        c11 = plt.Circle((pra-xx, pdec-yy), edgecolor="k", facecolor="none",
                         radius=sky_rad)
        c2 = plt.Circle((0, 0), edgecolor="orange", facecolor="none",
                        radius=aperture_rad, label="Adjusted centroid")
        c22 = plt.Circle((0, 0), edgecolor="orange", facecolor="none",
                         radius=sky_rad)
        plt.gca().add_artist(c1)
        plt.gca().add_artist(c11)
        plt.gca().add_artist(c2)
        plt.gca().add_artist(c22)
        plt.colorbar(im)
        
        myhandles = []
        markers = ["o", "o"]
        labels = ["Initial position", "Adjusted centroid"]
        cols = ["k", "orange"]
        for i in np.arange(len(markers)):
            myhandles.append(mlines.Line2D([], [], mec=cols[i], mfc="none",
                             marker=markers[i], ls="None", markersize=10,
                             label=labels[i]))
        plt.legend(handles=myhandles, loc="lower left", labelspacing=0.3,
                   fontsize=11, numpoints=1, frameon=False, ncol=5,
                   bbox_to_anchor=(0.0, 0.00), fancybox=False, shadow=True)

        plt.title("MIN: %.0f MAX: %.0f" %
                  (float(np.nanmin(impf[0].data.T[xx-50:xx+50, yy-50:yy+50])),
                   float(np.nanmax(impf[0].data.T[xx-50:xx+50, yy-50:yy+50]))))
        plt.savefig(os.path.join(plotdir, imname+"_zoom.png"))
        plt.clf()


def get_upper_limit(app_phot_file, sigma=3, zp=0):
    """
    Provides the upper X-sigma limit for that aperture.
    """
    ff = np.genfromtxt(app_phot_file, comments="#", dtype=[("id", "<f4"),
                                                           ("image", "|S20"),
                                                           ("X", "<f4"),
                                                           ("Y", "<f4"),
                                                           ("Xshift", "<f4"),
                                                           ("Yshift", "<f4"),
                                                           ("fwhm", "<f4"),
                                                           ("msky", "<f4"),
                                                           ("stdev", "<f4"),
                                                           ("flags", np.int),
                                                           ("rapert", "<f4"),
                                                           ("sum", "<f4"),
                                                           ("area", "<f4"),
                                                           ("nsky", "<f4"),
                                                           ("flux", "<f4"),
                                                           ("itime", "<f4"),
                                                           ("fit_mag", "<f4"),
                                                           ("fiterr", "<f4")])
    
    rmsflux = sigma * ff['stdev'] * np.sqrt(ff['area'])
    uplimmag = zp - 2.5*np.log10(rmsflux) + 2.5*np.log10(ff['itime'])
    
    return uplimmag


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""

        Aperture photometry package
            
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--reddir', type=str, dest="reduced",
                        help='Fits directory file with reduced images.',
                        default=None)

    args = parser.parse_args()
    
    reduced = args.reduced
    
    if reduced is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        reduced = os.path.join(_photpath, timestamp, "reduced")

    os.chdir(reduced)

    for f in glob.glob("*.fits"):
        if fitsutils.has_par(f, "IMGTYPE") and \
                fitsutils.get_par(f, "IMGTYPE") == "SCIENCE" or \
                fitsutils.get_par(f, "IMGTYPE") == "ACQUISITION":
            # print (f)
            get_app_phot_target(f, box=5)
         
    cmd = 'echo "FILE ONTARGET NAME FILTER JD LST APPMAG APPMAGER INSMAG ' \
          'INSMAGER ZEROPT ZEROPTU" > photometry/magnitudes.dat; gethead ' \
          'ONTARGET NAME FILTER JD LST APPMAG APPMAGER INSMAG INSMAGER ' \
          'ZEROPT ZEROPTU *fits | grep -E "fits[\ ]+1" >> ' \
          'photometry/magnitudes.dat'
    subprocess.call(cmd, shell=True)
