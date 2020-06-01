# :4-*- coding: utf-8 -*-
"""
Created on Sat May 21 10:26:37 2016

@author: nadiablago
"""
import datetime
import glob
import os

try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
try:
    import coordinates_conversor as cc
except ImportError:
    import drprc.coordinates_conversor as cc
import numpy as np
try:
    import sextractor
except ImportError:
    import drprc.sextractor as sextractor
from astropy.io import fits as pf
from matplotlib import pylab as plt
try:
    import time_utils
except ImportError:
    import drprc.time_utils as time_utils
import matplotlib.dates as md
import argparse

from configparser import ConfigParser
import codecs

parser = ConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.read_file(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')

plt.switch_backend('Agg')


def get_sextractor_stats(files):
    
    files.sort()
    sexfiles = [os.path.join(os.path.join(os.path.dirname(ff),
                                          "sextractor"),
                             os.path.basename(ff).replace(".fits", ".sex"))
                for ff in files]
    sexfiles.sort()

    if not os.path.isdir(os.path.join(os.path.dirname(files[0]), "stats")):
        os.makedirs(os.path.join(os.path.dirname(files[0]), "stats"))

    with open(os.path.join(os.path.dirname(files[0]),
                           "stats/stats.log"), "w") as out:
        for i, ff in enumerate(files):
            try:
                if fitsutils.has_par(ff, "IMGTYPE"):
                    imtype = fitsutils.get_par(ff, "IMGTYPE")
                    imtype = imtype.upper()
                else:
                    imtype = "NONE"

                if not ("ACQ" in imtype or imtype == "SCIENCE" or
                        imtype == "FOCUS" or imtype == "GUIDER"):
                    continue 

                sflist = []
                if not os.path.isfile(sexfiles[i]):
                    sflist = sextractor.run_sex([ff])
                if sflist and len(sflist) > 0:
                    sf = sflist[0]
                else:
                    sf = sexfiles[i]

                hd = pf.open(ff)[0].header
                try:
                    # get observation parameters
                    jd = hd["JD"]
                    obj = hd["OBJECT"]
                    airmass = hd["AIRMASS"]
                    in_temp = hd["IN_AIR"]
                    out_temp = hd["OUT_AIR"]
                    in_hum = hd["IN_HUM"]
                    # get weather string (includes imtype)
                    weather_string = ''
                    if in_temp < -99:
                        weather_string += 'nan,'
                    else:
                        weather_string += '%.1f,' % in_temp
                    weather_string += imtype
                    if out_temp < -99:
                        weather_string += ',nan,'
                    else:
                        weather_string += ',%.2f,' % out_temp
                    if in_hum < -99:
                        weather_string += ',nan'
                    else:
                        weather_string += ',%.2f' % in_hum
                    # get image stats
                    ns, fwhm, ellipticity, bkg = sextractor.analyse_image(sf)
                    out.write(
                        "%s,%s,%.3f,%d,%.2f,%.3f,%.3f,%.2f,%s\n"
                        % (os.path.abspath(ff), obj, jd, ns, fwhm, ellipticity,
                           bkg, airmass, weather_string))
                except Exception as e:
                    print("Error when retrieving the stats parameters from the "
                          "header of file %s.\n Error %s" % (ff, e))
            except IOError:
                print("Error when opening file %s" % ff)


def plot_stats(statfile):
    
    colors = {"ACQUISITION": "b", "SCIENCE": "r", "FOCUS": "g", "GUIDER": "k"}
    
    s = np.genfromtxt(statfile, delimiter=",", dtype=None, encoding=None)
    s.sort(order="f2")
    s = s[s["f3"] > 1]

    day_frac_diff = datetime.timedelta(
        np.ceil((datetime.datetime.now() -
                 datetime.datetime.utcnow()).total_seconds()) / 3600 / 24)
    datestat = np.array([time_utils.jd2utc(jd) for jd in s["f2"]])
    datestat = datestat + day_frac_diff
    
    # We add 5h to the UTC date, so it always keeps
    # the date of the end of the night.
    day = ("%s" % (datestat[-1]+datetime.timedelta(5./24))).split()[0]

    xfmt = md.DateFormatter('%H:%M')

    fg, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
    plt.suptitle("Statistics %s" % day)
    fg.set_figwidth(14)
    fg.set_figheight(11)
    ax1.plot(datestat, s["f3"], ".-")
    ax1.set_title('Number of bright sources extracted')
    
    for im in set(s["f9"]):
        mask = s["f9"] == im
        ax2.plot(datestat[mask], s["f4"][mask], ".", color=colors[im], label=im)
    ax2.set_title('FWHM [arcsec]')
    ax3.plot(datestat, s["f6"], ".-")
    ax3.set_title('Background')
    ax4.plot(datestat, s["f7"], ".-")
    ax4.set_title('Airmass')
    # get rid of bad values
    tin = s["f8"]
    tout = s["f10"]
    tin[tin < -100] = np.nan
    tout[tout < -100] = np.nan
    ax5.plot(datestat, tin, ".-", label="Inside")
    ax5.plot(datestat, tout, ".-", label="Outside")
    # ax5.plot(datestat, s["f11"], ".-")
    ax5.set_title('Temperature')
    ax6.plot(datestat, s["f5"], ".-")
    ax6.set_title('Ellipticity')
    
    ax1.xaxis.set_major_formatter(xfmt)
    ax2.xaxis.set_major_formatter(xfmt)
    ax3.xaxis.set_major_formatter(xfmt)
    ax4.xaxis.set_major_formatter(xfmt)
    ax5.xaxis.set_major_formatter(xfmt)
    ax6.xaxis.set_major_formatter(xfmt)

    labels = ax1.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)
    labels = ax2.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)
    labels = ax3.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)
    labels = ax4.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)
    labels = ax5.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)
    labels = ax6.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)

    ax2.legend(labelspacing=0.3, loc="upper right", fontsize=11, numpoints=1,
               frameon=True, ncol=1, fancybox=False, shadow=True,
               bbox_to_anchor=(1., 1.))

    ax5.legend(labelspacing=0.3, loc="upper left", fontsize=11, numpoints=1,
               frameon=True, ncol=1, fancybox=False, shadow=True,
               bbox_to_anchor=(0., 1.))
    
    plt.savefig(statfile.replace(".log", "%s.png" % day), bbox="tight")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""

        Runs sextractor on the rc images in the specified directory and
        plots image quality statistics
            
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--photdir', type=str, dest="photdir",
                        help='Fits directory file with tonight images.',
                        default=None)

    args = parser.parse_args()

    if args.photdir is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        photdir = os.path.join(_photpath, timestamp)
        print("Auto directory %s" % photdir)
    else:
        photdir = args.photdir
        timestamp = os.path.basename(os.path.abspath(photdir))
        print("Input directory %s" % photdir)

    rclist = glob.glob(os.path.join(os.path.abspath(photdir), "rc*[0-9].fits"))
    print("Running stats on %d rc images in %s" % (len(rclist), photdir))
    if len(rclist) > 0:
        get_sextractor_stats(rclist)
        plot_stats(os.path.join(os.path.abspath(photdir), "stats/stats.log"))
