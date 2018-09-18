#:4-*- coding: utf-8 -*-
"""
Created on Sat May 21 10:26:37 2016

@author: nadiablago
"""
import datetime
import matplotlib
matplotlib.use("Agg")
import glob, os
import recenter_ifu
import fitsutils
import coordinates_conversor as cc
import numpy as np
import sextractor 
from astropy.io import fits as pf
from matplotlib import pylab as plt
import time_utils
import matplotlib.dates as md
import argparse
import subprocess

from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')

def compile_stats_pointing():
    ra = 0
    dec = 0
    
    out = open("/tmp/pointing", "w")
    out.write("#f, imtype, obj, jd, filter, radeg, decdeg, dra, ddec\n")
    myfiles = glob.glob(_photpath + "/20160616/a_*[0-9].fits")
    myfiles.sort()
    for f in myfiles:
        #try:
        imtype = fitsutils.get_par(f, "IMGTYPE")
	imtype = imtype.upper()
        newra = fitsutils.get_par(f, "OBJRA")
        newdec = fitsutils.get_par(f, "OBJDEC")
        newra, newdec = cc.hour2deg(newra, newdec)        
        myfilter =  fitsutils.get_par(f, "FILTER")
        if ("ACQU" in imtype or imtype == "SCIENCE"):#: and np.round(ra, 2) != np.round(newra, 2) and np.round(dec, 2) != np.round(newdec, 2):
            obj  = fitsutils.get_par(f, "OBJECT")
            jd = fitsutils.get_par(f, "JD")
            status, dra, ddec = recenter_ifu.get_offset_center(f, plot=False, interactive=False)
        
            print f, ra, newra, imtype, jd, dra, ddec
            out.write("%s,%s,%s,%.2f,%s,%.5f,%.5f,%.2f,%.2f\n"%(f, imtype, obj, jd, myfilter, newra, newdec, dra, ddec))
        ra = newra
        dec = newdec
    
        #except:
        #    pass
    out.close()
    
def plot_stats_pointing():
    
    from matplotlib import pylab as plt
    
    d = np.genfromtxt("/tmp/pointing", dtype=None, delimiter=",", names=True)    
    
    
    plt.scatter(d["dra"], d["ddec"], c=d["jd"]-np.min(d["jd"]))
    cb = plt.colorbar()
    cb.set_label("JD- MIN(JD)")
    plt.xlabel("dRA [arcsec]")
    plt.ylabel("dDEC [arcsec]")
    plt.savefig("/tmp/pointing_errors.png")
    
    plt.xlim(-60,60)
    plt.ylim(-60,60)
    plt.savefig("/tmp/pointing_errors_60.png")
    plt.clf()
    
    d1 = d[d["imtype"]!="SCIENCE"]
    d = d[d["imtype"]=="SCIENCE"]
    dr = d[d["filter"]=="r"]
    dg = d[d["filter"]=="g"]
    di = d[d["filter"]=="i"]
    du = d[d["filter"]=="u"]
    
    filters = ["ACQ", "r", "g", "i", "u" ]    
    
    dar = np.array([dr, dg, di, du, d1])
    for i in np.arange(len(filters)):
        plt.figure(i+1)
        plt.scatter(dar[i]["dra"], dar[i]["ddec"], c=dar[i]["jd"]-np.min(d["jd"]))
        cb = plt.colorbar()
        plt.title(filters[i])
        cb.set_label("JD- MIN(JD)")
        plt.xlabel("dRA [arcsec]")
        plt.ylabel("dDEC [arcsec]")
        plt.savefig("/tmp/pointing_errors_%s.png"%filters[i])
        plt.clf()
    plt.close("all")


def get_sextractor_stats(files):
    
    files.sort()
    sexfiles = [os.path.join(os.path.join(os.path.dirname(f), "sextractor"), os.path.basename(f).replace(".fits", ".sex")) for f in files]    
    sexfiles.sort()
    

    if not os.path.isdir(os.path.join( os.path.dirname(files[0]), "stats")):
        os.makedirs(os.path.join(os.path.dirname(files[0]), "stats"))

    with open(os.path.join( os.path.dirname(files[0]), "stats/stats.log"), "w") as out:
        for i, f in enumerate(files):
            try:
                if (fitsutils.has_par(f, "IMGTYPE")):
                    imtype = fitsutils.get_par(f, "IMGTYPE")
                    imtype = imtype.upper()
                else:
                    imtype = "NONE"

                if not ("ACQ" in imtype or imtype == "SCIENCE" or imtype=="FOCUS" or imtype=="GUIDER"):
                    continue 

                if not os.path.isfile(sexfiles[i]):
                    sflist =  sextractor.run_sex([f])
		    if (not sflist is None and len(sflist) > 0):
			sf = sflist[0]
                else:
                    sf = sexfiles[i]

                hd = pf.open(f)[0].header
                try:
                    jd = hd["JD"]
                    obj = hd["OBJECT"]
                    airmass = hd["AIRMASS"]
                    in_temp = hd["IN_AIR"]
                    out_temp = hd["OUT_AIR"]
                    in_hum = hd["IN_HUM"]
                    ns, fwhm, ellipticity, bkg = sextractor.analyse_image(sf)
                    out.write("%s,%s,%.3f,%d,%.2f,%.3f,%.3f,%.2f,%.1f,%s,%.2f,%.2f\n"%(os.path.abspath(f),obj,jd,ns,fwhm,ellipticity,bkg,airmass,in_temp,imtype,out_temp,in_hum))
                except Exception as e:
                    print "Error when retrieving the stats parameters from the header of file %s.\n Error %s"%(f, e)
            except IOError:
                print "Error when opening file %s"%f
            
def plot_stats(statfile):
    
    colors = {"ACQUISITION":"b", "SCIENCE":"r", "FOCUS":"g", "GUIDER":"k"}
    
    s = np.genfromtxt(statfile, delimiter=",", dtype=None)
    s.sort(order="f2")
    s = s[s["f3"]>1]

    day_frac_diff = datetime.timedelta(np.ceil((datetime.datetime.now() - datetime.datetime.utcnow() ).total_seconds())/3600/24)
    datestat = np.array([ time_utils.jd2utc(jd) for jd in s["f2"]])
    datestat = datestat + day_frac_diff
    
    #We add 5h to the UTC date, so it alwasy keeps the date of the end of the night.
    day = ("%s"%(datestat[-1]+datetime.timedelta(5./24))).split()[0]

    xfmt = md.DateFormatter('%H:%M')

    f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
    plt.suptitle("Statistics %s"%day)
    f.set_figwidth(16)
    f.set_figheight(12)
    ax1.plot(datestat, s["f3"], ".-")
    ax1.set_title('Number of bright sources extracted')
    
    for im in set(s["f9"]):
        mask = s["f9"]==im
        ax2.plot(datestat[mask], s["f4"][mask], ".", color=colors[im], label=im)
    ax2.set_title('FWHM [arcsec]')
    ax3.plot(datestat, s["f6"], ".-")
    ax3.set_title('Background')
    ax4.plot(datestat, s["f7"], ".-")
    ax4.set_title('Airmass')
    ax5.plot(datestat, s["f8"], ".-", label="Inside")
    ax5.plot(datestat, s["f10"], ".-", label="Outside")
    #ax5.plot(datestat, s["f11"], ".-")
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

    ax2.legend(labelspacing=0.3, loc="upper right", fontsize=11, numpoints=1, frameon=False, ncol=1, fancybox=False, shadow=True, bbox_to_anchor=(1., 1.))

    ax5.legend(labelspacing=0.3, loc="upper left", fontsize=11, numpoints=1, frameon=False, ncol=1, fancybox=False, shadow=True, bbox_to_anchor=(0., 1.))
    
    plt.savefig(statfile.replace(".log", "%s.png"%(day)), bbox="tight")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.', default=None)

    args = parser.parse_args()
    
    photdir = args.photdir
    print "Parameter directory where stats are run : %s."%photdir
    
    if (photdir is None):
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-","")
        photdir = os.path.join(_photpath, timestamp)
	print "New directory %s"%photdir
    else:
        timestamp=os.path.basename(os.path.abspath(photdir))
    print "Running stats on", glob.glob(os.path.join(os.path.abspath(photdir), "rc*[0-9].fits"))
    get_sextractor_stats(glob.glob(os.path.join(os.path.abspath(photdir), "rc*[0-9].fits")))
    plot_stats(os.path.join(os.path.abspath(photdir), "stats/stats.log")) 

