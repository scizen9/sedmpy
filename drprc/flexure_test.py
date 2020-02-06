# -*- coding: utf-8 -*-
"""
Created on Mon May  9 12:25:51 2016
Updated on Thu Feb  6 2020 by neill

@author: nblago
"""

try:
    import sextractor
except ImportError:
    import drprc.sextractor as sextractor
try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils
import glob
import os
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib import pylab as plt
import datetime
import argparse
import logging

# Log into a file
log_format = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
root_dir = "/tmp/"
now = datetime.datetime.utcnow()
timestamp = datetime.datetime.isoformat(now)
timestamp = timestamp.split("T")[0]
logging.basicConfig(format=log_format,
                    filename=os.path.join(root_dir,
                                          "listener_{0}.log".format(timestamp)),
                    level=logging.INFO)
logger = logging.getLogger('flexure')


def run_flexure_test(sexfiles, plotdir):
    """
    Compares the positions for the extracted lines
    for all the files given as a parameter.
    """
    sexfiles.sort()
    posfiles = []

    logger.info("Running flexure test with %s files." % sexfiles)
    with open(os.path.join(plotdir, "flexure.log"), "w") as out:
        for sf in sexfiles:
            c = np.genfromtxt(sf)

            # We want bright sources.
            c = c[c[:, 4] < -9]

            # We want round-shaped sources.
            c = c[c[:, 8] < 0.6]

            # Save the positions in a separate file with only X Y.
            np.savetxt(sf.replace(".sex", ".pos"), c[:, 0:2])
            posfiles.append(sf.replace(".sex", ".pos"))

        c0 = np.genfromtxt(posfiles[0])
        i = 1
        while len(c0) == 0 and i < len(posfiles):
            c0 = np.genfromtxt(posfiles[i])
            i = i + 1

        flexure = []
        for f in posfiles[i:]:
            print(f)
            c1 = np.genfromtxt(f)

            c = SkyCoord(x=c0[:, 0], y=c0[:, 1], z=np.zeros(len(c0)), unit='m',
                         representation='cartesian')
            catalog = SkyCoord(x=c1[:, 0], y=c1[:, 1], z=np.zeros(len(c1)),
                               unit='m', representation='cartesian')
            # Now idx are indices into catalog that are the closest objects to
            # each of the coordinates in c, d2d are the
            # on-sky distances between them.
            idx, d2d, d3d = c.match_to_catalog_sky(catalog)

            d3d = d3d[d3d.value < 10]

            # Add the flexure for plotting later.
            flexure.append(np.median(d3d.value))

            out.write("%s,%s,%.3f\n" %
                      (os.path.basename(posfiles[0]), os.path.basename(f),
                       float(np.median(d3d.value))))

    plot_flexure(os.path.join(plotdir, "flexure.log"), plotdir)


def plot_flexure(flexfile, plotdir):
    import matplotlib.dates as md

    xfmt = md.DateFormatter('%d %b %H h')

    flexure = np.genfromtxt(flexfile, dtype=None, delimiter=",")

    dateflex = []
    for fl in flexure:
        # Format of the file: ifu20160630_00_08_38
        t = datetime.datetime.strptime(fl["f1"], 'ifu%Y%m%d_%H_%M_%S.pos')
        dateflex.append(t)

    plt.plot(dateflex, flexure["f2"], "o-")
    plt.xlabel("Date")
    plt.ylabel("Median (flexure) [pixels]")
    plt.gca().xaxis.set_major_formatter(xfmt)
    labels = plt.gca().get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)

    curdate = datetime.datetime.strftime(dateflex[-1], "%Y%m%d")

    plt.savefig(os.path.join(plotdir, "flexure_%s.png" % curdate),
                bbox_inches='tight')

    plt.clf()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object 
        coordinates in the reference pixel.
            
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', dest="raw", type=str,
                        help='Directory containing the raw fits for the night.',
                        default=None)

    args = parser.parse_args()

    timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
    timestamp = timestamp.split("T")[0].replace("-", "")

    if args.raw is None:
        raw = os.path.join("/scr2/sedm/raw/", timestamp)
    else:
        raw = args.raw

    phot = os.path.join("/scr2/sedmdrp/phot/", timestamp)

    # Get the files from the same day directory.
    files = glob.glob(os.path.join(raw, "ifu*fits"))

    # Get the files from one day before directory
    # to see the differences between days.
    daybefore = datetime.datetime.isoformat(
        datetime.datetime.utcnow() - datetime.timedelta(1))
    daybefore = daybefore.split("T")[0].replace("-", "")
    daybefore = os.path.join("/scr2/sedm/raw/", daybefore)

    if os.path.isdir(daybefore):
        filesold = glob.glob(os.path.join(daybefore, "ifu*fits"))
        files.extend(filesold)

    files_hg = [
        f for f in files if fitsutils.has_par(f, "OBJECT") and
        "Calib: Hg" in fitsutils.get_par(f, "OBJECT")]

    logger.info("Found the following Hg files: %s" % files_hg)
    if len(files_hg) > 1:
        files_hg.sort()
        sexlist = sextractor.run_sex(files_hg, mask=False)

        plot_dir = os.path.join(phot, "stats")
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        run_flexure_test(sexlist, plotdir=plot_dir)
