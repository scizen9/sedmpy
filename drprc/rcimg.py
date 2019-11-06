"""RCam Image functions

    Functions:
    * :func:`get_seeing`      Use sextractor to measure the FWHM seeing
"""
import os
import logging
import subprocess
import shutil
import numpy as np
from astropy.io.ascii import SExtractor
from matplotlib import pyplot as plt

logging.basicConfig(
    format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
    datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)


def get_seeing(imfile=None, destdir=None, save_fig=False, show_plot=False):
    """Measure the FWHM of stars in image"""
    seeing = -1.
    if imfile:
        # Guider image upon which to work
        fspec = os.path.join(destdir, 'guider_%s' % imfile)
        if os.path.exists(fspec):
            logging.info("Attempting to measure FWHM seeing in %s" % fspec)
            # run SExtractor
            cmd = ["sex", "-c", "%s" %
                   os.path.join(os.environ["SEDMPY"], "drprc", "config",
                                "sedmrc.sex"), fspec]
            logging.info(" ".join(cmd))
            retcode = subprocess.call(cmd)
            if retcode == 0:
                logging.info("SExtractor successful")
                sexfile = fspec.replace(".fits", ".sex")
                shutil.move("image.sex", sexfile)
                seeing = analyze_sex_cat(sexfile, save_fig=save_fig,
                                         show_plot=show_plot)
            else:
                logging.warning("SExtractor failed")
        else:
            logging.warning("%s not found, no seeing measured" % fspec)
    else:
        logging.warning("No image given so no seeing measured")

    return seeing
    # END: get_seeing


def analyze_sex_cat(catalog, show_plot=False, save_fig=False):
    """Get the FWHM in world coordinates from SExtractor"""
    seeing = -1.
    sex = SExtractor()
    data = sex.read(catalog)
    h = data.to_pandas()
    # Clean our list
    # No flags
    df = h.loc[h['FLAGS'] == 0]
    # Avoid crossbar shadows
    df = df[(df['X_IMAGE'] < 850) | (df['X_IMAGE'] > 1160)]
    df = df[(df['Y_IMAGE'] < 840) | (df['Y_IMAGE'] > 1170)]
    # Make sure we are round-ish
    df = df[(df['ELLIPTICITY'] < .3) & (df['FWHM_WORLD'] > 0.)]
    # Do we have enough good stars?
    if len(df) > 10:
        # Get stats on fwhms
        fwhms = df['FWHM_WORLD']*3600.  # convert from degrees to arcsecs
        seeing = np.nanmedian(fwhms)
        mean = np.nanmean(fwhms)
        # Plot or save?
        if show_plot or save_fig:
            if show_plot:
                plt.ion()
            else:
                plt.ioff()
            plt.hist(fwhms, bins=20)
            ylims = plt.gca().get_ylim()
            plt.title(catalog.split('/')[-1] + ' (%.2f ")' % seeing)
            plt.xlabel("Seeing FWHM (arcsec)")
            plt.ylabel("N stars")
            plt.plot([seeing, seeing], ylims, "-", label="Med")
            plt.plot([mean, mean], ylims, "-", label="Mean")
            plt.legend()
            if save_fig:
                plt.savefig(catalog.replace(".sex", "_seeing.png"))
            if show_plot:
                plt.show()
    else:
        logging.warning("Not enough good stars, seeing not measured")

    return seeing
