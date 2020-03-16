# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:18:46 2015
Updated on Wed Feb 12 2020: neill

@author: nadiablago
"""

from matplotlib import pylab as plt
import glob
from astropy.io import fits as pf
import os
import argparse
import numpy as np
from astropy.wcs import WCS
try:
    import zscale
except ImportError:
    import drprc.zscale as zscale
try:
    import time_utils
except ImportError:
    import drprc.time_utils as time_utils


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="""
        FITS file utilites
    """)
    parser.add_argument("-d", "--dir", type=str, default=".",
                        help="Directory to be plotted")
    parser.add_argument("-a", "--ra", type=float, default=0.0, help="ra")
    parser.add_argument("-b", "--dec", type=float, default=0.0, help="dec")
    parser.add_argument("-p", "--pattern", type=str, default="*.fits",
                        help="pattern")

    args = parser.parse_args()

    plot_dir = os.path.join(args.dir, "png")

    print(args.ra, args.dec)
    
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    for fn in glob.glob(args.dir + "/" + args.pattern):
        print("Plotting", os.path.basename(fn).replace(".fits", ".png"))
        hdul = pf.open(fn)
        if len(hdul) > 1:
            findices = np.arange(len(hdul)-1)+1
        else:
            findices = np.array([0])
        for fi in findices:
        
            hdr0 = hdul[fi].header
            image = hdul[fi].data * 1.
            
            im_nx, im_ny = image.shape
            if args.ra * args.dec != 0:
    
                # Get pixel coordinates of SN
                im_wcs = WCS(hdr0)
                try:
                    targ_pix = im_wcs.wcs_world2pix([(np.array([args.ra,
                                                                args.dec],
                                                               np.float_))],
                                                    1)[0]
                except:
                    print("ERROR when converting sky to wcs. Is astrometry in "
                          "place? Default coordinates assigned.")
                    targ_pix = [+im_nx/2., im_ny/2.]
    
                print(targ_pix)
            else:
                targ_pix = [+im_nx/2., im_ny/2.]
               
            sub_img = image - np.nanmin(image)
            av = np.median(sub_img.flatten())
            mi, ma = zscale.zscale(image)
            fg = plt.imshow(plt.log10(image), aspect="equal",
                            extent=(0, im_ny, 0, im_nx),
                            origin="lower", cmap=plt.get_cmap('gray_r'),
                            interpolation="none", vmin=np.log10(av),
                            vmax=np.log10(3*av))  # , interpolation="lanczos")
            plt.scatter(targ_pix[0], targ_pix[1], marker="x", s=10, c="red")
            plt.colorbar(fg)
            filename = os.path.basename(fn)
            plt.savefig(os.path.join(
                plot_dir, filename.replace("."+filename.split(".")[-1],
                                           "_{:}.png".format(fi))), dpi=200)
            plt.clf()


def get_par(myfits, par):
    """
    Returns the header parameter from the fits.
    """
    
    try:
        hdu = pf.open(myfits, ignore_missing_end=True)
        header = hdu[0].header
        if str.upper(par) in header.keys():
            return header[str.upper(par)]
        else:
            return None
    except IOError:
        return None       


def update_par(myfits, par, value):
    """
    Updates the fits files with the new parameter.
    """
    hdu = pf.open(myfits, ignore_missing_end=True)
    header = hdu[0].header
    header.set(par, value)
    hdu.writeto(myfits, overwrite=True)


def update_pars(myfits, pardic):
    """
    Updates the fits files with the new parameter.
    """
    hdu = pf.open(myfits, ignore_missing_end=True)
    header = hdu[0].header

    # TODO: do we really have to writeto for each keyword?
    for key in pardic:
        header.set(key, pardic[key])
        hdu.writeto(myfits, overwrite=True)


def has_par(myfits, par):
    """
    Updates the fits files with the new parameter.
    """
    hdu = pf.open(myfits, ignore_missing_end=True)
    header = hdu[0].header
    return par in header.keys()
