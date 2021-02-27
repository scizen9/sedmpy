# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:36:52 2021

@author: neill
"""
import math
from astropy.io import fits as pf
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
import numpy as np


def get_target_mag(imfile):
    """get target mag"""
    # initialize outputs
    targ_mag = None
    targ_magerr = None

    # Open reduced image file
    hdu = pf.open(imfile)
    on_target = hdu[0].header['ONTARGET']

    # Are we on target?
    if on_target:
        image = hdu[0].data.astype(float)
        targ_x = hdu[0].header['TARGXPX']
        targ_y = hdu[0].header['TARGYPX']

        # create source position
        positions = [(targ_x, targ_y)]

        # Do photometry
        apertures = CircularAperture(positions, r=4.)
        annulus_aperture = CircularAnnulus(positions, r_in=10, r_out=15)
        annulus_masks = annulus_aperture.to_mask(method='center')

        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(image)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        phot = aperture_photometry(image, apertures)
        phot['annulus_median'] = bkg_median
        phot['aper_bkg'] = bkg_median * apertures.area
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        ap_sum_bkgsub = phot['aper_sum_bkgsub'].data[0]
        if ap_sum_bkgsub > 0:
            targ_mag = 25.0 - math.log10(ap_sum_bkgsub)
            targ_magerr = targ_mag / 100.

    hdu.close()

    return targ_mag, targ_magerr
