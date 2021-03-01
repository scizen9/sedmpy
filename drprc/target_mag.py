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

sdss_r_band_extinction_coeff = 0.12

stds = {
    'sa95-42': 15.793,
    'hz4': 14.576,
    'lb227': 15.408,
    'hz2': 14.028,
    'bd+75d325': 9.811,
    'feige34': 11.449,
    'hz44': 11.917,
    'bd+33d2642': 11.014,
    'bd+28d4211': 10.776,
    'bd+25d4655': 9.923,
    'feige110': 12.082
    }    # r = V - 0.46 * (B-V) + 0.11


def get_target_mag(imfile, aper_r=4.0, sky_aper_in=10., sky_aper_out=15.,
                   zeropoint=None, verbose=False):
    """get target mag"""
    # initialize outputs
    targ_mag = None
    targ_magerr = None
    std_mag = None
    std_zeropoint = None

    # Open reduced image file
    hdu = pf.open(imfile)
    on_target = hdu[0].header['ONTARGET']

    # Are we on target?
    if on_target:
        image = hdu[0].data.astype(float)

        # Get header params
        targ_x = hdu[0].header['TARGXPX']
        targ_y = hdu[0].header['TARGYPX']
        targ_air = hdu[0].header['AIRMASS']
        targ_name = hdu[0].header['OBJECT']
        if verbose:
            print("x: %.2f, y: %.2f, air: %.3f, targ: %s" % (targ_x, targ_y, targ_air, targ_name))

        # Are we a standard star?
        if 'STD' in targ_name:
            std_name = targ_name.split()[0].split('-')[-1].lower()
            if std_name in stds:
                std_mag = stds[std_name]

        # create source position
        positions = [(targ_x, targ_y)]

        # Set up apertures
        apertures = CircularAperture(positions, r=aper_r)
        annulus_aperture = CircularAnnulus(positions, r_in=sky_aper_in, r_out=sky_aper_out)
        annulus_masks = annulus_aperture.to_mask(method='center')

        # Estimate background
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(image)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)

        # Calculate net aperture sum
        phot = aperture_photometry(image, apertures)
        phot['annulus_median'] = bkg_median
        phot['aper_bkg'] = bkg_median * apertures.area
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        ap_sum_bkgsub = phot['aper_sum_bkgsub'].data[0]

        # Is the net sum positive?
        if ap_sum_bkgsub > 0:

            # Calculate airmass-corrected instrumental magnitude
            int_mag = math.log10(ap_sum_bkgsub) - targ_air * sdss_r_band_extinction_coeff

            # Generate a new zeropoint, if we are a standard star
            if std_mag is not None:
                std_zeropoint = std_mag + int_mag
                zp = std_zeropoint
                if verbose:
                    print("std_mag: %2f, std_zp: %.2f" % (std_mag, zp))
            else:
                zp = 25.0

            # Use input zeropoint, if set
            if zeropoint is not None:
                zp = zeropoint

            # Calculated calibrated target magnitude
            targ_mag = zp - int_mag

            if verbose:
                print("cnts: %.1f, imag: %.2f, tmag: %.2f, zp: %.2f" % (ap_sum_bkgsub, int_mag, targ_mag, zp))

            targ_magerr = targ_mag / 100.
        else:
            print("Warning: negative net aperture sum, no mag calculated!")

    hdu.close()

    return targ_mag, targ_magerr, std_zeropoint
