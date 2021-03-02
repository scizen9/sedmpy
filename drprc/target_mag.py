# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:36:52 2021

@author: neill
"""
import math
from astropy.io import fits as pf
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry, CircularAperture, CircularAnnulus, centroid_sources, centroid_2dg
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


def get_target_mag(imfile, aper_r=10., sky_aper_in=15., sky_aper_out=25.,
                   zeropoint=None, verbose=False):
    """get target mag
    Inputs
    imfile: (str) - filename
    aper_r: (float) - aperture radius in pixels
    sky_aper_in: (float) - sky annulus inner radius
    sky_aper_out: (float) - sky annulus outer radius
    zeropoint: (float) - magnitude zeropoint
    verbose: (bool) - extra output?
    """
    # initialize outputs
    targ_mag = None
    targ_magerr = None
    std_mag = None
    std_zeropoint = None

    # default radii
    ap_r = aper_r
    sky_in = sky_aper_in
    sky_out = sky_aper_out

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
        targ_obj = hdu[0].header['OBJECT']
        targ_expt = hdu[0].header['EXPTIME']
        targ_gain = hdu[0].header['GAIN']
        targ_rnoise = hdu[0].header['RDNOISE']

        # Get target
        if 'STD' in targ_obj:
            targ_name = targ_obj.split()[0].split('STD-')[-1]
        elif 'ACQ' in targ_obj:
            targ_name = targ_obj.split()[0].split('ACQ-')[-1]
        else:
            targ_name = targ_obj.split()[0]

        if verbose:
            print("%s | x: %.2f, y: %.2f, expt: %.2f, air: %.3f, gain: %.3f, rn: %.1f" %
                  (targ_name, targ_x, targ_y, targ_expt, targ_air, targ_gain, targ_rnoise))

        # Are we a standard star?
        if targ_name.lower() in stds:
            std_name = targ_name.lower()
            if std_name in stds:
                std_mag = stds[std_name]

            # Adjust apertures
            ap_r = 40.
            sky_in = 50.
            sky_out = 65.

            # Centroid standard stars
            cent_x, cent_y = centroid_sources(image, targ_x, targ_y, box_size=int(ap_r),
                                              centroid_func=centroid_2dg)
            targ_x = cent_x[0]
            targ_y = cent_y[0]
            if verbose:
                print("Position updated | x: %.2f, y: %.2f" % (targ_x, targ_y))

        # Get error image
        erimg = np.sqrt(image * targ_gain + targ_rnoise ** 2)

        # Convert DN to electrons
        image *= targ_gain

        # create source position
        positions = [(targ_x, targ_y)]

        # Set up apertures
        apertures = CircularAperture(positions, r=ap_r)
        annulus_aperture = CircularAnnulus(positions, r_in=sky_in, r_out=sky_out)
        annulus_masks = annulus_aperture.to_mask(method='center')

        # Estimate background
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(image)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)

        # Calculate net aperture sums
        phot = aperture_photometry(image, apertures, error=erimg)
        phot['annulus_median'] = bkg_median
        phot['aper_bkg'] = bkg_median * apertures.area

        # Subtract background
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']

        # Extract data
        ap_sum_bkgsub = phot['aper_sum_bkgsub'].data[0]
        ap_sum_err = phot['aperture_sum_err'].data[0]

        # Get e-/s
        ap_sum_bkgsub_per_sec = ap_sum_bkgsub / targ_expt

        # Is the net sum positive?
        if ap_sum_bkgsub_per_sec > 0:

            # Calculate airmass-corrected instrumental magnitude
            air_cor = targ_air * sdss_r_band_extinction_coeff
            int_mag = -2.5 * math.log10(ap_sum_bkgsub_per_sec) - air_cor

            # Generate a new zeropoint, if we are a standard star
            if std_mag is not None:
                std_zeropoint = std_mag - int_mag
                zp = std_zeropoint
                if verbose:
                    print("std_mag: %2f, std_zp: %.2f" % (std_mag, zp))
            else:
                zp = 25.0

            # Use input zeropoint, if set
            if zeropoint is not None:
                zp = zeropoint

            # Calculated calibrated target magnitude
            targ_mag = int_mag + zp

            # Calculated magnitude error
            targ_magerr = 1.0857362 * ap_sum_err / ap_sum_bkgsub

            if verbose:
                print("e-/s: %.1f, err: %.1f, excor: %.2f, imag: %.2f, zp: %.2f, targ_mag: %.2f +- %.2f" %
                      (ap_sum_bkgsub_per_sec, ap_sum_err, air_cor, int_mag, zp, targ_mag, targ_magerr))

            # Update header
            hdu[0].header['RQMAG'] = (targ_mag, 'r-band quick magnitude')
            hdu[0].header['RQMGERR'] = (targ_magerr, 'r-band quick mag error')
            hdu[0].header['RQZP'] = (zp, 'r-band zeropoint used')
            if std_mag is not None:
                hdu[0].header['RQSTDZP'] = (std_zeropoint, 'calculated r-band zeropoint')
            hdu.writeto(imfile, overwrite=True)

        else:
            print("Warning: negative net aperture sum, no mag calculated!")

    hdu.close()

    return targ_mag, targ_magerr, std_zeropoint