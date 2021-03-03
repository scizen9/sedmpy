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

extinction_coeff = {'u': 0.594, 'g': 0.208, 'r': 0.120, 'i': 0.059}

# u = 1.28 * (U-B) + 1.14 + g
# g = V + 0.64 * (B-V) - 0.13
# r = V - 0.46 * (B-V) + 0.11
# i = -0.98 * (R-I) + 0.22 + r

stds = {
    'sa95-42':    {'u': 15.073, 'g': 15.340, 'r': 15.793, 'i': 16.138},
    'hz4':        {'u': 14.707, 'g': 14.431, 'r': 14.576, 'i': 14.855},
    'lb227':      {'u': 15.449, 'g': 15.228, 'r': 15.408, 'i': 15.734},
    'hz2':        {'u': 13.697, 'g': 13.689, 'r': 14.028, 'i': 14.357},
    'bd+75d325':  {'u':  8.793, 'g':  9.204, 'r':  9.811, 'i': 10.214},
    'feige34':    {'u': 10.403, 'g': 10.831, 'r': 11.449, 'i': 11.810},
    'hz44':       {'u': 10.966, 'g': 11.357, 'r': 11.917, 'i': 12.314},
    'bd+33d2642': {'u': 10.636, 'g': 10.592, 'r': 11.014, 'i': 11.308},
    'bd+28d4211': {'u':  9.706, 'g': 10.161, 'r': 10.776, 'i': 11.168},
    'bd+25d4655': {'u':  9.070, 'g':  9.364, 'r':  9.923, 'i': 10.3},   # i-band is not precise
    'feige110':   {'u': 11.153, 'g': 11.507, 'r': 12.082, 'i': 12.478}
    }


def get_target_mag(imfile, zeropoint=None, verbose=False):
    """get target mag
    Inputs
    imfile: (str) - filename
    zeropoint: (dict) - magnitude zeropoint dictionary with entries for ugri filters
    verbose: (bool) - extra output?
    """
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
        targ_obj = hdu[0].header['OBJECT']
        targ_expt = hdu[0].header['EXPTIME']
        targ_gain = hdu[0].header['GAIN']
        targ_rnoise = hdu[0].header['RDNOISE']
        targ_fwhm = hdu[0].header['FWHMPIX']

        # Get filter
        targ_filter = targ_obj.split()[-1]

        # Get target
        if 'STD' in targ_obj:
            targ_name = targ_obj.split()[0].split('STD-')[-1]
        elif 'ACQ' in targ_obj:
            targ_name = targ_obj.split()[0].split('ACQ-')[-1]
        else:
            targ_name = targ_obj.split()[0]

        if verbose:
            print("\n%s" % imfile)
            print("%s %s | x: %.2f, y: %.2f, fwhm: %.2f, expt: %.2f, air: %.3f, gain: %.3f, rn: %.1f" %
                  (targ_name, targ_filter, targ_x, targ_y, targ_fwhm, targ_expt, targ_air, targ_gain, targ_rnoise))

        # default radii
        if targ_fwhm > 0:
            ap_r = targ_fwhm
            sky_in = targ_fwhm + targ_fwhm * 0.2
            sky_out = targ_fwhm * 2.
        else:
            ap_r = 10.
            sky_in = 12.
            sky_out = 20.

        # Are we a standard star?
        if targ_name.lower() in stds:
            std_name = targ_name.lower()
            std_dict = stds[std_name]
            if targ_filter in std_dict:
                std_mag = std_dict[targ_filter]

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

        # Calculate aperture sums
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
            if targ_filter in extinction_coeff:
                air_cor = targ_air * extinction_coeff[targ_filter]
            else:
                air_cor = 0.
            int_mag = -2.5 * math.log10(ap_sum_bkgsub_per_sec) - air_cor

            # Generate a new zeropoint, if we are a standard star
            if std_mag is not None:
                std_zeropoint = std_mag - int_mag
                zp = std_zeropoint
                if verbose:
                    print("std_mag: %.2f, std_zp: %.2f" % (std_mag, zp))
            else:
                zp = 25.0

            # Use input zeropoint, if set
            if zeropoint is not None:
                if targ_filter in zeropoint:
                    if zeropoint[targ_filter] is not None:
                        zp = zeropoint[targ_filter]

            # Calculated calibrated target magnitude
            targ_mag = int_mag + zp

            # Calculated magnitude error
            targ_magerr = 1.0857362 * ap_sum_err / ap_sum_bkgsub

            if verbose:
                print("Ap results | e-/s: %.1f, err: %.1f, excor: %.2f, imag: %.2f, zp: %.2f, targ_mag: %.2f +- %.2f" %
                      (ap_sum_bkgsub_per_sec, ap_sum_err, air_cor, int_mag, zp, targ_mag, targ_magerr))

            # Update header
            key_stub = 'Q' + targ_filter.upper()
            hdu[0].header[key_stub + 'MAG'] = (targ_mag, targ_filter + '-band quick magnitude')
            hdu[0].header[key_stub + 'MGERR'] = (targ_magerr, targ_filter + '-band quick mag error')
            hdu[0].header[key_stub + 'ZP'] = (zp, targ_filter + '-band zeropoint used')
            if std_mag is not None:
                hdu[0].header[key_stub + 'STDZP'] = (std_zeropoint, targ_filter + '-band calculated zeropoint')
            hdu.writeto(imfile, overwrite=True)

        else:
            print("Warning: negative net aperture sum, no mag calculated!")

    hdu.close()

    return targ_mag, targ_magerr, std_zeropoint
