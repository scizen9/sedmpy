"""Bias, Overscan subtraction and Gain correction functions

Functions:
    * :func:`add_prefix` adds bias "b_" prefix to input name
    * :func:`full_frame` calculates smoothed overscan vector
    * :func:`remove` subtracts overscan and converts to electrons

Note:
    This is used as a python script as follows::

        Debias.py file [file ...]]

        positional arguments:
         file           An image file to be processed

"""
import time
import astropy.io.fits as pf
import numpy as np
from scipy.ndimage import filters
import sedmpy_version

drp_ver = sedmpy_version.__version__


def calculate_oscan(dat, header):
    """Subtract overscan from image"""
    # Get over/pre-scan regions
    try:
        # These should be in the Andor camera image headers
        ps_x0 = header['PSCANX0']
        ps_x1 = header['PSCANX1']
        os_x0 = header['OSCANX0']
        os_x1 = header['OSCANX1']
        andor = True
    except KeyError:
        # These are for the PIXIS camera images
        ps_x0 = 0
        ps_x1 = 0
        os_x0 = 2045
        os_x1 = 2048
        andor = False

    # Pre-scan value
    if ps_x1 > ps_x0:
        pscan = np.nanmedian(np.nanmedian(dat[:, :ps_x1], axis=1))
        header['PSCANVAL'] = (pscan, 'Pre-scan value')
    else:
        pscan = 0.
    # Over-scan value
    if os_x1 > os_x0:
        oscan = np.nanmedian(np.nanmedian(dat[:, os_x0:], axis=1))
        header['OSCANVAL'] = (oscan, 'Over-scan value')
    else:
        oscan = 0.

    if andor:
        if pscan < oscan:
            scan_value = pscan
        else:
            scan_value = oscan

        header['SCANVAL'] = (scan_value, 'Scan value subtracted')
        header['OSCANSUB'] = (True, 'Over-scan subtracted?')

    else:
        scan_value = 0.

    return scan_value, ps_x1, os_x0-1


def pixoscan(dat):
    """Calculate smoothed overscan vector

    Args:
        dat (numpy array): image frame

    Returns:
        numpy vector: median smoothed overscan image

    """

    oscan = np.nanmedian(dat[:, 2045:], axis=1)
    oscan = oscan.astype(np.float32)
    smooth = filters.median_filter(oscan, size=50)
    oscan_val = np.nanmedian(smooth)

    return np.tile(smooth, (2048, 1)), oscan_val


def remove(fits_obj):
    """Return the overscan-subtracted and gain-corrected version of input

    Args:
        fits_obj (fits object): fits science image to be oscan-subtracted

    Returns:
        data array: oscan-subtracted and gain corrected image

    """

    dat = fits_obj[0].data

    # Gain for electron conversion
    try:
        gain = fits_obj[0].header['GAIN']
    except KeyError:
        gain = 1.8  # Guess the gain

    # get overscan if correctly sized
    if dat.shape == (2048, 2048):
        oscan_img, val = pixoscan(dat)
        osub_img = (dat - oscan_img.T) * gain
    else:
        val = 0.
        osub_img = dat * gain

    # update header
    fits_obj[0].header['OSCANVAL'] = (val, 'Median Overscan Value')
    fits_obj[0].header['OSCANSUB'] = (True, 'Overscan subtracted?')

    return osub_img


def add_prefix(fname):
    """Adds bias prefix to input name

    Args:
        fname (str): file to be de-biased

    Returns:
        str: input filename with "b_" prepended

    """

    sp = fname.split("/")
    sp[-1] = 'b_' + sp[-1]

    return "/".join(sp)


if __name__ == '__main__':
    import sys

    files = sys.argv[1:]
    
    for ifile in files:
        try:
            if ifile[-5:] != '.fits':
                continue
        except:
            continue

        FF = pf.open(ifile)
        adcspeed = FF[0].header['ADCSPEED']

        bfname = "bias%1.1f.fits" % adcspeed
        try:
            mbias = pf.open(bfname)
        except FileNotFoundError:
            print("Master bias not found: %s" % bfname)
            continue

        # Are we an image from the Andor camera?
        do_andor = 'PSCANX0' in FF[0].header

        if do_andor:

            # Subtract overscan
            img = FF[0].data
            img = img.astype(np.float32)
            osval, x0, x1 = calculate_oscan(img, FF[0].header)
            if osval > 0.:
                img -= osval

            # Gain value
            try:
                and_gain = FF[0].header['GAIN']
            except KeyError:
                and_gain = 0.9

            # Bias frame subtraction and gain correction
            FF[0].data = (img - mbias[0].data) * and_gain

            # Trim overscan
            FF[0].data = FF[0].data[:, x0:x1]
            FF[0].header['OSCANTRM'] = (True, 'Overscan regions trimmed?')

        else:
            # Bias frame subtraction
            FF[0].data = FF[0].data - mbias[0].data

            # Overscan subtraction
            FF[0].data = remove(FF)

        outname = add_prefix(ifile)
        FF[0].header['BIASSUB'] = (True, 'Bias subtracted?')
        FF[0].header['MBIASFN'] = (bfname, 'Master Bias file used')
        try: 
            GAIN = FF[0].header['GAIN']
            FF[0].header['GAIN'] = (1.0, 'GAIN Adjusted (was %s)' % GAIN)
        except KeyError:
            GAIN = 1.8  # Guess the gain
            FF[0].header['GAIN'] = (1.0,
                                    'GAIN Adjusted (was guessed %s)' % GAIN)
        FF[0].header['BUNIT'] = 'electron'
        FF[0].header.add_history('drpifu.Debias run on %s' %
                                 time.strftime("%c"))
        FF[0].header['DRPVER'] = drp_ver
        FF.writeto(outname)
