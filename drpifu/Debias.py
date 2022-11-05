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
import numpy as np
import astropy.io.fits as pf
import scipy.ndimage.filters as FI
import sedmpy_version

drp_ver = sedmpy_version.__version__


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


def full_frame(dat):
    """Calculate smoothed overscan vector

    Args:
        dat (numpy array): image frame

    Returns:
        numpy vector: median smoothed overscan vector

    """

    bias = np.nanmedian(dat[:, 2045:], axis=1)
    bias = bias.astype(np.float32)
    smooth = FI.median_filter(bias, size=50)

    return np.tile(smooth, (2048, 1))


def remove(fits_obj):
    """Return the overscan-subtracted and gain-corrected version of the fits object

    Args:
        fits_obj (fits object): fits science image to be overscan-subtracted

    Returns:
        data array: bias-subtracted and gain corrected image

    """

    dat = fits_obj[0].data

    # Gain for electron conversion
    try:
        GAIN = fits_obj[0].header['GAIN']
    except:
        GAIN = 1.8  # Guess the gain

    # get overscan if correctly sized
    print(dat.shape)
    if dat.shape == (2048, 2048):
        bias_img = full_frame(dat)
        ret = (dat - bias_img.T) * GAIN
    else:
        ret = dat * GAIN
    # return overscan subtracted image
    return ret


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
        bias = pf.open(bfname)

        # Bias frame subtraction
        FF[0].data = FF[0].data - bias[0].data

        # Overscan subtraction
        FF[0].data = remove(FF)

        outname = add_prefix(ifile)
        FF[0].header['BIASSUB'] = ('Subtracted',
                                   'Ovrscn + bias handled by Debias.py')
        FF[0].header['BIASSUB2'] = (bfname, 'Bias file used')
        try: 
            GAIN = FF[0].header['GAIN']
            FF[0].header['GAIN'] = (1.0, 'GAIN Adjusted (was %s)' % GAIN)
        except: 
            GAIN = 1.8  # Guess the gain
            FF[0].header['GAIN'] = (1.0,
                                    'GAIN Adjusted (was guessed %s)' % GAIN)
        FF[0].header['BUNIT'] = 'electron'
        FF[0].header.add_history('drpifu.Debias run on %s' % time.strftime("%c"))
        FF[0].header['DRPVER'] = drp_ver
        FF.writeto(outname)

