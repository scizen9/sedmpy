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
import sedmpy_version
from Imcombine import subtract_oscan

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

        # Subtract overscan
        img = FF[0].data
        img = img.astype(np.float32)
        osval = subtract_oscan(img, FF[0].header)
        if osval > 0.:
            img -= osval

        # Gain value
        try:
            gain = FF[0].header['GAIN']
        except KeyError:
            gain = 1.8

        # Bias frame subtraction and gain correction
        FF[0].data = (img - mbias[0].data) * gain

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
