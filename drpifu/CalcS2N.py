# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 11:22:00 2021

@author: neill
"""

import argparse
import astropy.io.fits as pf
import numpy as np


def calc_s2n(spec_file=None, start_wave=4000., end_wave=8000.):
    """
    Calculates S/N from the input file.
    """

    s2nmed = None
    if spec_file is not None:
        # Update fits file
        ff = pf.open(spec_file.replace('.txt', '.fits'), mode='update')
        # get S/N ratio
        flux = ff[0].data
        noise = np.sqrt(ff[1].data)
        s2nspec = flux / noise
        # Get wavelength scale
        w0 = ff[0].header['CRVAL1']
        dw = ff[0].header['CDELT1']
        npx = ff[0].header['NAXIS1']
        w1 = w0 + dw * float(npx - 1)
        wave = np.arange(start=w0, end=w1, step=dw)
        # Get wavelength range
        if start_wave < w0:
            sw0 = w0
        else:
            sw0 = start_wave
        if end_wave > w1:
            sw1 = w1
        else:
            sw1 = end_wave
        sample_range = [i for i, w in enumerate(wave) if sw0 < w < sw1]
        s2nmed = np.nanmedian(s2nspec[sample_range])

        # Update ascii file
        with open(spec_file, "r") as specIn:
            spec_lines = specIn.readlines()

            comment_lines_count = 0
            for line in spec_lines:
                if line.split()[0] == '#':
                    comment_lines_count += 1
                else:
                    break

            ff[0].header['S2NMED'] = (s2nmed, 'Median S/N')
            ff[0].header['S2NWL0'] = (sw0, 'Start wave for S/N Calc')
            ff[0].header['S2NWHI'] = (sw1, 'End wave for S/N Calc')

            spec_lines.insert(comment_lines_count,
                              "# SN2_WL_END: " + str(sw1) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SN2_WL_START: " + str(sw0) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SN2_MEDIAN: " + str(s2nmed) + "\n")

    return s2nmed
    # END calc_s2n


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     """

        Calculates S/N for SEDM spectrum and records it in fits and ascii file.

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specfile', type=str, dest="specfile",
                        help='spec_*.txt ascii spectrum file.',
                        default=None)
    parser.add_argument('-s', '--start_wave', type=float, dest='s0',
                        help='starting wavelength in Angstroms',
                        default=4000.)
    parser.add_argument('-e', '--end_wave', type=float, dest='s1',
                        help='ending wavelength in Angstroms',
                        default=8000.)

    args = parser.parse_args()

    # Calculate S/N ratio from extracted spectra and record results in specfile
    print("Calcuating S/N from and recording results in %s" % args.specfile)
    s2n = calc_s2n(spec_file=args.specfile, start_wave=args.s0,
                   end_wave=args.s1)
    print("Median S/N is %f" % s2n)
