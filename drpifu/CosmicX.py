"""

.. module:: CosmicX

:synopsis: Remove cosmic rays from an astronomical image.

kludged from PyCosmic to run _lacosmicx instead, which is faster

"""
import time
import argparse
import numpy as np
import astropy.io.fits as pf

try:
    import _lacosmicx
except ImportError:
    print("Please install lacosmicx from github.com/cmccully/lacosmicx.")
    quit()

__version__ = "0.1"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Program to detect cosmics in single exposure CCD frames. Important: The image and the readout noise are assumed to be in units of electrons.
    The image also needs to be BIAS subtracted! The gain can be entered to convert the image from ADUs to electros, when this is down already set gain=1.0 as the default. A bad pixel mask of cosmics and a cleaned image will be provided by the routine as an output.""", formatter_class=argparse.ArgumentDefaultsHelpFormatter,prog='CosmicX.py')

    parser.add_argument("raw", type=str, help="""File name of the CCD raw frame FITS file from which to detect and reject cosmic ray hits. This frame is expected to be BIAS subtracted. If it is not provided in units of electrons please provide the appropriate gain factor to convert from ADUs to electrons with the --gain parameter.
    """)
    parser.add_argument("clean", type=str, help="""File name of the cosmics cleaned frame to be stored as a FITS file.
    """)
    parser.add_argument("mask", type=str, help="""File name of the cosmics mask  frame to be stored as a FITS file.
    """)
    parser.add_argument("--minexptime", type=float, default=60.0, help="""Minimum exposure time to run: below this don't run.""")
    parser.add_argument("--pssl", type=float, default=0.0, help="""Previously subtracted sky level.""")
    parser.add_argument("--gain", type=float, default=2.2, help="""CCD gain as float number when the bias-subtracted image was not yet converted to electrons.""")
    parser.add_argument("--readnoise", type=float, default=10.0, help="""CCD read-out noise in electrons.""")
    parser.add_argument("--sigclip", type=float, default=5.0, help="""Threshold value for the significance level of cosmics in units of the expected noise for the pixels.""")
    parser.add_argument("--sigfrac", type=float, default=0.3, help="""Fractional detection limit for neighbouring pixels.""")
    parser.add_argument("--objlim", type=float, default=5.0, help="""Minimum contrast between laplacian image and fine structure image. Use 5.0 if your image is undersampled, HST, ...""")
    parser.add_argument("--psffwhm", type=float, default=2.5, help="""Full Width Half Maximum of the PSF to use to generate the kernel.""")
    parser.add_argument("--fsmode", choices=['median', 'convolve'], default='median', help="""Method to build the fine structure image:
            'median': Use the median filter in the standard LA Cosmic algorithm
            'convolve': Convolve the image with the psf kernel to calculate the
            fine structure image.""")
    parser.add_argument("--psfmodel", choices=['gauss', 'gaussx', 'gaussy', 'moffat'], default='gauss', help="""Model to use to generate the psf kernel if --fsmode 'convolve' and
            psfk is None. The current choices are Gaussian and Moffat profiles.
            'gauss' and 'moffat' produce circular PSF kernels. The 'gaussx' and
            'gaussy' produce Gaussian kernels in the x and y directions
            respectively.""")
    parser.add_argument("--satlevel", type=float, default=50000.0, help="""If we find agglomerations of pixels above this level, we consider it to be a saturated star and do not try to correct and pixels around it. A negative satlevel skips this feature.""")
    parser.add_argument("--verbose", action="store_true", default=False, help="""Flag to print some progress information on the screen.""")
    parser.add_argument("--sepmed", action="store_true", default=False, help="""Flag to use separable median (faster).""")
    parser.add_argument("--niter", type=int, default=4, help="""Number of iteration to be performed by the algorithms. Usually 5-6 iterations are needed to converge to a stable solution.""")

    args = parser.parse_args()

    f = pf.open(args.raw)
    header = f[0].header
    array = np.array(f[0].data, dtype=np.float32)
    f.close()

    if header['EXPTIME'] >= args.minexptime:
        mask, clean = _lacosmicx.lacosmicx(array, gain=args.gain, readnoise=args.readnoise, psffwhm=args.psffwhm, sigclip=args.sigclip, sigfrac=args.sigfrac, objlim=args.objlim, fsmode=args.fsmode, psfmodel=args.psfmodel, verbose=args.verbose, sepmed=args.sepmed)

        header['history'] = "LA CosmicX: cleaned cosmic rays"
        header['history'] = "LA CosmicX params: sigclip=%5.2f sigfrac=%5.2f objlim=%5.2f" % (args.sigclip, args.sigfrac, args.objlim)
        header['history'] = "LA CosmicX params: fsmode=%s psfmodel=%s psffwhm=%5.2f" % (args.fsmode, args.psfmodel, args.psffwhm)
        header['history'] = "LA CosmicX params: sepmed=%s minexptime=%f" % (args.sepmed, args.minexptime)
        header['history'] = "LA CosmicX run on %s" % time.strftime("%c")

        pf.writeto(args.clean, clean, header)
        mask = np.cast["uint8"](mask)
        pf.writeto(args.mask, mask, header)
    else:
        header['history'] = "LA CosmicX: exptime < minexptime=%.1f" % args.minexptime
        pf.writeto(args.clean, array, header)
