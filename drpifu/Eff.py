#! /usr/bin/env python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    import sys
    import os
    import glob
    import subprocess
    import argparse
    import logging
    import pycalspec
    import numpy as np
    import pylab as pl
    import astropy.io.fits as pf

    logging.basicConfig(
        format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
        datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)

    # setup arguments parser
    parser = argparse.ArgumentParser(
        description="""Do an aperture extraction of a standard star.""",
        formatter_class=argparse.RawTextHelpFormatter)
    # setup arguments
    parser.add_argument('indate', type=str, default=None,
                        help='UT date as YYYYMMDD')
    parser.add_argument('--contains', type=str, default=None,
                        help='unique observation identifier')
    parser.add_argument('--area', type=float, default=18000.,
                        help='telescope area in cm (def: 18,000.)')
    parser.add_argument('--refl', type=float, default=0.82,
                        help='total reflectance (def: 0.82)')
    args = parser.parse_args()

    dd = args.indate
    area = args.area
    refl = args.refl
    ob_id = args.contains

    if not ob_id:
        logging.info("Usage - eff YYYYMMDD --contains <filespec>")
    else:
        # Check inputs and environment
        rd = '/'.join(os.getcwd().split('/')[:-1])
        reddir = os.environ['SEDMREDUXPATH']
        if rd not in reddir:
            logging.error("check SEDMREDUXPATH env var")
            sys.exit(1)
        else:
            # Calculate observed instrumental spectrum (uncalibrated)
            logging.info("Aperture extracting observation %s in %s" % (ob_id, dd))
            pars = ["extract_star.py", dd, "--aperture", ob_id, "--buffer",
                    "10", "--nofluxcal", "--maxpos"]
            logging.info("Running " + " ".join(pars))
            res = subprocess.run(pars)
            if res.returncode != 0:
                logging.error("Extraction failed.")
                sys.exit(1)

            # Get spec file
            specname = glob.glob("spec_aperture*_%s_*.fits" % ob_id)
            if not specname:
                logging.error("No files found for observation id %s" % ob_id)
                sys.exit(1)
            # Read in spectrum
            ff = pf.open(specname[0])
            # Observed spectrum in e-/Angstrom
            spec = ff[0].data
            # Aperture weight
            apwgt = ff[2].data
            # Original header
            ohdr = ff[0].header
            # Close spectrum
            ff.close()
            # Scale spectrum
            spec = spec * apwgt / ohdr['EXPTIME']
            # Calculate wavelengths
            lam0 = ohdr['CRVAL1']
            dlam = ohdr['CDELT1']
            lbda = (1 + np.arange(len(spec))) * dlam + lam0
            # Object name
            obname = specname[0].split('_')[-1].split('.')[0].split('-')[-1]
            # Get reference spectrum
            refspec = pycalspec.std_spectrum(obname)
            # Resample to our wavelengths
            rspec = refspec.reshape(lbda)
            # Convert to photons/s/cm^2/Angstrom
            rspho = 5.03411250e7 * rspec.data * lbda * dlam
            # Effective area (cm^2)
            earea = spec / rspho
            # Write out effective area spectrum
            areaname = specname[0].replace(".fits", "_ea.fits")
            ohdr['BUNIT'] = ('cm^2', 'Brightness units')
            hdu = pf.PrimaryHDU(earea, header=ohdr)
            fo = pf.HDUList([hdu])
            fo.writeto(areaname)
            # Calculate efficiency spectrum
            eff = 100. * earea/(area * refl)
            # Plot efficiency spectrum
            pl.figure(1)
            pl.plot(lbda, eff)
            pl.xlabel('Wave(A)')
            pl.ylabel('Eff(%)')
            pl.title("%s, Area = %.0f cm^2, Refl = %.0f %%" % (obname, area,
                                                               refl*100.))
            pl.grid(True)
            pl.ioff()
            # Save efficiency plot to file
            plotname = specname[0].replace(".fits", "_eff.png")
            pl.savefig(plotname)

