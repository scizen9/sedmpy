#! /usr/bin/env python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    import argparse
    import numpy as np
    import pylab as pl
    import astropy.io.fits as pf
    from scipy.interpolate import interp1d

    # setup arguments parser
    parser = argparse.ArgumentParser(
        description="""Compare two efficiencies.""",
        formatter_class=argparse.RawTextHelpFormatter)
    # setup arguments
    parser.add_argument('-a', '--andor_file', type=str, default=None,
                        help='Andor EA file')
    parser.add_argument('-p', '--pixis_file', type=str, default=None,
                        help='Pixis EA file')
    parser.add_argument('--area', type=float, default=18000.,
                        help='telescope area in cm (def: 18,000.)')
    parser.add_argument('--refl', type=float, default=0.82,
                        help='total reflectance (def: 0.82)')
    args = parser.parse_args()

    area = args.area
    refl = args.refl
    f1 = args.file1
    f2 = args.file2

    # Get first spec file
    # Read in spectrum
    ff = pf.open(f1)
    # Original header
    ohdr = ff[0].header
    # EA for #1
    ea1 = ff[0].data
    # Close spectrum
    ff.close()
    # Calculate wavelengths
    lam0 = ohdr['CRVAL1']
    dlam = ohdr['CDELT1']
    lbda1 = (1 + np.arange(len(ea1))) * dlam + lam0
    # Object name
    obname = f1.split('_')[-2].split('-')[-1]
    # Calculate efficiency spectrum
    eff1 = 100. * ea1/(area * refl)

    # Get second spec file
    # Read in spectrum
    ff = pf.open(f2)
    # Original header
    ohdr = ff[0].header
    # EA for #1
    ea2 = ff[0].data
    # Close spectrum
    ff.close()
    # Calculate wavelengths
    lam0 = ohdr['CRVAL1']
    dlam = ohdr['CDELT1']
    lbda2 = (1 + np.arange(len(ea2))) * dlam + lam0
    # Calculate efficiency spectrum
    eff2 = 100. * ea2/(area * refl)


    # Plot efficiency spectrum
    pl.figure(1)
    pl.plot(lbda1, eff1, label="Andor")
    pl.plot(lbda2, eff2, label="Pixis")
    pl.xlabel('Wave(A)')
    pl.ylabel('Eff(%)')
    pl.legend()
    pl.title("%s, Area = %.0f cm^2, Refl = %.0f %%" % (obname, area,
                                                       refl*100.))
    pl.grid(True)
    pl.ioff()
    # Save efficiency plot to file
    plotname = "compare_eff.png"
    pl.savefig(plotname)

    # get the ratio
    eff2int = interp1d(lbda2, eff2, kind='cubic', bounds_error=False,
                       fill_value='extrapolate')
    eff2rs = eff2int(lbda1)
    effrat = (eff1 / eff2rs - 1.0) * 100.

    pl.figure(2)
    pl.plot(lbda1, effrat)
    pl.xlabel('Wave(A)')
    pl.ylabel('% Improvement')
    pl.title("Andor vs. Pixis: %s" % obname)
    pl.grid(True)
    plotname = "eff_ratio.png"
    pl.savefig(plotname)


