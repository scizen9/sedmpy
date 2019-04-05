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
        description="""Re-do an extraction.""",
        formatter_class=argparse.RawTextHelpFormatter)
    # setup arguments
    parser.add_argument('obs_id', type=str, default=None,
                        help='observation timestamp as HH_MM_SS')
    args = parser.parse_args()

    if not args.obs_id:
        logging.info("Usage - eff <obs_id>")
    else:
        # Check inputs and environment
        ob_id = args.obs_id
        dd = os.getcwd().split('/')[-1]
        rd = '/'.join(os.getcwd().split('/')[:-1])
        reddir = os.environ['SEDMREDUXPATH']
        if rd not in reddir:
            logging.error("check SEDMREDUXPATH env var")
            sys.exit(1)
        # Recover a quality 5 observation that is already good
        else:
            logging.info("Aperture extracting observation %s in %s" % (ob_id, dd))
            pars = ["extract_star.py", dd, "--aperture", ob_id, "--buffer",
                    "10", "--nofluxcal"]
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
            spec = ff[0].data
            lam0 = ff[0].header['CRVAL1']
            dlam = ff[0].header['CDELT1']
            lbda = (1 + np.arange(len(spec))) * dlam + lam0
            # Object name
            obname = specname[0].split('_')[-1].split('.')[0].split('-')[-1]
            # Get reference spectrum
            refspec = pycalspec.std_spectrum(obname)
            rspec = refspec.reshape(lbda)
            ratio = spec / rspec.data
            rspho = 5.03411250e7 * rspec.data * lbda * dlam
            earea = (spec * dlam) / rspho
            eff = earea/18000.
            pl.plot(lbda, eff)
            pl.show()
            # Re-verify
            pars = ["verify", dd, "--contains",
                    "aperture_notfluxcal__crr_b_ifu%s_%s" % (dd, ob_id)]
            # logging.info("Running " + " ".join(pars))
            # ret = subprocess.call(pars)
            # if ret:
            #    logging.error("Verify failed!")
            #    sys.exit(1)
