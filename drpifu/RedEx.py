#! /usr/bin/env python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    import sys
    import os
    import glob
    import subprocess
    import argparse
    import logging
    import datetime

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
    parser.add_argument('new_x', type=str, default=None, nargs='?',
                        help='new x position (spaxels)')
    parser.add_argument('new_y', type=str, default=None, nargs='?',
                        help='new y position (spaxels)')
    parser.add_argument('--slack', action='store_true', default=False,
                        help='update slack pysedm-report channel')
    args = parser.parse_args()

    if not args.obs_id:
        logging.info("Usage - redex <obs_id> [<x> <y>] [--slack]")
    else:
        # Get tag id
        now = datetime.datetime.now()
        tagstr = "redo%02d%02d%02.0f" % (now.hour, now.minute, now.second)
        # Check inputs and environment
        ob_id = args.obs_id
        dd = os.getcwd().split('/')[-1]
        rd = '/'.join(os.getcwd().split('/')[:-1])
        reddir = os.environ['SEDMREDUXPATH']
        if rd not in reddir:
            logging.error("check SEDMREDUXPATH env var")
            sys.exit(1)

        logging.info("Re-extracting observation %s in %s" % (ob_id, dd))
        if args.new_x and args.new_y:
            xs = args.new_x
            ys = args.new_y
            pars = ["extract_star.py", dd, "--auto", ob_id, "--autobins", "6",
                    "--centroid", xs, ys, "--tag", tagstr]
            logging.info("Running " + " ".join(pars))
            res = subprocess.run(pars)
            if res.returncode != 0:
                logging.error("Extraction failed.")
                sys.exit(1)
        else:
            pars = ["extract_star.py", dd, "--auto", ob_id, "--autobins", "6",
                    "--display", "--tag", tagstr]
            logging.info("Running " + " ".join(pars))
            res = subprocess.run(pars)
            if res.returncode != 0:
                logging.error("Extraction failed.")
                sys.exit(1)

        # Object name
        obname = glob.glob("spec_*_%s*.fits" %
                           ob_id)[0].split('_')[-1].split('.')[0]
        # Re-classify
        flist = glob.glob("spec_*_%s_%s_*" % (ob_id, obname))
        for f in flist:
            logging.info("removing %s" % f)
            os.remove(f)
        logging.info("make classify")
        res = subprocess.run(["make", "classify"])
        if res.returncode != 0:
            home = os.environ['HOME']
            cmd = [os.path.join(home, "spy"),
                   os.path.join(home, "sedmpy/drpifu/Classify.py"), "-d",
                   os.getcwd()]
            logging.info(" ".join(cmd))
            res = subprocess.run(cmd)
            if res.returncode != 0:
                logging.error("SNID classification failed!")
                sys.exit(1)
        # Re-verify
        pars = ["verify", dd, "--contains", tagstr]
        logging.info("Running " + " ".join(pars))
        res = subprocess.run(pars)
        if res.returncode != 0:
            logging.error("Verify failed!")
            sys.exit(1)
        # Re-report
        if args.slack:
            pars = ["pysedm_report.py", dd, "--contains", tagstr, "--slack"]
        else:
            logging.info("Be sure to update slack manually")
            pars = ["pysedm_report.py", dd, "--contains", tagstr]
        logging.info("Running " + " ".join(pars))
        res = subprocess.run(pars)
        if res.returncode != 0:
            logging.error("pysedm_report.py failed!")
            sys.exit(1)
        # Prepare for upload
        logging.info("be sure to run make ztfupload when you are done.")
