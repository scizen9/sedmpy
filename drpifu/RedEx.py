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
    parser.add_argument('--recover', action='store_true', default=False,
                        help='recover Quality 5 extraction')
    args = parser.parse_args()

    if not args.obs_id:
        logging.info("Usage - redex <obs_id> [<x> <y>] [--slack] [--recover]")
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
        # Recover a quality 5 observation that is already good
        if args.recover:
            logging.info("Recovering a quality 5 spectrum %s in %s" % (ob_id,
                                                                       dd))
            import AutoReduce as ar
            import astropy.io.fits as pf
            # Update quality in fits file
            fits_file = glob.glob(os.path.join(rd, dd,
                                  "spec_auto_robot_lstep1__*_%s_*.fits" %
                                               ob_id))[0]
            ff = pf.open(fits_file, mode='update')
            # Make sure we are a quality 5 observation
            if ff[0].header['QUALITY'] != 5:
                logging.error("not a quality 5 observation!")
                sys.exit(1)
            ff[0].header['QUALITY'] = 1
            ff.flush()
            ff.close()
            # Update database entry
            ar.update_spec(fits_file, update=True)
            # Update quality in text file
            text_file = glob.glob(os.path.join(rd, dd,
                                  "spec_auto_robot_lstep1__*_%s_*.txt" %
                                               ob_id))[0]
            with open(text_file, "r") as textIn:
                spec_lines = textIn.readlines()
                index = [i for i, s in enumerate(spec_lines) if 'QUALITY' in s]
                if len(index) > 0:
                    spec_lines[index[0]] = "# QUALITY: 1\n"
                else:
                    comment_lines_count = 0
                    for line in spec_lines:
                        if line.split()[0] == '#':
                            comment_lines_count += 1
                        else:
                            break
                    spec_lines.insert(comment_lines_count, "# QUALITY: 1\n")
            with open(text_file, "w") as textOut:
                textOut.write("".join(spec_lines))
            # remove *.upl
            upl_file = glob.glob(os.path.join(rd, dd,
                                 "spec_auto_robot_lstep1__*_%s_*.upl" %
                                              ob_id))
            if upl_file:
                os.remove(upl_file[0])
            # now ready to re-upload to marshal
        # Re-extract a recoverable observation
        else:
            logging.info("Re-extracting observation %s in %s" % (ob_id, dd))
            if args.new_x and args.new_y:
                xs = args.new_x
                ys = args.new_y
                pars = ["extract_star.py", dd, "--auto", ob_id, "--autobins",
                        "6", "--centroid", xs, ys, "--tag", tagstr]
                logging.info("Running " + " ".join(pars))
                res = subprocess.run(pars)
                if res.returncode != 0:
                    logging.error("Extraction failed.")
                    sys.exit(1)
            else:
                pars = ["extract_star.py", dd, "--auto", ob_id, "--autobins",
                        "6", "--display", "--tag", tagstr]
                logging.info("Running " + " ".join(pars))
                res = subprocess.run(pars)
                if res.returncode != 0:
                    logging.error("Extraction failed.")
                    sys.exit(1)

            # Get spec file
            specname = glob.glob("spec_*_%s_*.fits" % ob_id)
            if not specname:
                logging.error("No files found for observation id %s" % ob_id)
                sys.exit(1)
            # Object name
            obname = specname[0].split('_')[-1].split('.')[0]
            # Re-classify
            flist = glob.glob("spec_*_%s_%s_*" % (ob_id, obname))
            for f in flist:
                logging.info("removing %s" % f)
                os.remove(f)
            logging.info("make classify")
            ret = subprocess.call(["make", "classify"])
            if ret:
                home = os.environ['HOME']
                cmd = [os.path.join(home, "spy"),
                       os.path.join(home, "sedmpy/drpifu/Classify.py"), "-d",
                       os.getcwd()]
                logging.info(" ".join(cmd))
                ret = subprocess.call(cmd)
                if ret:
                    logging.error("SNID classification failed!")
                    sys.exit(1)
            # Re-verify
            pars = ["verify", dd, "--contains", tagstr]
            logging.info("Running " + " ".join(pars))
            ret = subprocess.call(pars)
            if ret:
                logging.error("Verify failed!")
                sys.exit(1)
            else:
                # Display verify image and then prompt for quality
                verify_file = glob.glob("verify_auto_%s_*.png" % tagstr)[0]
                pars = ["display", verify_file]
                ret = subprocess.call(pars)
                good = input("Re-extraction good? (N/y): ")
                if good[0].capitalize() != "Y":
                    logging.error("Try re-extraction again.")
                    logging.info("Removing *_%s_*" % tagstr)
                    flist = glob.glob("*_%s_*" % tagstr)
                    for f in flist:
                        os.remove(f)
                    sys.exit(1)
            # Re-report
            if args.slack:
                pars = ["pysedm_report.py", dd, "--contains", tagstr, "--slack"]
            else:
                logging.info("Be sure to update slack manually if extraction good")
                logging.info("Use: pysedm_report.py %s --contains %s --justpush" %
                             (dd, tagstr))
                pars = ["pysedm_report.py", dd, "--contains", tagstr]
            logging.info("Running " + " ".join(pars))
            ret = subprocess.call(pars)
            if ret:
                logging.error("pysedm_report.py failed!")
                sys.exit(1)
        # Prepare for upload
        logging.info("be sure to run make ztfupload when you are done.")
