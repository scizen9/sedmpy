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
    import AutoReduce as ar
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
    parser.add_argument('new_x', type=str, default=None, nargs='?',
                        help='new x position (spaxels)')
    parser.add_argument('new_y', type=str, default=None, nargs='?',
                        help='new y position (spaxels)')
    parser.add_argument('--recover', action='store_true', default=False,
                        help='recover Quality 5 extraction')
    parser.add_argument('--local', action='store_true', default=False,
                        help='Process data locally only (no push to marshal or '
                             'slack or db')
    parser.add_argument('--lstep', type=str, default=None,
                        help='new lstep value (default is 1)')
    parser.add_argument('--oldext', action='store_true', default=False,
                        help='re-extract using extract_star.py instead of'
                             'extractstar.py')
    parser.add_argument('--testing', action="store_true", default=False,
                        help='just testing: extract spec, but no uploads')
    parser.add_argument('--doab', action="store_true", default=False,
                        help='re-extract an A/B pair manually')
    args = parser.parse_args()

    if args.lstep:
        lstep = args.lstep
    else:
        lstep = "1"

    if not args.obs_id:
        logging.info("Usage - redex <obs_id> [<x> <y>] [--recover] [--local]"
                     " [--lstep <n>] [--oldext] [--doab]")
    else:

        # Get tag id
        now = datetime.datetime.now()
        tagstr = "redo%02d%02d%02.0f" % (now.hour, now.minute, now.second)
        abtagstr = "ABredo%02d%02d%02.0f" % (now.hour, now.minute, now.second)
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
            # Update quality in fits file
            if args.doab:
                fspec = "spec_auto_AB*_%s_*.fits" % ob_id
            else:
                fspec = "spec_auto_robot_lstep1__*_%s_*.fits" % ob_id
            flist = glob.glob(os.path.join(rd, dd, fspec))
            if not flist:
                logging.error("No spec file with id %s" % ob_id)
                sys.exit(1)
            else:
                fits_file = flist[0]
            ff = pf.open(fits_file, mode='update')
            # Make sure we are a quality 5 observation
            if ff[0].header['QUALITY'] != 5:
                logging.error("not a quality 5 observation!")
                sys.exit(1)
            ff[0].header['QUALITY'] = 1
            try:
                obj = ff[0].header['OBJECT'].replace(" [A]", "").replace(" ",
                                                                         "-")
            except KeyError:
                logging.warning(
                    "Could not find OBJECT keyword, setting to Test")
                obj = 'Test'
            ff.flush()
            ff.close()
            # Update database entry
            ar.update_spec(fits_file, update_db=True)
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
            # make ready to re-upload to marshal
            upl_file = glob.glob(
                os.path.join(rd, dd, "spec_auto_robot_lstep1__*_%s_*.upl" %
                             ob_id))
            if upl_file:
                os.remove(upl_file[0])
            if args.local:
                pars = ["pysedm_report.py", dd, "--contains", ob_id]
            else:
                # Upload spectrum to marshal
                cmd = ("make", "ztfupload")
                retcode = subprocess.call(cmd)
                if retcode != 0:
                    logging.error("Error uploading spectra to marshal")
                else:
                    # e-mail user that recovery was made
                    ar.email_user(fits_file, dd, obj)
                # Re-report
                pars = ["pysedm_report.py", dd, "--contains", ob_id, "--slack"]
            logging.info("Running " + " ".join(pars))
            ret = subprocess.call(pars)
            if ret:
                logging.error("pysedm_report.py failed!")
        # Re-extract an A/B extraction manually
        elif args.doab:
            logging.info("Re-extracting an A/B observation %s in %s" % (ob_id,
                                                                        dd))
            # get reducer
            def_reducer = os.getenv("SEDM_USER", default='manual')
            reducer = input("Your name (<cr> - %s): " % def_reducer)
            if not reducer:
                reducer = def_reducer
            reducer = reducer.replace(" ", "_")
            # set up script
            pars = ["extractstar_ab.py", dd, "--auto", ob_id,
                    "--autobins", "6", "--display", "--lstep", lstep,
                    "--centroidA", "brightest", "--centroidB", "brightest",
                    "--tag", abtagstr, "--reducer", reducer]
            logging.info("Running " + " ".join(pars))
            res = subprocess.run(pars)
            if res.returncode != 0:
                logging.error("Extraction failed.")
                sys.exit(1)
            # Re-classify
            logging.info("make classify")
            ret = subprocess.call(["make", "classify"])
            if ret:
                logging.error("SNID classification failed!")
                sys.exit(1)
            # Re-verify
            pars = ["verify", dd, "--contains", abtagstr]
            logging.info("Running " + " ".join(pars))
            ret = subprocess.call(pars)
            if ret:
                logging.error("Verify failed!")
                sys.exit(1)
            else:
                # Display verify image and then prompt for quality
                verify_file = glob.glob("verify_auto_%s_*.png" % abtagstr)[0]
                pars = ["display", verify_file]
                ret = subprocess.call(pars)
                good = input("Re-extraction good? (N/y): ")
                if not good:
                    good = 'N'
                if good[0].capitalize() != "Y":
                    logging.error("Try re-extraction again.")
                    logging.info("Removing *_%s_*" % tagstr)
                    flist = glob.glob("*_%s_*" % tagstr)
                    for f in flist:
                        os.remove(f)
                    sys.exit(1)
            # Re-report
            pars = ["pysedm_report.py", dd, "--contains", abtagstr]
            logging.info("Running " + " ".join(pars))
            ret = subprocess.call(pars)
            if ret:
                logging.error("pysedm_report.py failed!")
                sys.exit(1)
        # Re-extract a recoverable observation
        else:
            logging.info("Re-extracting observation %s in %s" % (ob_id, dd))
            # get reducer
            def_reducer = os.getenv("SEDM_USER", default='manual')
            reducer = input("Your name (<cr> - %s): " % def_reducer)
            if not reducer:
                reducer = def_reducer
            reducer = reducer.replace(" ", "_")

            if args.new_x and args.new_y:
                xs = args.new_x
                ys = args.new_y
                if args.oldext:
                    pars = ["extract_star.py", dd, "--auto", ob_id,
                            "--autobins", "6", "--centroid", xs, ys,
                            "--lstep", lstep, "--tag", tagstr,
                            "--reducer", reducer]
                else:
                    pars = ["extractstar.py", dd, "--auto", ob_id,
                            "--autobins", "6", "--centroid", xs, ys,
                            "--lstep", lstep, "--tag", tagstr,
                            "--reducer", reducer]
                logging.info("Running " + " ".join(pars))
                res = subprocess.run(pars)
                if res.returncode != 0:
                    logging.error("Extraction failed.")
                    sys.exit(1)
            else:
                if args.oldext:
                    pars = ["extract_star.py", dd, "--auto", ob_id,
                            "--autobins", "6", "--display", "--lstep", lstep,
                            "--tag", tagstr, "--reducer", reducer]
                else:
                    pars = ["extractstar.py", dd, "--auto", ob_id,
                            "--autobins", "6", "--display", "--lstep", lstep,
                            "--centroid", "auto",
                            "--tag", tagstr, "--reducer", reducer]
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
                if not good:
                    good = 'N'
                if good[0].capitalize() != "Y":
                    logging.error("Try re-extraction again.")
                    logging.info("Removing *_%s_*" % tagstr)
                    flist = glob.glob("*_%s_*" % tagstr)
                    for f in flist:
                        os.remove(f)
                    sys.exit(1)
            # Re-report
            if args.local or args.testing:
                pars = ["pysedm_report.py", dd, "--contains", tagstr]
            else:
                pars = ["pysedm_report.py", dd, "--contains", tagstr, "--slack"]
            logging.info("Running " + " ".join(pars))
            ret = subprocess.call(pars)
            if ret:
                logging.error("pysedm_report.py failed!")
                sys.exit(1)
            # Upload spectrum to marshal
            if not args.local:
                home = os.environ['HOME']
                if args.testing:
                    cmd = [os.path.join(home, "spy"),
                           os.path.join(home, "sedmpy/growth/growth.py"), dd,
                           "--reducedby", reducer, "--testing"]
                else:
                    cmd = [os.path.join(home, "spy"),
                           os.path.join(home, "sedmpy/growth/growth.py"), dd,
                           "--reducedby", reducer]
                logging.info(" ".join(cmd))
                retcode = subprocess.call(cmd)
                if retcode != 0:
                    logging.error("Error uploading spectra to marshal")
            # Get the resulting fits file
            fits_file = glob.glob(
                os.path.join(rd, dd, "spec_auto_%s_*_%s.fits" %
                             (tagstr, obname)))
            if len(fits_file) == 1 and not args.testing:
                if not args.local:
                    # Update database entry
                    ar.update_spec(fits_file[0])
                    # e-mail user that re-extraction was made
                    ar.email_user(fits_file[0], dd, obname)
            else:
                if args.testing:
                    logging.info("TESTING RedEx.py script")
                else:
                    logging.error("Error finding fits file")
        ret = subprocess.call(["make", "report"])
