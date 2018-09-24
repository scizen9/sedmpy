#! /usr/bin/env python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    import sys
    import os
    import glob
    import subprocess

    if len(sys.argv) < 2:
        print("Usage - redex <obs_id> [<x> <y>]")
    else:
        # Check inputs and environment
        ob_id = sys.argv[1]
        dd = os.getcwd().split('/')[-1]
        rd = '/'.join(os.getcwd().split('/')[:-1])
        reddir = os.environ['SEDMREDUXPATH']
        if rd not in reddir:
            print("Environment error - check SEDMREDUXPATH env var")
            sys.exit(1)

        print("Re-extracting observation %s in %s" % (ob_id, dd))
        if len(sys.argv) == 4:
            xs = sys.argv[2]
            ys = sys.argv[3]
            pars = ["extract_star.py", dd, "--auto", ob_id, "--autobins", "6",
                    "--centroid", xs, ys]
            print("Running" + " ".join(pars))
            res = subprocess.run(pars)
            if res.returncode != 0:
                print("Extraction failed.")
                sys.exit(1)
        else:
            pars = ["extract_star.py", dd, "--auto", ob_id, "--autobins", "6",
                    "--display"]
            print("Running" + " ".join(pars))
            res = subprocess.run(pars)
            if res.returncode != 0:
                print("Extraction failed.")
                sys.exit(1)

        # Object name
        obname = glob.glob("spec_*_%s*.fits" %
                           ob_id)[0].split('_')[-1].split('.')[0]
        # Re-classify
        flist = glob.glob("spec_*_%s_%s_*" % (ob_id, obname))
        for f in flist:
            print("removing %s" % f)
            os.remove(f)
        print("make classify")
        res = subprocess.run(["make", "classify"])
        if res.returncode != 0:
            print("make classify failed!")
            sys.exit(1)
        # Re-verify
        cfile = glob.glob("crr_b_ifu%s_%s.fits" % (dd, ob_id))[0].split('.')[0]
        pars = ["verify", dd, "--contains", cfile]
        print("Running" + " ".join(pars))
        res = subprocess.run(pars)
        if res.returncode != 0:
            print("Verify failed!")
            sys.exit(1)
        # Re-report
        pars = ["pysedm_report.py", dd, "--contains", ob_id, "--slack"]
        print("Running" + " ".join(pars))
        res = subprocess.run(pars)
        if res.returncode != 0:
            print("pysedm_report.py failed!")
            sys.exit(1)
        # Prepare for upload
        upf = glob.glob("spec_*_%s_%s.upl" % (ob_id, obname))
        for uf in upf:
            print("removing %s" % uf)
            os.remove(uf)
        print("be sure to run make ztfupload when you are done.")
