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
        obid = sys.argv[1]
        dd = os.getcwd().split('/')[-1]
        print("Re-extracting observation %s in %s" % (obid, dd))
        if len(sys.argv) == 4:
            xs = sys.argv[2]
            ys = sys.argv[3]
            print("extract_star.py %s --auto %s --autobins 6 --centroid %s %s" %
                  (dd, obid, xs, ys))
            res = subprocess.run(["extract_star.py", dd, "--auto", obid,
                                  "--autobins", "6", "--centroid", xs, ys])
            if res.returncode != 0:
                print("Extraction failed.")
                sys.exit(1)
        else:
            print("extract_star.py %s --auto %s --autobins 6 --display" %
                  (dd, obid))
            res = subprocess.run(["extract_star.py", dd, "--auto", obid,
                                  "--autobins", "6", "--display"])
            if res.returncode != 0:
                print("Extraction failed.")
                sys.exit(1)
        # Object name
        obname = glob.glob("spec_*_%s*.fits" %
                           obid)[0].split('_')[-1].split('.')[0]
        # Re-classify
        flist = glob.glob("spec_*_%s_%s_*" % (obid, obname))
        for f in flist:
            print("removing %s" % f)
            os.remove(f)
        print("make classify")
        res = subprocess.run(["make", "classify"])
        if res.returncode != 0:
            print("make classify failed!")
            sys.exit(1)
        # Re-verify
        cfile = glob.glob("crr_b_ifu%s_%s.fits" % (dd, obid))[0].split('.')[0]
        print("verify %s --contains %s" % (dd, cfile))
        res = subprocess.run(["verify", dd, "--contains", cfile])
        if res.returncode != 0:
            print("Verify failed!")
            sys.exit(1)
        # Re-report
        print("pysedm_report.py %s --contains %s --slack" % (dd, obid))
        res = subprocess.run(["pysedm_report.py", dd, "--contains", obid,
                              "--slack"])
        if res.returncode != 0:
            print("pysedm_report.py failed!")
            sys.exit(1)
        # Prepare for upload
        upf = glob.glob("spec_*_%s_%s.upl" % (obid, obname))[0]
        print("removing %s" % upf)
        os.remove(upf)
        print("be sure to run make ztfupload when you are done.")
