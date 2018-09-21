#! /usr/bin/env python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    import sys
    import os
    import glob

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
        else:
            print("extract_star.py %s --auto %s --autobins 6 --display" %
                  (dd, obid))
        # Object name
        obname = glob.glob("spec_*_%s*.fits" %
                           obid)[0].split('_')[-1].split('.')[0]
        # Re-classify
        flist = glob.glob("spec_*_%s_%s_*" % (obid, obname))
        for f in flist:
            print(f)
        print("make classify")
        # Re-verify
        cfile = glob.glob("crr_b_ifu%s_%s.fits" % (dd, obid))[0].split('.')[0]
        print("~/sedmpy/drpifu/Verify.py %s --contains %s" % (dd, cfile))
        # Re-report
        print("pysedm_report.py %s --contains %s --slack" % (dd, obid))
        # Prepare for upload
        upf = glob.glob("spec_*_%s_%s.upl" % (obid, obname))[0]
        print("rm %s" % upf)
        print("be sure to run make ztfupload when you are done.")
