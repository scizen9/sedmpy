#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import subprocess
import argparse
rdir = '/scr2/sedm/raw/'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="update old headers",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--start_date', type=int, default=None,
            help='starting date for fixing headers')
    parser.add_argument('--end_date', type=int, default=None,
            help='ending date for fixing headers')
    parser.add_argument('--rezip', action="store_true", default=False,
            help='Re-gzip image files')

    args = parser.parse_args()

    if not args.start_date or not args.end_date:
        print("must supply both --start_date and --end_date")
    elif args.start_date < 20151115:
        print("start_date must be > 20151114")
    elif args.end_date > 20201231:
        print("end_date must be < 20210101")
    else:
        print("start: %d, end: %d" % (args.start_date, args.end_date))
        dir_list = glob.glob(os.path.join(rdir, '20??????'))
        for rd in dir_list:
            fdate = int(rd.split('/')[-1])
            if args.start_date <= fdate <= args.end_date:
                print(rd)
                os.chdir(rd)
                uzlist = os.listdir(rd)
                print("gunzipping %d files" % len(uzlist))
                for fl in uzlist:
                    if '.fits.gz' in fl:
                        subprocess.run(["gunzip", os.path.join(rd, fl)])
                subprocess.run(["/scr2/sedmdrp/spy", "/scr2/sedmdrp/sedmpy/drpifu/HdrFix.py", "--rawdir", rdir, "--date", rd.split('/')[-1]])
                if args.rezip:
                    zlist = os.listdir(rd)
                    print("gzipping %d files" % len(zlist))
                    for fl in zlist:
                        if '.fits' in fl:
                            subprocess.run(["gzip", os.path.join(rd, fl)])
                else:
                    print("skipping re-zipping")

