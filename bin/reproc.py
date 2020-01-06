#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import subprocess
import argparse
rdir = '/scr2/sedm/raw/'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="re-process old data",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--start_date', type=str, default='20160601',
            help='starting date for fixing headers')
    parser.add_argument('--end_date', type=str, default='20161231',
            help='ending date for fixing headers')
    parser.add_argument('--extract', action="store_true", default=False,
            help='Perform extraction as well')

    args = parser.parse_args()

    if not args.start_date or not args.end_date:
        print("must supply both --start_date and --end_date")
    else:
        sdate = int(args.start_date)
        edate = int(args.end_date)
        print("start: %d, end: %d" % (sdate, edate))
        dir_list = glob.glob(os.path.join(rdir, '20??????'))
        for rd in dir_list:
            fdate = rd.split('/')[-1]
            if sdate <= int(fdate) <= edate:
                print(rd)
                # reduce
                cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --reduce --date %s >& rp_reduce.log" % fdate
                print(cmd)
                os.system(cmd)
                os.system("mv rp_reduce.log %s" % fdate)
                # calibrate
                cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --calibrate --date %s >& rp_calibrate.log" % fdate
                print(cmd)
                os.system(cmd)
                os.system("mv rp_calibrate.log %s" % fdate)
                # cube
                cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --cube --date %s >& rp_cube.log" % fdate
                print(cmd)
                os.system(cmd)
                os.system("mv rp_cube.log %s" % fdate)
                # extract?
                if args.extract:
                    cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --extract --date %s >& rp_extract.log" % fdate
                    print(cmd)
                    os.system(cmd)
                    os.system("mv rp_extract.log %s" % fdate)
