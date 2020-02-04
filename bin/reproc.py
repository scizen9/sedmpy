#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import argparse
rdir = '/scr2/sedm/raw/'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="re-process old data",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--start_date', type=int, default=None,
            help='starting date for fixing headers')
    parser.add_argument('--end_date', type=int, default=None,
            help='ending date for fixing headers')
    parser.add_argument('--extract', action="store_true", default=False,
            help='Perform extraction as well')
    parser.add_argument('--manual', action="store_true", default=False,
                        help='Allow manual extraction')

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
            fdate = rd.split('/')[-1]
            if args.start_date <= int(fdate) <= args.end_date:
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
                    if args.manual:
                        cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --extract --manual --date %s >& rp_extract.log" % fdate
                    else:
                        cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --extract --date %s >& rp_extract.log" % fdate
                    print(cmd)
                    os.system(cmd)
                    os.system("mv rp_extract.log %s" % fdate)
