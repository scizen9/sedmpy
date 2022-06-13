#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import sedmpy_version
import json

# Get pipeline configuration
# Find config file: default is sedmpy/config/sedmconfig.json
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(sedmpy_version.CONFIG_DIR, 'sedmconfig.json')
with open(configfile) as config_file:
    sedm_cfg = json.load(config_file)

# Get paths
_rawpath = sedm_cfg['paths']['rawpath']
_reduxpath = sedm_cfg['paths']['reduxpath']


if __name__ == '__main__':

    # Make sure we are in the redux dir
    print("Changing directory to %s" % _reduxpath)
    os.chdir(_reduxpath)

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="re-process old data",
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
        dir_list = glob.glob(os.path.join(_rawpath, '20??????'))
        for rd in dir_list:
            fdate = rd.split('/')[-1]
            if args.start_date <= int(fdate) <= args.end_date:
                print(rd)
                # reduce
                cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --reduce " \
                      "--date %s >& rp_reduce.log" % fdate
                print(cmd)
                os.system(cmd)
                os.system("mv rp_reduce.log %s" % fdate)
                # calibrate
                cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --calibrate " \
                      "--date %s >& rp_calibrate.log" % fdate
                print(cmd)
                os.system(cmd)
                os.system("mv rp_calibrate.log %s" % fdate)
                # cube
                cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb --cube " \
                      "--date %s >& rp_cube.log" % fdate
                print(cmd)
                os.system(cmd)
                os.system("mv rp_cube.log %s" % fdate)
                # extract?
                if args.extract:
                    if args.manual:
                        cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb " \
                              "--extract --manual --date %s >& rp_extract.log" \
                              % fdate
                    else:
                        cmd = "~/spy ~/sedmpy/drpifu/ReProcess.py --nodb " \
                              "--extract --date %s >& rp_extract.log" % fdate
                    print(cmd)
                    os.system(cmd)
                    os.system("mv rp_extract.log %s" % fdate)
