# -*- coding: utf-8 -*-
"""
links raw rc images into appropriate /scr2/sedmdrp/phot directory

Created on Fri Mar  4 23:31:27 2016
Updated on Thu Feb 13 by neill

@author: nadiablago
"""
import datetime
import subprocess
import os
import glob
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""
        link rc files into phot directory
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--utdate', type=str, default=None,
                        help='UT date string: YYYYMMDD')

    args = parser.parse_args()

    if args.utdate:
        directory = args.utdate
    else:
        directory = datetime.datetime.strftime(datetime.datetime.utcnow(),
                                               "%Y%m%d")

    if not os.path.isdir("/scr2/sedmdrp/redux/phot/%s" % directory):
        os.makedirs("/scr2/sedmdrp/redux/phot/%s" % directory)
    fls = glob.glob("/scr2/sedm/raw/%s/rc*" % directory)
    flphot = glob.glob("/scr2/sedmdrp/redux/phot/%s/rc*" % directory)
    if len(fls) > 0 and len(fls) > len(flphot):
        cmd = "ln -s /scr2/sedm/raw/%s/rc* /scr2/sedmdrp/redux/phot/%s/"\
              % (directory, directory)
        subprocess.call(cmd, shell=True)
        os.chdir("/scr2/sedmdrp/redux/phot/%s" % directory)
        cmd = "/scr2/sedmdrp/miniconda3/bin/python " \
              "/scr2/sedmdrp/sedmpy/drpifu/What.py rc%s_*.fits > rcwhat.list" \
              % directory
        subprocess.call(cmd, shell=True)
    elif len(flphot) > len(fls):
        # remove bad links
        for fl in flphot:
            if os.path.islink(fl):
                if not os.path.exists(os.readlink(fl)):
                    os.unlink(fl)
