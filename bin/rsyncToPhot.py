# -*- coding: utf-8 -*-
"""
links raw rc images into appropriate phot directory

Created on Fri Mar  4 23:31:27 2016
Updated on 29 May, 2022 by neill

@author: nadiablago
"""
import datetime
import subprocess
import os
import glob
import json
import argparse
import version

cfg_path = os.path.join(version.CONFIG_DIR, 'sedmconfig.json')
with open(cfg_path) as cfg_file:
    sedm_cfg = json.load(cfg_file)

_phot_dir = sedm_cfg['paths']['photpath']
_raw_dir = sedm_cfg['paths']['rawpath']

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

    source = os.path.join(_raw_dir, directory)
    dest = os.path.join(_phot_dir, directory)

    if not os.path.isdir(dest):
        os.makedirs(dest)
    fls = glob.glob(os.path.join(source, "rc2*"))
    flphot = glob.glob(os.path.join(dest, "rc2*"))
    if len(fls) > 0 and len(fls) > len(flphot):
        cmd = "ln -s %s/rc* %s/" % (source, dest)
        subprocess.call(cmd, shell=True)
        os.chdir(dest)
        cmd = "~/spy what rc%s_*.fits > rcwhat.list " % directory
        subprocess.call(cmd, shell=True)
    elif len(flphot) > len(fls):
        # remove bad links
        for fl in flphot:
            if 'rcwhat' in fl or 'Bias' in fl or 'Flat' in fl:
                continue
            if os.path.islink(fl):
                if not os.path.exists(os.readlink(fl)):
                    os.unlink(fl)
