# -*- coding: utf-8 -*-
"""
links raw rc images into appropriate /scr2/sedm/phot directory

Created on Fri Mar  4 23:31:27 2016
Updated on Thu Feb  6 by neill

@author: nadiablago
"""
import datetime
import subprocess
import os
import glob

if __name__ == '__main__':

    directory = datetime.datetime.strftime(datetime.datetime.utcnow(), "%Y%m%d")
    if not os.path.isdir("/scr2/sedm/phot/%s" % directory):
        os.makedirs("/scr2/sedm/phot/%s" % directory)
    fls = glob.glob("/scr2/sedm/raw/%s/rc*" % directory)
    flphot = glob.glob("/scr2/sedm/phot/%s/rc*" % directory)
    if len(fls) > 0 and len(fls) > len(flphot):
        cmd = "ln -s /scr2/sedm/raw/%s/rc* /scr2/sedm/phot/%s/" % (directory,
                                                                   directory)
        subprocess.call(cmd, shell=True)


