import argparse
import datetime
import os
import subprocess
import logging

logging.basicConfig()
logFormatter = '%(asctime)s - %(levelname)s - %(message)s'
logging.basicConfig(format=logFormatter, level=logging.DEBUG)
logger = logging.getLogger(__name__)

computer = os.uname()[1]  # a quick fix to switch between development computers

if computer == 'pele':
    _logpath = '/scr7/rsw/sedm/log/'
    _photpath = '/scr7/rsw/sedm/phot/'
    pypath = '/scr7/rsw/anaconda3/bin/python'

elif computer == 'pharos':
    raw_dir = '/scr2/sedm/raw/'
    phot_dir = '/scr2/sedm/phot/'
    redux_dir = '/scr2/sedmdrp/redux/'
    host = 'localhost'


def reduce_all_dir(photdir, overwrite=False):
    # Reduce the data that is already in the directory.
    cmd = "%s rcred_v3.py -d %s" % (pypath, photdir)

    if (overwrite):
        cmd = cmd + " -o"
    subprocess.call(cmd, shell=True)
    logger.info("Reduce all dir: %s" % cmd)

    # Copy the content of the reduced directory into a new directory with the date of the observations.
    dayname = os.path.basename(photdir)
    reducedname = os.path.join(photdir, "reduced")

    # Reduce the data that is already in the directory.
    cmd = "%s zeropoint.py  %s" % (pypath, reducedname)
    subprocess.call(cmd, shell=True)
    logger.info("zeropoint for all dir: %s" % cmd)

    if (os.path.isdir(reducedname)):
        cmd = "rcp -r %s grbuser@transient.caltech.edu:/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s" % (
        reducedname, dayname)
        subprocess.call(cmd, shell=True)
    else:
        os.makedirs(reducedname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.',
                        default=None)
    parser.add_argument('-f', '--fullred', action='store_true', dest="fullred", default=False,
                        help='Whether we should do a full reduction.')
    parser.add_argument('-o', '--overwrite', action="store_true", help='re-reduce and overwrite the reduced images?',
                        default=False)

    args = parser.parse_args()

    photdir = args.photdir
    fullred = args.fullred
    overwrite = args.overwrite
    fullred = True
    photdir = '/scr7/rsw/sedm/phot/20181018/'

    if (photdir is None):
        photdir = os.path.join(_photpath, datetime.datetime.utcnow().strftime("%Y%m%d"))
    if (fullred):
        reduce_all_dir(os.path.abspath(photdir), overwrite=overwrite)
    #reduce_on_the_fly(os.path.abspath(photdir))

    # After 12h, invoke the zeropoint calibration.
    #zeropoint.main(os.path.join(os.path.abspath(photdir), "reduced"))
