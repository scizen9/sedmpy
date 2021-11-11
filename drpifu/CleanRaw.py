import os
import glob
import argparse
from configparser import ConfigParser
import codecs

# Get pipeline configuration
cfg_parser = ConfigParser()
# Find config file: default is sedmpy/config/sedmconfig.cfg
try:
    configfile = os.environ["SEDMCONFIG"]
except KeyError:
    configfile = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                              "../config/sedmconfig.cfg")
# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    cfg_parser.read_file(f)
# Get paths
_rawpath = cfg_parser.get('paths', 'rawpath')


def clean_post_raw(outdir, utdstr):
    """Remove uncompressed raw files and add record to backup file"""
    ndel = 0
    # Remove raw uncompressed file
    flist = glob.glob(os.path.join(outdir, "*%s_*.fits.gz" % utdstr))
    for fl in flist:
        rawf = fl.split('.gz')[0]
        if os.path.exists(rawf):
            os.remove(rawf)
            ndel += 1
    # append dir to backup file
    back_file = cfg_parser.get('backup', 'raw_backup_file')
    if os.path.exists(back_file):
        with open(back_file, 'a') as bf:
            bf.writelines(utdstr + "\n")
        print("%s appended to %s, ready for rsync" % (utdstr, back_file))
    else:
        print("Cannot open backup file for update: %s" % back_file)

    return ndel


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Clean uncompressed raw files

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--date', type=str, default=None,
                        help='Select date to process')

    args = parser.parse_args()

    if args.date is not None:
        odir = os.path.join(_rawpath, args.date)
        print("Cleaning raw dir %s" % odir)
        nrm = clean_post_raw(odir, args.date)
        print("%d raw files removed" % nrm)
    else:
        print("Error: Must provide a UTDate to raw_clean with --date")