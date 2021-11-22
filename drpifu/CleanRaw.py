import os
import glob
import argparse
from configparser import ConfigParser
import codecs
import datetime

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
    if not os.path.exists(outdir):
        print("Error: outdir not found: %s" % outdir)
        return 0
    flist = glob.glob(os.path.join(outdir, "*%s_*.fits.gz" % utdstr))
    for fl in flist:
        rawf = fl.split('.gz')[0]
        if os.path.exists(rawf):
            os.remove(rawf)
            ndel += 1
    # append dir to backup file
    if ndel > 0:
        back_file = cfg_parser.get('backup', 'raw_backup_file')
        try:
            with open(back_file, 'w') as bf:
                bf.writelines(utdstr + "\n")
            print("%s written to %s, ready for rsync" % (utdstr, back_file))
        except OSError:
            print("Cannot open backup file for update: %s" % back_file)

    return ndel


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Clean uncompressed raw files

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--date', type=str, default=None,
                        help='Select date to process (default: today)')

    args = parser.parse_args()

    if args.date is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        odir = os.path.join(_rawpath, timestamp)
    else:
        timestamp = args.date
        odir = os.path.join(_rawpath, args.date)
    print("Cleaning raw dir %s" % odir)
    nrm = clean_post_raw(odir, timestamp)
    print("%d raw files removed" % nrm)
