import glob
import os
import sys
import datetime
import subprocess
import argparse
import db.SedmDb


def update_status(target_dir, upfil=None, dbase=None):
    """Generate DRP report using output spec_*.txt files

    :param target_dir:
    :param upfil:
    :param dbase:

    """

    if target_dir[-1] != '/':
        target_dir += '/'

    flist = glob.glob("%sspec_*.txt" % target_dir)
    flist.sort()

    for f in flist:
        # Have we already updated this one?
        if os.path.exists(f.split('.')[0] + ".upd"):
            print("Already updated: %s" % f)
            continue
        # Get object name
        tname = f.split('ifu')[-1].split('_')[4:]
        if len(tname) > 1:
            objname = '_'.join(tname).split('.txt')[0]
        else:
            objname = tname[0].split('.txt')[0]

        # check the ascii spectrum file for status data
        with open(f, "r") as sfl:
            lines = sfl.readlines()

            # check for request id
            req_id = [li for li in lines if "REQ_ID" in li]
            if req_id:
                req_id = req_id[0].split()
                if len(req_id) > 2:
                    req_id = req_id[2]
                else:
                    print("no request ID found for %s" % f)
                    continue
            else:
                print("no request ID found for %s" % f)
                continue
            # Get observation id
            res = dbase.get_from_observation(["id"], {"request_id": req_id})[0]

            # check for Quality
            quality = [li for li in lines if "QUALITY" in li]
            if quality:
                quality = quality[0].split()
                if len(quality) > 2:
                    quality = int(quality[2])
                else:
                    quality = None
            else:
                quality = None

            # check for e-mail
            email = [li for li in lines if "EMAIL" in li]
            if email:
                email = email[0].split()
                if len(email) > 2:
                    email = email[2]
                else:
                    email = None
            else:
                email = None

            # close spectrum file
            sfl.close()

        # Do we make an update?
        if req_id:
            if quality < 3:
                stat = "no change"
            else:
                stat = "to pending"
        else:
            stat = "no request"

    # Check for failed extractions
    flist = glob.glob("spec_*failed.fits")
    if len(flist) > 0:
        flist.sort()
        print("\nThere were/was %d failed extraction(s)" % len(flist))
        for f in flist:
            tstr = ':'.join(f.split('ifu')[-1].split('_')[1:4])
            tok = f.split("_failed")[0].split("ifu")[-1]
            out = subprocess.check_output(('grep', tok, 'what.list'),
                                          universal_newlines=True)
            print("%8s %-25s FAILED" % (tstr, out.split()[3]))
    else:
        print("\nThere were no failed extractions")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""

    Updates database status.

    """,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('indate', type=str, default=None,
                        help='input directory date (UT date as YYYYMMDD)')
    parser.add_argument('--data_file', type=str, default=None,
                        help='Data file to upload.')

    args = parser.parse_args()

    # Check environment
    try:
        reddir = os.environ["SEDMREDUXPATH"]
    except KeyError:
        print("please set environment variable SEDMREDUXPATH")
        sys.exit(1)

    # Get source dir
    if args.indate:
        utc = args.indate
    else:
        utc = datetime.datetime.utcnow().strftime("%Y%m%d")
    srcdir = reddir + '/' + utc + '/'

    # Check source dir
    if not os.path.exists(srcdir):
        print("Dir not found: %s" % srcdir)
    else:
        print("Getting spec status from %s" % srcdir)

        sedmdb = db.SedmDb.SedmDB()
        update_status(srcdir, upfil=args.data_file, dbase=sedmdb)
