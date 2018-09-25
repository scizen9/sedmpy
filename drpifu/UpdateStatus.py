import glob
import os
import time
import subprocess


def update_status():
    """Generate DRP report using output spec_*.txt files"""

    flist = glob.glob("spec_*.txt")
    flist.sort()
    print("\nStatus update report generated on %s" % time.strftime("%c"))
    print("\nSEDM DRP run in %s\nFound %d spec_*.txt files" %
          (os.getcwd(), len(flist)))

    print("UTStart  Object                    Q Status")
    recs = []
    for f in flist:
        # Get object name
        tname = f.split('ifu')[-1].split('_')[4:]
        if len(tname) > 1:
            objname = '_'.join(tname).split('.txt')[0]
        else:
            objname = tname[0].split('.txt')[0]

        # Get time string
        tstr = ':'.join(f.split('ifu')[-1].split('_')[1:4])

        # check the ascii spectrum file for quality data
        with open(f, "r") as sfl:
            lines = sfl.readlines()

            # check for request id
            req_id = [li for li in lines if "REQ_ID" in li]
            if req_id:
                req_id = req_id.split()
                if len(req_id) > 2:
                    req_id = req_id[2]
                else:
                    req_id = None
            else:
                req_id = None

            # check for Quality
            quality = [li for li in lines if "QUALITY" in li]
            if len(quality) > 0:
                quality = int(quality.split()[-1])
            sfl.close()

        # Do we make an update?
        if req_id:
            if quality < 3:
                stat = "no change"
            else:
                stat = "to pending"
        else:
            stat = "no request"

        recs.append("%8s %-25s %d %7s" % (tstr, objname, quality, stat))
    recs.sort()
    for r in recs:
        print(r)
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
    update_status()
