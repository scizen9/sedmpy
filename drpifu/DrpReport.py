import glob
import os
import time
import subprocess


def report():
    """Generate DRP report using output spec_*.txt files"""

    flist = glob.glob("spec_*.txt")
    flist.sort()
    print("\nReport generated on %s" % time.strftime("%c"))
    print("\nSEDM DRP run in %s\nFound %d spec_*.txt files" %
          (os.getcwd(), len(flist)))
    print("See http://pharos.caltech.edu/data_access/ifu?obsdate=%s\n" %
          os.getcwd().split('/')[-1])

    print("UTStart  Object                    Exptime Air    Qual"
          "                           method  Allocation                     "
          "Type Subtype  z           Rlap")
    recs = []
    for f in flist:
        # Get object name
        tname = f.split('_ifu')[-1].split('_')[4:]
        if len(tname) > 1:
            objname = '_'.join(tname).split('.txt')[0]
        else:
            objname = tname[0].split('.txt')[0]
        # Get time string
        tstr = ':'.join(f.split('_ifu')[-1].split('_')[1:4])
        # check the ascii spectrum file for SNID data
        with open(f, "r") as sfl:
            lines = sfl.readlines()
            # check for SNID classification
            clas = [li for li in lines if "SNIDMATCHTYPE" in li]
            ctype = ""
            if len(clas) > 0:
                for cl in clas:
                    ctype += (" %s" % cl.split()[-1])
            # check for SNID classification
            clas = [li for li in lines if "SNIDMATCHSUBTYPE" in li]
            stype = ""
            if len(clas) > 0:
                for cl in clas:
                    stype += (" %s" % cl.split()[-1])
            # get redshift
            zmch = [li for li in lines if "REDSHIFT" in li]
            if len(zmch) > 0:
                zmch = ("%.4f" % float(zmch[0].split()[-1]))
            else:
                zmch = ""
            # get rlap
            rlap = [li for li in lines if "RLAP" in li]
            if len(rlap) > 0:
                rlap = ("%.2f" % float(rlap[0].split()[-1]))
            else:
                rlap = ""
            # get flux calibration
            flxcal = [li for li in lines if "FLUXCAL" in li]
            if len(flxcal) > 0:
                flxcal = flxcal[0].split()[-1]
            else:
                flxcal = "False"
            # get method
            meth = f.split('__crr')[0].split('spec_')[-1]
            # get exposure time
            expt = [li for li in lines if "EXPTIME" in li]
            if len(expt) > 0:
                expt = ("%.1f" % float(expt[0].split()[-1]))
            else:
                expt = "-"
            # get airmass
            air = [li for li in lines if "AIRMASS" in li]
            if len(air) > 0:
                air = ("%.3f" % float(air[0].split()[-1]))
            else:
                air = "-"
            # get program ID
            prid = [li for li in lines if "P60PRID" in li]
            if len(prid) > 0:
                prid = prid[0].split()[-1]
            else:
                prid = "-"
            # get Quality
            quality = [li for li in lines if "QUALITY" in li]
            if len(quality) > 0:
                quality = int(quality[0].split()[-1])
            else:
                quality = 9
            sfl.close()
        if ctype == "" or quality == 5:
            if "STD" in f:
                ctype = " STD"
            else:
                ctype = " QUALITY_%d" % quality
                stype = ""
                zmch = ""
                rlap = ""

        recs.append("%8s %-25s %7s %5s     %d %32s  %-21s  %12s  "
                    "%6s %-9s  %6s" %
                    (tstr, objname, expt, air, quality, meth, prid, ctype,
                     stype, zmch, rlap))
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
            tok = f.split("_failed")[0].split("_ifu")[-1]
            try:
                out = subprocess.check_output(('grep', tok, 'what.list'),
                                              universal_newlines=True)
                out = out.split()[3]
            except subprocess.CalledProcessError:
                out = tok
            print("%8s %-25s FAILED" % (tstr, out))
    else:
        print("\nThere were no failed extractions")


if __name__ == '__main__':
    report()
