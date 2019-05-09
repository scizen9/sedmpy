# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""
import datetime
import glob
import os
import argparse
import re
import drpifu.RunSnid as RunSnid


def classify(spec_dir='./', overwrite=False):
    """
    Runs snid in batch mode on all the *_SEDM.txt files found in the given
    directory.  If a given file was already classified, it skips it, unless
    overwrite is requested.
    """

    summary = []
    flb = glob.glob(os.path.join(spec_dir, "spec_*.txt"))
    for fl in flb:
        print(fl)
        # don't classify standard stars
        if "STD" in fl or "BD" in fl or "Feige" in fl or "HZ" in fl:
            print("standard star")
            continue
        if "Hiltner" in fl or "Kopff" in fl or "LTT" in fl:
            print("standard star")
            continue
        # don't classify galaxies
        if "NGC" in fl or "PGC" in fl or "MCG" in fl or "2MASX" in fl:
            print("galaxy")
            continue
        # don't classify stars
        if "TYC" in fl or "SAO" in fl or "HD" in fl or "Tycho" in fl:
            print("star")
            continue
        # skip uncalibrated objects
        if "notfluxcal" in fl:
            print("uncalibrated")
            continue
        # retrieve the quality and tracking rate of the spectra.
        with open(fl, "r") as sfl:
            lines = sfl.readlines()

            q = [li for li in lines if "QUALITY" in li]

            if len(q) > 0:
                token = re.search(r'\(?([0-9]+)\)?', q[0])
                q = int(token.group(1))
            else:
                if "crr_b_ifu" in fl:
                    q = 1
                else:
                    q = 5

            t = [li for li in lines if "RA_RATE" in li]

            if len(t) > 0:
                ra_rate = float(t[0].split()[2])
            else:
                ra_rate = 0.

            t = [li for li in lines if "DEC_RATE" in li]

            if len(t) > 0:
                dec_rate = float(t[0].split()[-1])
            else:
                dec_rate = 0.

            # If quality is good, check for previous classification
            if q < 3 or q == 5:
                clas = [li for li in lines if "SNID" in li]
                # If the file has been classified, move to the next
                if len(clas) > 0 and not overwrite:
                    print("already classified")
                    continue
            else:
                print("bad quality")
                continue

            # non-sidereal, then skip
            if ra_rate != 0. or dec_rate != 0.:
                print("non-sidereal")
                continue
            
        # If we are here, we run the classification with snid
        specfl = fl.split('/')[-1]
        good = RunSnid.run_snid(spec_file=specfl, overwrite=overwrite)
        # If we actually ran, record the results
        if good:
            res = specfl + " " + RunSnid.record_snid(spec_file=specfl)
            summary.append(res)
    # END loop over each file matching *_SEDM.txt
    for res in summary:
        print(res)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """

        Classifies the output ascii spectra

        Specify input directory with -d, or --specdir parameters.
        If none specified, use current date directory in /scr2/sedmdrp/redux/

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specdir', type=str, dest="specdir",
                        help='Directory with output spec_*.txt files.',
                        default=None)
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='Overwrite existing classification')

    args = parser.parse_args()
    
    specdir = args.specdir
    
    if specdir is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        try:
            reddir = os.environ['SEDMREDUXPATH']
        except KeyError:
            reddir = "/scr2/sedmdrp/redux"
        specdir = os.path.join(reddir, timestamp)
    else:
        specdir = os.path.abspath(specdir)
        timestamp = os.path.basename(specdir)

    os.chdir(specdir)
    print(os.getcwd())

    # Run snid on extracted spectra
    classify(spec_dir=specdir, overwrite=args.overwrite)
