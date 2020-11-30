# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 2020

@author: neill
"""

import os
import argparse
import subprocess
import re
import astropy.io.fits as pf


def parse_and_fill(spec, sniascore_output):
    """
    Receive the SEDM file and the output from snid. It parses the output
    from snid, fills a dictionary with the results and fills some of the
    comments in the original file, so that these could be used in the marshal.
    """
    
    pars = {"SNIascore": -1, "SNIascore_err": -1,
            "SNIa_z": -1, "SNIa_z_err": -1}

    if os.path.exists(sniascore_output):
        with open(sniascore_output, "r") as snia:

            line = snia.readline()

            results = line.split(',')

            pars["SNIascore"] = float(results[0])
            pars["SNIascore_err"] = float(results[1])
            pars["SNIa_z"] = float(results[2])
            pars["SNIa_z_err"] = float(results[3])

    print("SNIaScore RESULTS: Score=%(SNIascore).4f,"
          " Err=%(SNIascore_err).4f, "
          "z=%(SNIa_z).4f+-%(SNIa_z_err).4f" % pars)

    # Update fits file
    ff = pf.open(spec.replace('.txt', '.fits'), mode='update')
    # Update ascii file
    with open(spec, "r") as specIn:
        spec_lines = specIn.readlines()
        
        comment_lines_count = 0
        for line in spec_lines:
            if line.split()[0] == '#':
                comment_lines_count += 1
            else:
                break

        if os.path.exists(sniascore_output):
            ff[0].header['SNIASCOR'] = (pars["SNIascore"], 'SNIa score')
            ff[0].header['SNIASCER'] = (pars["SNIascore_err"], 'SNIa score err')
            ff[0].header['SNIASCZ'] = (pars["SNIa_z"], 'SNIa redshift')
            ff[0].header['SNIASCZE'] = (pars["SNIa_z_err"], 'SNIa redshift err')

            spec_lines.insert(comment_lines_count,
                              "# SNIASCORE_ZERR: " + str(pars["SNIa_z_err"])
                              + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIASCORE_Z: " + str(pars["SNIa_z"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIASCORE_ERR: " + str(pars["SNIascore_err"])
                              + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIASCORE: " + str(pars["SNIascore"]) + "\n")
        else:
            spec_lines.insert(comment_lines_count,
                              "# SNIASCORE: -1\n")
            ff[0].header['SNIASCOR'] = (-1, 'SNIa score')

    with open(spec, "w") as specOut:
        specOut.write("".join(spec_lines))
    ff.close()

    return pars


def run_sniascore(spec_file=None, overwrite=False):
    """
    Runs SNIascore in batch mode on the input file.  If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    ran = False
    if spec_file is not None:
        fl = spec_file
        # retrieve the quality and classification of the spectra.
        with open(fl, "r") as sfl:
            ll = sfl.readlines()

            q = [li for li in ll if "QUALITY" in li]
            n = [li for li in ll if "nan" in li and "#" not in li]

            if len(n) > 0:
                print("Found nan's in spectrum, skipping")
                return False

            if len(q) > 0:
                token = re.search(r'\(?([0-9]+)\)?', q[0])
                q = int(token.group(1))
            else:
                if "crr_b_ifu" in fl:
                    q = 1
                else:
                    q = 5
            print("quality = %d" % q)

            # check for previous classification
            clas = [li for li in ll if "SNIASCOR" in li]
            if clas:
                print("SNIascore: ", clas)

        if (q < 3 or q == 5) and (len(clas) <= 0 or overwrite):
            # If we are here, we run SNIascore
            cm = "matlab -nodisplay -r \"addpath('/home/cfremling/SEDM_ML/');" \
                 "SNIascore('%s');quit\"" % fl.replace(" ", "\\ ")
            print(cm)
            try:
                subprocess.call(cm, shell=True)
                ran = True
            except:
                print("Error running SNIascore")
        else:
            if q >= 3:
                print("low quality spectrum")
            if len(clas) > 0:
                print("already classified")
    return ran
    # END run_sniascore


def record_sniascore(spec_file=None):
    """
    Records snid results in specfile. If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    if spec_file is not None:
        fl = spec_file
        try:
            sniascore_output = fl.replace(".txt", ".txt.SNIascore")
            pars = parse_and_fill(fl, sniascore_output)
            res = ("SNIascore RESULTS: Score=%(SNIascore).3f,"
                   "Score_err=%(SNIascore_err).3f, "
                   "z=%(SNIa_z).4f+-%(SNIa_z_err).4f" % pars)
        except:
            print("Error recording SNIascore")
            res = ""
        return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""

        Classifies an input SEDM ascii spectrum with SNIascore.

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specfile', type=str, dest="specfile",
                        help='spec_*.txt ascii spectrum file.',
                        default=None)
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='Overwrite existing classification')
    parser.add_argument('--update', action="store_true", default=False,
                        help='Update existing classification')

    args = parser.parse_args()

    # Run SNIascore on extracted spectra and record results in specfile
    if args.update:
        print("Updating SNIascore results in %s" % args.specfile)
        record_sniascore(spec_file=args.specfile)
    else:
        print("Running SNIascore on and recording results in %s"
              % args.specfile)
        good = run_sniascore(spec_file=args.specfile, overwrite=args.overwrite)
        if good:
            record_sniascore(spec_file=args.specfile)
