# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""

import os
import glob
import argparse
import subprocess
import re
import astropy.io.fits as pf


def find_line_match(lines, match_dict):
    """
    Returns the line number of first line matching the position-word pairs
    provided within matchDict, if present within lines.
    If no match is found, -1 is returned.
    """
    line_number = -1

    for index in range(0, len(lines)):
        for position_key in match_dict:
            if (len(lines[index].split()) == 0 or 
                    len(lines[index].split()) <= position_key or
                    lines[index].split()[position_key] !=
                    match_dict[position_key]):
                pass
            else:
                line_number = index
                break

    return line_number


def parse_and_fill(spec, snidoutput):
    """
    Receive the SEDM file and the output from snid. It parses the output
    from snid, fills a dictionary with the results and fills some of the
    comments in the original file, so that these could be used in the marshal.
    """
    
    pars = {"zmed": -1, "zmederr": -1, "agemed": -1, "agemederr": -1,
            "templ": "", "bestMatchAge": -1,
            "Ia": 0, "Ib": 0, "Ic": 0, "II": 0, "NotSN": 0, "rlap": 0,
            "bestMatchType": "", "bestMatchSubtype": "", "bestMatchRedshift": 0}

    if os.path.exists(snidoutput):
        with open(snidoutput, "r") as snid:

            lines = snid.readlines()

            zmed_line = find_line_match(lines, {0: 'zmed'})
            agem_line = find_line_match(lines, {0: 'agem'})
            type_ia_line = find_line_match(lines, {0: 'Ia'})
            type_ib_line = find_line_match(lines, {0: 'Ib'})
            type_ic_line = find_line_match(lines, {0: 'Ic'})
            type_ii_line = find_line_match(lines, {0: 'II'})
            type_not_sn_line = find_line_match(lines, {0: 'NotSN'})
            best_match_line = find_line_match(lines, {0: '1'})

            pars["zmed"] = float(lines[zmed_line].split()[1])
            pars["zmederr"] = float(lines[zmed_line].split()[2])
            pars["agemed"] = float(lines[agem_line].split()[1])
            pars["agemederr"] = float(lines[agem_line].split()[2])
            pars["Ia"] = float(lines[type_ia_line].split()[2])
            pars["Ib"] = float(lines[type_ib_line].split()[2])
            pars["Ic"] = float(lines[type_ic_line].split()[2])
            pars["II"] = float(lines[type_ii_line].split()[2])
            pars["NotSN"] = float(lines[type_not_sn_line].split()[2])
            pars["rlap"] = float(lines[best_match_line].split()[4])
            pars["templ"] = lines[best_match_line].split()[1]
            pars["bestMatchAge"] = float(lines[best_match_line].split()[7])
            # Test type
            typstr = lines[best_match_line].split()[2]
            sntype = ''.join(c for c in typstr if ord(c) >= 32)
            if len(sntype) > 0:
                pars["bestMatchType"] = sntype.split("-")[0]
                if len(sntype.split("-")) > 1:
                    pars["bestMatchSubtype"] = sntype.split("-")[1]
                else:
                    pars["bestMatchSubtype"] = '-'
            else:
                pars["bestMatchType"] = "None"
                pars["bestMatchSubtype"] = '-'
                pars["templ"] = '-'
            pars["bestMatchRedshift"] = float(lines[best_match_line].split()[5])

    print("SNID RESULTS: Type=%(bestMatchType)s, Rlap=%(rlap).2f, "
          "Age=%(agemed).2f+-%(agemederr)s day, "
          "z=%(zmed).4f+-%(zmederr).4f" % pars)

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

        if os.path.exists(snidoutput):
            ff[0].header['SNIDTYPE'] = (pars["bestMatchType"],
                                        'SNID mtach type')
            ff[0].header['SNIDSUBT'] = (pars["bestMatchSubtype"],
                                        'SNID match subtype')
            ff[0].header['SNIDRLAP'] = (pars["rlap"], 'SNID match rlap')
            ff[0].header['SNIDZMAT'] = (pars["bestMatchRedshift"],
                                        'SNID match z')
            ff[0].header['SNIDAMAT'] = (pars["bestMatchAge"],
                                        'SNID match age (days)')
            ff[0].header['SNIDTEMP'] = (pars["templ"], 'SNID match template')
            ff[0].header['SNIDZMED'] = (pars["zmed"], 'SNID z med')
            ff[0].header['SNIDZERR'] = (pars["zmederr"], 'SNID z med err')
            ff[0].header['SNIDAMED'] = (pars["agemed"], 'SNID Age med (days)')
            ff[0].header['SNIDAERR'] = (pars["agemederr"],
                                        'SNID Age med err (days)')
            ff[0].header['SNIDFIA'] = (pars["Ia"], 'SNID Fraction SN Ia')
            ff[0].header['SNIDFIB'] = (pars["Ib"], 'SNID Fraction SN Ib')
            ff[0].header['SNIDFIC'] = (pars["Ic"], 'SNID Fraction SN IC')
            ff[0].header['SNIDFII'] = (pars["II"], 'SNID Fraction SN II')
            ff[0].header['SNIDFNOT'] = (pars["NotSN"], 'SNID Fraction Not SN')

            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_NOTSN: " + str(pars["NotSN"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_II: " + str(pars["II"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_IC: " + str(pars["Ic"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_IB: " + str(pars["Ib"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_IA: " + str(pars["Ia"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDAGEMEDERR: " + str(pars["agemederr"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDAGEMED: " + str(pars["agemed"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDZMEDERR: " + str(pars["zmederr"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDZMED: " + str(pars["zmed"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHTEMPLATE: " + pars["templ"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHAGE: " +
                              str(pars["bestMatchAge"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHREDSHIFT: " +
                              str(pars["bestMatchRedshift"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHRLAP: " + str(pars["rlap"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHSUBTYPE: " +
                              pars["bestMatchSubtype"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHTYPE: " +
                              pars["bestMatchType"] + "\n")
        else:
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHTYPE: NONE\n")
            ff[0].header['SNIDTYPE'] = ("NONE", 'SNID match type')

    with open(spec, "w") as specOut:
        specOut.write("".join(spec_lines))
    ff.close()

    return pars["bestMatchType"], pars


def run_snid(spec_file=None, overwrite=False):
    """
    Runs snid in batch mode on the input file.  If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    ran = False
    if spec_file is not None:
        fl = spec_file
        # retrieve the quality and classification of the spectra.
        with open(fl, "r") as sfl:
            l = sfl.readlines()

            q = [li for li in l if "QUALITY" in li]

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
            clas = [li for li in l if "SNID" in li]
            if clas:
                print("classification: ", clas)

        if (q < 3 or q == 5) and (len(clas) <= 0 or overwrite):
            # If we are here, we run the classification with snid
            cm = "snid wmin=4000 wmax=9500 skyclip=1 medlen=20 aband=1" \
                 " rlapmin=4 inter=0 plot=2 %s" % fl.replace(" ", "\\ ")
            print(cm)
            try:
                subprocess.call(cm, shell=True)
                ran = True
            except:
                print("Error running snid")
        else:
            if q >= 3:
                print("low quality spectrum")
            if len(clas) > 0:
                print("already classified")
    return ran
    # END run_snid


def record_snid(spec_file=None):
    """
    Records snid results in specfile. If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    if spec_file is not None:
        fl = spec_file
        try:
            snidoutput = fl.replace(".txt", "_snid.output")
            snid_type, pars = parse_and_fill(fl, snidoutput)
            psfspec = fl.replace(".txt", "_comp0001_snidflux*")
            psoutput = glob.glob(psfspec)
            if len(psoutput) > 0:
                psoutput = psoutput[0]
            else:
                psoutput = 'snid_output_file_not_found'
            if os.path.exists(psoutput):
                # handle spaces, if they exist
                fl = fl.replace(" ", "\\ ")
                psoutput = psoutput.replace(" ", "\\ ")
                pngfile = fl.replace(".txt", "_" + snid_type + ".png")
                cm = "convert -flatten -rotate 90 " + psoutput + " " + pngfile
                subprocess.call(cm, shell=True)
            res = ("SNID RESULTS: Type=%(bestMatchType)s, Rlap=%(rlap).2f, "
                   "Age=%(agemed).2f+-%(agemederr)s day, "
                   "z=%(zmed).4f+-%(zmederr).4f" % pars)
        except:
            print("Error recording snid")
            res = ""
        return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """

        Classifies an input SEDM ascii spectrum with SNID.

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specfile', type=str, dest="specfile",
                        help='spec_*.txt ascii spectrum file.',
                        default=None)
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='Overwrite existing classification')
    parser.add_argument('--update', action="store_true", default=False,
                        help='Update existing classification')

    args = parser.parse_args()

    # Run snid on extracted spectra and record results in specfile
    if args.update:
        print("Updating snid results in %s" % args.specfile)
        record_snid(spec_file=args.specfile)
    else:
        print("Running snid on and recording results in %s" % args.specfile)
        good = run_snid(spec_file=args.specfile, overwrite=args.overwrite)
        if good:
            record_snid(spec_file=args.specfile)
