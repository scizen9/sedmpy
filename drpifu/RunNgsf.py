# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""

import os
import csv
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


def parse_and_fill(spec, ngsfoutput):
    """
    Receive the SEDM file and the output from ngsf. It parses the output
    from ngsf, fills a dictionary with the results and fills some of the
    comments in the original file, so that these could be used in the marshal.
    """

    # pre-process pars
    try:
        z = float(ngsfoutput[5])
    except ValueError:
        z = -100.
    try:
        a_v = float(ngsfoutput[6])
    except ValueError:
        a_v = -100.
    try:
        phase = float(ngsfoutput[7])
    except ValueError:
        phase = -100.
    try:
        frac_sn = float(ngsfoutput[9])
    except ValueError:
        frac_sn = -1.
    try:
        frac_gal = float(ngsfoutput[10])
    except ValueError:
        frac_gal = -1.
    try:
        chi2_dof = float(ngsfoutput[11])
    except ValueError:
        chi2_dof = -1.
    try:
        chi2_dof2 = float(ngsfoutput[12])
    except ValueError:
        chi2_dof2 = -1.
    pars = {"HostType": ngsfoutput[1],
            "z": z,
            "A_v": a_v,
            "Phase": phase,
            "Band": ngsfoutput[8],
            "Frac(SN)": frac_sn,
            "Frac(gal)": frac_gal,
            "CHI2/dof": chi2_dof,
            "CHI2/dof2": chi2_dof2,
            "Templ": "", "bestMatchType": "", "bestMatchSubtype": "",
            "Instr": ""}

    sn = ngsfoutput[2]
    segs = sn.split('/')
    pars["Templ"] = segs[1]
    pars["bestMatchType"] = segs[0].split()[0]
    if " " in segs[0]:
        pars["bestMatchSubtype"] = segs[0].split()[-1]
    pars["Instr"] = segs[2].split()[0]

    # print("NGSF RESULTS: Type=%(bestMatchType)s, Host=%(HostType)s, "
    #       "CHI2=%(CHI2/dof).2f, Age=%(Phase).2f day, z=%(z).4f" % pars)

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

        if ngsfoutput is not None:
            ff[0].header['NGSFTYPE'] = (pars["bestMatchType"],
                                        'NGSF mtach type')
            ff[0].header['NGSFSUBT'] = (pars["bestMatchSubtype"],
                                        'NGSF match subtype')
            ff[0].header['NGSFCHI2'] = (pars["CHI2/dof"], 'NGSF fit CHI2/dof')
            ff[0].header['NGSFCH22'] = (pars["CHI2/dof2"], 'NGSF fit CHI2/dof2')
            ff[0].header['NGSFZ'] = (pars["z"], 'NGSF z')
            ff[0].header['NGSFAV'] = (pars["A_v"], 'NGSF A_v')
            ff[0].header['NGSFPHAS'] = (pars["Phase"], 'NGSF Phase (days)')
            ff[0].header['NGSFTEMP'] = (pars["Templ"], 'NGSF fit template')
            ff[0].header['NGSFFSN'] = (pars["Frac(SN)"], 'NGSF Fraction SN')
            ff[0].header['NGSFFGAL'] = (pars["Frac(gal)"], 'NGSF Fraction Gal')
            ff[0].header['NGSFHTYP'] = (pars["HostType"], 'NGSF Host Gal Type')
            ff[0].header['NGFSINST'] = (pars["Instr"], 'NGSF Templ Instrume')

            spec_lines.insert(comment_lines_count,
                              "# NGSFTEMPLINSTR: " + pars["Instr"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFHOSTTYPE: " + pars["HostType"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFFRAC_GAL: " + str(
                                  pars["Frac(gal)"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFFRAC_SN: " + str(pars["Frac(SN)"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFMATCHTEMPLATE: " + pars["Templ"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFPHASE: " + str(pars["Phase"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFAV: " + str(pars["A_v"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFREDSHIFT: " + str(pars["z"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFFITCHI2_DOF2: " + str(
                                  pars["CHI2/dof2"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFFITCHI2_DOF: " + str(
                                  pars["CHI2/dof"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFSUBTYPE: " +
                              pars["bestMatchSubtype"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# NGSFTYPE: " +
                              pars["bestMatchType"] + "\n")
        else:
            spec_lines.insert(comment_lines_count,
                              "# NGSFTYPE: NONE\n")
            ff[0].header['NGSFTYPE'] = ("NONE", 'NGSF match type')

    with open(spec, "w") as specOut:
        specOut.write("".join(spec_lines))
    ff.close()

    return pars["bestMatchType"], pars


def run_ngsf(spec_file=None, overwrite=False):
    """
    Runs ngsf in batch mode on the input file.  If a given file was already
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
            clas = [li for li in l if "NGSF" in li]
            if clas:
                print("classification: ", clas)

        if (q < 3 or q == 5) and (len(clas) <= 0 or overwrite):
            # If we are here, we run the classification with ngsf
            cm = "ngsf_run %s" % fl.replace(" ", "\\ ")
            print(cm)
            try:
                subprocess.call(cm, shell=True)
                ran = True
            except:
                print("Error running ngsf")
        else:
            if q >= 3:
                print("low quality spectrum")
            if len(clas) > 0:
                print("already classified")
    return ran
    # END run_ngsf


def record_ngsf(spec_file=None):
    """
    Records ngsf results in specfile. If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    if spec_file is not None:
        fl = spec_file
        try:
            ngsf_results = []
            ngsfoutput = fl.replace(".txt", ".csv")
            with open(ngsfoutput, mode='r') as inp:
                reader = csv.reader(inp)
                for row in reader:
                    ngsf_results.append(row)
            ngsf_type, pars = parse_and_fill(fl, ngsf_results[1])

            res = ("NGSF RESULTS: Type=%(bestMatchType)s, Host=%(HostType)s, "
                   "CHI2=%(CHI2/dof).2f, Age=%(Phase).2f day, z=%(z).4f" % pars)
        except:
            print("Error recording ngsf")
            res = ""
        return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """

        Classifies an input SEDM ascii spectrum with NGSF.

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specfile', type=str, dest="specfile",
                        help='spec_*.txt ascii spectrum file.',
                        default=None)
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='Overwrite existing classification')
    parser.add_argument('--update', action="store_true", default=False,
                        help='Update existing classification')

    args = parser.parse_args()

    # Run ngsf on extracted spectra and record results in specfile
    if args.update:
        print("Updating ngsf results in %s" % args.specfile)
        record_ngsf(spec_file=args.specfile)
    else:
        print("Running ngsf on and recording results in %s" % args.specfile)
        good = run_ngsf(spec_file=args.specfile, overwrite=args.overwrite)
        if good:
            record_ngsf(spec_file=args.specfile)
