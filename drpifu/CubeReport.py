#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import pil
except ImportError:
    import drpifu.pil as pil
import datetime
import glob
import os


def build_image_report(indir=None):
    """ """
    now = datetime.datetime.now()
    timestring = "made on %d-%02d-%02d at %02d:%02d:%02d" % (now.year,
                                                             now.month,
                                                             now.day,
                                                             now.hour,
                                                             now.minute,
                                                             now.second)
        
    footer = pil.get_buffer([20, 1], timestring, fontsize=15,
                            textprop=dict(color="0.5"), barcolor="k")

    # Missing plot format
    prop_missing = dict(fontsize=20, textprop=dict(color="C1"))

    # ---------
    # Wavelength Cube
    # ---------
    try:
        wl_file = glob.glob(indir + "_wavesolution_dispersionmap.png")[0]
        img_wl = pil.Image.open(wl_file)
    except IndexError:
        print("Cannot find wavelength cal cube image")
        img_wl = pil.get_buffer([8, 4], "WL Cal Cube image missing",
                                **prop_missing)

    # ---------
    # Flat Cube
    # ---------
    try:
        fl_file = glob.glob(indir + "_flat3d.png")[0]
        img_fl = pil.Image.open(fl_file)
    except IndexError:
        print("Cannot find wavelength cal cube image")
        img_fl = pil.get_buffer([8, 4], "Flat Cal Cube image missing",
                                **prop_missing)

    # ============== #
    #  Combination   #
    # ============== #    
    # title
    crtitle = "%s" % indir

    # Auto
    title_img = pil.get_buffer([8, 1], crtitle, fontsize=16, hline=[0.3, 0.7],
                               barprop=dict(lw=1))

    img_upper = pil.get_image_row([title_img])

    img_lower = pil.get_image_row([img_wl, img_fl])

    img_combined = pil.get_image_column([img_upper, img_lower, footer])
            
    return img_combined, "calcube_report_"+indir+".png"
        

#################################
#
#   MAIN 
#
#################################
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=""" verify calibration cubes
            """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('infile', type=str, default=None,
                        help='cube filepath (UT date as YYYYMMDD)')

    parser.add_argument('--slack', action="store_true", default=False,
                        help='Submit the report to slack')

    # ================ #
    # END of Option    #
    # ================ #
    args = parser.parse_args()

    print("Verifying cal cubes in directory %s" % args.infile)

    img_report, report_filename = build_image_report(indir=args.infile)
    if img_report is not None and report_filename is not None:
        img_report.save(report_filename, dpi=(1000, 1000))

    if args.slack:
        try:
            from pysedmpush import slack
        except ImportError:
            raise ImportError("you need to install pysedmpush"
                              " to be able to push on slack")
    # ================= #
    #   The Scripts     #
    # ================= #
    SLACK_CHANNEL = "pysedm-report"

    # Slack push report
    if args.slack:
        if not os.path.isfile(report_filename):
            print("No file-image created by the "
                  "CubeReport.build_image_report(). Nothing to push on slack")
        else:
            print("pushing the report to %s" % SLACK_CHANNEL)
            # Title & caption
            title = "pysedm-report: %s | Cal Cubes" % args.infile
            caption = ""
            # Push
            slack.push_image(report_filename, caption=caption, title=title,
                             channel=SLACK_CHANNEL)
