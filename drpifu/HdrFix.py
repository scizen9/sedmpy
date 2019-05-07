#! /usr/bin/env python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    import os
    import glob
    import argparse
    import logging
    import astropy.io.fits as pf

    logging.basicConfig(
        format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
        datefmt='%Y%m%d %H:%M:%S', level=logging.INFO)

    # setup arguments parser
    parser = argparse.ArgumentParser(
        description="""Do an aperture extraction of a standard star.""",
        formatter_class=argparse.RawTextHelpFormatter)
    # setup arguments
    parser.add_argument('indate', type=str, default=None,
                        help='UT date as YYYYMMDD')
    parser.add_argument('--fspec', type=str, default=None,
                        help='unique observation identifier')
    args = parser.parse_args()

    dd = args.indate
    fspec = args.fspec

    if not fspec:
        logging.info("Usage - hdrfix YYYYMMDD --fspec <filespec>")
    else:
        # Check inputs and environment
        reddir = os.environ['SEDMREDUXPATH']
        # get a list of files
        flist = glob.glob(os.path.join(reddir, dd, fspec))
        # loop over files
        for fl in flist:
            logging.info(fl)
            ff = pf.open(fl, 'update')
            # image type
            if 'OBJECT' in ff[0].header:
                obj = ff[0].header['OBJECT']
                if 'STD-' in obj:
                    ff[0].header['IMGTYPE'] = 'Standard'
                elif 'Calib' in obj:
                    if 'dome' in obj:
                        ff[0].header['IMGTYPE'] = 'dome'
                    elif 'bias' in obj:
                        ff[0].header['IMGTYPE'] = 'bias'
                    else:
                        ff[0].header['IMGTYPE'] = 'lamp'
                else:
                    ff[0].header['IMGTYPE'] = 'Science'
            # Humidity
            if 'Inside_Rel_Hum' in ff[0].header:
                rel_hum = ff[0].header['Inside_Rel_Hum']
                ff[0].header['IN_HUM'] = rel_hum
            else:
                ff[0].header['IN_HUM'] = -1.
                logging.warning("No relative humidity")
            # Temperature
            if 'Inside_Air_Temp' in ff[0].header:
                in_temp = ff[0].header['Inside_Air_Temp']
                ff[0].header['IN_AIR'] = in_temp
            else:
                ff[0].header['IN_AIR'] = -1.
                logging.warning("No inside temperature")
            # Parallactic Angle
            if 'PRLLTC' in ff[0].header:
                tel_pa = ff[0].header['PRLLTC']
                ff[0].header['TEL_PA'] = tel_pa
            else:
                logging.warning("No telescope PA")
                ff[0].header['TEL_PA'] = -1.
            # Equinox
            ff[0].header['EQUINOX'] = 2000
            ff[0].header['DOMEST'] = 'Open'
            # Close
            ff.close()
