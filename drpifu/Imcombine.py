import argparse
from stsci.image.numcombine import numCombine
import astropy.io.fits as pf
import numpy as np
import time
import sedmpy_version

drp_ver = sedmpy_version.__version__


def subtract_oscan(dat, header):
    """Subtract overscan from image"""
    # Get over/pre-scan regions
    try:
        # These should be in the Andor camera image headers
        ps_x0 = header['PSCANX0']
        ps_x1 = header['PSCANX1']
        os_x0 = header['OSCANX0']
        os_x1 = header['OSCANX1']
        andor = True
    except KeyError:
        # These are for the PIXIS camera images
        ps_x0 = 0
        ps_x1 = 0
        os_x0 = 2045
        os_x1 = 2048
        andor = False

    # Pre-scan value
    if ps_x1 > ps_x0:
        pscan = np.nanmedian(np.nanmedian(dat[:, :ps_x1], axis=1))
        header['PSCANVAL'] = (pscan, 'Pre-scan value')
    else:
        pscan = 0.
    # Over-scan value
    if os_x1 > os_x0:
        oscan = np.nanmedian(np.nanmedian(dat[:, os_x0:], axis=1))
        header['OSCANVAL'] = (oscan, 'Over-scan value')
    else:
        oscan = 0.

    if andor:
        if pscan < oscan:
            scan_value = pscan
        else:
            scan_value = oscan

        header['SCANVAL'] = (scan_value, 'Scan value subtracted')
        header['OSCANSUB'] = (True, 'Over-scan subtracted?')

    else:
        scan_value = 0.

    return scan_value


def imcombine(flist, fout, listfile=None, combtype="mean",
              nlow=0, nhigh=0, sub_oscan=False):

    """Convenience wrapper around STSCI python task numCombine

    Args:
        flist (list of str): The list of files to imcombine
        listfile (string): The file to store the list of images in
        fout (str): The full path to the output file
        combtype (str): median, mean, sum, minimum
        nlow (int): Number of low pixels to throw out in median calculation
        nhigh (int): Number of high pixels to throw out in median calculation
        sub_oscan (bool): If True, subtract overscan before combining
    
    Returns:
        None

    Side effects:
        Creates the imcombined file at location `out`

    """

    imstack = []
    hdr = {}
    for fl in flist:
        inhdu = pf.open(fl)
        img = inhdu[0].data
        img = img.astype(np.float32)
        if sub_oscan:
            scan_val = subtract_oscan(img, inhdu[0].header)
            if scan_val > 0.:
                img -= scan_val
        imstack.append(img)

        hdr = inhdu[0].header

    result = numCombine(imstack, combinationType=combtype,
                        nlow=nlow, nhigh=nhigh)

    oimg = result.combArrObj

    ncom = 1
    for fl in flist:
        key = "IMCMB%03d" % ncom
        ncom += 1
        hdr[key] = fl

    hdr['NCOMBINE'] = (len(flist), 'number of images combined')
    hdr['COMBTYPE'] = (combtype, 'type of combine')
    hdr['COMBNLO'] = (nlow, 'number of low pixels to reject')
    hdr['COMBNHI'] = (nhigh, 'number of high pixels to reject')
    hdr.add_history('SEDMr.Imcombine run on %s' % time.strftime("%c"))
    hdr['DRPVER'] = drp_ver

    try:
        pf.writeto(fout, oimg, hdr)
    except OSError:
        print("%s already exists!" % fout)

    if listfile is None:
        path = "imcombine.lst"
    else:
        path = listfile
    f = open(path, "w")
    for file in flist:
        f.write(file + "\n")
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Imcombine.py performs:

        1) Median combination
        2) Mean combine
        3) Mean combine w/ sigma clipping

    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--files', type=str, nargs='*', default=[])
    parser.add_argument('--listfile', type=str, default=None)
    parser.add_argument('--Nhi', type=int, default=0)
    parser.add_argument('--Nlo', type=int, default=0)
    parser.add_argument('--combtype', type=str, default='mean')
    parser.add_argument('--reject', type=str, default='none')
    parser.add_argument('--outname', type=str, default=None)
    parser.add_argument('--sub_oscan', action="store_true", default=False)
    args = parser.parse_args()

    filelist = args.files
    out = args.outname
    if args.outname is None:
        print("Set --outname")

    imcombine(filelist, out, listfile=args.listfile, combtype=args.combtype,
              nlow=args.Nlo, nhigh=args.Nhi, sub_oscan=args.sub_oscan)
