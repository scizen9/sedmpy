import argparse
from stsci.image.numcombine import numCombine as nc
import astropy.io.fits as pf
import numpy as np
import time
import Version as Version

drp_ver = Version.ifu_drp_version()


def imcombine(flist, fout, listfile=None, combtype="mean",
              nlow=0, nhigh=0):

    """Convenience wrapper around STSCI python task numCombine

    Args:
        flist (list of str): The list of files to imcombine
        listfile (string): The file to store the list of images in
        fout (str): The full path to the output file
        combtype (str): median, mean, sum, minimum
        nlow (int): Number of low pixels to throw out in median calculation
        nhigh (int): Number of high pixels to throw out in median calculation
    
    Returns:
        None

    Side effects:
        Creates the imcombined file at location `out`

    """

    imstack = []
    for fl in flist:
        inhdu = pf.open(fl)
        img = inhdu[0].data
        img = img.astype(np.float32)
        imstack.append(img)

    hdr = inhdu[0].header

    result = nc(imstack, combinationType=combtype, nlow=nlow, nhigh=nhigh)

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
    parser = argparse.ArgumentParser(description=\
        """Imcombine.py performs:

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
    args = parser.parse_args()

    filelist = args.files
    out = args.outname
    if args.outname is None:
        print("Set --outname")

    imcombine(filelist, out, listfile=args.listfile, combtype=args.combtype,
              nlow=args.Nlo, nhigh=args.Nhi)


