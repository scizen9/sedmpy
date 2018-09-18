import astropy.io.fits as pf


def go(imfiles):
    for ifile in imfiles:
        FF = pf.open(ifile)
        header = FF[0].header
        header['FNAME'] = ifile
        try:
            print("%(FNAME)28s (%(AIRMASS)1.3f/%(ADCSPEED)1.1f/%(EXPTIME)s s):"
                  " %(OBJECT)-s" % header)
        except:
            print("%28s : ?" % ifile)


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        raise Exception("not enough arguments")

    go(sys.argv[2:])
