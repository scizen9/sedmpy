import astropy.io.fits as pf


def go(imfiles):
    for ifile in imfiles:
        try:
            FF = pf.open(ifile)
            header = FF[0].header
            header['FNAME'] = ifile.split('.gz')[0]
            if 'AIRMASS' not in header:
                header['AIRMASS'] = -1.
            if 'ADCSPEED' not in header:
                header['ADCSPEED'] = -1.
            if 'EXPTIME' not in header:
                header['EXPTIME'] = -1.
            if 'OBJECT' not in header:
                header['OBJECT'] = 'ERR!'
            try:
                print("%(FNAME)28s (%(AIRMASS)1.3f/%(ADCSPEED)1.1f/%(EXPTIME)s"
                      " s): %(OBJECT)-s" % header)
            except:
                print("%28s : ?" % ifile.split('.gz')[0])
        except FileNotFoundError:
            pass
        except OSError:
            pass


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        raise Exception("not enough arguments")

    go(sys.argv[2:])
