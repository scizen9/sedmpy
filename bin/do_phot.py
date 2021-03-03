import glob
try:
    from target_mag import get_target_mag
except ImportError:
    from drprc.target_mag import get_target_mag

try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils

flist = glob.glob('rc*.fits')
flist.sort()

phot_zp = {'u': None, 'g': None, 'r': None, 'i': None}

for rf in flist:
    if fitsutils.get_par(rf, "ONTARGET"):
        targ_obj = fitsutils.get_par(rf, "OBJECT")
        targ_filt = targ_obj.split()[-1]
        targ_name = targ_obj.split()[0]
        print("\nGetting quick %s-band mag for %s in %s"
              % (targ_filt, targ_name, rf))
        target_mag, target_magerr, std_zp = get_target_mag(rf, zeropoint=phot_zp, verbose=True)
        print("Quick MAG = %.3f +- %.3f" % (target_mag, target_magerr))
        if std_zp is not None:
            print("Quick MAG_ZP: %.3f" % std_zp)
            if phot_zp[targ_filt] is None:
                phot_zp[targ_filt] = std_zp
