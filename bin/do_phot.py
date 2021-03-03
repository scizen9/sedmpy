import glob
try:
    from target_mag import get_target_mag
except ImportError:
    from drprc.target_mag import get_target_mag

try:
    import fitsutils
except ImportError:
    import drprc.fitsutils as fitsutils

phot_zp = None
flist = glob.glob('rc*.fits')
for rf in flist:
    if fitsutils.get_par(rf, "ONTARGET"):
        print("\nGetting target mag for %s" % rf)
        target_mag, target_magerr, std_zp = get_target_mag(rf, r_zeropoint=phot_zp, verbose=True)
        print("Quick MAG = %.3f +- %.3f" % (target_mag, target_magerr))
        if std_zp is not None:
            print("Quick MAG_ZP: %.3f" % std_zp)
            if phot_zp is None:
                phot_zp = std_zp