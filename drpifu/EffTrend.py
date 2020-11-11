"""Plot efficiency trend for SEDM

"""
import glob
import os
import csv

import pylab as pl
import numpy as np
from astropy.time import Time
import astropy.io.fits as pf

sdir = '/scr2/sedmdrp/redux'
fspec = os.path.join(sdir, '20??????')
dlist = sorted([d for d in glob.glob(fspec) if os.path.isdir(d)])[1:]

recent = False  # do the whole data set

area = 18000.0  # P60 area in cm^2
refl = 0.82     # P60 reflectance fraction
da = []
ef1 = []
ef2 = []
ef3 = []
ef4 = []
ef5 = []
efs = []
for d in dlist:
    ddate = d.split('/')[-1]
    dtime = Time(ddate[0:4]+'-'+ddate[4:6]+'-'+ddate[6:])
    ddate_int = int(ddate)

    # Recent data
    if recent and (ddate_int < 20190101):
        continue

    print(d)

    fspec = os.path.join(d, 'spec_aperture_*_STD-*_ea.fits')
    slist = glob.glob(fspec)

    for s in slist:
        ff = pf.open(s)
        hdr = ff[0].header
        ea = ff[0].data
        # Avoid bad quality
        if hdr['QUALITY'] != 0:
            print("Bad quality: %s" % s)
            continue
        # Avoid edge objects
        xx = ff[0].header['XPOS']
        yy = ff[0].header['YPOS']
        if abs(xx) > 10 or abs(yy) > 10:
            print("Close to the edge: %s" % s)
            continue

        # Calculate efficiency
        ef = 100. * ea / (18000.0 * 0.82)
        # Calculate wavelengths
        lam0 = hdr['CRVAL1']
        dlam = hdr['CDELT1']
        wl = (1 + np.arange(len(ea))) * dlam + lam0
        wl = wl / 10.   # Convert to nm

        # Check each 100 nm bin

        # 400 - 500 nm
        vec = ef[(wl > 400) * (wl < 500)]
        if len(vec) < 1:
            continue
        e1 = np.nanmean(vec)
        if e1 > 100 or e1 < 0:
            continue

        # 500 - 600 nm
        vec = ef[(wl > 500) * (wl < 600)]
        if len(vec) < 1:
            continue
        e2 = np.nanmean(vec)
        if e2 > 100 or e2 < 0:
            continue

        # 600 - 700 nm
        vec = ef[(wl > 600) * (wl < 700)]
        if len(vec) < 1:
            continue
        e3 = np.nanmean(vec)
        if e3 > 100 or e3 < 0:
            continue

        # 700 - 800 nm
        vec = ef[(wl > 700) * (wl < 800)]
        if len(vec) < 1:
            continue
        e4 = np.nanmean(vec)
        if e4 > 100 or e4 < 0:
            continue

        # 800 - 900 nm
        vec = ef[(wl > 800) * (wl < 900)]
        if len(vec) < 1:
            continue
        e5 = np.nanmean(vec)
        if e5 > 100 or e5 < 0:
            continue

        # Fill efficiency vectors
        ef1.append(e1)
        ef2.append(e2)
        ef3.append(e3)
        ef4.append(e4)
        ef5.append(e5)
        efs.append(s)

        da.append(ddate[0:4]+'-'+ddate[4:6]+'-'+ddate[6:])

t = Time(da)
pl.plot_date(t.plot_date, ef1, '^', linestyle='None', markersize=2.0,
             label='400-500 nm')
pl.plot_date(t.plot_date, ef2, 'v', linestyle='None', markersize=2.0,
             label='500-600 nm')
pl.plot_date(t.plot_date, ef3, 'x', linestyle='None', markersize=2.0,
             label='600-700 nm')
pl.plot_date(t.plot_date, ef4, 'D', linestyle='None', markersize=2.0,
             label='700-800 nm')
pl.plot_date(t.plot_date, ef5, 'o', linestyle='None', markersize=2.0,
             label='800-900 nm')
pl.gcf().autofmt_xdate()
pl.xlabel('Date (UTC)')
pl.ylabel('Efficiency (%)')
pl.title('Efficiency Trend')
pl.legend(loc=2)
pl.grid(True)
pl.ylim(-1, 45)
if recent:
    ofil = os.path.join(sdir, 'SEDM_eff_recent_pysedm.pdf')
    tfil = os.path.join(sdir, 'SEDM_eff_recent_pysedm.txt')
else:
    ofil = os.path.join(sdir, 'SEDM_eff_trend_pysedm.pdf')
    tfil = os.path.join(sdir, 'SEDM_eff_trend_pysedm.txt')
pl.savefig(ofil)
with open(tfil, 'w') as dfil:
    dat_writer = csv.writer(dfil, delimiter=" ", quoting=csv.QUOTE_MINIMAL)
    dat_writer.writerows(zip(da, ef1, ef2, ef3, ef4, ef5, efs))


