import os
import subprocess
import shutil
from astropy.io.ascii import SExtractor
from astropy.io import fits as pf
import numpy as np
from matplotlib import pylab as plt
import datetime

SITE_ROOT = os.path.realpath(os.path.dirname(__file__))
config_dir = os.path.join(SITE_ROOT, 'sexConfig')

root_dir = '/data2/sedm/'
print(config_dir)


def run_sex(image, outdir="", default_to_image_dir=True, overwrite=True):
    """

    :param image:
    :param mask:
    :param cosmics:
    :param overwrite:
    :return:
    """

    if not os.path.exists(image):
        return False

    if not outdir and default_to_image_dir:
        outdir = os.path.dirname(image)
        outdir = os.path.join(outdir, 'sextractor')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    sexcat = os.path.join(outdir,
                          os.path.basename(image).replace(".fits", ".sex"))

    if not overwrite and os.path.exists(sexcat):
        print('Sextractor already run on %s' % image)
        return sexcat

    cmd = "sex -c %s/daofind.sex %s" % (config_dir, image)
    subprocess.call(cmd, stdout=subprocess.DEVNULL, shell=True)

    shutil.move("image.sex", sexcat)

    return sexcat


def analyze_sex_by_pandas(catalog, filter_data=True, return_type='dict'):
    """

    :param catalog:
    :param return_type:
    :return:
    """

    sex = SExtractor()
    data = sex.read(catalog)
    h = data.to_pandas()
    if filter_data:
        df = h.loc[h['FLAGS'] == 0]
        df = df[(df['X_IMAGE'] < 885) | (df['X_IMAGE'] > 1540)]
        df = df[(df['Y_IMAGE'] < 850) | (df['Y_IMAGE'] > 1125)]
        filterd_df = df[(df['ELLIPTICITY'] < .3)]

        return filterd_df


def analyse_sex(sexfileslist, plot=True, interactive=False):
    '''
    Analyses the sextractor filelist to determine the best focus for the RC camera.

	#   1 X_IMAGE                Object position along x                                    [pixel]
	#   2 Y_IMAGE                Object position along y                                    [pixel]
	#   3 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
	#   4 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
	#   5 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
	#   6 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
	#   7 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
	#   8 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
	#   9 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
	#  10 BACKGROUND             Background at centroid position                            [count]
	#  11 FLAGS                  Extraction flags
	#  12 A_IMAGE                Isophotal image mejor axis
	#  13 B_IMAGE                Isophotal image minor axis
   	#  14 THETA_IMAGE            Isophotal image position angle
   	#  15 PETRO_RADIUS           Petrosian radius

    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    '''
    sexfileslist = list(sexfileslist)
    sexfileslist.sort()

    focpos = []
    fwhms = []
    std_fwhm = []
    for i, f in enumerate(sexfileslist):
        fits = f.replace("sextractor/", "").replace(".sex", ".fits")
        FF = pf.open(fits)
        pos = float(FF[0].header['focpos'])

        s = np.genfromtxt(f, comments="#")
        print(f, "Initial number of sources", len(s))

        s = s[s[:, 1] < 2000]
        s = s[s[:, 10] != 4]

        # Only objects with FWHM less than 40 pixels... but larger than 2
        s = s[s[:, 7] < 60]
        s = s[s[:, 7] > 1]

        # Select bright magnitudes
        s = s[s[:, 4] < np.percentile(s[:, 4], 30)]
        # Select round sources (ellipticity is 1-axis_ratio)
        s = s[s[:, 8] < np.percentile(s[:, 8], 30)]
        print(f, "number of sources", len(s))

        focpos.append(pos)
        fwhms.append(np.nanmean(s[:, 7] * 0.394))
        mad = np.median(np.abs(s[:, 7] - np.nanmean(s[:, 7]))) / 0.67448975019608171 * 0.394
        # std_fwhm.append(np.std(s[:,6]*0.394))
        std_fwhm.append(mad)

    focpos = np.array(focpos)
    fwhms = np.array(fwhms)
    std_fwhm = np.array(std_fwhm)

    n = len(fwhms)

    best_seeing_id = np.nanargmin(fwhms)
    # We will take 4 datapoints on the left and right of the best value.
    selected_ids = np.arange(-4, 5, 1)
    selected_ids = selected_ids + best_seeing_id
    selected_ids = np.minimum(selected_ids, n - 1)
    selected_ids = np.maximum(selected_ids, 0)
    print("FWHMS: %s, focpos: %s, Best seeing id: %d. Selected ids %s" % (fwhms, focpos, best_seeing_id, selected_ids))
    selected_ids = np.array(list(set(selected_ids)))

    focpos = focpos[selected_ids]
    fwhms = fwhms[selected_ids]
    std_fwhm = std_fwhm[selected_ids]

    std_fwhm = np.maximum(1e-5, np.array(std_fwhm))

    coefs = np.polyfit(focpos, fwhms, w=1 / std_fwhm, deg=2)

    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    print("Best focus:%.2f" % x[np.argmin(p(x))], coefs, std_fwhm)

    if (plot == True):
        plt.title("Best focus:%.2f" % x[np.argmin(p(x))])
        with open(os.path.join("focus"), "w") as f:
            f.write(str(focpos))
            f.write(str(fwhms))
        plt.errorbar(focpos, fwhms, yerr=std_fwhm, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("FWHM (arcsec)")
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_%s.png" % (datetime.datetime.utcnow()).strftime("%Y%m%d-%H:%M:%S")))
            plt.clf()
    return x[np.argmin(p(x))], coefs[0]


def analyse_sex_ifu(sexfileslist, plot=True, interactive=False, debug=False, ifu=False):
    '''
    Analyses the sextractor filelist to determine the best focus for the IFU camera stage controller.

	#   1 X_IMAGE                Object position along x                                    [pixel]
	#   2 Y_IMAGE                Object position along y                                    [pixel]
	#   3 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
	#   4 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
	#   5 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
	#   6 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
	#   7 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
	#   8 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
	#   9 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
	#  10 BACKGROUND             Background at centroid position                            [count]
	#  11 FLAGS                  Extraction flags
	#  12 A_IMAGE                Isophotal image mejor axis
	#  13 B_IMAGE                Isophotal image minor axis
   	#  14 THETA_IMAGE            Isophotal image position angle
   	#  15 PETRO_RADIUS           Petrosian radius

    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    '''
    sexfileslist = list(sexfileslist)
    sexfileslist.sort()

    focpos = []
    mags = []
    std_mags = []
    mags98perc = []

    for i, f in enumerate(sexfileslist):
        print(f)
        fits = f.replace("sextractor/", "").replace(".sex", ".fits")
        FF = pf.open(fits)

        pos = float(FF[0].header['ifufoc2'])

        s = np.genfromtxt(f, comments="#", dtype=[("x", np.float), ("y", np.float), ("ra", np.float), ("dec", np.float), \
                                                  ("mag", np.float), ("magerr", np.float), ("fwhm_world", np.float),
                                                  ("fwhm_image", np.float), ("ellipticity", np.float), \
                                                  ("background", np.float), ("flags", np.float), ("a_image", np.float),
                                                  ("b_image", np.float), ("theta_image", np.float), \
                                                  ("petro_radius", np.float)])

        # for the focus purpose, we only want traces, meaning that their ellipticity needs to be large.
        sb = s[(s["x"] > 200) * (s["x"] < 1848) * (s["y"] > 200) * (s["y"] < 1848)]
        # sb = sb[sb["a_image"]> 30]

        if (debug):
            plt.figure(figsize=(10, 7))
            plt.suptitle('Focus %.2f' % pos, fontsize=10, horizontalalignment="right")
            ax = plt.subplot2grid((2, 2), (0, 0))
            magmap = ax.scatter(sb["x"], sb["y"], c=10 ** (-0.4 * sb["mag"]), s=10, edgecolors='face')
            plt.title("flux / spaxel")
            plt.colorbar(magmap, label="flux", format="%.1e")

            ax = plt.subplot2grid((2, 2), (0, 1))
            ax.hist(sb["a_image"], bins=30, range=(0, 60))
            plt.title("A image")

            ax = plt.subplot2grid((2, 2), (1, 0))
            ax.hist(sb["b_image"], bins=20, range=(0, 5))
            plt.title("B image")

            ax = plt.subplot2grid((2, 2), (1, 1))
            ax.hist(sb["mag"], bins=25, range=(-15, -4))
            plt.title("mags image")

            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]), os.path.basename(f).replace(".sex", ".png")))
            plt.clf()

        focpos.append(pos)
        avex = np.average(sb["x"])
        avey = np.average(sb["y"])
        mags.append(np.std(np.sqrt((sb["x"] - avex) ** 2 + (sb["y"] - avey) ** 2)))
        # mags.append(np.percentile(sb["mag"], 25))
        print()
        std_mags.append(np.std(sb["mag"] < np.percentile(sb["mag"], 25)))

    focpos = np.array(focpos)
    mags = np.array(mags)
    std_mags = np.maximum(1e-5, np.array(std_mags))

    coefs = np.polyfit(focpos, mags, w=1 / std_mags, deg=2)

    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    print("Best focus:%.2f" % x[np.argmax(p(x))], coefs, std_mags)

    if (plot == True):
        plt.figure()
        plt.title("Best focus IFU:%.2f" % x[np.argmax(p(x))])
        with open(os.path.join("focus"), "w") as f:
            f.write(str(focpos))
            f.write(str(mags))
        plt.errorbar(focpos, mags, yerr=std_mags, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("Mag brightest spaxel")
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_ifu_%s.png" % (datetime.datetime.utcnow()).strftime("%Y%m%d-%H:%M:%S")))
            plt.clf()
        plt.close("all")

    return x[np.argmax(p(x))], coefs[0]


def analyse_sex_ifu_spectrograph(sexfileslist, plot=True, interactive=False, debug=False):
    '''
    Analyses the sextractor filelist to determine the best focus for the IFU spectrograph instrument.

	#   1 X_IMAGE                Object position along x                                    [pixel]
	#   2 Y_IMAGE                Object position along y                                    [pixel]
	#   3 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
	#   4 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
	#   5 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
	#   6 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
	#   7 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
	#   8 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
	#   9 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
	#  10 BACKGROUND             Background at centroid position                            [count]
	#  11 FLAGS                  Extraction flags
	#  12 A_IMAGE                Isophotal image mejor axis
	#  13 B_IMAGE                Isophotal image minor axis
   	#  14 THETA_IMAGE            Isophotal image position angle
   	#  15 PETRO_RADIUS           Petrosian radius

    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    '''
    sexfileslist = list(sexfileslist)
    sexfileslist.sort()

    focpos = []
    fwhms = []
    std_fwhm = []
    mags = []
    mags_std = []
    trace_width = []
    trace_width_std = []

    for i, f in enumerate(sexfileslist):
        print(f)
        fits = f.replace("sextractor/", "").replace(".sex", ".fits")
        FF = pf.open(fits)

        pos = float(FF[0].header['ifufoc2'])

        s = np.genfromtxt(f, comments="#", dtype=[("x", np.float), ("y", np.float), ("ra", np.float), ("dec", np.float), \
                                                  ("mag", np.float), ("magerr", np.float), ("fwhm_world", np.float),
                                                  ("fwhm_image", np.float), ("ellipticity", np.float), \
                                                  ("background", np.float), ("flags", np.float), ("a_image", np.float),
                                                  ("b_image", np.float), ("theta_image", np.float) \
                                                  ])

        sb = s[(s["x"] > 200) * (s["x"] < 1848) * (s["y"] > 200) * (s["y"] < 1848)]

        sb = sb[np.abs(sb["theta_image"]) < 10]
        sb = sb[np.abs(sb["a_image"]) > 20]

        if (debug):
            plt.figure(figsize=(10, 7))
            plt.suptitle('Focus %.2f' % pos, fontsize=10, horizontalalignment="right")
            ax = plt.subplot2grid((2, 3), (0, 0))
            magmap = ax.scatter(sb["x"], sb["y"], c=10 ** (-0.4 * sb["mag"]), s=5, edgecolors='face')
            plt.title("flux / spaxel")
            plt.colorbar(magmap, label="flux", format="%.2e")

            ax = plt.subplot2grid((2, 3), (0, 1))
            ax.hist(sb["a_image"], bins=20, range=(20, 80))
            plt.title("A image")

            ax = plt.subplot2grid((2, 3), (0, 2))
            ax.hist(sb["ellipticity"], bins=50, range=(0.9, 1))
            plt.title("Ellipticity")

            ax = plt.subplot2grid((2, 3), (1, 0))
            ax.hist(sb["b_image"], bins=20, range=(0.8, 2))
            plt.title("B image")

            ax = plt.subplot2grid((2, 3), (1, 1))
            ax.hist(sb["theta_image"], bins=20, range=(-1, 1))
            plt.title("Theta (A/B)")

            ax = plt.subplot2grid((2, 3), (1, 2))
            ax.hist(sb["mag"], bins=20)
            plt.title("mag")

            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]), os.path.basename(f).replace(".sex", ".png")))
            plt.clf()

        focpos.append(pos)
        fwhms.append(np.median(sb["a_image"]))
        mags.append(np.average(sb["mag"]))
        mags_std.append(np.std(sb["mag"]))
        trace_width.append(np.median(sb["b_image"]))
        trace_width_std.append(np.minimum(0.3, np.std(sb["b_image"])))

        # std_fwhm.append(np.std(sb["mag"] < np.percentile(sb["mag"], 15)))
        std_fwhm.append(np.std(sb["a_image"]))  # < np.percentile(sb["b_image"], 15)))

    focpos = np.array(focpos)
    fwhms = np.array(fwhms)
    std_fwhm = np.minimum(6.5, np.array(std_fwhm))
    print(focpos)
    trace_width = np.array(trace_width)
    trace_width_std = np.array(trace_width_std)

    mags_std = np.array(mags_std)

    coefs = np.polyfit(focpos, fwhms, w=1 / std_fwhm, deg=2)
    coefs_mag = np.polyfit(focpos, mags, w=1 / mags_std, deg=2)
    coefs_width = np.polyfit(focpos, trace_width, w=1 / trace_width_std, deg=2)

    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    pmag = np.poly1d(coefs_mag)
    pwidth = np.poly1d(coefs_width)

    print("Best focus:%.2f" % x[np.argmax(p(x))], coefs, std_fwhm)

    if (plot == True):
        plt.title("Best focus:%.2f" % x[np.argmax(p(x))])
        with open(os.path.join("focus"), "w") as f:
            f.write(str(focpos))
            f.write(str(fwhms))
        plt.errorbar(focpos, fwhms, yerr=std_fwhm, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("Median trace longitude [pixels]")
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_ifu_tracelength_%s.png" % (datetime.datetime.utcnow()).strftime(
                                         "%Y%m%d-%H:%M:%S")))
            plt.clf()

    if (plot == True):
        plt.title("Best focus:%.2f" % x[np.argmin(pmag(x))])
        plt.errorbar(focpos, mags, yerr=mags_std, fmt="o")
        plt.plot(x, pmag(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("Mag")
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),
                                     "focus_ifu_mag_%s.png" % (datetime.datetime.utcnow()).strftime("%Y%m%d-%H:%M:%S")))
            plt.clf()

    return x[np.argmin(pmag(x))], np.abs(coefs_mag[0])


def get_focus(lfiles, get_path_by_basename=True, plot=True, interactive=False,
              basedir='/data2/sedm/', obsdate="", foctype='rc'):
    '''
    Receives a list of focus files and returns the best focus value.

    '''
    out_files = []

    if get_path_by_basename:
        basename = os.path.basename(lfiles[0])
        basename = basename.split('_')
        obsdate = basename[0][-8:]
    print(basename)
    for l in lfiles:
        out_files.append(run_sex(os.path.join(basedir, obsdate, os.path.basename(l))))
    #focus, sigma = analyse_sex(out_files, plot=plot, interactive=interactive)
    focus, sigma = analyse_sex_ifu_spectrograph(out_files) #analyse_sex(out_files, plot=plot, interactive=interactive)
    return focus, sigma


if __name__ == "__main__":
    files = """ifu20190424_04_47_25.fits
ifu20190424_04_48_00.fits
ifu20190424_04_48_35.fits
ifu20190424_04_49_10.fits
ifu20190424_04_49_46.fits
ifu20190424_04_50_21.fits
ifu20190424_04_50_56.fits
ifu20190424_04_51_31.fits"""
    files = files.split('\n')
    print(get_focus(files))