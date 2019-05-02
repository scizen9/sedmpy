from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from matplotlib import pylab as plt
from astropy.io import fits
import os
import time
import numpy as np
import datetime
from astropy.stats import sigma_clipped_stats



def get_stats(image, data_extension=0, fwhm=5.0, threshold=2, sigma=3.0,
              header_extension=0, header_keyword='FOCPOS'):
    """
    
    :param image: 
    :return: 
    """
    # 1. Start by making sure the image exists
    start = time.time()
    if not os.path.exists(image):
        return -1*(time.time() - start), "Image not found"

    # 2. Try extracting data from the image
    try:
        hdu = fits.open(image)
        data = hdu[data_extension].data
        #data = data[100:900, 1200:2400]
        #print(data)
        focus_stage = hdu[header_extension].header[header_keyword]
        hdu.close()
    except Exception as e:
        return -1 * (time.time() - start), "Image not found"

    print("Time to open images:%ss" % float(time.time()-start))
    start=time.time()
    # 3. Now get the initial stats from the image

    mean, median, std = sigma_clipped_stats(data, sigma=sigma)
    print("Time to get stats:%ss" % float(time.time() - start))
    start = time.time()
    # 4. Find the sources
    daofind = IRAFStarFinder(fwhm=fwhm, threshold=threshold * std, sigma_radius=2.5, minsep_fwhm=5,
                             sharplo=0.5, sharphi=5.0, roundlo=0.0, roundhi=0.5, sky=None,
                             exclude_border=True)

    print("Time to setup star finder:%ss" % float(time.time() - start))
    start = time.time()
    sources = daofind(data-median)
    print(sources)
    print("Time to get sources:%ss" % float(time.time() - start))
    start = time.time()
    focus_list = []
    #print(sources)
    if sources:

    #    #print(np.mean(sources['sharpness']),np.mean(sources['roundness1']), np.mean(sources['roundness2']),np.mean(sources['npix']), np.mean(sources['sky']),np.mean(sources['peak']),np.mean(sources['flux']), np.mean(sources['mag']), focus_stage, image)
        #print(sources[0]['fwhm', 'xcentroid', 'ycentroid'], focus_stage)
        print(np.mean(sources['fwhm']), np.mean(sources['roundness']), np.mean(sources['pa']), focus_stage, image)
        return [np.mean(sources['fwhm']), np.std(sources['fwhm']), focus_stage]


def analyse_sex(fwhm_list, focpos_list, std_list, plot=True, interactive=False):
    """

    :param fwhm_list:
    :param focpos_list:
    :param plot:
    :param interactive:
    :return:
    """
    focpos = np.array(focpos_list)
    fwhms = np.array(fwhm_list)
    std_fwhm = np.array(std_list)

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
        with open("focus", "w") as f:
            f.write(str(focpos))
            f.write(str(fwhms))
        plt.errorbar(focpos, fwhms, yerr=std_fwhm, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("FWHM (arcsec)")
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join("focus_%s.png" % (datetime.datetime.utcnow()).strftime("%Y%m%d-%H:%M:%S")))
            plt.clf()
    return x[np.argmin(p(x))], coefs[0]


if __name__ == "__main__":
    import glob
    data_list = sorted(glob.glob("/scr/rsw/sedm/rc20190318_02*.fits"))
    fwhm = []
    focus = []
    error = []
    for i in data_list:
        y = get_stats(i)
        if y:
            fwhm.append(y[0])
            error.append(y[1])
            focus.append(y[2])

    print(analyse_sex(fwhm, focus, error))