import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy import units as u
import photutils
import numpy as np
import subprocess
import time


def solve_astrometry(img, radius=2.5, with_pix=True,
                     first_call=False, tweak=3, make_plots=False):
    """
    
    :param img: 
    :param radius: 
    :param with_pix: 
    :param first_call: 
    :param tweak: 
    :return: 
    """
    s = time.time()
    # 1. Get the needed header information for the solve field command
    image_header = fits.getheader(img)
    try:
        ra, dec = image_header['OBJRA'], image_header['OBJDEC']
    except Exception as e:
        print(str(e))
        ra, dec = image_header['RA'], image_header['DEC']

    # 2. Create the solved astronomy base name
    astro = os.path.join(os.path.dirname(img), "a_" + os.path.basename(img))

    print("Solving astrometry on field with (ra,dec)=",
          ra, dec, "Image", img, "New image", astro)

    # 3. Create the base solve-field command
    if make_plots:
        cmd = (" solve-field --no-fits2fits --ra %s --dec %s --radius "
               "%.4f -t %d --overwrite %s "
               "" % (ra, dec, radius, tweak, img))
    else:
        cmd = (" solve-field --no-fits2fits --ra %s --dec %s --radius "
           "%.4f -p --new-fits %s -W none -B none -M none "
           "-R none -S none -t %d --overwrite %s --parity neg "
           "" % (ra, dec, radius, astro, tweak, img))

    if with_pix:
        cmd = cmd + " --scale-units arcsecperpix --scale-low 0.355 --scale-high 0.400 --"

    print(cmd)

    cmd = cmd + " > /tmp/astrometry_fail  2>/tmp/astrometry_fail"
    try:
        subprocess.call(cmd, shell=True, timeout=120)
    except Exception as e:
        print(str(e))
        print(time.time()-s)
        print("astrometry failed")

    # Cleaning after astrometry.net
    """if os.path.isfile(img.replace(".fits", ".axy")):
        os.remove(img.replace(".fits", ".axy"))
    if os.path.isfile(img.replace(".fits", "-indx.xyls")):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if os.path.isfile("none"):
        os.remove("none")"""

    # This is legacy as it was decided that it is very unlikely that the
    # pointing will be so far out the we need to redo the search with an
    # 8degree radius
    if not os.path.isfile(astro) and first_call:
        print("Astrometry failed on file %s! Trying with a larger radius..."
              % img)

        #solve_astrometry(img, radius=8, with_pix=True, first_call=False)

        if not os.path.isfile(astro):
            print("Astrometry FAILED!")

    # We return the basename of the file even though there is a chance it
    # doesn't exist.
    print(time.time()-s)
    return astro

def get_offset_to_brightest2(image):
    """

    :param image:
    :return:
    """

    image_data = fits.getdata(image)
    data = image_data[1293 - 500:1293 + 500, 1280 - 500:1280 + 500]
    mean, median, std = sigma_clipped_stats(data, sigma=2.0, iters=5)
    daofind = photutils.DAOStarFinder(fwhm=3.0, threshold=5. * std)
    sources = daofind(data - median)
    for i in sources:
        print(i)
    print(sources)


def make_cutouts(image):
    """

    :param image:
    :return:
    """


    image_data = fits.getdata(image)
    image_header = fits.getheader(image)
    wcs = WCS(image_header)

    rdata_center = wcs.all_pix2world(1561.25, 1538.45, 0)
    print(rdata_center)
#    new_center = wcs.all_world2pix()
#    ref_ra = rdata_center[0]
#    ref_dec = rdata_center[1]

    objCoords = SkyCoord(169.90527008, 61.920226159, unit=(u.deg, u.deg), frame='icrs')

    refCoords = SkyCoord(170, 62, unit=(u.deg, u.deg), frame='icrs')

    dra, ddec = objCoords.spherical_offsets_to(refCoords)

    print(rdata_center)
    print(dra.arcsec, ddec.arcsec)

    rdata = image_data[1079:2048, 1030:2048]

    hdu = fits.PrimaryHDU(rdata, uint=False)
    hdu.scale('int16', bzero=32768)

    hdu.writeto('rdata.fits', output_verify="fix", )





def get_offset_to_brightest(image):
    """
    Get the offset to the brightest object near the reference pixel
    
    :param image: 
    :return: 
    """

    image_header = fits.getheader(image)
    image_data = fits.getdata(image)
    obj_ra, obj_dec = image_header['OBJRA'], image_header['OBJDEC']

    if not isinstance(obj_ra, float) and not isinstance(obj_dec, float):
        objCoords = SkyCoord(obj_ra, obj_dec, unit=(u.hour, u.deg), frame='icrs')
    else:
        objCoords = SkyCoord(obj_ra, obj_dec, unit=(u.deg, u.deg), frame='icrs')

    wcs = WCS(image_header)
    pra, pdec = wcs.all_world2pix(objCoords.ra.degree, objCoords.dec.degree, 0)

    imageloc = image_data[1293 - 50:1293 + 50, 1280 - 50:1280 + 50]

    nx = 50
    ny = 50

    def_x = np.argmax(np.sum(imageloc, axis=0))
    def_y = np.argmax(np.sum(imageloc, axis=1))
    print(def_x, def_y)
    newx = pra - nx / 2. + def_x
    newy = pdec - ny / 2. + def_y

    pra, pdec = wcs.all_pix2world(newx, newy, 0)
    refCoords = SkyCoord(pra, pdec, unit='deg', frame='icrs')
    dra, ddec = refCoords.spherical_offsets_to(objCoords)
    print(newx, newy)
    return 1, dra.arcsec/2, ddec.arcsec/2


def get_offset_to_reference(image, get_ref_from_header=False,
                            reference=(1293, 1280)):
    """
    Given an image with solved WCS. Compute the offset to the 
    reference pixel 
    :param image: 
    :param get_ref_from_header: 
    :param reference: 
    :return: 
    """

    # 1. Start by checking if the image exists
    if not os.path.exists(image):
        print("File %s does not exist! Returning Zero offsets..." % image)
        return -1, 0, 0

    # 2. Get the object coordinates from header and convert to degrees when
    # needed.
    image_header = fits.getheader(image)
    obj_ra, obj_dec = image_header['OBJRA'], image_header['OBJDEC']

    if not isinstance(obj_ra, float) and not isinstance(obj_dec, float):
        objCoords = SkyCoord(obj_ra, obj_dec, unit=(u.hour, u.deg), frame='icrs')
    else:
        objCoords = SkyCoord(obj_ra, obj_dec, unit=(u.deg, u.deg), frame='icrs')

    # 3. Get the WCS reference pixel position
    wcs = WCS(fits.getheader(image))

    if get_ref_from_header:
        x, y = image_header['crpix1'], image_header['crpix2']
    else:
        x, y = reference

    ref_ra, ref_dec = wcs.all_pix2world(x, y, 0)

    # 4. Get the ra and dec separation between the object and reference pixel
    refCoords = SkyCoord(ref_ra, ref_dec, unit='deg', frame='icrs')

    dra, ddec = objCoords.spherical_offsets_to(refCoords)

    return 0, -1*dra.arcsec, -1*ddec.arcsec


def find_closest_offsets(az, el):
    """
    Find the closest offset value of previously solved astrometry
    :param az: 
    :param el: 
    :return: 
    """

    return 3, 0, 0


def find_average_offsets(obsdate):
    """
    
    :param obsdate: 
    :return: 
    """

    return 4, 0, 0


def calculate_offset(raw_image, overwrite=True, use_brightest=True,
                     use_closest=False, use_average=False,
                     parse_directory_from_file=True,
                     base_dir="/data2/sedm/", make_plots=False):
    """
    
    :param raw_image: 
    :param overwrite: 
    :param use_brightest: 
    :param use_closest: 
    :param use_average: 
    :return: 
    """

    if parse_directory_from_file:
        raw_image = raw_image.split('/')[-2:]
        raw_image = os.path.join(base_dir, raw_image[0], raw_image[1])

    # 1. Check if the solved output file exist for the raw file, we
    # expect the files to be in the same directory.
    astro = os.path.join(os.path.dirname(raw_image),
                         "a_" + os.path.basename(raw_image))

    if overwrite:
        astro = solve_astrometry(raw_image, make_plots=make_plots)

    # 2. Now check if the solved file exist and solve the offset
    if os.path.exists(astro):
        return get_offset_to_reference(astro)

    # 3. If it doesn't exist then try the next method of return
    if use_brightest:
        return get_offset_to_brightest(raw_image)

    if use_closest:
        return find_closest_offsets(1, 1)

    if use_average:
        return find_average_offsets("20180928")


if __name__ == "__main__":

    t = calculate_offset('/data2/sedm/20190427/rc20190427_05_48_48.fits', parse_directory_from_file=False, make_plots=True)

    """import glob
    files = glob.glob('/data1/sedm/astrom/rc*.fits')
    data = open('astrometry_solve.txt', 'w')
    for i in files:
        t = calculate_offset(i, parse_directory_from_file=False, make_plots=True)
        print(t)
        try:
            st = ''
            for k in t:
                st += str(k)+','
        except:
            t = '-999,999,999'
        head = fits.getheader(i)
        airmass = head['AIRMASS']
        az = head['TEL_AZ']
        el = head['TEL_EL']
        ha = head['TEL_HA']
        dec = head['TEL_DEC']
        data.write('%s,%s,%s,%s,%s,%s,%s\n' % (i, ha, dec, az, el, airmass, st))
    data.close()
    #print(ret)
    #print(get_offset_to_reference('/scr/rsw/sedm/data/raw/20181204/rc20181204_09_59_04.new'))
    #print(calculate_offset(image2, False))"""

