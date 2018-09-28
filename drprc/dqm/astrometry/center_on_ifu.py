import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import subprocess
import time


def solve_astrometry(img, radius=3.0, with_pix=True,
                     first_call=False, tweak=3):
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
    cmd = (" solve-field --no-fits2fits --ra %s --dec %s --radius "
           "%.4f -p --new-fits %s -W none -B none -P none -M none "
           "-R none -S none -t %d --overwrite %s " % (ra, dec, radius,
                                                      astro, tweak, img))

    if with_pix:
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.425 --"

    print(cmd)

    cmd = cmd + " > /tmp/astrometry_fail  2>/tmp/astrometry_fail"
    try:
        subprocess.call(cmd, shell=True, timeout=60)
    except Exception as e:
        print(str(e))
        print(time.time()-s)
        print("astrometry failed")

    # Cleaning after astrometry.net
    if os.path.isfile(img.replace(".fits", ".axy")):
        os.remove(img.replace(".fits", ".axy"))
    if os.path.isfile(img.replace(".fits", "-indx.xyls")):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if os.path.isfile("none"):
        os.remove("none")

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
                     use_closest=False, use_average=False):
    """
    
    :param raw_image: 
    :param overwrite: 
    :param use_brightest: 
    :param use_closest: 
    :param use_average: 
    :return: 
    """

    # 1. Check if the solved output file exist for the raw file, we
    # expect the files to be in the same directory.
    astro = os.path.join(os.path.dirname(raw_image),
                         "a_" + os.path.basename(raw_image))

    if overwrite:
        astro = solve_astrometry(raw_image)

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
    image2 = "/scr7/rsw/sedm/phot/20180813/rc20180813_07_24_09.fits"
    # image = "/data2/sedm/20180914/a_rc20180914_12_36_56.fits"
    #ret = solve_astrometry(image2)
    #print(ret)
    #print(get_offset_to_reference(ret))
    print(calculate_offset(image2, False))
    #print(get_offset_to_brightest(image2))
