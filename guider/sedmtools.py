# IMPORTS
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.optimize import fmin
import os
from calendar import timegm

# This method uses images in 'focus_dir' to find the best focus for the
# Telescope The directory 'focus_dir' should contain a series of images at
# different secondary focus positions The code then runs through, gets the FWHM
# of each, fits a quadratic function and minimizes FWHM for the best focus
# position


def focusfind(focus_dir):

    # Get image list
    images = getImgList(focus_dir)
    focuses = []
    metrics = []

    # Get star coordinates
    c = []
    for x in images:
        c1 = findstar(x)
        if list(c1) != [-1, -1]:
            c.append(c1)
            break

    if len(c) == 0:
        print("No stars.")
        return
    print(c)
    s = 20

    # Run through images and add fwhm & sec focus to lists
    for im in images:
        fwhm_avg = 0
        count = 0
        for x, y in c:

            # Extract star
            starBox = im[0].data[x-s:x+s, y-s:y+s]

            # Background subtract
            median = np.median(starBox)
            for a in range(starBox.shape[0]):
                for b in range(starBox.shape[1]):
                    starBox[a, b] -= median

            # Fit 2D gaussian and deduce FWHM
            try:
                sig_x, sig_y = fitgaussian(starBox)[3:5]
            except:
                print("Could not fit Gaussian")
                continue

            fwhm_x = 2*np.sqrt(-2*(abs(sig_x**2))*np.log(0.5))
            fwhm_y = 2*np.sqrt(-2*(abs(sig_y**2))*np.log(0.5))

            fwhm_avg += 0.395*(fwhm_x+fwhm_y)/2
            count += 1

        if count == 0:
            print("Image analysis failed.")
            continue

        fwhm_avg /= count
        metrics.append(fwhm_avg)
        focuses.append(im[0].header['FOCPOS'])

    # Fit Gaussian to lists
    metrics = np.array(metrics)
    focuses = np.array(focuses)
    a, b, c = curve_fit(quadratic, focuses, metrics)[0]

    # Minimize Gaussian
    foc = focuses[np.argmin(metrics)]
    minFoc = fmin(quadratic, foc, [a, b, c])[0]

    # Plot data
    smooth = np.arange(focuses[np.argmin(focuses)],
                       focuses[np.argmax(focuses)],
                       0.005)
    plt.figure()
    plt.plot(focuses, metrics, 'bo')
    plt.plot(minFoc, quadratic(minFoc, a, b, c), 'ro')
    plt.plot(smooth, quadratic(smooth, a, b, c), 'g-')
    plt.annotate("Optimised Focus: %.4fmm" % (minFoc),
                 xy=(minFoc, quadratic(minFoc, a, b, c)),
                 xytext=(minFoc, quadratic(minFoc, a, b, c)+3),
                 arrowprops=dict(arrowstyle='->', shrinkA=0))
    plt.ylabel('FWHM(px)')
    plt.xlabel('Secondary Focus Position (mm)')

    return minFoc


# This method convolves two 2D matrices and uses
# the correlation matrix to return the offset
def getshift(img1, img2, n=1):

    print("These are the images being compared")
    print(img1, img2)
    print("******************")
    if type(img1) == pf.hdu.hdulist.HDUList:
        img1 = img1[0].data
    if type(img2) == pf.hdu.hdulist.HDUList:
        img2 = img2[0].data
    x = img1.shape[0]
    y = img1.shape[1]

    conv = convolve(img1[0:x/n, 0:y/n], img2[0:x/n, 0:y/n])
    peak = conv[np.unravel_index(np.nanargmax(conv), conv.shape)]
    for a in range(conv.shape[0]):
        for b in range(conv.shape[1]):
            if conv[a, b] < peak/2:
                conv[a, b] = 0
    params = fitgaussian(conv)

    # Convert to RA and Dec
    return -params[2]*0.395, -params[1]*0.395


# Returns the Full-width-half-max of a star (given or auto-found) in an image
# Works by fitting a 2D Gaussian function to the data and using
# the beam width to calculate FWHM
def getfwhm(img, c=(-1, -1)):

    # Check type and make sure to convert to np.array
    # if type(img) == pf.hdu.hdulist.HDUList: img = img[0].data
    # elif type(img) == str: img = pf.open(img, ignore_missing_end=True)[0].data
    img = img[0].data
    # No star coords given, try find star, return -1 if failed
    if c == (-1, -1):
        c = findstar(img)
        if c == (-1, -1):
            return -1

    # If star coords given, extract star & fit 2d gaussian
    else:
        print("These are the coords that are being used for getfwhm")
        print(c)
        print("*******************")
        s = 10
        x, y = c
        starBox = img[x-s:x+s+1, y-s:y+s+1]

        # Background subtract
        median = np.median(starBox)
        for a in range(starBox.shape[0]):
            for b in range(starBox.shape[1]):
                starBox[a, b] -= median

        # Fit 2D gaussian and deduce FWHM
        try:
            sig_x, sig_y = fitgaussian(starBox)[3:5]
        except:
            print("Could not fit Gaussian")
            return -1
        fwhm_x = 2*np.sqrt(-2*(abs(sig_x)**2)*np.log(0.5))
        fwhm_y = 2*np.sqrt(-2*(abs(sig_y)**2)*np.log(0.5))

        li = [fwhm_x, fwhm_y, fwhm_x/fwhm_y, (fwhm_x+fwhm_y)/2]
        for i in range(len(li)):
            if i == 2:
                continue
            li[i] *= 0.395

        return li


# This method takes in an image and finds a useable star
# It works by roughly the following algorithm:
# Select bright but not saturated pixel
def findstar(img, m=25):

    cap = 40000
    print(type(img))
    data = img[0].data
    # if type(img) == pf.hdu.hdulist.HDUList: data = img[0].data
    # elif type(img) == np.ndarray: data = img #If array
    # else:
    #    print("Unrecognized input type for findStar(img)")
    #    return (-1,-1)

    # Crop data within a certain margin, use deep copy
    dataCropped = np.empty_like(data[m:data.shape[0]-m, m:data.shape[1]-m])
    dataCropped[:] = data[m:data.shape[0]-m, m:data.shape[1]-m]

    # Keep taking brightest until px < 90% of Full Well Value
    found = False

    while not found:

        i, j = np.unravel_index(np.nanargmax(dataCropped), dataCropped.shape)

        import sys
        sys.stdout.flush()

        # Mask entire saturated area so that we don't select
        # any more of these saturated pixels
        if dataCropped[i, j] > cap:

            # Create a box around the original saturated pixel
            u, d, r, l = 1, 1, 1, 1
            box = None
            done = False
            count = 0

            # Keep expanding box till all connected
            # saturated pixels are enclosed
            while not done:
                count += 1
                if count > 100:
                    done = True

                # Expand box
                box = dataCropped[max(0, i-l):min(i+r+1, data.shape[0]-1),
                                  max(0, j-d):min(j+u+1, data.shape[1]-1)]

                # By default, we do not want to expand
                expand_up, expand_down, expand_left, expand_right = \
                    False, False, False, False

                # Run vertically through right- and left-most columns of box
                x_0 = 0
                x_e = box.shape[0]-1
                for y in range(box.shape[1]):
                    # If there is a saturated pixel in the left-most column,
                    # we need to keep expanding left
                    if box[x_0, y] > cap:
                        expand_left = True
                        l += 1
                    # If there is a saturated pixel in the right-most column,
                    # we need to keep expanding right
                    if box[x_e, y] > cap:
                        expand_right = True
                        r += 1

                # Run horizontally across top and bottom rows of box
                y_0 = 0
                y_e = box.shape[1]-1
                for x in range(box.shape[0]):
                    # If there is a saturated pixel in the top row,
                    # we need to keep expanding up
                    if box[x, y_0] > cap:
                        expand_down = True
                        d += 1
                    # If there is a saturated pixel in the bottom row,
                    # we need to keep expanding down
                    if box[x, y_e] > cap:
                        expand_up = True
                        u += 1

                # Check if box has reached edge of image on any side and
                # stop expanding if so
                if i-l <= 0:
                    expand_left = False
                if j-d <= 0:
                    expand_down = False
                if i+r+1 >= data.shape[0]-1:
                    expand_right = False
                if j+u+1 >= data.shape[1]-1:
                    expand_up = False

                # If we have finally enclosed all saturated pixels,
                # exit the loop
                if not(expand_up or expand_down or expand_left or expand_right):
                    done = True

            border = 20
            u += border
            d += border
            l += border
            r += border

            # Expand box again by 'border' amount, so we mask
            # the surrounding area too
            box = dataCropped[max(0, i-l):min(i+r+1, data.shape[0]-1),
                              max(0, j-d):min(j+u+1, data.shape[1]-1)]

            # Set all pixels to zero; this is masking the saturated pixels
            # so we don't select them again
            for a in range(box.shape[0]):
                for b in range(box.shape[1]):
                    dataCropped[a+max(0, i-l), b+max(0, j-d)] = 0

        else:
            found = True  # If our selection is not saturated, it may be used

    # Add back on margin values to convert to main coord-system
    i += m
    j += m

    minVal = 2000

    # Set a lower threshold on how bright star must be
    if data[i, j] < minVal:
        print("No star brighter than %s found" % minVal)
        return -1, -1

    # Return position of star
    return i, j


# Simple quadratic function
def quadratic(x, a, b, c):
    return a*x**2 + b*x + c


# Take in directory containing FITS images and
# return python list of pyfits objects
def getImgList(path):

    try:
        fileList = os.listdir(path)  # Parse filenames to list

        imageList = []  # Create empty list to hold FITS objects

        for infile in fileList:

            imageList.append(pf.open(path+"\\"+infile),
                             ignore_missing_end=True)  # Build FITS list

        return imageList

    except Exception:
        print("Error in getting image list")


# Extract r,i,g,B bands from Rainbow Camera fits image
def extractBands(img, runType='all'):

    if type(img) == pf.hdu.hdulist.HDUList:
        data = img[0].data
    elif type(img) == str:
        data = (pf.open(img, ignore_missing_end=True))[0].data
    else:
        data = img

    try:
        # Get 1/8 of dimensions of matrix
        x = data.shape[0]
        y = data.shape[1]

        # Parse central box from each quadrant using deep copies
        data1 = np.empty_like(data[0: x/2, 0: y/2])
        data1[:] = data[0: x/2, 0: y/2]

        data2 = np.empty_like(data1)
        data2[:] = np.array(data[x/2: x, 0: y/2])

        data3 = np.empty_like(data1)
        data3[:] = np.array(data[0: x/2, y/2: y])

        data4 = np.empty_like(data1)
        data4[:] = np.array(data[x/2: x, y/2: y])

        tuple1 = [data1, (0, 0)]
        tuple2 = [data2, (x/2, 0)]
        tuple3 = [data3, (0, y/2)]
        tuple4 = [data4, (x/2, y/2)]
        allTuples = [tuple1, tuple2, tuple3, tuple4]

        # Return appropriate band data
        if runType == 'B':
            return [tuple2]
        elif runType == 'g':
            return [tuple4]
        elif runType == 'r':
            return [tuple1]
        elif runType == 'i':
            return [tuple3]
        else:
            return allTuples  # Default

    except Exception:

        print("Error while parsing bands from FITS")


# Parse time from the 'utc' value in a fits header
def getTime(fits):

    utc = fits[0].header['utc']
    utc = utc.split(':')
    for i in range(len(utc)):
        utc[i] = int(float(utc[i]))
    utc = [utc[0], int(utc[1]/30), utc[1] % 30]+utc[2:]
    return timegm(utc)


# I DID NOT WRITE THE FOLLOWING METHODS: convolve,fitgaussian,gaussian,moments
# THEY ARE 'OFF-THE-SHELF'... OR 'OFF-THE-GOOGLE', I should say.

# Convolution function
def convolve(image1, image2, MinPad=True, pad=True):
    """ Not so simple convolution """
    try:
        # Just for comfort:
        FFt = np.fft.fft2
        iFFt = np.fft.ifft2

        # The size of the images:
        r1, c1 = image1.shape
        r2, c2 = image2.shape

        # MinPad results simpler padding,smaller images:
        if MinPad:
            r = r1+r2
            c = c1+c2
        else:
            # if the Numerical Recipies says so:
            r = 2*max(r1, r2)
            c = 2*max(c1, c2)

        # For nice FFT, we need the power of 2:
        if pad:
            pr2 = int(np.log(r)/np.log(2.0) + 1.0)
            pc2 = int(np.log(c)/np.log(2.0) + 1.0)
            rOrig = r
            cOrig = c
            r = 2**pr2
            c = 2**pc2
        # end of if pad

        # numpy fft has the padding built in, which can save us some steps
        # here. The thing is the s(hape) parameter:
        fftimage = FFt(image1, s=(r, c)) * FFt(image2, s=(r, c))

        if pad:
            # Crop to ignore correlations of less than 1/2 overlap
            x, y = ((iFFt(fftimage))[:rOrig, :cOrig]).real.shape
            return ((iFFt(fftimage))[:rOrig, :cOrig]).real[x/4:3*x/4, y/4:3*y/4]
        else:
            # Crop to ignore correlations of less than 1/2 overlap
            x, y = (iFFt(fftimage)).real.shape
            return (iFFt(fftimage)).real[x/4:3*x/4, y/4:3*y/4]

    except Exception:
        print("Error in convolution")
    return


# Gaussian fitting function
def fitgaussian(data):
    try:

        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        params = moments(data)
        errorfunction = lambda p:\
            np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
        p, success = leastsq(errorfunction, params)
        return p

    except Exception:

        print("Error in fitgaussian")
        return []


# Gaussian function
def gaussian(height, center_x, center_y, width_x, width_y):

    """Returns a gaussian function with the given parameters"""
    try:
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x, y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
    except Exception:

        print("Error in gaussian method")
        return


def moments(data):

    try:
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution by calculating its
        moments """
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
        row = data[int(x), :]
        width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
        height = data.max()
        return height, x, y, width_x, width_y

    except Exception:

        print("Error in moments()")
        return
