# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:57:05 2015

@author: nadiablago
"""
from __future__ import print_function

import matplotlib
matplotlib.use("Agg")
import zscale
from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import math, glob
from astropy.wcs import WCS
from astropy.io.votable import parse_single_table
import coordinates_conversor
import app_phot
from numpy.lib import recfunctions as rfn
import scipy.optimize as opt
import fitsutils
import os, shutil
from astropy import stats
import argparse
import logging
import datetime
import QueryCatalogue
import rcred
import time_utils
import matplotlib.dates as md

from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')


FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
now = datetime.datetime.utcnow()
timestamp=datetime.datetime.isoformat(now)
creationdate = timestamp
timestamp=timestamp.split("T")[0]

try:
    #Log into a file
    root_dir = _logpath
    logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "rcred_{0}.log".format(timestamp)), level=logging.INFO)
    logger = logging.getLogger('zeropoint')
except:
    logging.basicConfig(format=FORMAT, filename=os.path.join("/tmp", "rcred_{0}.log".format(timestamp)), level=logging.INFO)
    logger= logging.getLogger("zeropoint")
    




def are_isolated(rav, decv, r):
    '''
    Returns a mask saying whether the stars with coordinates rav, decv are isolated in a radius of r arcsec.
    '''
    
    index = np.arange(len(rav))
    mask = []
    
    for i in np.arange(len(rav)):
        d = coordinates_conversor.get_distance(rav[i], decv[i], rav[index!=i], decv[index!=i])
        mask.append(~ np.any(d*3600<r))
    
    return np.array(mask)

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''
    Produces a 2D gaussian centered in xo, yo with the parameters specified.
    xdata_tuple: coordinates of the points where the 2D Gaussian is computed.
    
    '''
    (x, y) = xdata_tuple                                                        
    xo = float(xo)                                                              
    yo = float(yo)                                                              
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))                                   
    return g.ravel()
    
def twoD_Gauss_test(theta=0):
    '''
    Generates a test Gaussian and fits it using the scisoft optimization software.
    '''
    # Create x and y indices
    x = np.linspace(0, 200, 201)
    y = np.linspace(0, 200, 201)
    x, y = np.meshgrid(x, y)
    
    #create data
    data = twoD_Gaussian((x, y), 3, 100, 100, 20, 40, theta, 10)
    
    # plot twoD_Gaussian data generated above
    plt.figure()
    plt.imshow(data.reshape(201, 201), origin="bottom", extent=(x.min(), x.max(), y.min(), y.max()))
    plt.colorbar()
    
    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = (3,100,100,20,40,0,10)
    
    data_noisy = data + 0.2*np.random.normal(size=data.shape)
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_noisy, p0=initial_guess)
    
    data_fitted = twoD_Gaussian((x, y), *popt)
    
    fig, ax = plt.subplots(1, 1)
    ax.hold(True)
    ax.imshow(data_noisy.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
        extent=(x.min(), x.max(), y.min(), y.max()))
    ax.contour(x, y, data_fitted.reshape(201, 201), 8, colors='w')
    plt.show()

def find_fwhm(imfile, xpos, ypos, plot=True):
    '''
    Finds and returns the best parameters for the FWHM in arcsec for the stars marked with X, Y
    '''
    
    f = pf.open(imfile)
    img = f[0].data
    pix2ang = 0.394
    def_fwhm = 2./pix2ang
    rad = math.ceil(30./pix2ang)
    
    out = np.zeros(len(xpos), dtype=[('detected', np.bool), ('fwhm', np.float), ('e', np.float)])
    
    for i, (x_i,y_i) in enumerate(zip(xpos, ypos)):

        x_i = int(x_i)
	y_i = int(y_i)
        hrad = int(math.ceil(rad/2.))
        sub = img[x_i-hrad:x_i+hrad, y_i-hrad:y_i+hrad]
        x = np.linspace(0, len(sub), len(sub))
        y = np.linspace(0, len(sub), len(sub))
        X, Y = np.meshgrid(x, y)
    
        #(xdata_tuple, amplitude, xo, yo, def_fwhm, def_fwhm, theta, offset):
        def_x = np.argmax(np.sum(sub, axis=0))
        def_y = np.argmax(np.sum(sub, axis=1))

        initial_guess = (100, def_x, def_y, def_fwhm, def_fwhm, 0, np.percentile(sub, 40))
        detected = True
        try:
            popt, pcov = opt.curve_fit(twoD_Gaussian, (X, Y), sub.flatten(), p0=initial_guess)
            fwhm_x = np.abs(popt[3])*2*np.sqrt(2*np.log(2))
            fwhm_y = np.abs(popt[4])*2*np.sqrt(2*np.log(2))
            amplitude=popt[0]
            background=popt[-1]
        except:
            detected = False
            fwhm_x = 0
            fwhm_y = 0
            amplitude = 0
            background = 0
        
        
        detected = detected * (amplitude> 0)
        
        if (not detected):
            fwhm_x = 0
            fwhm_y = 0     
        
        
        logger.info("%s %s Amplitude %.3f\t BG %.3f\t BG_stats %.3f\t  FWHM_x,FWHM_y=(%.3f, %.3f)"%(i, detected, amplitude, background, np.percentile(sub, 50), fwhm_x, fwhm_y))
        
        out[i] = (detected, np.average([fwhm_x, fwhm_y]), np.minimum(fwhm_x, fwhm_y) / np.maximum(fwhm_x, fwhm_y))
        
        if (detected & plot):
            data_fitted = twoD_Gaussian((X, Y), *popt)
            
            fig, (ax, ax2) = plt.subplots(1, 2)
            ax.hold(True)
            ax.imshow(sub, cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(X, Y, data_fitted.reshape(sub.shape[0], sub.shape[1]), 5, colors='w')
            ax2.imshow(sub-data_fitted.reshape(sub.shape[0], sub.shape[1]), cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            ax2.contour(X, Y, data_fitted.reshape(sub.shape[0], sub.shape[1]), 5, colors='w')
            ax.scatter(def_x, def_y, marker="*", s=100, color="yellow")
            plt.title("DETECTED for Position X,Y = %d,%d"%(x_i,y_i))
            plt.savefig(os.path.join(os.path.dirname(imfile), "gauss_%d"%i))
            plt.clf()
        if ((not detected) & plot):           
            fig, ax = plt.subplots(1)
            ax.hold(True)
            ax.imshow(sub, cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            plt.title("NOT DETECTED for Position X,Y = %d,%d"%(x_i,y_i))
            plt.savefig(os.path.join(os.path.dirname(imfile), "gauss_%d"%i))
            plt.clf()
            #plt.show()
            
    return out
    
def clean_tmp_files():
    
    files = ["/tmp/tmp_sdss_%s.cat"%creationdate, '/tmp/tmp_%s.cat'%creationdate, '/tmp/tmp_apass_%s.cat'%creationdate, '/tmp/tmp_sdss_%s.cat'%creationdate, '/tmp/sdss_cat_det_%s.txt'%creationdate]
    
    for f in files:
        if os.path.isfile(f):
            os.remove(f)
    
def extract_star_sequence(imfile, band, plot=True, survey='ps1', debug=False, refstars=None, plotdir=".", pix2ang = 0.394):
    '''
    Given a fits image: imfile and a the name of the band which we want to extract the sources from,
    it saves the extracted sources into  '/tmp/sdss_cat_det.txt' file.
    If the band does not match the bands in the survey, a change is performed to adapt to the new band.
    
    If plotting activated, plots the USNOB1 field of stars on top of the star field image.
    Red circles are stars identified from the catalogue in the correct magnitude range.
    Yellow circles are stars that are isolated.
    
    '''
    
    survey = str.lower(survey)
    minmag = 13.5
    maxmag = 20.0
        
    f = pf.open(imfile)
    wcs = WCS(f[0].header)

    img = f[0].data
    img[img<0] = 0
    
    ra, dec = wcs.wcs_pix2world(np.array([img.shape[0]/2, img.shape[1]/2], ndmin=2), 1)[0]
    ra0, dec0 = wcs.wcs_pix2world(np.array([img.shape[0], img.shape[1]], ndmin=2), 1)[0]


    sr = 2*np.abs(dec-dec0)
    logger.info("%.4f %.4f %.4f"%( ra,dec, sr))
    
    qc = QueryCatalogue.QueryCatalogue(ra, dec, sr/1.8, minmag, maxmag, logger)
    
    if not refstars is None:
        shutil.copy(refstars, "/tmp/tmp_sdss_%s.cat"%creationdate)
        catalog = np.genfromtxt("/tmp/tmp_sdss_%s.cat"%creationdate, names=True, dtype=None, delimiter=",")
        cat_ra = catalog["ra"]
        cat_dec = catalog["dec"]
        try:
            mag = catalog["R"]
        except:
            mag = catalog["r"]
        
    elif (survey=='sdss'):

        catalog = qc.query_sdss()

        if (np.ndim(catalog)==0 or len(catalog)<2):
            return False
            
        try:
            cat_ra = np.array(catalog['ra'], ndmin=1)
            cat_dec = np.array(catalog['dec'], ndmin=1)
            if (band in catalog.dtype.names):
                print ( "SDSS filter detected")
                mag = np.array(catalog[band], ndmin=1)
            else:
                print ("Unknown band!! %s"%band)
        except IOError:
            logger.error( "Problems with SDSS image %s"% band)
            return False
        except ValueError:
            logger.error( "Problems with the catalogue for the image")
            return False
            
    elif (survey=='ps1'):

        catalog = qc.query_catalogue(filtered=True)

        if (np.ndim(catalog)==0 or catalog is None):
            return False
            
        try:
            cat_ra = np.array(catalog['ra'], ndmin=1)
            cat_dec = np.array(catalog['dec'], ndmin=1)
            if (band in catalog.dtype.names):
                mag = np.array(catalog[band], ndmin=1)
                print ( "SDSS filter detected %s"%band)
            else:
                print ("Unknown band!! %s"%band)
        except IOError:
            logger.error( "Problems with SDSS image %s"% band)
            return False
        except ValueError:
            logger.error( "Problems with the catalogue for the image")
            return False


    print ("Catalogue has %d entries"%len(cat_ra))
    #Convert ra, dec position of all stars to pixels.
    star_pix = np.array([0,0])
    for i in range(len(cat_ra)):
        # Get pixel coordinates of reference stars
        ra, dec = rcred.get_xy_coords(imfile, cat_ra[i], cat_dec[i])
        #s = wcs.wcs_sky2pix(np.array([cat_ra[i], cat_dec[i]], ndmin=2), 1)[0]
        star_pix = np.row_stack((star_pix, np.array([ra, dec])))
    star_pix = star_pix[1:]

    rad = math.ceil(25./pix2ang)
    
    #Select only the stars within the image.
    mask = (star_pix[:,0]>-rad) * (star_pix[:,0]<img.shape[1]+rad)*(star_pix[:,1]>-rad) * (star_pix[:,1]<img.shape[0]+rad)
    if (band == 'u'):    
        mask = mask * (mag < 20)
        
    #Select only stars isolated in a radius of ~12 arcsec.
    mask2 = np.array(are_isolated(cat_ra[mask], cat_dec[mask], 10.))
    if (len(mask2)==0):
	logger.error("No good stars left")
	return False   
 
    #Select only stars that are within the proper magnitude range
    mask3 = (mag[mask][mask2] < maxmag) * (mag[mask][mask2] > minmag) 
    
    
    #Combine all masks
    mask3 = mask3 * (star_pix[:,0][mask][mask2]>rad) * (star_pix[:,0][mask][mask2]<img.shape[1]-rad)*(star_pix[:,1][mask][mask2]>rad) * (star_pix[:,1][mask][mask2]<img.shape[0]-rad)

        
    if (not np.any(mask) and not np.any(mask2) and not np.any(mask3)) or len(catalog[mask][mask2][mask3])==0:
        logger.debug (star_pix)
        print ( "No stars left...", mask, mask2, mask3)
        return False
    else:
        catalog = catalog[mask][mask2][mask3]
        s = star_pix[mask][mask2][mask3]

        print ("left %d stars"%(len(catalog)), catalog.dtype.names)

        z = np.zeros(len(s), dtype=[('x','f8'), ('y', 'f8')])
        z['x'] = s[:,0]
        z['y'] = s[:,1]
                
        header='x y '
        for n  in catalog.dtype.names:
            if n in ["objid", "ra", "dec", "u", "g", "r", "i", "z", "Err_u", "Err_g", "Err_r", "Err_i", "Err_z"] or\
                n in ['id', 'ra', 'dec', 'U', 'B', 'V', 'R', 'I', 'dU', 'dB', 'dV', 'dR', 'dI']:
                z = rfn.append_fields(z, names=n, data=catalog[n], usemask=False)
                header += n.replace("Err_", "d") + " "
           
        fmt = "%.5f"
        for i in range(len(z[0])-1):
            fmt +=  " %.5f"

        np.savetxt('/tmp/sdss_cat_%s.txt'%creationdate, z, fmt=fmt, header = header)
        logger.info( "Saved catalogue stars to %s"% ('/tmp/sdss_cat_%s.txt'%creationdate))
        
        #Find FWHM for this image            
        out = find_fwhm(imfile, star_pix[:,1][mask][mask2][mask3], star_pix[:,0][mask][mask2][mask3], plot=debug)
        mask_valid_fwhm = (out['detected']) * (out['e']>0.6) * ~np.isnan(out['fwhm']* (out['fwhm'] < 30))            

        if ((np.count_nonzero(mask_valid_fwhm) < 3) and (fitsutils.get_par(imfile, "FILTER")!="u")) or ( (np.count_nonzero(mask_valid_fwhm) < 2) and (fitsutils.get_par(imfile, "FILTER")=="u")):
            logger.error( "ERROR with FWHM!! Too few points for a valid estimation. %d"% np.count_nonzero(mask_valid_fwhm)+ ") points")
            logger.error( "%s %s"%(out["detected"], out["fwhm"]))
            return False

        outd = out[mask_valid_fwhm]

        logger.info( 'Average FWHM %.3f arcsec, %.3f pixels'%(np.median(outd['fwhm']),  np.median(outd['fwhm'])*pix2ang))
        


        if (band in 'ugriz' and survey=='sdss'):
            header='x y objid ra dec u g r i z du dg dr di dz'
        if (band in 'grizy' and survey=='ps1'):
            header='x y objid ra dec g r i z dg dr di dz'
        elif band in 'UBVRI':
            header='x y objid ra dec U B V R I dU dB dV dR dI'
        np.savetxt('/tmp/sdss_cat_det_%s.txt'%creationdate, z[mask_valid_fwhm], fmt=fmt, \
        header=header)
        #print ("Saved to /tmp/sdss_cat_det_%s.txt"%creationdate)
        
            
        
    #Plot results
    img = img - np.nanmin(img)
    zmin, zmax = zscale.zscale(img)
        
    logger.info( "Found %d stars in %s. "%(len(cat_dec), survey)+\
	     "%d of them within the FoV. "%len(cat_ra[mask]) +\
            "%d of them are isolated."%len(cat_ra[mask][mask2])+\
            "%d of them with suitable magnitudes. "%len(cat_ra[mask][mask2][mask3]) +\
            "%d of them with detected stars."%np.count_nonzero(mask_valid_fwhm)) 
    
    
    if (plot):
        im = plt.imshow(img, aspect="equal", origin="lower", cmap=matplotlib.cm.gray_r, interpolation="none", vmin=zmin, vmax=zmax)

        
        if (len(star_pix[:,0][mask]) >0):        
            plt.scatter(star_pix[:,0][mask], star_pix[:,1][mask], marker="o", s=np.minimum(150, 10000*(10./mag[mask][mask2])**9), edgecolor="red", facecolor="none", label="catalogue")
        
        if (len(star_pix[:,0][mask][mask2]) >0):        
            plt.scatter(star_pix[:,0][mask][mask2], star_pix[:,1][mask][mask2], marker="o", s=20, edgecolor="yellow", facecolor="none", label="isolated")

        if (len(star_pix[:,0][mask][mask2][mask3]) >0):        
            plt.scatter(star_pix[:,0][mask][mask2][mask3], star_pix[:,1][mask][mask2][mask3], marker="o", s=200, edgecolor="green", facecolor="none", label="wihtin farme and mag")


        selected = star_pix[:,:][mask][mask2][mask3][mask_valid_fwhm]
        if (len(selected) >0):        
            plt.scatter(selected[:,0], selected[:,1], marker="o", \
                s=400, edgecolor="blue", facecolor="none")
            for i in np.arange(len(selected)):
                plt.text(selected[i,0]+10, selected[i,1]+10, i+1)
        
        plt.legend(loc="best", frameon=False, framealpha=0.9)
        figname = os.path.join( plotdir, os.path.basename(imfile).replace('.fits', '.seqstars.png').replace('.new', '.seqstars.png'))
        plt.savefig( figname)
        logger.info( "Saved stars to %s"%figname)        
        plt.clf()
        
    return True
        
def add_to_zp_cal(ref_stars, image, logname):
    '''
    Records the parameters and instrumental and standard star magnitudes needed for zeropoint calibration.
    ref_stars: name of the file that contains the reference stars that are present in the image.
    image: fits image with the sources in them.
    logname: name of the file where the information on reference and measurements are logged.
    
    minst = M + c0 + c1(AIRMASS) + c2(color) + c3(UT)    
    '''
   
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'r', 'z':'i', 'U':'B', 'B':'V', 'V':'R', 'R':'I', 'I':'R'}

    r = np.genfromtxt(ref_stars, delimiter=" ", dtype=None, names=True)
    imapp = os.path.join(os.path.join(os.path.dirname(image), "zeropoint"), os.path.basename(image) + ".app.mag")
    #imapp = os.path.join(os.path.dirname(image), os.path.basename(image) + ".app.mag")

    my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"),  ("image","|S20"), ("X","<f4"), ("Y","<f4"), ("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("msky","<f4"), ("stdev","<f4"),\
        ("flags", np.int), ("rapert", "<f4"), ("sum", "<f4"), ("area", "<f4"), ("nsky","<f4") , ("flux", "<f4"), ("itime", "<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])

    if (my.size < 2):
        my = np.array([my])
    if (r.size < 2):
        r = np.array([r])
    
    band = fitsutils.get_par(image, 'filter')
    mask_valid1 = np.array(my["fwhm"]<9000) * np.array(~ np.isnan(r[band])) * np.array(~ np.isnan(my["fit_mag"])) * np.array(my["fit_mag"]<9000)
    
    r = r[mask_valid1]
    my = my[mask_valid1]
    N = len(r)
    
    logger.info("Adding %d stars to the ZP calibration file %s"%(N, logname)) 

    my["fiterr"][np.isnan( my["fiterr"])] = 100

    col_band = coldic[band]
    exptime = fitsutils.get_par(image, "exptime")
    fwhm = fitsutils.get_par(image, "FWHM")
    if fwhm is None:
	fwhm = 0
    airmass = 1.3
    name = "object"
    if (fitsutils.has_par(image, "NAME")):
        name = fitsutils.get_par(image, "NAME")
    if (fitsutils.has_par(image, "AIRMASS")):
        airmass = fitsutils.get_par(image, "AIRMASS")
    if (fitsutils.has_par(image, "JD")):
        date = fitsutils.get_par(image, "JD")
    elif (fitsutils.has_par(image, "MJD")):
        date = fitsutils.get_par(image, "MJD")
    elif (fitsutils.has_par(image, "MJD-OBS")):
        date = fitsutils.get_par(image, "MJD-OBS")
    else:
        date = 0
    
    if (not os.path.isfile(logname)):
            with open( logname, "a") as f:
                f.write("#object,filename,filter,std,stderr,inst,insterr,jd,airmass,color,exptime,x,y,dx,dy,fwhm\n")
    with open( logname, "a") as f:
        for i in range(N):
            f.write("%s,%s,%s,%.3f,%.3f,%3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f\n"%
            (name, image, band, r[i][band], r[i]["d"+band], my[i]['fit_mag'], my[i]['fiterr'], date, airmass, r[i][band]-r[i][col_band],exptime,my[i]['X'],my[i]['Y'],my[i]['Xshift'],my[i]['Yshift'], fwhm))

def lsq_test():

    x = np.random.normal(0, 1, 100)
    y = np.random.normal(0, 1, 100)
    
    X, Y = np.meshgrid(x, y)
    
    D = 23 + 2*X + 4*Y + np.random.rand(100, 100)*1
    O = np.zeros_like(D)
    
    M = np.zeros((len(X.flatten()), 3))
    M[:,0] = 1 
    M[:,1] = X.flatten()
    M[:,2] = Y.flatten()
    
    depend = D.flatten() - O.flatten()
        
    lsq_result = np.linalg.lstsq(M, depend)
    coef = lsq_result[0]
    res = lsq_result[1]
    
    f, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(X.flatten(), depend -coef[0] -Y.flatten()*coef[2] , ".")
    ax1.plot(X.flatten(),  X.flatten()*coef[1])
    ax1.set_xlabel("X")

    ax2.plot(Y.flatten(), depend -coef[0] -X.flatten()*coef[1] , ".")
    ax2.plot(Y.flatten(),  Y.flatten()*coef[2])
    ax2.set_xlabel("Y")
    
    plt.show()
       
       
def calibrate_zp_fourshot(logfile, plot=True):
    '''
    The field of view is quite small, therefore, all rc shots are used to calibrate the zeropoint.
    This routine retrieves all the stars taken with the same filter when pointing to the science object.
    '''

    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    a.sort(order=['jd'], axis=0)
    a = a[a['inst']!=0]

    plotdir = os.path.join(os.path.dirname(os.path.abspath(logfile)), "zeropoint")

    for name in set(a['object']):    
        for b in set(a['filter']):
            aib = a[(a["filter"]==b)*(a["object"]==name)*(np.abs(a["color"])<1)]
            
            if len(aib) < 3:
                print ("Less than 3 stars found with %s %s"%(name,b))
                continue
            
            #First fit for the linear and detect the deviations
            coefs, residuals, rank, singular_values, rcond = np.polyfit(aib["std"], aib["inst"], w=1./np.maximum(0.3, np.sqrt(aib["stderr"]**2 + aib["insterr"]**2)), deg=1, full=True)
            p = np.poly1d(coefs)
            
            if (plot):
                plt.figure()
                plt.title("%s Filter %s"%(name, b))
                plt.errorbar(aib["std"], aib["inst"], yerr=np.sqrt(aib["stderr"]**2 + aib["insterr"]**2), fmt="o")
                plt.plot(aib["std"], p(aib["std"]))
                
            diff = np.abs(aib["inst"] - p(aib["std"]))
            mad = stats.funcs.median_absolute_deviation(diff)
            aib = aib[diff<mad*5]
            
            if (plot):
                plt.errorbar(aib["std"], aib["inst"], yerr=np.sqrt(aib["stderr"]**2 + aib["insterr"]**2), fmt="o")
                plt.plot(aib["std"], p(aib["std"]))
                plt.savefig(os.path.join(plotdir, "zp_mag_mag_%s_%s.png"%(name, b)))
                plt.close()
            

            #Then fit for the colour
            coefs, residuals, rank, singular_values, rcond = np.polyfit(aib["color"], aib["std"] - aib["inst"], w=1./np.maximum(0.15, np.sqrt(aib["stderr"]**2 + aib["insterr"]**2)), deg=1, full=True)
            p = np.poly1d(coefs)
            
            color, zp = coefs
            
            mad = stats.funcs.median_absolute_deviation(p(aib["color"]) - (aib["std"] - aib["inst"]))
             
            
            for f in aib["filename"]:
                #Add these values to the header.
                pardic = {"IQZEROPT" : 1,\
                        "ZPCAT" : "SDSS4shot",\
                        "ZEROPTU" : float("%.3f"%mad),\
                        "ZEROPT" : float("%.3f"%zp),\
                        "ZP": float("%.3f"%zp),\
                        "ZPERR": float("%.3f"%mad)}
                fitsutils.update_pars(f, pardic)

            if (plot):
                plt.figure()
                plt.title("%s Filter %s"%(name, b))
                plt.errorbar(aib["color"], aib["std"] - aib["inst"], yerr=np.sqrt(aib["stderr"]**2 + aib["insterr"]**2), fmt="o")
                plt.plot(aib["color"], p(aib["color"]))
                plt.savefig(os.path.join(plotdir, "zp_col_mag_%s_%s.png"%(name, b)))
                plt.close()
    return
    
def lsq_zeropoint(logfile, plotdir=".", plot=True):
    '''
    Uses least squares approach to compute the optimum coefficients for ZP, colour term, airmass and time.
    Check this one:
    https://docs.scipy.org/doc/scipy/reference/odr.html
    
    '''
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    a.sort(order=['jd'], axis=0)
    '''a = a[a['inst']!=0]
    a = a[a['std']<20]
    a = a[a['insterr']<0.1]'''
    a = a[a['fwhm']<3.5]
    a = a[(a['std']-a['inst']<24) * (a['std']-a['inst']>20.0)]
        
    #Select only sources that have not been deviating too far fromt he predicted position.
    a = a[ (np.abs(a['dx'])<10)*(np.abs(a['dy'])<10)]
    
    a['jd'] = a['jd'] - np.min(a['jd'])
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}



    for b in ['g','r','i','u']:#set(a['filter']):
        ab = a[a['filter']==b]
        
        #Filter extreme colours which may bias the coefficient.
        '''mincol = np.median(ab["color"]) - 2*np.std(ab["color"])
        maxcol = np.median(ab["color"]) + 2*np.std(ab["color"])
        
        ab = ab[ (ab["color"]>mincol) * (ab["color"]<maxcol) ] '''       
        
        #Remove detections which are too far
        try:
            coefs, residuals, rank, singular_values, rcond = np.polyfit(ab["std"], ab["inst"], w=1./np.maximum(0.3, np.sqrt(ab["stderr"]**2 + ab["insterr"]**2)), deg=1, full=True)
            p = np.poly1d(coefs)
            
            if (plot):
                plt.figure()
                plt.title("Filter %s"%(b))
                plt.errorbar(ab["std"], ab["inst"], yerr=np.sqrt(ab["stderr"]**2 + ab["insterr"]**2), fmt="bo", alpha=0.5, ms=2)
                plt.plot(ab["std"], p(ab["std"]), "b-")
         
            #Removing outliers that deviate 5 times from the general trend.       
            mad = stats.funcs.median_absolute_deviation(ab["inst"] - p(ab["std"]))
            ab = ab[np.abs(ab["inst"] - p(ab["std"]))<mad*5]
            
                    
            if (plot):
                plt.errorbar(ab["std"], ab["inst"], yerr=np.sqrt(ab["stderr"]**2 + ab["insterr"]**2), fmt="ro", alpha=0.5, ms=2)
                plt.plot(ab["std"], p(ab["std"]), "r-")
                plt.savefig(os.path.join(plotdir, "filter_%s.png"%b))        
        except TypeError as e:
            print ("Could not eliminate outliers. Probably too few observations")
            
        #Find the coefficients.
        M = np.zeros((len(ab), 10))
        M[:,0] = 1      
        M[:,1] = ab['color'] 
        M[:,2] = ab['airmass'] - 1.3
        M[:,3] = ab['jd'] 
        M[:,4] = ab['jd']**2
        M[:,5] = ab['x']
        M[:,6] = ab['y']
        M[:,7] = ab['x']**2
        M[:,8] = ab['y']**2
        M[:,9] = ab['fwhm'] - 1.5
        #M[:,7] = ab['dx']
        #M[:,8] = ab['dy']
        
        #print M
        
        depend = ab['std']-ab['inst']
                
        lsq_result = np.linalg.lstsq(M, depend)
        coef = lsq_result[0]
        res = lsq_result[1]
                

        #Empirical and predicted values
        est_zp = coef[0] +ab['color']*coef[1] +(ab['airmass']-1.3)*coef[2] + coef[3]*ab['jd'] + coef[4]*ab['jd']**2 + coef[5]*ab['x'] + coef[6]*ab['y'] + coef[7]*ab['x']**2 + coef[8]*ab['y']**2 + coef[9]*(ab['fwhm']-1.5)#+ coef[5]*ab['jd']**3 


        emp_col = depend - (est_zp - ab['color']*coef[1])
        pred_col = ab['color']*coef[1]

        emp_airmass = depend - (est_zp - (ab['airmass']-1.3)*coef[2])
        pred_airmass = (ab['airmass']-1.3)*coef[2]
        
        emp_jd = depend - (est_zp - (coef[3]*ab['jd'] + coef[4]*ab['jd']**2) )
        pred_jd =  coef[3]*ab['jd'] + coef[4]*ab['jd']**2# + coef[5]*ab['jd']**3 
                
        rms =  np.sqrt(np.sum((depend-est_zp)**2)/(len(depend)-1))
        mad = stats.funcs.median_absolute_deviation(depend-est_zp)
        
        #print ("Filter %s ZP %.2f Color index %.2f RMS %.2f MAD %.2f"%(b, coef[0], coef[1], rms, mad))
        
        #Remove outliers
        mask = np.abs(depend-est_zp) < 5*mad
        
        ab = ab[mask]
        
        M = np.zeros((len(ab), 13))
        M[:,0] = 1      
        #M[:,1] = ab['inst']
        M[:,1] = ab['color'] 
        M[:,2] = ab['airmass'] - 1.3
        M[:,3] = ab['jd'] 
        M[:,4] = ab['jd']**2
        M[:,5] = ab['jd']**3
        M[:,6] = ab['x']
        M[:,7] = ab['y']
        M[:,8] = ab['x']**2
        M[:,9] = ab['y']**2
        M[:,10] = ab['x']**3
        M[:,11] = ab['y']**3
        #M[:,12] = ab['x']**4
        #M[:,13] = ab['y']**4
        M[:,12] = ab['fwhm']-1.5
        #M[:,6] = ab['fwhm']-1.5
        
        
        #print M
        
        depend = ab['std']-ab['inst']
        
        try:
            lsq_result = np.linalg.lstsq(M, depend)
            coef = lsq_result[0]
            res = lsq_result[1]
        except:
            print ("Failed for filter %s"%b)
            continue
        
        coefs_map = {'color': coef[1], 'airmass': coef[2], 'jd': coef[3], 'jd2': coef[4], 'jd3': coef[5], 'fwhm': coef[12],\
        'x': coef[6], 'y':coef[7], 'x2':coef[8], 'y2':coef[9], 'x3':coef[10], 'y3':coef[11]}#, 'x4':coef[12], 'y4':coef[13], 'fwhm': coef[14]}
        #Save the coefficients 
        np.savetxt("coefs_%s.txt"%b, coef)
        
        #Empirical and predicted values
        #estimated zeropoint for each observation
        #+ coefs_map['fwhm']*(ab['fwhm']-1.5) \
        est_zp = coef[0] + coefs_map['color']*ab['color'] + coefs_map['airmass']*(ab['airmass']-1.3) + coefs_map['jd']*ab['jd'] + coefs_map['jd2']*ab['jd']**2 + coefs_map['jd3']*ab['jd']**3 + \
            + coefs_map['x']*ab['x'] + coefs_map['y']*ab['y'] + coefs_map['x2']*ab['x']**2 + coefs_map['y2']*ab['y']**2 + coefs_map['x3']*ab['x']**3 + coefs_map['y3']*ab['y']**3 + \
             + coefs_map['fwhm']*(ab['fwhm']-1.5)#+ coef[5]*ab['jd']**3 
             #+ coefs_map['x4']*ab['x']**4 + coefs_map['y4']*ab['y']**4 


        emp_col = depend - (est_zp - coefs_map['color']*ab['color'])
        pred_col = ab['color']*coef[1]

        emp_airmass = depend - (est_zp - coefs_map['airmass']*(ab['airmass']-1.3))
        pred_airmass = (ab['airmass']-1.3)*coef[2]
        
        emp_jd = depend - (est_zp - (coefs_map['jd']*ab['jd'] + coefs_map['jd2']*ab['jd']**2 + coefs_map['jd3']*ab['jd']**3) )
        pred_jd =  coefs_map['jd']*ab['jd'] + coefs_map['jd2']*ab['jd']**2 + coefs_map['jd3']*ab['jd']**3 
                
        rms =  np.sqrt(np.sum((depend-est_zp)**2)/(len(depend)-1))
        mad = stats.funcs.median_absolute_deviation(depend-est_zp)

        np.savetxt("rms_%s.txt"%b, np.array([rms]))
        print ("Filter %s ZP %.2f Color index %.2f RMS %.2f MAD %.2f"%(b, coef[0], coef[1], rms, mad))

        if (plot):
            plt.close("all")
            f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
            f.set_figheight(12)
            f.set_figwidth(10)
            ax1.plot(ab['color'], emp_col, "o", color=cols[b], ms=4, alpha=0.4)
            #ax1.plot(ab['color'], depend - est_zp, "o", color=cols[b], ms=4, alpha=0.4)
            ax1.plot(ab['color'],  pred_col, color=cols[b])
            ax1.set_title("colour")

            ax2.plot(ab['airmass'], emp_airmass, "o", color=cols[b], ms=4, alpha=0.4)
            #ax2.plot(ab['airmass'], depend - est_zp, "o", color=cols[b], ms=4, alpha=0.4)

            ax2.plot(ab['airmass'], pred_airmass , color=cols[b])
            ax2.set_title("airmass")
            
            ax3.errorbar(ab['fwhm'], depend - (est_zp -coefs_map['fwhm']*(ab['fwhm']-1.5)), yerr=np.minimum(1, ab['insterr']), fmt="o", c=cols[b], alpha=0.4, ms=4)
            #ax3.errorbar(ab['fwhm'], depend - est_zp, yerr=np.minimum(1, ab['insterr']), fmt="o", c=cols[b], alpha=0.4, ms=4)
            ax3.plot(ab['fwhm'], coefs_map['fwhm']*(ab['fwhm']-1.5), color=cols[b])
            ax3.set_xlabel("FWHM")  
            ax3.set_title("FWHM")

            
            ax4.errorbar(ab['x'], depend - (est_zp -(coefs_map['x']*ab['x'] + coefs_map['x2']*ab['x']**2 + coefs_map['x3']*ab['x']**3 )), yerr=np.minimum(1, ab['insterr']), fmt="o", c=cols[b], alpha=0.4, ms=4)
            #ax4.plot(ab['x'], depend - est_zp, color=cols[b])
            ax4.errorbar(ab['y'], depend - (est_zp -(coefs_map['y']*ab['y'] + coefs_map['y2']*ab['y']**2 + coefs_map['y3']*ab['y']**3)), yerr=np.minimum(1, ab['insterr']), fmt="o", c="blue", alpha=0.4, ms=4)
            #ax4.plot(ab['x'], depend - est_zp, color=cols[b])
            ax4.set_xlabel("Y")  
            ax4.set_title("X and Y")

            
            ax5.plot(ab['jd']*24, emp_jd, "o", color=cols[b], ms=4, alpha=0.4)
            ax5.plot(ab['jd']*24, pred_jd, color=cols[b])
            ax5.set_xlabel("Elapsed hours")  
            ax5.set_title("Obs. Time")

                
            #ax4.plot(ab['jd'], depend, "*", color=cols[b], alpha=0.4)
            #ax4.errorbar(ab['jd'], est_zp, yerr=ab['insterr'], marker="o", c=cols[b], ls="none")
            ax6.errorbar(ab['jd']*24, depend - est_zp, yerr=np.minimum(0, ab['insterr']), fmt="o", c=cols[b], alpha=0.4, ms=3)
            ax6.set_ylim(-0.2, 0.2)
            ax6.set_xlabel("Elapsed hours")  
            ax6.set_title("Residuals (Measured_ZP - Estimated_ZP)")
            ax6.invert_yaxis()
            
            
            plt.tight_layout()
            
            if (not plotdir is None):
                plt.savefig(os.path.join(plotdir, "allstars_%s.png"%b))
            else:	
                plt.show()


            #Second set of figures.
            medians = []
            perc_10 = []
            perc_90 = []         
            magdiff = ab['std']-(ab['inst'] + est_zp)
            magbins = np.arange(14,20.5,0.5)
            for mi in magbins :
                mask = np.abs(ab['std']-mi)<0.25
                if np.any(mask):
                    medians.append(np.median(magdiff[mask]))
                    perc_10.append(np.percentile(magdiff[mask], 16.666667))
                    perc_90.append(np.percentile(magdiff[mask], 83.333333))
                else:
                    medians.append(np.nan)
                    perc_10.append(np.nan)
                    perc_90.append(np.nan)
                
            f  = plt.figure()


            plt.errorbar(ab['std'], ab['std']-(ab['inst'] + est_zp), yerr=np.sqrt(ab['insterr']**2 + ab['insterr']**2), fmt="o", c=cols[b], alpha=0.3, ms=3)
            #plt.hist2d(ab['std'], ab['std']-(ab['inst'] + est_zp), bins=(35,35), range=((13,20.5), (-0.2,0.2)), cmap=matplotlib.cm.gray_r)
            plt.plot(magbins, medians, "k-")
            plt.plot(magbins, perc_10, "k:")
            plt.plot(magbins, perc_90, "k:")
            xmin, xmax = plt.xlim()
            plt.hlines(0, xmin, xmax, linestyles="--", color="b")
            plt.hlines(0.05, xmin, xmax, linestyles="--", color="b", linewidth=0.3)
            plt.hlines(-0.05, xmin, xmax, linestyles="--", color="b", linewidth=0.3)
            plt.hlines(0.025, xmin, xmax, linestyles="--", color="red", linewidth=0.3)
            plt.hlines(-0.025, xmin, xmax, linestyles="--", color="red", linewidth=0.3)
            plt.ylim(-0.2,0.2)
            plt.xlabel("magnitude")  
            plt.ylabel("Measured - Estimated [mag]")           
                
            if (not plotdir is None):
                plt.savefig(os.path.join(plotdir, "zpfit_%s.png"%b))
                plt.close("all")
            else:	
                plt.show()
                
            #Third set of figures.
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            diff = depend - est_zp
            im = ax1.scatter(ab['x'], ab['y'], c=depend - (est_zp), s=4)# -(coefs_map['x']*ab['x'] + coefs_map['x2']*ab['x']**2 + coefs_map['x3']*ab['x']**3 ) \
                #-(coefs_map['y']*ab['y'] + coefs_map['y2']*ab['y']**2 + coefs_map['y3']*ab['y']**3) ))
            ax1.set_title("Filter %s"%b)
            ax1.set_xlabel("X")
            ax1.set_ylabel("Y")
            plt.colorbar(im)
                     
            myx = np.linspace(np.min(ab['x']), np.max(ab['x']), 15)
            myy = np.linspace(np.min(ab['y']), np.max(ab['y']), 15)

            diffsx = []
            diffsy = []

            for i in range(1,len(myx)):
                diffsx.append(np.average(diff[ (ab['x']<myx[i])*(ab['x']>myx[i-1]) ]))
                diffsy.append(np.average(diff[ (ab['x']<myy[i])*(ab['x']>myy[i-1]) ]))
                
            ax2.plot(myy[1:], diffsy)
            ax3.plot(myx[1:], diffsx)
            
            if (not plotdir is None):
                plt.savefig(os.path.join(plotdir, "x_y_depend_%s.png"%b))
                plt.close("all")
            else:	
                plt.show() 

def interpolate_zp(reduced, logfile):
    '''
    Uses the zeropoint coefficients derived from SDSS fields to interpolate the 
    zeropoint for the images that are outside of SDSS field.
    '''        
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    a.sort(order=['jd'], axis=0)
    a = a[a['inst']!=0]
    a = a[a['insterr']>0]
    
    jdmin =  np.min(a['jd'])
    
    jdmax = np.max(a['jd'])
    
    zpfiles = glob.glob(os.path.join(reduced, "*fits"))
    
    zpfiles = [zf for zf in zpfiles if fitsutils.has_par(zf, "IQZEROPT") and \
    (fitsutils.get_par(zf, "IQZEROPT")==0 or fitsutils.get_par(zf, "ZEROPT")==0 or fitsutils.get_par(zf, "ZPCAT") == "SDSSinterpolated") ]
    
    #Load the coefficients.
    coefs = {}    
    for fi in ["u", "g", "r", "i"]:
        coeffile = os.path.join(reduced, "coefs_%s.txt"%fi)
        if os.path.isfile(coeffile):
            coefs[fi] = np.genfromtxt(coeffile)

    #Load the rms.
    rms = {}    
    for fi in ["u", "g", "r", "i"]:
        rmsfile = os.path.join(reduced, "rms_%s.txt"%fi)
        if os.path.isfile(rmsfile):
            rms[fi] = np.genfromtxt(os.path.join(reduced, "rms_%s.txt"%fi))
        
    for image in zpfiles:
        filt = fitsutils.get_par(image, "FILTER")
        
        #To not extrapolate outside of the valid interval.
        jd = np.maximum(np.percentile(a['jd']-jdmin, 10), fitsutils.get_par(image, "JD") - jdmin)
        jd = np.minimum(np.percentile(a['jd']-jdmin, 90), fitsutils.get_par(image, "JD") - jdmin)
        airmass = fitsutils.get_par(image, "AIRMASS")
        x= a['x']
        y = a['y']
        fwhm = a['fwhm']
        #If there are coefficients for that filter, load them and interpolate.
        #Otherwise, skip this file.
        if not coefs.has_key(filt):
            continue
        
        #est_zp = coef[0] +ab['color']*coef[1] +(ab['airmass']-1.3)*coef[2] + coef[3]*ab['jd'] + coef[4]*ab['jd']**2 + coef[5]*ab['jd']**3 + coef[6]*ab['jd']**4 + coef[7]*ab['jd']**5

        
        values = np.array([1, 0, airmass-1.3, jd, jd**2, jd**3, x, y, x**2, y**2, x**3, y**3, fwhm-1.5])
        est_zp = np.sum(coefs[filt][0:-1]*values)
        
        #Update the header with the computed zeropoint.
        pardic = {
                "IQZEROPT" : 1,\
                "ZPCAT" : "SDSSinterpolated",\
                "ZEROPTU" : float(rms[filt]),\
                "ZEROPT" : est_zp}
        fitsutils.update_pars(image, pardic)
    
    
def lsq_zeropoint_partial(logfile, plot=True):
    '''
    Uses least squares approach to compute the optimum coefficients for ZP, colour term, airmass and time.
    
    '''
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    
    a = a[a['insterr']<0.1]
    a['jd'] = a['jd'] - np.min(a['jd'])
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}


    for b in set(a['filter']):
        ab = a[a['filter']==b]
        M = np.zeros((len(ab), 2))
        M[:,0] = 1      
        #M[:,1] = ab['inst']
        M[:,1] = ab['color'] 

        
        #print M
        
        depend = ab['std']-ab['inst']
        
        lsq_result = np.linalg.lstsq(M, depend)
        coef = lsq_result[0]
        res = lsq_result[1]
                
        print ("Band", b, "zp %.2f, col %.2f"%(coef[0], coef[1]), "Residuals", res)
        
        if (plot):
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            ax1.plot(ab['color'], depend -coef[0], "o", color=cols[b])
            ax1.plot(ab['color'],  ab['color']*coef[1] , color=cols[b])
            ax1.set_xlabel("color")
     
    if (plot):
        plt.show()
        
def find_zeropoint_noid(ref_stars, image, plot=True, plotdir="."):
    '''
    Finds the zeropoint by comparig the magnitude of the stars measured in the image,
    vs. the magnitude of the stars from the reference catalogue.
    band: band for which we want to measure the zeropoint.
    col_band: band for the colour term we want to use.
    
    returns the zeropoint, the colour term and the standard deviation in the 
    zeropoint from all the measured stars.
    
    '''
    
    logger.info( 'Finding the optimum ZP fit...')
    
    logger.info("Reference stars used: %s"% ref_stars)
    
    r = np.genfromtxt(ref_stars, delimiter=" ", dtype=None, names=True)
    imapp = os.path.join(os.path.join(os.path.dirname(image), "zeropoint"), os.path.basename(image) + ".app.mag")

    my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"),  ("image","|S20"), ("X","<f4"), ("Y","<f4"), ("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("msky","<f4"), ("stdev","<f4"),\
        ("flags", np.int), ("rapert", "<f4"), ("sum", "<f4"), ("area", "<f4"), ("nsky","<f4") , ("flux", "<f4"), ("itime", "<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])
    logger.info("Reading aperture magnitude from %s file. %d stars found here.\n The reference stars are in file %s, of length %d"%(imapp, len(my), ref_stars, len(r)))
    '''try:
    except:
        my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"), ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])
    '''
    
    if (my.size < 2):
        my = np.array([my])
    if (r.size < 2):
        r = np.array([r])
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'z', 'z':'i', 'U':'B', 'B':'V', 'V':'R', 'R':'I', 'I':'R'}
    band = fitsutils.get_par(image, 'filter')
    col_band = coldic[band]
    
    mask_valid1 = np.array(my["fwhm"]<9000) * np.array(my["fit_mag"]<9000) * np.array(~ np.isnan(r[band])) * np.array(~ np.isnan(my["fit_mag"]))
    
    N = len(r)
    r = r[mask_valid1]
    my = my[mask_valid1]
    
    if len(my) == 0:
        logger.warn( "Warning, no reliable measurements for file %s. Returning 0,0,0."% image)
        return 0, 0, 0
        
    my["fiterr"][np.isnan( my["fiterr"])] = 100
    
   
    mask_color = np.abs(r[band]-r[col_band]) < 0.8

    ids = np.arange(N)+1
    ids = ids[mask_valid1]

    if np.count_nonzero(mask_color) > 2:
    	r = r[mask_color]
    	my = my[mask_color]
        ids = ids[mask_color]
 

    coefs, residuals, rank, singular_values, rcond = np.polyfit(r[band]-r[col_band], r[band] - my["fit_mag"], w=1./np.maximum(0.1, np.sqrt(my["fiterr"]**2 + r['d'+band]**2)), deg=1, full=True)
    p = np.poly1d(coefs)
    
    logger.info("Coefficients for 1 deg polynomial fit to the zeropoint: %s"% coefs)
    logger.info("%s - %s = %.3f, %.3f"%(band, col_band, p[0], p[1]))


    pred = p(r[band]-r[col_band])
    measured = r[band]- my["fit_mag"]
    mad = stats.funcs.median_absolute_deviation(pred-measured)
    
    #print ("MAD1",  mad)
    
    if len(r) > 4:
        mask = np.abs(pred-measured)/mad < 3
            
        r = r[mask]
        my = my[mask]
        ids = ids[mask]
        
        try:
		coefs, residuals, rank, singular_values, rcond = np.polyfit(r[band]-r[col_band], r[band] - my["fit_mag"], w=1./np.maximum(0.2, np.sqrt(my["fiterr"]**2 + r['d'+band]**2)), deg=1, full=True)
        except:
		logger.error("Critical failure for polyfit for object %s"%image)
		return 0,0,0
        
        logger.info("Coefficients for 1 deg polynomial fit to the zeropoint: %s. After outlier rejection."% coefs)
        p = np.poly1d(coefs)
    
        pred = p(r[band]-r[col_band])
        measured = r[band]- my["fit_mag"]
        mad = stats.funcs.median_absolute_deviation(pred-measured)
        #print ("MAD2",  mad)

        
    if (plot):
        plt.errorbar(r[band]-r[col_band], r[band] - my["fit_mag"] , yerr=np.sqrt(my["fiterr"]**2 + r['d'+band]**2), marker="o", ls="None")
        for i, myid in enumerate(ids):
            plt.text(r[band][i]-r[col_band][i] + 0.01, r[band][i] - my["fit_mag"][i]+ 0.01, str(myid))
        x = np.linspace(np.min(r[band]-r[col_band]), np.max(r[band]-r[col_band]), 100)
        plt.plot(x, p(x))
        plt.title("Best fit ZP: %.2f colour term: %.2f MAD: %.2f"%(p[0], p[1], mad))
        plt.xlabel("{:} - {:}".format(band, col_band))
        plt.ylabel("ZP")
        plotname = os.path.join( plotdir, os.path.basename(image).replace('.fits', ".zp.png").replace('.new', ".zp.png"))
	plt.tight_layout()
        plt.savefig( plotname)

        print ("Plotting into %s."%plotname)
        plt.clf()
    
    logger.info("%s - %s = %.3f, %.3f"%(band, col_band, p[0], p[1]))
    
    pred = p(r[band]-r[col_band])
    measured = r[band]- my["fit_mag"]

    fitsutils.update_par(image, "ZP", np.round(p[0], 3))
    fitsutils.update_par(image, "COLTERM", np.round(p[1], 3))
    fitsutils.update_par(image, "ZPERR", np.round(mad, 3)) #np.std(pred - measured))


    return np.round(p[0], 3), np.round(p[1], 3), np.round(mad, 3)


def calibrate_zeropoint(image, plot=True, plotdir=None, astro=False, debug=False, refstars=None, minexp=1):
    '''
    Calibrates the zeropoint using SDSS catalogue.    
    '''
    
    
    if plot and plotdir is None:
        plotdir = os.path.join(os.path.dirname(image), "zeropoint")
        if not os.path.isdir(plotdir):
            print ("Creating diretory", plotdir)
            os.makedirs(plotdir)
        
    filt = fitsutils.get_par(image, 'filter')
    exptime = fitsutils.get_par(image, "exptime")
    if fitsutils.has_par(image, "JD"):
        date = fitsutils.get_par(image, "JD")
    elif fitsutils.has_par(image, "MJD"):
        date = fitsutils.get_par(image, "MJD")
    elif fitsutils.has_par(image, "MJD-OBS"):
        date = fitsutils.get_par(image, "MJD-OBS")
    else:
        date=0
        
    if fitsutils.has_par(image, "AIRMASS"): 
        airmass = fitsutils.get_par(image, "AIRMASS")
    else:
        airmass = 1.3
    objname = fitsutils.get_par(image, "OBJECT")
    band = fitsutils.get_par(image, "FILTER")
        
    
    logger.info( "Starting calibration of ZP for image %s for object %s with filter %s."%(image, objname, band))

    if (exptime < minexp):
        logger.error( "ERROR. Exposure time too short for image (%s) to see anything..."%image)
        return 


    #repeat astrometry
    if(astro):
        rcred.solve_astrometry(image, radius=3, with_pix=True, overwrite=True, tweak=3)

    calcat = ""
    if (filt == 'u'):
        calcat = "SDSS"
        extracted = extract_star_sequence(os.path.abspath(image), filt, plot=plot, survey='sdss', debug=debug, refstars=refstars, plotdir=plotdir)
    else:
        calcat = "PS1"
        extracted = extract_star_sequence(os.path.abspath(image), filt, plot=plot, survey='ps1', debug=debug, refstars=refstars, plotdir=plotdir)
        
    if (not extracted):
        logger.warn( "Field+filter not in SDSS or error when retrieving the catalogue... Skipping. Image %s not zeropoint calibrated."%image)
        #Add these values to the header.
        pardic = {"IQZEROPT" : 0,\
            "ZPCAT" : "None",\
            "ZEROPTU" : 0,\
            "ZEROPT" : 0, \
            "ZP":0,\
            "ZPERR":0,\
            "%s_coef"%filt : 0}
        fitsutils.update_pars(image, pardic)
        return 

    #If extraction worked, we can get the FWHM        
    fwhm = fitsutils.get_par(image, "fwhm")
    if (not fwhm is None):
	fwhm_pix = fwhm / 0.394
    else:
	fwhm = 0
	fwhm_pix = 0
	logger.error("No FWHM value has been computed for soruce %s."%image)

    #Run aperture photometry on the positions of the stars.
    app_phot.get_app_phot("/tmp/sdss_cat_det_%s.txt"%creationdate, image, wcsin='logic', plotdir=plotdir, box=15)
    
    #Compute the zeropoint for the specific image.
    z, c, err = find_zeropoint_noid("/tmp/sdss_cat_det_%s.txt"%creationdate, image, plot=plot, plotdir=plotdir)
    
    #Add these values to the header.
    pardic = {"IQZEROPT" : 1,\
            "ZPCAT" : calcat,\
            "ZEROPTU" : np.round(err, 3),\
            "ZEROPT" : np.round(z, 3),
            "%s_coef"%filt : c}
    fitsutils.update_pars(image, pardic)
            
    #Log the current zeropoint for this image
    logname = os.path.join(os.path.dirname(image), "zeropoint.log")
    
    #Add the data to a later stage zeropoint calibrtion with all-sky data.
    zplogname = os.path.join(os.path.dirname(image), "allstars_zp.log")
    
    add_to_zp_cal("/tmp/sdss_cat_det_%s.txt"%creationdate, image, zplogname)


    if (not os.path.isfile(logname)):
            with open( logname, "a") as f:
                f.write("#filename,exptime,filter,date,airmass,fwhm_pix,fwhm_as,zeropoint,color,err\n")
    with open( logname, "a") as f:
        f.write("%s,%.1f,%s,%3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n"%(image,exptime,filt,date,airmass,fwhm_pix,fwhm,z,c,err))
        
    clean_tmp_files()
    
        
def plot_zp(zpfile, plotdir=None):
    import datetime
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}
    
    def mkdate(text):
        return datetime.datetime.strptime(text, '%Y-%m-%d %H:%M:%S.%f') 
    
    day_frac_diff = datetime.timedelta(np.ceil((datetime.datetime.now() - datetime.datetime.utcnow() ).total_seconds())/3600/24)


    xfmt = md.DateFormatter('%H:%M')
    

    plt.clf()
    a = np.genfromtxt(zpfile, names=True, dtype=None, delimiter=',')
    
    #We add 5h to the UTC date, so it alwasy keeps the date of the end of the night.
    day = time_utils.jd2utc(a['date'][-1], string=True).split()[0]

    
    for fi in ['u', 'g', 'r', 'i']:
        print ("Found %d points for filter %s"%(len(a[a['filter']==fi]), fi))
        for i in range(len(a[a['filter']==fi])):
            jd = a[a['filter']==fi]['date'][i]
            datestat = time_utils.jd2utc(jd) + day_frac_diff
    
            plt.errorbar( datestat, a[a['filter']==fi]['zeropoint'][i], yerr=a[a['filter']==fi]['err'][i], marker='o', mfc=cols[fi], \
                mec='k', ecolor=cols[fi], ls='none', ms=a[a['filter']==fi]['fwhm_pix'][i]*2)
        logger.info( "Median zeropoint for filter %s: %.2f mag"%(fi, np.median(a[a['filter']==fi]['zeropoint'])))

    plt.gca().xaxis.set_major_formatter(xfmt)

    plt.gca().invert_yaxis()
    plt.title("ZP for day: %s"%day)
    plt.xlabel("Time")
    plt.ylabel("ZP [mag]")
    plt.ylim(24,20.5)
    if (plotdir is None):
        plt.show()
    else:
	plt.savefig(os.path.join(plotdir, "zeropoint_per_exposure.png"))
	plt.clf()
    
def plot_zp_airmass(zpfile):
    import datetime
    
    def mkdate(text):
        return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S') 
    
    a = np.genfromtxt(zpfile, names=True, dtype=None, delimiter=',')
    
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}
    
    for fi in ['u', 'g', 'r', 'i']:
        for i in range(len(a[a['filter']==fi])):
            plt.errorbar( a[a['filter']==fi]['airmass'][i], a[a['filter']==fi]['zeropoint'][i], yerr=a[a['filter']==fi]['err'][i], marker='o', mfc=cols[fi], mec='k', ecolor=cols[fi], ls='none', ms=a[a['filter']==fi]['fwhm_pix'][i])
        logger.info( "Median zeropoint for filter %s: %.2f mag"%(fi, np.median(a[a['filter']==fi]['zeropoint'])))

    plt.gca().invert_yaxis()
    plt.xlabel("Airmass")
    plt.ylabel("ZP [mag]")
    plt.show()


def reset_zp(directory):
    '''
    Resets the zeropoint keyword to if the result is interpolated.
    '''    
    
    files = glob.glob(directory+"/*fits")
    
    for f in files:
        if (fitsutils.get_par(f, "ZPCAT")=="SDSSinterpolated"):
            fitsutils.update_par(f, "IQZEROPT", 0)
            fitsutils.update_par(f, "ZEROPT", 0)
            fitsutils.update_par(f, "ZEROPTU", 0)
    
    
def main(reduced):
    '''
    Performs the main zeropoint calculations for the folder and plots the results.
    
    '''
    os.chdir(reduced)
    
    plotdir = "zeropoint"
    if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)
    

    for f in glob.glob("*_a_*.fits"):
        logger.info("Starting calibration of zeropoint for %s"% f)
        if (not fitsutils.has_par(f, "IMGTYPE") or fitsutils.get_par(f, "IMGTYPE").upper() == "SCIENCE"):
            calibrate_zeropoint(f, plotdir=plotdir)
    if (os.path.isfile("zeropoint.log")):
        plot_zp("zeropoint.log", plotdir)
    if (os.path.isfile("allstars_zp.log")):
        lsq_zeropoint("allstars_zp.log", plotdir)
        #interpolate_zp(reduced, "allstars_zp.log")
        
    clean_tmp_files()
     
if __name__ == '__main__':
    

    parser = argparse.ArgumentParser(description=\
        '''

        Computes the zeropoints for all the images in the folder.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.', default=None)

    args = parser.parse_args()
    
    photdir = args.photdir
    
    if (photdir is None):
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-","")
        photdir = os.path.join(_photpath, timestamp)

    main(os.path.join(os.path.abspath(photdir), "reduced"))



