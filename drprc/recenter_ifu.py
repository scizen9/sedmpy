# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:30:32 2016

@author: nadiablago
"""
import fitsutils
import subprocess, os, sys
from astropy.io import fits as pf
from astropy.wcs import WCS
import coordinates_conversor as cc
import numpy as np
import argparse
from matplotlib import pylab as plt
import matplotlib.patches as patches
import zscale
import matplotlib
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import fit_utils
import datetime
import logging


from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_logpath = parser.get('paths', 'logpath')
vim  = parser.get('paths', 'photpath')


#Log into a file
FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
root_dir = _logpath
timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
timestamp=timestamp.split("T")[0]
logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "listener_{0}.log".format(timestamp)), level=logging.INFO)
logger = logging.getLogger('listener')



def solve_astrometry(img, radius=3.0, with_pix=True, tweak=3):
    '''
    img: fits image where astrometry should be solved.
    radius: radius of uncertainty on astrometric position in image.
    '''

    try:
        ra = fitsutils.get_par(img, 'OBJRA')
        dec = fitsutils.get_par(img, 'OBJDEC')
    except:
        ra = fitsutils.get_par(img, 'RA')
        dec = fitsutils.get_par(img, 'DEC')

    mydir = os.path.dirname(img)
    if mydir=="":
        mydir = "."
    os.chdir(mydir)
    
    astro = os.path.join( os.path.dirname(img), "a_" + os.path.basename(img))
    
    print "Solving astrometry on field with (ra,dec)=", ra, dec, "Image",img, "New image", astro
    
    cmd = " solve-field --ra %s --dec %s --radius %.4f -p --new-fits %s --cpulimit 45 -W none -B none -P none -M none -R none -S none -t %d --overwrite %s "%(ra, dec, radius, astro, tweak, img)
    if (with_pix):
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.425"
    print cmd
    logger.info(cmd)
    
    cmd = cmd + " > /tmp/astrometry_fail  2>/tmp/astrometry_fail"

    subprocess.call(cmd, shell=True)
    
    #Cleaning after astrometry.net
    if (os.path.isfile(img.replace(".fits", ".axy"))):
        os.remove(img.replace(".fits", ".axy"))
    if (os.path.isfile(img.replace(".fits", "-indx.xyls"))):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if (os.path.isfile("none")):
        os.remove("none")
        
    try:
        os.remove("/tmp/tmp.*")
    except:
        pass
        
        
    if(not os.path.isfile(astro)):
        print "Astrometry FAILED!"
        logger.error("Astrometry FAILED!")
        
    return astro
    
def get_offset_center_failed_astro(f, plot=False, interactive=True):
    '''
    For fields where astrometry is challenging, there is a simple solution.
    Find the brightest peak within the pointing error of the telescope.
    As this fields will usually be centered in Standard stars and very short exposure time,
    the fit intends to 
    
    '''

    image = pf.open(f)
    data = image[0].data
    wcs = WCS(image[0].header)
    ra, dec = cc.hour2deg(image[0].header['OBJRA'], image[0].header['OBJDEC'] )

    pra, pdec = wcs.wcs_sky2pix(ra, dec, 0)
    #Get a local image
    #xl, yl = np.array(wcs.wcs_sky2pix(ra+(60./3600)*np.cos(np.deg2rad(dec)), dec-60./3600, 0), dtype=np.int)
    #xu, yu = np.array(wcs.wcs_sky2pix(ra-(60./3600)*np.cos(np.deg2rad(dec)), dec+60./3600, 0), dtype=np.int)
    imageloc = image[0].data.T[1293-150:1293+150,1280-150:1280+150]
    
    nx = 300
    ny=300
            
    def_x = np.argmax(np.sum(imageloc, axis=0))
    def_y = np.argmax(np.sum(imageloc, axis=1))
    
    newx = pra-nx/2.+def_x
    newy = pdec-ny/2.+def_y
    
    
    pra, pdec = wcs.wcs_pix2sky(np.array([[newx[0], newy[0]]] , np.float_), 0)[0]
    dra, ddec = cc.get_offset(ra, dec, pra, pdec)
    
    print "Offset", dra, ddec, "Position RA,DEC", pra, pdec

    x,y, fwhmx, fwhmy, bkg, amp = fit_utils.fit_gauss(imageloc)   
    
    if (plot):
        plt.figure(figsize=(8,8))
        obj = fitsutils.get_par(f, "OBJECT")

        plt.suptitle(obj, fontsize=20)
        zmin, zmax = zscale.zscale(imageloc)
        plt.imshow(imageloc, aspect="auto", interpolation="none", origin="lower", vmin=zmin, vmax=zmax, extent=(0,+300,0,+300))
        plt.plot(x, y, "go", ms=20, label="Centroid using gaussiuan fit.")    
        plt.plot(def_x, def_y, "b*", ms=20, label="Centroid using max/min.")
        plt.plot(150,150,"wo", ms=20, label="Initial pointing")
        plt.legend()
            
        '''zmin, zmax = zscale.zscale(data)
        plt.imshow(data, aspect="auto", interpolation="none", origin="lower", vmin=zmin, vmax=zmax)
        plt.plot(newx, newy, "go")    
        plt.show()'''
        
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(f).replace("raw", "phot"), os.path.basename(f).replace(".fits", "_std.png")))

        plt.clf()
        
    
    return 1, ddec, dra
            
    
def get_offset_center(f, plot=False, interactive=False):
    '''
    Given a fits image, returns the offset in Ra, DEC, that needs to be applied for the telescope tp go
    from the current pointing position, to the coodinates of the object specified in the fits file.
    '''
    
    if(not os.path.isfile(f)):
        print "File %s does not exist! Returning Zero offsets..."%f
        return -1, 0,0
    else:
        image = pf.open(f)
        wcs = WCS(image[0].header)
        rra, rdec = cc.hour2deg(image[0].header['OBJRA'],image[0].header['OBJDEC'] )
        x, y = np.round(wcs.wcs_sky2pix(rra, rdec, 0), 0)
        pra, pdec = wcs.wcs_pix2sky(np.array([[1293., 1280.]] , np.float_), 0)[0]
        dra, ddec = cc.get_offset(pra, pdec, rra, rdec)
            
        xl, yu = np.round(wcs.wcs_sky2pix(rra+90./3600, rdec-90./3600, 0), 0)
        xu, yl = np.round(wcs.wcs_sky2pix(rra-90./3600, rdec+90./3600, 0), 0)

        imageloc = image[0].data.T[xl:xu,yl:yu]

        if imageloc.shape[0]==0 or imageloc.shape[1]==0:
            logger.warn( "Astrometry has FAILED on this! The object is outside the frame! Resending to the numb astrometric solution")
            logger.error("Astrometry has FAILED on this! The object is outside the frame! Resending to the numb astrometric solution")
            print "Pixels are", xl, xu, yl, yu
            try:
                code, dra, ddec = get_offset_center_failed_astro(f, plot=plot, interactive=interactive)
                return 2, dra, ddec
            except:
                return -1,0,0
        if(plot):
            plt.figure(figsize=(8,8))
            
            zmin, zmax = zscale.zscale(imageloc)

            #print zmin, zmax, imageloc, (xl,xu,yl,yu)
    
            obj = fitsutils.get_par(f, "OBJECT")
            plt.suptitle(obj, fontsize=20)
            plt.imshow(imageloc.T,  extent=(xl[0],xu[0],yl[0],yu[0]), aspect="equal", interpolation="none", origin="lower", vmin=zmin, vmax=zmax)
            plt.plot(1293., 1280., "ws", ms=7, label="Current pointing")
            plt.plot(x, y, "b*", ms=10, label="Target pointing")
            plt.gca().invert_xaxis()
            plt.legend()
            if (interactive):
                plt.show()
            else:
                plt.savefig(os.path.join(os.path.dirname(f).replace("raw", "phot"), os.path.basename(f).replace(".fits", "_a.png")))
            plt.clf()


        return 0, dra, ddec

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))
    
def get_offsets_A_B(f, plot=False, interactive=False):
    '''
    Returns the offsets for A and B, so that when offseting, the images do not overlap.  
    Example fits:/scr2/nblago/Projects/SEDM/data/finders/f_b_a_rPTF15fks_r.fits
    '''
    from scipy import stats
    
    image = pf.open(f)
    data = image[0].data.T
    wcs = WCS(image[0].header)
    ra, dec = cc.hour2deg(image[0].header['OBJRA'], image[0].header['OBJDEC'] )
    obj = fitsutils.get_par(f, "OBJECT")
    pra, pdec = wcs.wcs_sky2pix(ra, dec, 0)
    
   
    #Get a local image
    xl, yl = np.round(wcs.wcs_sky2pix(ra+30./3600, dec-30./3600, 0), 0)
    xu, yu = np.round(wcs.wcs_sky2pix(ra-30./3600, dec+30./3600, 0), 0)
    imageloc = image[0].data[xl:xu,yu:yl]
    
    bkg = np.median(imageloc)
    perc10 = np.percentile(imageloc, 15)
    perc90 = np.percentile(imageloc, 85)
    mask = (image[0].data > perc10) * (image[0].data < perc90)
    bkg_std = np.std(image[0].data[mask])# mad(image[0].data)
    
        
    linestyles = ['solid' , 'dashed' , 'dashdot' , 'dotted', 'solid']
    
    offsets = np.array([[+3, -4], [+3, +4], [+2, +4], [+2, +2]])

    Noff = len(offsets)
    pvalues = np.zeros((Noff, 2))
    
    if (plot):
        fig, axarr = plt.subplots(2, Noff, sharey='row', figsize=(6*len(offsets), 12))
        axarr = np.array(axarr, ndmin=1)
        

    for i, off in enumerate(offsets):
        ra_off = ra - off[0]/3600.
        dec_off = dec - off[1]/3600.

        ra_off_op = ra + off[0]/3600.
        dec_off_op = dec + off[1]/3600.
        
        prao, pdeco = wcs.wcs_sky2pix(ra+ off[0]/3600., dec + off[1]/3600., 0)
        praosym, pdecosym = wcs.wcs_sky2pix(ra- 2*off[0]/3600., dec - 2*off[1]/3600., 0)

        #Extract sample window to check wether the background matches well.
        #number on samples on each side
        ns = 7
        sample = data[prao-ns:prao+ns, pdeco-ns:pdeco+ns]
        sm = np.median(sample)
        sstd = np.std(sample[(sample<6000)*(sample>0)])
        bkg_prob = np.minimum(1 - stats.norm.cdf(sm, bkg, bkg_std), stats.norm.cdf(sm, bkg, bkg_std))


        #Extract sample window to check wether the background matches well.        
        sample = data[praosym-7:praosym+7, pdecosym-7:pdecosym+7]
        smo = np.median(sample)
        sstdo = np.std(sample)
        bkg_probo = np.minimum(1 - stats.norm.cdf(smo, bkg, bkg_std), stats.norm.cdf(smo, bkg, bkg_std))
        
        pvalues[i] = np.array([bkg_prob, bkg_probo])
        

        if(plot):
            
            #Retrieve the image of the object
            xl, yl = wcs.wcs_sky2pix(ra+20./3600, dec-20./3600, 0)
            xu, yu = wcs.wcs_sky2pix(ra-20./3600, dec+20./3600, 0)
            ifuwin = data[xl:xu,yu:yl]
            
            
            #Retrieve the offset A image of the object
            x0, y0 = wcs.wcs_sky2pix(ra_off+20./3600, dec_off-20./3600, 0)
            x1, y1 = wcs.wcs_sky2pix(ra_off-20./3600, dec_off+20./3600, 0)
            ifuwin1 = data[x0:x1,y1:y0]

            nx, ny = ifuwin1.shape
            
            #print nx,ny, ifuwin1.shape
            
            #Retrieve the offset A image of the object
            x0, y0 = wcs.wcs_sky2pix(ra_off_op+20./3600, dec_off_op-20./3600, 0)
            x1, y1 = wcs.wcs_sky2pix(ra_off_op-20./3600, dec_off_op+20./3600, 0)
            ifuwin2 = data[x0:x0+nx,y1:y1+ny]
            
            #Plot the A and B
            zmin, zmax = zscale.zscale(ifuwin)
            zmin1, zmax1 = zscale.zscale(ifuwin1-ifuwin2)

                       

            axarr[0, i].imshow(ifuwin.T, aspect="auto", vmin=zmin, vmax=zmax, extent=(20,-20,-20,20), alpha=0.5)
            axarr[1, i].imshow(ifuwin1.T-ifuwin2.T, aspect="auto", vmin=zmin1, vmax=zmax1, extent=(20,-20,-20,20), alpha=0.5)#, cmap=matplotlib.cm.RdBu_r)  
            #axarr[1, i].imshow(-1 * ifuwin2.T, aspect="auto", vmin=zmin2, vmax=zmax2, extent=(20,-20,-20,20), alpha=0.5)#, cmap=matplotlib.cm.RdBu_r) 

            axarr[0, i].scatter(0, 0, marker="x", s=20)
            axarr[1, i].scatter(0, 0, marker="x", s=20)

            axarr[0, i].set_xlabel("OBJRA")
            axarr[0, i].set_ylabel("OBJDEC")
            axarr[1, i].set_xlabel("OBJRA")
            axarr[1, i].set_ylabel("OBJDEC")
            axarr[0, i].text(19, 17, "Red stats: $\mu=$%.2f, $\sigma=$%.2f, p-value=%.5f"%(sm, sstd, bkg_prob),  bbox=dict(facecolor='white', alpha=0.5))
            axarr[0, i].text(19, 14, "Blue stats: $\mu=$%.2f, $\sigma=$%.2f, p-value=%.5f"%(smo, sstdo, bkg_probo),  bbox=dict(facecolor='white', alpha=0.5))
            axarr[0, i].text(19, 11, "Background stats: $\mu=$%.2f, $\sigma=$%.2f"%( bkg, bkg_std),  bbox=dict(facecolor='white', alpha=0.5))
        
            #r = plt.Circle((prao-pra, pdeco-pdec), 5, facecolor="none", edgecolor="red", lw=2)
            #b = plt.Circle((praosym-pra, pdecosym-pdec), 5, facecolor="none", edgecolor="blue", lw=2)
            r = plt.Circle((2*off[0], + 2*off[1]), 2, facecolor="none", edgecolor="red", lw=3, ls=linestyles[i])
            b = plt.Circle((-2*off[0], -2*off[1]), 2, facecolor="none", edgecolor="blue", lw=3, ls=linestyles[i])
            
            #Plot the true location of the IFU
            patches = []
            Path = mpath.Path
            path_data = [
                (Path.MOVETO, [20, -12]),
                (Path.LINETO, [20, 20]),
                (Path.LINETO, [-13,  13]),
                (Path.LINETO, [-20 , -18]),
                (Path.CLOSEPOLY, [20, 10])
                ]
            codes, verts = zip(*path_data)
            path = mpath.Path(verts, codes)
            patch = mpatches.PathPatch(path)
            patches.append(patch)
            collection = PatchCollection(patches, cmap=plt.cm.YlGn, alpha=0.2)
            colors = np.linspace(0, 1, len(patches))
            collection.set_array(np.array(colors))

            axarr[0, i].add_artist(r)
            axarr[0, i].add_artist(b)
            axarr[1, i].add_collection(collection)


            plt.suptitle(obj, fontsize=20)
    
    prod = np.prod(pvalues, axis=1)
        
    if (plot and interactive):
        axarr[0, np.argmax(prod)].text(0,0,"WINNER")
        plt.show()
    elif(plot):
        axarr[0, np.argmax(prod)].text(0,0,"WINNER")
        plt.savefig(os.path.join(os.path.dirname(f).replace("raw", "phot"), os.path.basename(f).replace(".fits", "_ab.png")))
        plt.clf()
       

    return 0, offsets[np.argmax(prod)], -2*offsets[np.argmax(prod)]
    
def main(infile, isAB, astro=True, plot=False):
    '''
    Computes the coffsets for the acquisition image.
    Returns:
    - Code to indicate if the operation has succeeded.
    - A offset only
    - AB offsets
    
    
    Code meaning:
        -  0: Success. All offsets are ok.
        -  1: Offsets computed with Standard Star approach.
        -  2: Offsets were found to be out of the image.
        - -1: Something got very wrong. The astrumetry file does not exist.
    '''
    offset_file = "/tmp/%s_dra_ddec.txt"%(os.path.basename(infile).replace(".fits", ""))

    print "Found image %s as first acquisition image after the slew. Computing offset for IFU..."%infile
    '''if (os.path.isfile(offset_file)):
        #print "Offset file %s already exists for the image."%offset_file
        return
    else:
	print "Found image %s as first acquisition image after the slew. Computing offset for IFU..."%infile'''
        
    newfile = ""
    retcode = 0
    dra = 0
    ddec = 0
    
    #Comment whenever we have the new astrometry file.
    if (astro):
        try:
            newfile = solve_astrometry(infile)
            retcode, dra, ddec = get_offset_center(newfile, plot=True, interactive=False)
        except Exception as e:
            logger.error(str(sys.exc_info()[0]))
            logger.error(e)
            logger.error("Astrometry failed on file %s. Computing the \"Failed solution option\""%infile)
            newfile = infile.replace("rc_", "a_rc_")
            
        if (not os.path.isfile(newfile)):
            retcode, dra, ddec = get_offset_center_failed_astro(infile, plot=True, interactive=False)
    else:
        retcode, dra, ddec = get_offset_center(infile, plot=True, interactive=False)
        newfile = infile

    
    if ( isAB and os.path.isfile(newfile) ):
        retcode, aoff, boff = get_offsets_A_B(newfile, plot=plot, interactive=False)
    
        aoff[0] += dra
        aoff[1] += ddec

        np.savetxt(offset_file, np.array([("A", "%.2f"%(aoff[0]), "%.2f"%aoff[1]), ("B", "%.2f"%(boff[0]), "%.2f"%(boff[1]))]), fmt="%s")
        logger.info( "Offsets computed for A+B: \n A %.4f %.4f \n B %.4f %.4f"%(aoff[0], aoff[1], boff[0], boff[1]))

        return retcode, aoff[0], aoff[1], boff[0], boff[1]
    elif (isAB and astro and not os.path.isfile(newfile)):
        
        np.savetxt(offset_file, np.array([("A", "%.2f"%(dra), "%.2f"%ddec), ("B", "%.2f"%(3.), "%.2f"%(3.))]), fmt="%s")
        logger.warn("Astrometry file not found. Using default offsets (3,3).")        
        logger.info( "Offsets computed for AB: \n AB %.4f %.4f %.4f %.4f"%(dra, ddec, 3, 3 ))
        return retcode,dra,ddec,3,3 
    else:
        np.savetxt(offset_file, np.array([("CENTER", "%.2f"%dra, "%.2f"%ddec)]), fmt="%s")
        logger.info( "Offsets computed for A: \n A %.4f %.4f"%(dra, ddec))

        return retcode, dra, ddec

    
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
           
        -i image
        -b It is AB shot.
        -a Astrometry is needed.
        -p plot results.
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-i', '--image', type=str, dest="image", help='Fits file with acquisition image.')
    parser.add_argument('-b', '--isAB', action='store_true', dest="isAB", default=False, help='Whether we need A, B pair pointing or only A.')
    parser.add_argument('-a', '--astro', action='store_true', dest="astro", default=False, help='Whether we need to solve the astrometry for the image.')
    parser.add_argument('-p', '--plot', action='store_true', dest="plot", default=False, help='Whether we plot the astrometry for the image.')

    
    args = parser.parse_args()
    
    infile = args.image
    isAB = args.isAB
    astro = args.astro
    plot = args.plot
    
    ret = main(infile, isAB, astro, plot)
    print ret
    logger.info("Returned from offseets with code %s"% ret[0])
    #retcode, offsets 
