# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 13:26:06 2017

@author: nblago
"""
from __future__ import print_function

import urllib
import datetime
from astropy.io import votable
import numpy as np
import os, sys
import logging
import warnings
from numpy.lib import recfunctions as rfn
from astropy import units as u
from astropy.coordinates import SkyCoord

class QueryCatalogue:

            
    def __init__(self, ra=0, dec=0, radius=0, minmag=5, maxmag=23, logger=None):
        
        if (type(ra)==str or type(ra)==str ):
            c = SkyCoord('%s %s'%(ra, dec), unit=(u.hourangle, u.deg),  frame='icrs')
            self.ra = c.ra.value
            self.dec = c.dec.value
        else:
            self.ra = ra
            self.dec= dec
        self.rad = float(radius)
        self.minmag = minmag
        self.maxmag = maxmag
        self.logger = logger
        
        if (logger is None):
            FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
            '''now = datetime.datetime.utcnow()
            timestamp=datetime.datetime.isoformat(now)
            creationdate = timestamp
            timestamp=timestamp.split("T")[0]'''

            console = logging.StreamHandler()
            
            formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            console.setFormatter(formatter)
            
            # add the handler to the root logger
            logging.getLogger('').addHandler(console)
            
            #logging.basicConfig(format=FORMAT, level=logging.INFO)
            self.logger = logging.getLogger('QueryCatalogue')
    
            '''ch = logging.StreamHandler(sys.stdout)
            ch.setFormatter(format)
            self.logger.addHandler(ch)'''
            
            #logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
  



    def query_usnob1(self):
    
        #ra, dec = coordinates_conversor.hour2deg(f[0].header['RA'], f[0].header['DEC'])
        #SEDM FoV is 6.5 arcmin, due to uncertainties in the position, 4 arcmin radius assumed.
        # Download USNO-B1 catalog for the position
    
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
    
    
        catalog_url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-B1&RA=%.5f&DEC=%.5f&SR=%.4f&VERB=1' % (self.ra, self.dec, self.rad)
        self.logger.info( "Downloading USNO-B1 catalog...")
        self.logger.info(catalog_url)
        
        tmp_file = '/tmp/tmp_usnob1_%s.cat'%timestamp
        
        urllib.urlretrieve(catalog_url, tmp_file)
        
        # Read RA, Dec and magnitude from XML format USNO catalog
        catalog = votable.parse_single_table(tmp_file).to_table()
        
        #Clean temporary file.
        if (os.path.isfile(tmp_file)):
            os.remove(tmp_file)

        return catalog.as_array().data
        
        '''cat_ra = catalog.array['RA'].data
        cat_dec = catalog.array['DEC'].data
        cat_R1mag = catalog.array['R1'].data
        cat_R2mag = catalog.array['R2'].data
        cat_B1mag = catalog.array['B1'].data
        cat_B2mag = catalog.array['B2'].data
        cat_I2mag = catalog.array['I2'].data
        
        cat_R1mag[cat_R1mag==0] = np.nan
        cat_R2mag[cat_R2mag==0] = np.nan
        cat_B1mag[cat_B1mag==0] = np.nan
        cat_B2mag[cat_B2mag==0] = np.nan
        cat_I2mag[cat_I2mag==0] = np.nan
        
        Bmag = np.nanmean(np.array([cat_B1mag, cat_B2mag]), axis=0)    
        Rmag = np.nanmean(np.array([cat_R1mag, cat_R2mag]), axis=0)    
        Imag = cat_I2mag
            
        mag = Rmag'''
        
    
            
    def query_apass(self):
        '''
        Queries the APASS catalogue
        '''
        
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
    
        
        catalog_url = 'http://www.aavso.org/cgi-bin/apass_download.pl?ra=%.5f&dec=%.5f&radius=%.4f8&outtype=1' % (self.ra, self.dec, self.rad)
        self.logger.info( "Downloading APASS catalog...")
        self.logger.info(catalog_url)
        
        tmp_file = '/tmp/tmp_apass_%s.cat'%timestamp
        
        urllib.urlretrieve(catalog_url, tmp_file)
        catalog = np.genfromtxt(tmp_file, delimiter=",", names=True)
            
        #Clean temporary file.
        if (os.path.isfile(tmp_file)):
            os.remove(tmp_file)
            
        return catalog
            
    
    def query_sdss(self):
        '''
        Queries the SDSS catalogue. The minmag and maxmag apply to the r-band.
        If there is no SDSS, an empty array will be returned.
        '''
        
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
    
        
        catalog_url='http://skyserver.sdss.org/dr9/en/tools/search/x_radial.asp?ra=%.5f&dec=%.5f&check_type=type&type=6&radius=%.4f&check_u=u&min_u=%.2f&max_u=%.2f&check_g=g&min_g=%.2f&max_g=%.2f&check_r=r&min_r=%.2f&max_r=%.2f&check_i=i&min_i=%.2f&max_i=%.2f&check_z=z&min_z=%.2f&max_z=%.2f&entries=top&topnum=500&format=csv'%\
            (self.ra, self.dec, self.rad*60,self.minmag,self.maxmag,self.minmag,self.maxmag,self.minmag,self.maxmag,self.minmag,self.maxmag,self.minmag,self.maxmag)
        self.logger.info( "Downloading SDSS catalog...")
        self.logger.info( "%s"%catalog_url )
        
        tmp_file = '/tmp/tmp_sdss_%s.cat'%timestamp
        
        urllib.urlretrieve(catalog_url, tmp_file)
        catalog = np.genfromtxt(tmp_file, delimiter=",", names=True)
        
        #Clean temporary file.
        if (os.path.isfile(tmp_file)):
            os.remove(tmp_file)
        
        if len(catalog.dtype) ==1:
            catalog = np.array([], dtype=[('objid', '<f8'), ('run', '<f8'), ('rerun', '<f8'), ('camcol', '<f8'), ('field', '<f8'), ('obj', '<f8'), \
            ('type', '<f8'), ('ra', '<f8'), ('dec', '<f8'), ('u', '<f8'), ('g', '<f8'), ('r', '<f8'), ('i', '<f8'), ('z', '<f8'), \
            ('Err_u', '<f8'), ('Err_g', '<f8'), ('Err_r', '<f8'), ('Err_i', '<f8'), ('Err_z', '<f8')])

        return catalog
            
    def query_catalogue(self, catalogue="PS1V3OBJECTS", filtered=True):
        '''
        Sends a VO query to the PS1 catalogue.
        Filters the result by mangitude.
        
        From: http://gsss.stsci.edu/Software/WebServices.htm
        
        General Catalog Access : http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?Parameters...

        Required Parameter List
        1 of the following 3 queries - VO ConeSearch, BoxSearch, IDsearch
        
            RA=ra(deg) &DEC=dec(deg) &SR=search radius(deg)
            BBOX=raMin(deg),decMin(deg),raMax(deg),decMax(deg)
            ID=catID
        
        Optional Parameters
        
            FORMAT= VOTABLE(default) | HTML | KML | CSV | TSV | JSON | TEXT(limited set of catalogs)
            CATALOG=GSC23(default) | GSC11 | GSC12 | USNOB | SDSS | FIRST | 2MASS | IRAS | GALEX | GAIA | TGAS | WISE
            | CAOM_OBSCORE | CAOM_OBSPOINTING | PS1V3OBJECTS | PS1V3DETECTIONS
            FILENAME=outputname (directs output to file)
            MAXOBJ=n (limits number of entries returned by brightest magnitude)
            MAGRANGE=bright,faint (limits number of entries returned by limits)
            MINDET=n (minimum numbr of detections PanSTARRS only)


        '''
        
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
    

            
        url = "http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?CAT=%s&RA=%.5f&DEC=%.5f&SR=%.5f&MAGRANGE=%.3f,%.3f"%(catalogue, self.ra, self.dec, self.rad, self.minmag, self.maxmag)
        
        u = urllib.urlopen(url)
        
        tmp_file = '/tmp/ps1_cat_%s.xml'%timestamp
        f = open(tmp_file, "w")
        f.writelines(u.readlines())
        f.close()
        
        # Read RA, Dec and magnitude from XML format USNO catalog
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            catalog = votable.parse_single_table(tmp_file).to_table()

        if catalog.as_array() is None:
        	#Clean temporary file.\
        	if (os.path.isfile(tmp_file)):
            		os.remove(tmp_file)
		return None

        catalog = catalog.as_array().data

        print (catalog)

        #If it is PS1, we know what fields we want. 
        #Otherwise, we just return everything.
        if (catalogue == "PS1V3OBJECTS"):
            
            if (filtered):
                #Filter spurious sources/ Objects where the majority of pixels where not masked (QFperfect >=0.9) and likely stars (rmeanpsfmag - rmeankronmag < 0.5)
                catalog = catalog[ (catalog["ng"]>3)*(catalog["nr"]>3)* (catalog["ni"]>3)\
                *(catalog["gQfPerfect"]>=0.95) *(catalog["rQfPerfect"]>=0.95)*(catalog["iQfPerfect"]>=0.95) * (catalog["rMeanPSFMag"] - catalog["rMeanKronMag"] < 0.5)]        
            
            newcat = np.zeros(len(catalog), dtype=[("ra", np.double), ("dec", np.double), ("objid", np.long), ("mag", np.float), \
                ("g", np.float), ("r", np.float), ("i", np.float), ("z", np.float), ("y", np.float), \
                ("Err_g", np.float), ("Err_r", np.float), ("Err_i", np.float), ("Err_z", np.float), ("Err_y", np.float), ("distance", np.double)])
            newcat["objid"] = catalog["objID"]
            newcat["ra"] = catalog["RAmean"]
            newcat["dec"] = catalog["DECmean"]
            newcat["mag"] = catalog["rMeanPSFMag"]
            newcat["g"] = catalog["gMeanPSFMag"]
            newcat["r"] = catalog["rMeanPSFMag"]
            newcat["i"] = catalog["iMeanPSFMag"]
            newcat["z"] = catalog["zMeanPSFMag"]
            newcat["y"] = catalog["yMeanPSFMag"]
            newcat["Err_g"] = catalog["gMeanPSFMagErr"]
            newcat["Err_r"] = catalog["rMeanPSFMagErr"]
            newcat["Err_i"] = catalog["iMeanPSFMagErr"]
            newcat["Err_z"] = catalog["zMeanPSFMagErr"]
            newcat["Err_y"] = catalog["yMeanPSFMagErr"]  
            newcat["distance"] = catalog["distance"]  
        else:
            newcat = catalog
            
        #Clean temporary file.\
        if (os.path.isfile(tmp_file)):
            os.remove(tmp_file)
        
        return newcat
    
    
    
