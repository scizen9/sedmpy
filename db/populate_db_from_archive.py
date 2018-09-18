# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:50:10 2017

@author: nadiablago
"""

from __future__ import print_function 
import os, glob, shutil
import numpy as np
import argparse
from SEDMrph import fitsutils
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import datetime
from astropy.io import fits
import archive_mod_sedmdb  # this file requires being able to submit observations with pre-determined id values, which isn't allowed in the original


#cur_dir = os.getcwd()
#os.chdir('/home/sedm/kpy/SEDMDb')
sdb = archive_mod_sedmdb.SedmDB(dbname='sedmdbtest', host='localhost')
#os.chdir(cur_dir)

program_dict = {}

def ra_to_decimal(ra):
    hms = ra.split(':')
    return float(hms[0])*15+float(hms[1])/4+float(hms[2])/240


def dec_to_decimal(dec):
    dms = dec.split(':')
    if dms[0][0] == '-':
        return -(abs(float(dms[0]))+float(dms[1])/60+float(dms[2])/3600)
    else:
        return float(dms[0])+float(dms[1])/60+float(dms[2])/3600


def fill_par_dic_obs(fitsfile):
    '''
    Parses the fits file to fill the needed parameters for the database
    '''
    pardic_obs = {
    'id': 0,
    'object_id': 0, 
    'request_id': 0, 
    'mjd': 0.0, 
    'airmass': 0.0,
    'exptime': 0.0, 
    'fitsfile': '', 
    'lst': '', 
    'ra': 0, 
    'dec': 0, 
    'tel_ra': '',
    'tel_dec': '', 
    'tel_az': 0, 
    'tel_el': 0, 
    'tel_pa': 0, 
    'ra_off': 0,
    'dec_off': 0, 
    'imtype': '', 
    'camera': ''
            }
            
    f = fits.open(fitsfile)
    
    for k in pardic_obs.keys():
        if k in f[0].header:
            pardic_obs[k] = f[0].header[k]
    pardic_obs['mjd'] = f[0].header['JD'] - 2400000.5
    if fitsutils.has_par(fitsfile, 'ra'):
        pardic_obs['ra'] = ra_to_decimal(f[0].header['RA'])
        pardic_obs['dec'] = dec_to_decimal(f[0].header['DEC'])
    else:
        pardic_obs['ra'] = ra_to_decimal(f[0].header['TEL_RA'])
        pardic_obs['dec'] = dec_to_decimal(f[0].header['TEL_DEC'])
    pardic_obs['camera'] = f[0].header['CAM_NAME']   
    
        
    #Rtrieve the old obs_id from the header
    if fitsutils.has_par(fitsfile, 'obs_id'):
        obs_id = f[0].header['obs_id']
        pardic_obs['id'] = obs_id  # if there is no obs_id, add_observation will add its own

    return pardic_obs


def create_user(f):
    '''
    Reads the email encoded in the header and assigns it to a user.
    It tries to locate a user in te DB with that email.
    If not, a new user is created.
    '''
    email = fitsutils.get_par(f, "EMAIL")
    res = sdb.get_from_users(["id"], {"email":email})
    if len(res) == 1:
        return res[0]
    else:
        username = email.split("@")[0]
        sdb.add_user({'username':username,
                      'name':'Unknown',
                      'email':email,
                      'password':username})
        res = sdb.get_from_users(["id"], {"email":email})
        
        sdb.add_group({'designator': ''})
        return res[0]
    
def create_allocation(f, allocation_desig):
    '''
    Creates an allocation registry in the program table.
    
    '''
    #program_name = fitsutils.get_par(f, "P60PRNM")
    program_pi = fitsutils.get_par(f, "P60PRPI")
    #target_id = fitsutils.get_par(f, "TARGID")
    
    if allocation_desig:
        program_ident = allocation_desig.split('-')[1]
    else:
        return -1
    
    # based on the content after ...- assign the allocation to the correct program
    if program_ident[:2] == 'C0':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002171957334})
    elif program_ident[:2] == 'T0':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002171926014})
    elif program_ident[:2] == 'I0':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002165908206})
    elif program_ident[:2] == 'J0':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002171903862})
    elif program_ident[:2] == 'EN':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 3})
    elif program_ident[:2] == 'CA':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 3})
    elif program_ident[:2] == 'GU':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002163457445})
    elif program_ident[:2] == 'D0':
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002163457445})
    elif program_ident[:2] == "SE":
        al = sdb.add_allocation({'designator': allocation_desig, 'program_id': 20171002163457445})
    return al[0]

def create_sci_request(files, inidate, enddate, scitype, object_id, allocation):
    filt = fitsutils.get_par(files[0], 'FILTER')
    if filt not in ['u','g','r','i']:
        nex = '{0, %s, 0, 0, 0}' % len(files)
    elif filt == 'u':
        nex = '{0, %s, 0, 0, 0}' % len(files)
    elif filt == 'g':
        nex = '{0, 0, %s, 0, 0}' % len(files)
    elif filt == 'r':
        nex = '{0, 0, 0, %s, 0}' % len(files)
    elif filt == 'i':
        nex = '{0, 0, 0, 0, %s}' % len(files)
    
    exptime = fitsutils.get_par(files[0], "EXPTIME")
    pardic = {'object_id':object_id,
              'user_id':2,  #user_id 2 is the admin user
              'allocation_id':allocation,
              'exptime':'{0, %d}'%exptime,
              'priority':.1,
              'inidate': inidate, #('year-month-day') (start of observing window),
              'enddate': enddate, #('year-month-day') (end of observing window),
              'nexposures': nex}
    req = sdb.add_request(pardic)
    if req >= 0:
        return req[0]
    else:
        print(req[1])

                            
def create_cal_request(files, inidate, enddate, caltype, obj_id):
    '''
    Creates a default calibration request assigned to sedmcal, the calibration user.
    It also looks at each file and creates an atomic request associated to each observation.
    '''    
    #First check that there is no request for that night
    res = sdb.get_from_request(['object_id'],
    {'object_id': obj_id, 'user_id':32, 'allocation_id':1, 'inidate':inidate, 'enddate':enddate},
    {'inidate':'>=', 'enddate':'<='})
    
    #If there is no request, we add it.
    if len(res) == 0:
        exptime = fitsutils.get_par(files[0], "EXPTIME")
        pardic = {'object_id':obj_id,
                  'user_id':32,
                  'allocation_id':1,
                  'exptime':'{0, %d}'%exptime,
                  'priority':.1,
                  'inidate': inidate, #('year-month-day') (start of observing window),
                  'enddate': enddate, #('year-month-day') (end of observing window),
                  'nexposures':'{0, %s, 0, 0, 0}' % len(files)}
        req = sdb.add_request(pardic)
        if req[0] == -1:
            print(req[1], files)
        return req[0]
        
def create_obs(reqid, files, inidate, enddate, caltype, obj_id):
    '''
    For each file, creates the observation associated with the file.
    
    '''
    for f in files:
                
        jd_init = fitsutils.get_par(f, "JD")
        jd_end = jd_init + fitsutils.get_par(f, "EXPTIME")/(3600*24.)

        inidate = Time(jd_init, format='jd').iso
        enddate = Time(jd_end, format='jd').iso
                        
        pardic_obs = fill_par_dic_obs(f)
        pardic_obs["object_id"] = obj_id
        pardic_obs["request_id"] = reqid
        pardic_obs['fitsfile'] = f
        pardic_obs['imtype'] = caltype
        
        #Also add the observation

        obs = sdb.add_observation(pardic_obs)
        if obs[0] < 0:
            print(obs[1])


    
def log_calibrations(lfiles, caltype="test"):
    '''
    Logs the calibrations.
    
    These are the default objects of calibration.
    
    (id, marshal_id, name, ra, dec) (3, 0, 'bias', 0,0);
    (id, marshal_id, name, ra, dec) values (4, 1, 'twilight', 0,0);
    (id, marshal_id, name, ra, dec) values (5, 2, 'dome', 0,0);
    (id, marshal_id, name, ra, dec) values (6, 3, 'focus', 0,0);
    (id, marshal_id, name, ra, dec) values (7, 4, 'test', 0,0);
    (id, marshal_id, name, ra, dec) values (8, NULL, 'test2', 0,0);
    (id, marshal_id, name, ra, dec) values (9, NULL, 'test3', 0,0);

    '''
    
    #Retrieve the date of the first file
    lfiles.sort()
    
    jd_init = fitsutils.get_par(lfiles[0], "JD")
    jd_end = fitsutils.get_par(lfiles[-1], "JD")
    
    inidate = Time(jd_init, format='jd').iso
    enddate = Time(jd_end, format='jd').iso
    
    #Select the object_id
    caltype_obj = {'bias': 3, 'twilight':4, 'dome':5, 'focus':6, 'test':7, 'test2':8, 'test3':9}
    obj_id = caltype_obj[caltype]
    
    #First make sure that we have a Calibration Request to which assign all the registers.
    reqid = create_cal_request(lfiles, inidate, enddate, caltype, obj_id)
    
    #Then associate all the files to that request id.
    create_obs(reqid, lfiles, inidate, enddate, caltype, obj_id)

    sdb.update_request({'id': reqid, 'STATUS': 'COMPLETED'})
    
    
    
def log_science(lfiles, scitype):
    '''
    Logs the science files.
    '''
    
    for f in lfiles:
        #Initial time of the request is the starting JD.
        jd_init = fitsutils.get_par(f, "JD")
        #The end time of the request is the starting point plus the exposure plus 60s overhead.
        jd_end = fitsutils.get_par(f, "JD") + (fitsutils.get_par(f, "EXPTIME") + 60)/(24*3600.) 
        
        inidate = Time(jd_init, format='jd').iso
        enddate = Time(jd_end, format='jd').iso
    
        #First make sure that its allocation exists
        if fitsutils.has_par(f, "P60PRID"):
            allocation = fitsutils.get_par(f, "P60PRID").upper()
        else:
            continue  #skip if there is no allocation
        allocations = sdb.get_from_allocation(['designator', 'id'])
        if allocation not in [al[0] for al in allocations]:
            all_id = create_allocation(f, allocation)
        else:
            all_id = [al[1] for al in allocations if al[0] == allocation][0]
        if all_id == -1:  # if there was an issue creat_allocation will return -1
            print("Error finding or creating the allocation for file %s!" % (f,))
            continue  # skip to next file if we can't get a valid allocation
        

        #Next make sure that its object exists
        name = fitsutils.get_par(f, "OBJNAME").replace('"', '').lower()
        RA = ra_to_decimal(fitsutils.get_par(f, "RA"))
        DEC = dec_to_decimal(fitsutils.get_par(f, "DEC"))
        exists = sdb.get_from_object(['id', 'ra', 'dec'],{'name': name})

        # check if there was an object found with the same name and coordinates
        if exists and SkyCoord(ra=exists[0][1]*u.deg, dec=exists[0][2]*u.deg, frame='icrs').separation(SkyCoord(ra=RA*u.deg, dec=DEC*u.deg, frame='icrs')) < .001*u.deg:
            obj_id = exists[0][0]
        else:
            obj_id = sdb.add_object({'name': name, 'typedesig': 'f', 'ra': RA, 'dec': DEC})[0]
        
        if obj_id == -1:
            print("Error adding the object for file %s!" % (f,))
            continue  # skip to next file if we can't get a valid object
        
        # create individual requests for each file
        req_id = create_sci_request([f], inidate, enddate, scitype, obj_id, all_id)
        if req_id == -1:
            continue
        # create observations linked to the requests
        create_obs(req_id, [f], inidate, enddate, scitype, obj_id)
    
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=\
        '''

        Parses a directory with nights of photometric data in sub_directories.
        Obtains the relevant values from the headers and populates the DB.
        
        %run populate_db_from_archive.py -d DIR 
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--dir', type=str, help='Directory containing the nightly science fits directories.', default=None)
    
    args = parser.parse_args()
    directory = args.dir

    """
    parser = argparse.ArgumentParser(description=\
        '''

        Parses a directory with one night photometric data.
        Obtains the relevant values from the headers and populates the DB.
        
        %run populate_db_from_archive.py -d PHOTDIR 
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--photdir', type=str, help='Directory containing the science fits for the night.', default=None)

    args = parser.parse_args()
    
    photdir = args.photdir"""
        
    for photdir in os.listdir(directory):
        print(directory+photdir)
        if os.path.isdir(directory+photdir):
            myfiles = {
                "ACQUISITION":[],
                "BIAS":[],
                "DOME":[],
                "FOCUS":[],
                "GUIDER":[],
                "SCIENCE":[],
                "TWILIGHT":[]}
                
            mydir = os.path.abspath(directory+photdir)
            #Gather all RC fits files in the folder with the keyword IMGTYPE=SCIENCE
            for f in glob.glob(os.path.join(mydir, "rc*fits")):
                try:
                    if (fitsutils.has_par(f, "IMGTYPE")):
                        imgtype = fitsutils.get_par(f, "IMGTYPE").upper()
                        if imgtype == "ACQUISTION":
                            myfiles["ACQUISITION"].append(f)
                        myfiles[imgtype].append(f)
                except:
                    print("problems opening file %s" %f)

            if myfiles["BIAS"]:
                log_calibrations(myfiles["BIAS"], caltype="bias")
            if myfiles["DOME"]:
                log_calibrations(myfiles["DOME"], caltype="dome")
            if myfiles["FOCUS"]:
                log_calibrations(myfiles["FOCUS"], caltype="focus")
            if myfiles["TWILIGHT"]:
                log_calibrations(myfiles["TWILIGHT"], caltype="twilight")
            
            if myfiles["SCIENCE"]:
                log_science(myfiles["SCIENCE"], "science")
            if myfiles["GUIDER"]:
                log_science(myfiles["GUIDER"], "guider")
            if myfiles["ACQUISITION"]:
                log_science(myfiles["ACQUISITION"], "acquisition")
