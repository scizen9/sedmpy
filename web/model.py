from werkzeug.security import generate_password_hash, check_password_hash
from db.SedmDb import SedmDB
from db.SedmDb_tools import DbTools
import datetime
import os
import re
import pandas as pd
import numpy as np
import requests
import glob
from decimal import Decimal
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, Label
from bokeh.models import CDSView, GroupFilter
from bokeh.models.tools import HoverTool
from bokeh.models.ranges import Range1d
from bokeh.models.axes import LinearAxis
from bokeh.models.annotations import BoxAnnotation
from bokeh.plotting import figure
from bokeh.embed import components
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun, get_moon
from scheduler.scheduler import ScheduleNight

superuser_list = ['rsw', 'SEDm_admin', 189, 2, 20180523190352189]

pd.set_option('display.max_colwidth', -1)

request_values = ['id', 'object_id', 'marshal_id',
                  'user_id', 'allocation_id',
                  'exptime', 'priority', 'inidate',
                  'enddate', 'maxairmass',
                  'cadence', 'phasesamples',
                  'sampletolerance', 'filters',
                  'nexposures', 'obs_seq', 'status',
                  'creationdate', 'lastmodified',
                  'max_fwhm', 'min_moon_dist',
                  'max_moon_illum', 'max_cloud_cover',
                  'seq_repeats', 'seq_completed',
                  'last_obs_jd']

object_values = ['id', 'marshal_id', 'name', 'iauname',
                 'ra', 'dec', 'epoch', 'typedesig', 'magnitude']

request_form_values = ['request_id', 'object_id', 'marshal_id',
                       'user_id', 'allocation', 'ifu', 'ifu_use_mag',
                       'ab', 'rc', 'do_r', 'do_g', 'do_i', 'do_u',
                       'r_exptime', 'g_exptime', 'i_exptime', 'u_exptime',
                       'r_repeats', 'g_repeats', 'i_repeats', 'u_repeats',
                       'ifu_exptime', 'priority', 'inidate', 'rc_use_mag',
                       'enddate', 'maxairmass', 'status', 'max_fwhm',
                       'min_moon_dist', 'max_moon_illum', 'max_cloud_cover',
                       'seq_repeats', 'seq_completed']

rc_filter_list = ['r', 'g', 'i', 'u']
schedule = ScheduleNight()

# this all needs to go in some sort of config file instead of changing the source code constantly
computer = os.uname()[1] # a quick fix

if computer == 'pele':
    raw_dir = '/scr7/rsw/sedm/raw/'
    phot_dir = '/scr7/rsw/sedm/phot/'
    redux_dir = '/scr7/rsw/sedm/redux/'
    host = 'pharos.caltech.edu'

elif computer == 'pharos':
    raw_dir = '/scr2/sedm/raw/'
    phot_dir = '/scr2/sedm/phot/'
    redux_dir = '/scr2/sedmdrp/redux/'
    host = 'localhost'

db = SedmDB(host=host, dbname='sedmdb')


def get_db():
    return SedmDB(host=host, dbname='sedmdb')


def tools():
    return DbTools(db)


def get_from_users(user_id):
    return db.get_from_users(['id', 'username'], {'id': user_id})


def get_object_values(objid=None):
    try:
        objid = int(objid)
    except Exception as e:
        return {'error': str(e)}

    ret = db.get_from_object(values=object_values,
                             where_dict={'id': int(objid)})
    if not ret:
        return {'error': "No object found with that id number"}

    return make_dict_from_dbget(object_values, ret[0])


def get_object_info(name=None, ra=None, dec=None, radius=5, out_type='html'):
    """

    :param radius: 
    :param out_type: 
    :param name: 
    :param ra: 
    :param dec: 
    :return: 
    """

    # 1. Start by looking for objects by name or coordinates
    if name:
        ids = db.get_object_id_from_name(name)
    elif ra and dec:
        ids = db.get_objects_near(ra=ra, dec=dec, radius=radius)
    else:
        ids = None

    # 2. If no results return empty list and message
    if not ids:
        return {'message': 'No objects found with that name or coordinates',
                'objects': False}

    # 3. If there are objects get a list of values for each match
    obj_list = []
    for obj in ids:
        ret = db.get_from_object(values=object_values,
                                 where_dict={'id': obj[0]})
        obj_list.append(ret[0])

    # 4. Convert the list into the appropriate output
    if out_type == 'html':
        df = pd.DataFrame(obj_list, columns=object_values)

        df['Select'] = df['id'].apply(add_link)
        return {'message': df.to_html(escape=False, classes='table', index=False)}
    else:
        return obj_list


def get_homepage(userid, username):
    sedm_dict = {'enddate': datetime.datetime.utcnow() + datetime.timedelta(days=1),
                 'inidate': datetime.datetime.utcnow() - datetime.timedelta(days=7, hours=8)}

    # 1. Get a dataframe of all requests for the current user
    requests = get_requests_for_user(userid, sedm_dict['inidate'], sedm_dict['enddate'])

    # organize requests into dataframes by whether they are completed or not
    complete = requests[(requests['status'] == 'COMPLETED') | (requests['status'] == 'REDUCED')]
    active = requests[(requests['status'] == 'PENDING') | (requests['status'] == 'ACTIVE')]
    expired = requests[(requests['status'] == 'EXPIRED')]

    # retrieve information about the user's allocations
    ac = get_allocations_user(userid)

    # Create html tables
    sedm_dict['active'] = {'table': active.to_html(escape=False,
                                                   classes='table',
                                                   index=False),
                           'title': 'Active Requests for the last 7 days'}

    sedm_dict['complete'] = {'table': complete.to_html(escape=False,
                                                       classes='table',
                                                       index=False),
                             'title': 'Completed Requests in the last 7 days'}

    sedm_dict['expired'] = {'table': expired.to_html(escape=False,
                                                     classes='table',
                                                     index=False),
                            'title': 'Expired in the last 7 days'}

    sedm_dict['allocations'] = {'table': ac.to_html(escape=False,
                                                    classes='table table-striped',
                                                    index=False,
                                                    col_space=10),
                                'title': 'Your Active Allocations'}

    sedm_dict['visibility'] = {'title': 'Visibilities for active requests',
                               'url':   '/visibility'}

    # Make a greeting statement
    sedm_dict['greeting'] = 'Hello %s!' % username
    return sedm_dict


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE SCHEDULER PAGE.              #
# KEYWORD:SCHEDULER                                                           #
###############################################################################
def get_schedule(start_time="", end_time="", return_type='table'):
    """

    :param return_type: 
    :param start_time:
    :param end_time:
    :return:
    """

    current_time = datetime.datetime.utcnow()

    if start_time:
        start_time = Time(start_time)
    else:
        start_time = schedule.obs_times['evening_nautical']

    if current_time > start_time.to_datetime():
        start_time = Time(current_time)

    if end_time:
        end_time = Time(end_time)
    else:
        end_time = schedule.obs_times['morning_nautical']

    out = schedule.simulate_night(start_time=start_time, end_time=end_time, return_type=return_type)

    return out


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE REQUEST PAGE.                #
# KEYWORD:REQUEST                                                             #
###############################################################################
def add_link(objid):
    return """<input type = "button" onclick = 'addValues("%s")' 
    value = "Use" />""" % objid


def get_request_page(userid, form1, content=None):
    req_dict = {}

    # Start by getting all the allocations the user can add targets under
    # Create a list for a select field on the request page

    alloc = get_allocations_user(userid)

    if alloc is None or len(alloc) == 0:
        choices = [(0, "You have none active!")]
    else:
        choices = [z for z in zip(alloc['id'], alloc['allocation'])]
        choices.insert(0, (0, '------'))

    form1.allocation.choices = choices

    # if an object id is given or a request id, prepopulate the field
    if content:
        req_dict.update(populate_form(content, form1))

    # Set the start date if one is not given
    if not form1.inidate.data:
        form1.inidate.data = datetime.datetime.today()

    if not form1.enddate.data:
        form1.enddate.data = datetime.datetime.today() + datetime.timedelta(2)

    #
    return req_dict, form1


def populate_form(content, form):
    """

    :param content: 
    :param form: 
    :return: 
    """
    # We look to see if there is a request id to begin with because if so then
    # the object info will automatically be tied into that.  If it is not in
    # there then we should start looking at other keywords.
    if 'request_id' in content:

        data = db.get_from_request(values=request_values,
                                   where_dict={'id': content['request_id'][0]})

        ret_dict = make_dict_from_dbget(headers=request_values, data=data[0])
        obj_dict = get_object_values(objid=ret_dict['object_id'])

        # Remove the identical id so as not to get confused
        del obj_dict['id']
        del obj_dict['marshal_id']

        ret_dict.update(obj_dict)
        ret_dict['request_id'] = content['request_id'][0]

        # I am setting a status_id because when I do
        # form.status.data = x['status']
        # it doesn't set the proper select option.
        ret_dict['status_id'] = ret_dict['status']

        ret_dict.update(parse_db_target_filters(ret_dict['obs_seq'], ret_dict['exptime']))
    elif 'object_id' in content:
        ret_dict = get_object_values(objid=content['object_id'][0])

        if 'error' in ret_dict:
            ret_dict['message'] = ret_dict['error']

    else:
        ret_dict = {'message': "There was nothing to process"}

    # TODO: There has to be a better way to do this...
    if 'request_id' in ret_dict:
        form.request_id.data = ret_dict['request_id']
    if 'object_id' in ret_dict:
        form.object_id.data = ret_dict['object_id']
    if 'marshal_id' in ret_dict:
        form.marshal_id.data = ret_dict['marshal_id']
    if 'allocation_id' in ret_dict:
        form.allocation_id.data = ret_dict['allocation_id']
        form.allocation.data = str(ret_dict['allocation_id'])
    if 'ra' in ret_dict:
        form.obj_ra.data = ret_dict['ra']
    if 'dec' in ret_dict:
        form.obj_dec.data = ret_dict['dec']
    if 'epoch' in ret_dict:
        form.obj_epoch.data = ret_dict['epoch']
    if 'magnitude' in ret_dict:
        form.obj_mag.data = ret_dict['magnitude']
    if 'name' in ret_dict:
        form.obj_name.data = ret_dict['name']
    if 'inidate' in ret_dict:
        form.inidate.data = ret_dict['inidate']
    if 'enddate' in ret_dict:
        form.enddate.data = ret_dict['enddate']
    if 'min_moon_distance' in ret_dict:
        form.min_moon_distance = ret_dict['min_moon_distance']
    if 'priority' in ret_dict:
        form.priority.data = ret_dict['priority']
    if 'maxairmass' in ret_dict:
        form.maxairmass.data = ret_dict['maxairmass']
    if 'cadence' in ret_dict:
        form.cadence.data = ret_dict['cadence']
    if 'phasesamples' in ret_dict:
        form.phasesamples.data = ret_dict['phasesamples']
    if 'sampletolerance' in ret_dict:
        form.sampletolerance.data = ret_dict['sampletolerance']
    if 'status_id' in ret_dict:
        form.status.data = ret_dict['status_id']
    if 'max_fwhm' in ret_dict:
        form.max_fwhm.data = ret_dict['max_fwhm']
    if 'seq_repeats' in ret_dict:
        if not ret_dict['seq_repeats']:
            form.seq_repeats.data = 1
        else:
            form.seq_repeats.data = ret_dict['seq_repeats']
    if 'seq_completed' in ret_dict:
        form.seq_repeats.data = ret_dict['seq_completed']
    if 'max_moon_illum' in ret_dict:
        form.max_moon_illum.data = ret_dict['max_moon_illum']
    if 'max_cloud_cover' in ret_dict:
        form.max_cloud_cover.data = ret_dict['max_cloud_cover']
    if 'creationdate' in ret_dict:
        form.creationdate.data = ret_dict['creationdate']
    if 'lastmodified' in ret_dict:
        form.lastmodified.data = ret_dict['lastmodified']
    if 'last_obs_jd' in ret_dict:
        form.last_obs_jd.data = ret_dict['last_obs_jd']
    if 'do_ifu' in ret_dict:
        form.ifu.data = ret_dict['do_ifu']
    if 'ifu_exptime' in ret_dict:
        form.ifu_exptime.data = ret_dict['ifu_exptime']
    if 'ab' in ret_dict:
        form.ab.data = ret_dict['ab']
    if 'do_rc' in ret_dict:
        form.rc.data = ret_dict['do_rc']
    if 'do_r' in ret_dict:
        form.do_r.data = ret_dict['do_r']
    if 'do_i' in ret_dict:
        form.do_i.data = ret_dict['do_i']
    if 'do_g' in ret_dict:
        form.do_g.data = ret_dict['do_g']
    if 'do_u' in ret_dict:
        form.do_u.data = ret_dict['do_u']
    if 'r_exptime' in ret_dict:
        form.r_exptime.data = ret_dict['r_exptime']
    if 'g_exptime' in ret_dict:
        form.g_exptime.data = ret_dict['g_exptime']
    if 'i_exptime' in ret_dict:
        form.i_exptime.data = ret_dict['i_exptime']
    if 'u_exptime' in ret_dict:
        form.u_exptime.data = ret_dict['u_exptime']
    if 'r_repeats' in ret_dict:
        form.r_repeats.data = ret_dict['r_repeats']
    if 'r_repeats' in ret_dict:
        form.r_repeats.data = ret_dict['r_repeats']
    if 'g_repeats' in ret_dict:
        form.g_repeats.data = ret_dict['g_repeats']
    if 'i_repeats' in ret_dict:
        form.i_repeats.data = ret_dict['i_repeats']
    if 'u_repeats' in ret_dict:
        form.u_repeats.data = ret_dict['u_repeats']
    if 'seq_completed' in ret_dict:
        form.seq_completed.data = ret_dict['seq_completed']
    return ret_dict


def parse_db_target_filters(obs_seq, exptime):
    """
    Parse database target scheme

    :param obs_seq: 
    :param exptime: 
    :return:
    """
    # Prep the variables
    rc_filter_list = ['r', 'g', 'i', 'u']

    return_dict = {
        'do_ifu': False, 'ifu_exptime': 0,
        'do_rc': False,
        'do_r': False, 'r_exptime': 0, 'r_repeats': 1,
        'do_g': False, 'g_exptime': 0, 'g_repeats': 1,
        'do_i': False, 'i_exptime': 0, 'i_repeats': 1,
        'do_u': False, 'u_exptime': 0, 'u_repeats': 1,

    }

    # 1. First we extract the filter sequence
    seq = list(obs_seq)
    exptime = list(exptime)

    # 2. Remove ifu observations first if they exist
    index = [i for i, s in enumerate(seq) if 'ifu' in s]

    if index:
        for j in index:
            seq.pop(j)
            return_dict['ifu_exptime'] = int(exptime.pop(j))
            return_dict['do_ifu'] = True
            if return_dict['ifu_exptime'] == 0:
                return_dict['do_ifu'] = False

    # 3. If the seq list is empty then there is no photmetry follow-up
    # and we should exit
    if not seq:
        return return_dict

    # 4. If we are still here then we need to get the photometry sequence
    return_dict['do_rc'] = True

    for i in range(len(seq)):
        flt = seq[i][-1]
        flt_exptime = int(exptime[i])
        flt_repeat = int(seq[i][:-1])

        # 4a. After parsing the indivual elements we need to check that they are
        # valid values
        if flt in rc_filter_list:
            if 0 < flt_exptime < 600:
                if 0 < flt_repeat < 100:
                    return_dict['do_%s' % flt] = True
                    return_dict['%s_exptime' % flt] = flt_exptime
                    return_dict['%s_repeats' % flt] = flt_repeat
        else:
            continue

    return return_dict


def add_object_to_db(content):
    """

    :param content: 
    :return: 
    """

    return_dict = {}

    # Add the object to the database
    if content['obj_ra'] and content['obj_dec']:
        ra = content['obj_ra']
        dec = content['obj_dec']
        if ":" in ra or ":" in dec:
            return_dict['message'] += ("This is embarrassing, I know I said "
                                       "I would take coordinates in this "
                                       "format but I haven't actually "
                                       "implemented it yet...  Please "
                                       "resubmit using degrees and then "
                                       "bug me to fix this...--")
            objid = False
        else:
            if content['obj_epoch']:
                epoch = content['obj_epoch']
            else:
                epoch = 2000

            objdict = {
                'name': content['obj_name'],
                'ra': ra,
                'dec': dec,
                'typedesig': 'f',
                'epoch': epoch,
            }

            if content['obj_mag']:
                objdict['magnitude'] = content['obj_mag']

            objid, msg = db.add_object(objdict)

            if objid == -1:
                return_dict['message'] += msg + ('--For now I am going to '
                                                 'assume that you want to '
                                                 'use this object and will '
                                                 'go ahead with the rest '
                                                 'of the request--')
                objid = msg.split()[-1]

    else:
        return_dict['message'] = ("How am I suppose to add your request "
                                  "if you don't give me any coordinates--")

        objid = False

    return return_dict, objid


def process_request_form(content, form, userid):
    """

    :param content: 
    :param userid: 
    :param form: 
    :return: 
    """
    request_dict = {}
    process_dict = {'message': ''}
    obs_seq_dict = {}

    alloc = get_allocations_user(userid)

    if alloc is None or len(alloc) == 0:
        choices = [(0, "You have none active!")]
    else:
        choices = [z for z in zip(alloc['id'], alloc['allocation'])]
        choices.insert(0, (0, '------'))

    form.allocation.choices = choices

    # 1. Let's start by making sure we have all the information needed for the
    #    object id

    if content['object_id']:
        objid = content['object_id']
    else:
        message, objid = add_object_to_db(content)
        if not objid:
            return {**process_dict, **message}, form
        else:
            if 'message' in message:
                process_dict['meesage'] += message['message']

    request_dict['object_id'] = int(objid)

    # 2. Now let's put together the request by getting all the values into a dictionary
    obs_seq_key_list = ['ifu', 'rc', 'ab', 'do_r', 'do_i', 'do_u', 'do_g',
                        'r_repeats', 'g_repeats', 'i_repeats', 'u_repeats',
                        'r_exptime', 'g_exptime', 'i_exptime', 'u_exptime',
                        'ifu_use_mag', 'rc_use_mag', 'ifu_exptime']
    for key in request_form_values:
        try:
            # This should handle both the case when an object id has already
            # been added and when we had to generate a new one
            if key == 'object_id':
                pass
            # This case will handle new requests when a user id is not given
            elif key == 'user_id' and not content['user_id']:
                request_dict['user_id'] = userid
            elif key == 'allocation':
                request_dict['allocation_id'] = content[key]
            # This section should handle all the observation data such as if
            # we want ifu/rc follow-up and exposure times.  Note that because
            # of the database format we must handle this data outside the request
            # dictionary.
            elif key in obs_seq_key_list:
                if key in content:
                    obs_seq_dict[key] = content[key]
                else:
                    obs_seq_dict[key] = False
            else:
                request_dict[key] = content[key]
        except Exception as e:
            print(str(e), key)

    # 3. Now we need to create the obs_seq and exptime entries
    #    We need to also make sure and add the object magnitude
    #    to calculate exposure times
    if content['obj_mag']:
        obs_seq_dict['obj_mag'] = content['obj_mag']
    else:
        obs_seq_dict['obj_mag'] = 17.5

    filter_dict = (make_obs_seq(obs_seq_dict))

    if 'ERROR' in filter_dict:
        process_dict['message'] += filter_dict['ERROR']
    else:
        process_dict['message'] += filter_dict.pop("proc_message")

        request_dict = {**filter_dict, **request_dict}

        if content['request_id']:

            request_dict['id'] = int(content['request_id'])
            request_dict.pop('request_id')
            for k, v in request_dict.items():
                if not v:
                    request_dict[k] = -1
            ret = db.update_request(request_dict)
        else:

            request_dict.pop('request_id')
            print(request_dict)
            ret = db.add_request(request_dict)

    return process_dict, form


def get_add_csv(user_id, form, content):
    """

    :param user_id: 
    :param form: 
    :param content: 
    :return: 
    """

    return {'test': 'test'}, form


def process_add_csv(content, form, user_id):
    """

    :param content: 
    :param form: 
    :param user_id: 
    :return: 
    """

    return {'test': 'test'}, form


def make_obs_seq(obs_seq_dict):
    """

    :param obs_seq_dict: 
    :return: 
    """
    filters_list = []
    exptime_list = []
    ret_dict = {"proc_message": ""}

    if obs_seq_dict['ifu']:
        # There may be case in the future where people want more than one IFU
        # at a time.  In which case this code will need to be changed.
        if obs_seq_dict['ifu_use_mag']:
            if obs_seq_dict['ifu_exptime'] and int(obs_seq_dict['ifu_exptime']) > 0:
                ret_dict["proc_message"] += ("You should know that you "
                                             "supplied a non-zero value "
                                             "in the ifu exposure time "
                                             "field.  However because you "
                                             "checked the use magnitude box "
                                             "I will be ignoring the supplied "
                                             "value.--")
            try:
                mag = float(obs_seq_dict['obj_mag'])
                if mag == 0:
                    ret_dict['proc_message'] += ("I find it hard to believe "
                                                 "that you really wanted to "
                                                 "observe something zero "
                                                 "magnitude.  So I can't let "
                                                 "this go through. Feel free"
                                                 "to contact me and dispute "
                                                 "this.--")
                    ifu_exptime = False
                else:
                    ifu_exptime = get_filter_exptime('ifu', mag)

            except Exception as e:
                ret_dict['proc_message'] += ("For some reason I couldn't "
                                             "process your magnitude.  If you "
                                             "didn't add one then that is on "
                                             "you. Otherwise there is something "
                                             "wrong with this '%s' value.  For "
                                             "the record here is the error "
                                             "message %s--" % (obs_seq_dict['obj_mag'],
                                                               str(e)))

                ifu_exptime = False
        else:
            try:
                ifu_exptime = int(obs_seq_dict['ifu_exptime'])
                if 0 <= ifu_exptime <= 7200:
                    pass
                else:

                    ret_dict['proc_message'] += ("I don't know what you are "
                                                 "trying to but %s is not an "
                                                 "acceptable IFU exposure time. "
                                                 "It's either less than 0 or "
                                                 "more than two hours.--" %
                                                 str(ifu_exptime))
                    ifu_exptime = False
            except Exception as e:
                ret_dict['proc_message'] += ("There is something wrong with "
                                             "your exposure time value.  '%s' "
                                             "is not a proper value.  Here is "
                                             "the error message return: %s--" %
                                             (obs_seq_dict['ifu_exptime'], str(e)))

                ifu_exptime = False

        if ifu_exptime:
            filters_list.append("1ifu")
            exptime_list.append(str(ifu_exptime))

    if obs_seq_dict['rc']:
        for flt in rc_filter_list:
            if obs_seq_dict['do_%s' % flt]:
                repeats = obs_seq_dict['%s_repeats' % flt]
                if 1 <= int(repeats) <= 100:
                    pass
                else:
                    ret_dict['proc_message'] += ("There is something wrong "
                                                 "with the number of "
                                                 "repeats you have "
                                                 "requested. Forcing it to 1"
                                                 "--")
                    repeats = 1
                if obs_seq_dict['rc_use_mag']:
                    mag = obs_seq_dict['obj_mag']
                    exptime = get_filter_exptime(flt, mag)

                    filters_list.append("%s%s" % (str(repeats), flt))
                    exptime_list.append(str(exptime))
                else:
                    exptime = int(obs_seq_dict['%s_exptime' % flt])
                    if 0 <= exptime <= 600:
                        pass
                    else:
                        ret_dict['proc_message'] += ("The exposure time (%s) "
                                                     "you entered for filter "
                                                     "(%s) makes no sense. If "
                                                     "you entered something "
                                                     "more than 10mins it is "
                                                     "wasting time.  Feel free "
                                                     "to contact me to disput "
                                                     "this--" % (str(exptime),
                                                                 flt))
                        exptime = False

                    if exptime:
                        filters_list.append("%s%s" % (str(repeats), flt))
                        exptime_list.append(str(exptime))

    if not filters_list:
        ret_dict["ERROR"] = "NO FILTERS COULD BE DETERMINE"
        return ret_dict
    else:
        if len(filters_list) == len(exptime_list):
            ret_dict['obs_seq'] = '{%s}' % ','.join(filters_list)
            ret_dict['exptime'] = '{%s}' % ','.join(exptime_list)
        else:
            ret_dict["ERROR"] = ("Filter and exposure time list don't match "
                                 "%s : %s" % (','.join(filters_list),
                                              ','.join(exptime_list)))
    return ret_dict


def get_allocations_user(user_id, return_type=''):
    res = db.execute_sql(""" SELECT a.id, a.designator, p.designator, g.designator, a.time_allocated, a.time_spent
                            FROM allocation a, program p, groups g, usergroups ug
                            WHERE a.program_id = p.id AND p.group_id = g.id 
                            AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d""" % user_id)

    # create the dataframe and set the allocation names to be linked
    if return_type == 'list':
        data = []
        for i in res:
            data.append(i[0])
    else:
        data = pd.DataFrame(res, columns=['id', 'allocation', 'program', 'group', 'time allocated', 'time spent'])

    return data


def get_requests_for_user(user_id, inidate=None, enddate=None):
    """

    :param user_id: 
    :param inidate: 
    :param enddate: 
    :return: 
    """
    if not inidate:
        inidate = datetime.datetime.utcnow() - datetime.timedelta(days=7, hours=8)
    if not enddate:
        enddate = datetime.datetime.utcnow() + datetime.timedelta(days=1)

    request_query = ("""SELECT a.designator, o.name, o.ra, o.dec, r.inidate, r.enddate, r.priority, r.status, r.lastmodified, r.obs_seq, r.exptime, r.id 
                        FROM request r, object o, allocation a 
                        WHERE o.id = r.object_id AND a.id = r.allocation_id  
                            AND ( r.lastmodified >= DATE('%s') AND r.lastmodified <= DATE('%s') )
                            AND r.allocation_id IN
                           (SELECT a.id
                            FROM allocation a, groups g, usergroups ug, users u, program p
                            WHERE ug.user_id = u.id AND ug.group_id = g.id AND u.id = %d AND p.group_id = g.id AND a.program_id = p.id
                            ) ORDER BY r.lastmodified DESC;""" % (inidate, enddate, user_id))

    data = db.execute_sql(request_query)
    data = pd.DataFrame(data,
                        columns=['allocation', 'object', 'RA', 'DEC', 'start date', 'end date', 'priority', 'status',
                                 'lastmodified', 'obs_seq', 'exptime', 'UPDATE'])

    if user_id in superuser_list:
        data['UPDATE'] = data['UPDATE'].apply(convert_to_link)
    else:
        data.drop(columns=['UPDATE', 'RA', 'DEC'])

    return data


def convert_to_link(reqid):
    return """<a href='request?request_id=%s'>+</a>""" % reqid


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE LOGIN PAGE.                  #
# KEYWORD:LOGIN                                                               #
###############################################################################
def check_login(username, password):
    """

    :param username: 
    :param password: 
    :return: 
    """

    user_pass = db.get_from_users(['username', 'password', 'id'], {'username': username})

    if not user_pass:
        return False, 'Incorrect username or password!'

    if check_password_hash(user_pass[0][1], password=password):
        return True, user_pass[0][2]
    else:
        return False, 'Incorrect username or password!!'


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE STATS PAGE.                  #
# KEYWORD:STATS                                                               #
###############################################################################
def get_project_stats(content, user_id=""):
    """

    :param content: 
    :param user_id: 
    :return: 
    """

    # Start by getting all the allocations for a user
    if 'inidate' not in content:
        inidate = None
    else:
        inidate = content['inidate']

    if 'enddate' not in content:
        enddate = None
    else:
        enddate = content['enddate']

    data = get_allocation_stats(user_id, inidate=inidate, enddate=enddate)
    plots = plot_stats_allocation(data)

    script, div = components(plots)
    return {'script': script, 'div': div}


def get_allocation_stats(user_id, inidate=None, enddate=None):
    """
    Obtains a list of allocations that belong to the user and 
    query the total allocated name and time spent for that allocation.

    If no user_id is provided, all active allocations are returned.
    """
    if (user_id is None):
        res = db.get_from_allocation(["designator", "time_allocated", "time_spent"], {"active": True})
        df = pd.DataFrame(res, columns=["designator", "time_allocated", "time_spent"])

        alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
        spent_hours = np.array([ts.total_seconds() / 3600. for ts in df["time_spent"]])
        free_hours = alloc_hours - spent_hours

        df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)

    else:
        if (inidate is None or enddate is None):
            res = db.execute_sql(""" SELECT a.designator, a.time_allocated, a.time_spent
                                    FROM allocation a, program p, groups g, usergroups ug
                                    WHERE a.program_id = p.id AND p.group_id = g.id 
                                    AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d""" % (user_id))

            df = pd.DataFrame(res, columns=["designator", "time_allocated", "time_spent"])

            alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
            spent_hours = np.array([ts.total_seconds() / 3600. for ts in df["time_spent"]])
            free_hours = alloc_hours - spent_hours

            df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)


        else:
            res = db.execute_sql(""" SELECT DISTINCT a.id, a.designator, a.time_allocated
                                    FROM allocation a, program p, groups g, usergroups ug
                                    WHERE a.program_id = p.id AND p.group_id = g.id 
                                    AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d;""" % (user_id))
            allocdes = []
            spent_hours = []
            alloc = []
            for ais in res:
                spent = db.get_allocation_spent_time(ais[0], inidate, enddate)
                allocdes.append(ais[1])
                spent_hours.append(int(spent) / 3600.)
                alloc.append(ais[2])
            res = np.array([allocdes, alloc, spent_hours])

            df = pd.DataFrame(res.T, columns=["designator", "time_allocated", "time_spent"])
            alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
            free_hours = alloc_hours - spent_hours
            df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)

    df = df.sort_values(by=["alloc_hours"], ascending=False)

    alloc_names = df["designator"].values
    category = ["alloc_hours", "spent_hours", "free_hours"]

    data = {'allocations': alloc_names}

    for cat in category:
        data[cat] = df[cat]

    return data


def plot_stats_allocation(data):
    """
    Plots in the shape of bars the time available and spent for each active allocation.
    """

    data = {key:np.nan_to_num(data[key]) for key in data}

    # Create the first plot with the allocation hours
    alloc_names = data['allocations']
    categories = ["spent_hours", "free_hours"]
    colors = ["#e84d60", "darkgreen"]  # "#c9d9d3"

    N = len(alloc_names)

    source = ColumnDataSource(data=data)
    p = figure(x_range=alloc_names, plot_height=420, plot_width=80 * 8,
               title="Time spent/available for SEDM allocations this term",
               toolbar_location=None, tools="")

    p.vbar_stack(categories, x='allocations', width=0.9, color=colors, source=source, legend=["Spent", "Available"])
    p.y_range.start = 0
    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    p.outline_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "horizontal"
    p.yaxis.axis_label = 'Hours'
    p.xaxis.major_label_orientation = 0.3

    # Create the second plot with the % spent
    alloc_names = data['allocations']
    percentage = (data["spent_hours"] / data["alloc_hours"]) * 100

    colors = N * ['#084594']
    '''for i, p in enumerate(percentage):
        if p<50: colors[i] = '#22A784'
        elif p>50 and p<75: colors[i] = '#FD9F6C'
        else: colors[i] = '#DD4968'''

    source = ColumnDataSource(data=dict(alloc_names=alloc_names, percentage=percentage, color=colors))

    p2 = figure(x_range=alloc_names, y_range=(0, 100), plot_height=420, plot_width=80 * 8,
                title="Percentage of time spent",
                toolbar_location=None, tools="")

    p2.vbar(x='alloc_names', top='percentage', width=0.9, color='color', source=source)

    p2.xgrid.grid_line_color = None
    p2.legend.orientation = "horizontal"
    p2.legend.location = "top_center"
    p2.yaxis.axis_label = '% time spent'
    p2.xaxis.major_label_orientation = 0.3

    # Create the pie charts
    pieColors = 10 * ["red", "green", "blue", "orange", "yellow", 'lime', 'brown', 'cyan',
                      'magenta', 'olive', 'black', 'teal', 'gold', 'crimson', 'moccasin', 'greenyellow', 'navy',
                      'ivory', 'lightpink']

    # First one with the time spent

    # define starts/ends for wedges from percentages of a circle
    percents_only = np.round(np.array(list(data["spent_hours"] / np.sum(data["spent_hours"]))) * 100, 1)
    percents = np.cumsum([0] + list(data["spent_hours"] / np.sum(data["spent_hours"])))
    starts = [per * 2 * np.pi for per in percents[:-1]]
    ends = [per * 2 * np.pi for per in percents[1:]]

    p3 = figure(x_range=(-1, 2.5), y_range=(-1.1, 1.1), plot_height=420, plot_width=600, title="% spent")

    # Add individual wedges:
    for i in range(N):
        p3.wedge(x=0, y=0, radius=.9, start_angle=starts[i], end_angle=ends[i], color=pieColors[i],
                 legend="[{0}%] {1}".format(percents_only[i], alloc_names[i]))

    p3.xgrid.grid_line_color = None
    p3.ygrid.grid_line_color = None
    p3.legend.orientation = "vertical"
    p3.legend.location = "top_right"
    p3.legend.border_line_alpha = 0
    p3.legend.background_fill_color = None
    p3.xaxis.visible = False
    p3.yaxis.visible = False

    # Second one with the time allocated

    # define starts/ends for wedges from percentages of a circle
    percents_only = np.round(np.array(list(data["alloc_hours"] / np.sum(data["alloc_hours"]))) * 100, 1)
    percents = np.cumsum([0] + list(data["alloc_hours"] / np.sum(data["alloc_hours"])))
    starts = [per * 2 * np.pi for per in percents[:-1]]
    ends = [per * 2 * np.pi for per in percents[1:]]

    p4 = figure(x_range=(-1, 2.5), y_range=(-1.1, 1.1), plot_height=420, plot_width=600,
                title="% time allocated to each program")
    # Add individual wedges:
    for i in range(N):
        p4.wedge(x=0, y=0, radius=.9, start_angle=starts[i], end_angle=ends[i], color=pieColors[i],
                 legend="[{0}%] {1}".format(percents_only[i], alloc_names[i]))

    p4.xgrid.grid_line_color = None
    p4.ygrid.grid_line_color = None
    p4.legend.orientation = "vertical"
    p4.legend.location = "top_right"
    p4.legend.border_line_alpha = 0
    p4.legend.background_fill_color = None
    p4.xaxis.visible = False
    p4.yaxis.visible = False

    layout = row(column(p, p2), column(p4, p3))

    curdoc().add_root(layout)
    curdoc().title = "Allocation stats"

    return layout


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE VIEW_DATA PAGE.              #
# KEYWORD:VIEW_DATA                                                           #
###############################################################################
def get_science_products(user_id="", obsdate="", camera_type=""):
    """

    :param user_id: 
    :param obsdate: 
    :param camera_type: 
    :return: 
    """
    if not obsdate:
        obsdate = datetime.datetime.utcnow().strftime("%Y%m%d")
    else:
        obsdate = obsdate[0]

    if not camera_type:
        camera_type = 'ifu'

    data_dir = '%s%s/' % (redux_dir, obsdate)
    if camera_type == 'ifu':
        return get_ifu_products(data_dir, user_id, obsdate)
    else:
        return {'message': "Need to add RC case"}


def get_ifu_products(obsdir, user_id, obsdate="", show_finder=True,
                     product_type='all'):
    """

    :param obsdir: 
    :param user_id: 
    :param obsdate: 
    :param product_type: 
    :return: 
    """
    ifu_dict = {}

    # Look first to make sure there is a data directory.
    if not os.path.exists(obsdir):
        return {'message': 'No data directory could be located for %s UT' %
                           os.path.basename(os.path.normpath(obsdir)),
                'obsdate': obsdate}

    if not obsdate:
        obsdate = os.path.basename(os.path.normpath(obsdir))
    sedm_dict = {'obsdate': obsdate,
                 'sci_data': ''}

    # Now lets get the non-science products (i.e. calibrations)
    calib_dict = {'flat3d': os.path.join(obsdir, '%s_flat3d.png' % obsdate),
                  'wavesolution': os.path.join(obsdir,
                                               '%s_wavesolution'
                                               '_dispersionmap.png' % obsdate),
                  'cube_lambdarms': os.path.join(obsdir, 'cube_lambdarms.png'),
                  'cube_trace_sigma': os.path.join(obsdir,
                                                   'cube_trace_sigma.png')}
    # If a calibration frame doesn't exist then pop it out to avoid bad links
    # on the page
    remove_list = []
    div_str = ''
    for k, v in calib_dict.items():
        if not os.path.exists(v):
            remove_list.append(k)

    if remove_list:
        for i in remove_list:
            calib_dict.pop(i)
    print(calib_dict, 'calib products')

    div_str += """<div class="row">"""
    div_str += """<h4>Calibrations</h4>"""
    for k, v in calib_dict.items():
        impath = "/data/%s/%s" % (obsdate, os.path.basename(v))
        impathlink = "/data/%s/%s" % (obsdate,
                                      os.path.basename(v.replace('.png', '.pdf')))
        if not os.path.exists(impathlink):
            impathlink = impath
        div_str += """<div class="col-md-{0}">
          <div class="thumbnail">
            <a href="{1}">
              <img src="{2}" width="{3}px" height="{4}px">
            </a>
          </div>
        </div>""".format(2, impathlink, impath, 400, 400)
    div_str += "</div>"
    sedm_dict['sci_data'] += div_str
    # To get ifu products we first look to see if a what.list file has been
    # created. This way we will know which files to add to our dict and
    # whether the user has permissions to see the file
    if not os.path.exists(os.path.join(obsdir, 'what.list')):
        return {'message': 'Could not find summary file (what.list) for %s UT' %
                           os.path.basename(os.path.normpath(obsdir))}

    # Go throught the what list and return all non-calibration entries
    with open(os.path.join(obsdir, 'what.list')) as f:
        what_list = f.read().splitlines()

    science_list = []
    standard_list = []
    for targ in what_list:
        if 'Calib' in targ:
            pass
        elif '[A]' in targ or '[B]' in targ or 'STD' in targ:
            science_list.append(targ)
        elif 'STD' in targ:
            pass
            #standard_list.append(targ)
        else:
            # There shouldn't be anything here but should put something in
            # later to verify this is the case
            pass

    # Now we go through and make sure the user is allowed to see this target
    show_list = []
    if len(science_list) >= 1:
        allocation_id_list = get_allocations_user(user_id=user_id,
                                                  return_type='list')

        for sci_targ in science_list:

            # Start by pulling up all request that match the science target
            targ_name = sci_targ.split(':')[1].split()[0]
            if 'STD' not in targ_name:
                # 1. Get the object id
                object_ids = db.get_object_id_from_name(targ_name)

                if len(object_ids) == 1:
                    object_id = object_ids[0][0]
                else:
                    # TODO       what really needs to happen here is that we need to
                    # TODO cont: find the id that is closest to the obsdate.
                    # TODO cont: For now I am just going to use last added
                    object_id = object_ids[-1][0]

                target_requests = db.get_from_request(values=['allocation_id'],
                                                      where_dict={'object_id':
                                                                      object_id,
                                                                  'status':
                                                                      'COMPLETED'})

                # Right now I am only seeing if there exists a match between
                # allocations of all request.  It's possible the request could
                # have been made by another group as another follow-up and thus
                # the user shouldn't be able to see it.  This should be able to
                # be fixed once all request are listed in the headers of the
                # science images.
                for req in target_requests:
                    if req[0] in allocation_id_list:
                        show_list.append((sci_targ, targ_name))
                    else:
                        print("You can't see this")
            else:
                targ_name = sci_targ.split(':')[1].split()[0].replace('STD-', '')
                show_list.append((sci_targ, targ_name))

    if len(standard_list) >= 1:
        for std_targ in standard_list:
            targ_name = std_targ.split(':')[1].split()[0].replace('STD-', '')
            show_list.append((std_targ, targ_name))

    # We have our list of targets that we can be shown, now lets actually find
    # the files that we will show on the web page.  To make this backwards
    # compatible I have to look for two types of files
    if len(show_list) >= 1:
        science_dict = {}
        count = 0
        div_str = ''
        for targ in show_list:
            targ_params = targ[0].split()
            fits_file = targ_params[0].replace('.fits', '')
            name = targ[1]

            image_list = (glob.glob('%sifu_spaxels_*%s*.png' % (obsdir,
                                                                fits_file)) +
                          glob.glob('%simage_%s*.png' % (obsdir, name)))

            spec_list = (glob.glob('%s%s_SEDM.png' % (obsdir, name)) +
                         glob.glob('%sspec_forcepsf*%s*.png' % (obsdir,
                                                                fits_file)))

            e3d_list = (glob.glob('%se3d*%s*.fits' % (obsdir, fits_file)))
            if name not in science_dict:
                science_dict[name] = {'image_list': image_list,
                                      'spec_list': spec_list,
                                      'e3d_list': e3d_list}
            else:
                # We do this to handle cases where there are two or more of
                # the same object name
                science_dict[name+'_xRx_%s' % str(count)] = {'image_list': image_list,
                                                             'spec_list': spec_list,
                                                             'e3d_list': e3d_list}
            count += 1
        # Alright now we build the table that will show the spectra, image file
        # and classification.

        count = 0

        for obj, obj_data in science_dict.items():
            if '_xRx_' in obj:
                obj = obj.split('_xRx_')[0]

            if 'ZTF' in obj:
                obj_link = ('<a href="http://skipper.caltech.edu:8080/'
                            'cgi-bin/growth/view_source.cgi?name=%s">%s</a>' %
                            (obj, obj))

                div_str += """<div class="row">"""
                div_str += """<h4>%s</h4>""" % obj_link
            else:
                div_str += """<div class="row">"""
                div_str += """<h4>%s</h4>""" % obj

            if obj_data['e3d_list']:
                for j in obj_data['e3d_list']:
                    impath = "/data/%s/%s" % (obsdate, os.path.basename(j))
                    div_str += ('<div class="col-md-{2}">'
                                '<a href="%s">E3D File</a>'
                                '</div>' % impath)

            # ToDO: Grab data from somewhere to put in the meta data column
            if obj_data['image_list']:
                for i in obj_data['image_list']:

                    impath = "/data/%s/%s" % (obsdate, os.path.basename(i))
                    impathlink = "/data/%s/%s" % (obsdate,
                                                  os.path.basename(i.replace('.png', '.pdf')))
                    if not os.path.exists(impathlink):
                        impathlink = impath

                    div_str += """<div class="col-md-{0}">
                    <div class="thumbnail">
                      <a href="{1}">
                        <img src="{2}" width="{3}px" height="{4}px">
                      </a>
                    </div>
                  </div>""".format(2, impathlink, impath, 400, 400)
            if show_finder:
                finder_path = os.path.join(phot_dir, obsdate, 'finders')
                if os.path.exists(finder_path):
                    finder_img = glob.glob(finder_path + '/*%s*.png' % obj)
                    if finder_img:
                        impathlink = "/data/%s/%s" % (obsdate, os.path.basename(finder_img[-1]))
                        div_str += """<div class="col-md-{0}">
                                              <div class="thumbnail">
                                                <a href="{1}">
                                                  <img src="{2}" width="{3}px" height="{4}px">
                                                </a>
                                              </div>
                                            </div>""".format(4, impathlink, impathlink, 250, 250)
            if obj_data['spec_list']:
                for i in obj_data['spec_list']:
                    impath = "/data/%s/%s" % (obsdate, os.path.basename(i))
                    impathlink = "/data/%s/%s" % (obsdate,
                                                  os.path.basename(i.replace('.png', '.pdf')))
                    if not os.path.exists(impathlink):
                        impathlink = impath
                    div_str += """<div class="col-lg-{0}">
                      <div class="thumbnail">
                        <a href="{1}">
                          <img src="{2}" width="{3}px" height="{4}px">
                        </a>
                      </div>
                    </div>""".format(4, impathlink, impath, 400, 400)

            div_str += "</div>"

        sedm_dict['sci_data'] += div_str

    return sedm_dict

###############################################################################
# THIS SECTION HANDLES THE ACTIVE_VISIBILITIES PAGE.                          #
# KEYWORD:VISIBILITIES #???                                                   #
###############################################################################
def get_active_visibility(userid):
    sedm_dict = {'enddate': datetime.datetime.utcnow() + datetime.timedelta(days=1),
                 'inidate': datetime.datetime.utcnow() - datetime.timedelta(days=7, hours=8)}

    # 1. Get a dataframe of all requests for the current user
    requests = get_requests_for_user(userid, sedm_dict['inidate'], sedm_dict['enddate'])

    # organize requests into dataframes by whether they are completed or not
    active = requests[(requests['status'] == 'PENDING') | (requests['status'] == 'ACTIVE')]

    # retrieve information about the user's allocations
    ac = get_allocations_user(userid)

    # Create html tables
    sedm_dict['active'] = {'table': active.to_html(escape=False,
                                                   classes='table',
                                                   index=False),
                           'title': 'Active Requests for the last 7 days'}

    sedm_dict['script'], sedm_dict['div'] = plot_visibility(userid, sedm_dict)
    return sedm_dict

def plot_visibility(userid, sedm_dict, obsdate=None):
    '''
     plots visibilities for active requests at the current date. Will be adapted to plot previous observations and arbitrary objects
    userid: user whose allocations will be shown in color with details. Others will be greyed out
    sedm_dict: should have ['active']['table'] and ['enddate'] and ['inidate']
    obsdate: <str> YYYYMMDD. if "None", will use current date

    returns: components of a bokeh figure with the appropriate plot
    '''

    allocpalette = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a', '#b15928', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']
    requests = get_requests_for_user(2, sedm_dict['inidate'], sedm_dict['enddate']) #admin
    active = requests[(requests['status'] == 'PENDING') | (requests['status'] == 'ACTIVE')]
    # ['allocation', 'object', 'RA', 'DEC', 'start date', 'end date', 'priority', 'status', 'lastmodified', 'obs_seq', 'exptime', 'UPDATE']
    
    allowed_allocs = get_allocations_user(userid)
    active['allocation'].mask(~np.in1d(active['allocation'], allowed_allocs['allocation']), other='other', inplace=True)

    programs = {i['allocation']:i['program'] for _, i in allowed_allocs.iterrows()}
    programs['other'] = ['other']
    active.sort_values('allocation') # this needs to be alphabetical for the legend to look correct
    
    p = figure(plot_width=700, plot_height=500, toolbar_location='above',
               y_range=(0, 90), y_axis_location="right")
    
    ### setup with axes, sun/moon, frames, background 
    # TODO Dima says to never ever use SkyCoord in production code
    palomar_mountain = EarthLocation(lon=243.1361*u.deg, lat=33.3558*u.deg, height=1712*u.m)
    utcoffset = -7 * u.hour  # Pacific Daylight Time

    
    if obsdate is None: # plotting a single object, or the pending objects in future
        time = (Time.now() - utcoffset).datetime # date is based on local time
        time = Time(datetime.datetime(time.year, time.month, time.day))
    else: # past observations on a particular night
        time = Time(datetime.datetime(int(obsdate[:4]), int(obsdate[4:6]), int(obsdate[6:8])))
        all_requests = all_requests[all_requests['status'] == 'COMPLETED']
        all_requests = all_requests[time - 12 * u.hour <= all_requests['lastmodified'] < time + 12 * u.hour]
    midnight = time - utcoffset # 7am local time of correct date, midnight UTC

    delta_midnight = np.linspace(-8, 8, 500) * u.hour
    t = midnight + delta_midnight
    abstimes = np.asarray([i.datetime.strftime('%I:%M %p') for i in t + utcoffset])

    frame = AltAz(obstime=t, location=palomar_mountain)
    sun_alt  =  get_sun(t).transform_to(frame).alt
    moon_alt = get_moon(t).transform_to(frame).alt
    
    # shading for nighttime and twilight
    dark_times    = delta_midnight[sun_alt < 0].value
    twilit_times  = delta_midnight[sun_alt < -18 * u.deg].value
    plotted_times = delta_midnight[sun_alt <   5 * u.deg].value
    
    twilight = BoxAnnotation(left=min(twilit_times), right=max(twilit_times), bottom=0, 
                             fill_alpha=0.15, fill_color='black', level='underlay')
    night    = BoxAnnotation(left=min(dark_times),    right=max(dark_times),    bottom=0, 
                             fill_alpha=0.25, fill_color='black', level='underlay')
    earth    = BoxAnnotation(top=0, fill_alpha=0.8, fill_color='sienna')
    
    p.add_layout(night)
    p.add_layout(twilight)
    p.add_layout(earth)
    
    # sun and moon
    sun  = p.line(delta_midnight, sun_alt,  line_color='red', name="Sun", legend='Sun', line_dash='dashed')
    moon = p.line(delta_midnight, moon_alt, line_color='yellow', line_dash='dashed', 
                                                   name="Moon", legend='Moon')
    # labels and axes
    p.title.text = "Visibility for %s UTC" %midnight
    p.xaxis.axis_label = "Hours from PDT Midnight"
    p.x_range.start = min(plotted_times)
    p.x_range.end   = max(plotted_times)
    p.yaxis.axis_label = "Airmass"
    
    # primary airmass label on right
    airmasses = (1.01, 1.1, 1.25, 1.5, 2., 3., 6.)
    ticker = [90 - np.arccos(1./i) * 180/np.pi for i in airmasses]
    p.yaxis.ticker = ticker
    p.yaxis.major_label_overrides = {tick: str(airmasses[i]) for i, tick in enumerate(ticker)}
    
    # add supplementary alt label on left
    p.extra_y_ranges = {"altitude": Range1d(0, 90)}
    p.add_layout(LinearAxis(y_range_name="altitude", axis_label='Altitude [deg]'), 'left')

    ##########################################################################
    ### adding data from the actual objects
    #objs = SkyCoord(np.array(ras,  dtype=np.float), 
    #                np.array(decs, dtype=np.float), unit="deg")
    
    approx_midnight = int(Time.now().jd - .5) + .5 - utcoffset.value/24.
    palo_sin_lat = 0.549836545
    palo_cos_lat = 0.835272275
    palo_long = 243.1362

    alloc_color = {}
    for i, val in allowed_allocs.iterrows():
        alloc_color[val['allocation']] = allocpalette[i % len(allocpalette)]
    alloc_color['other'] = 'lightgray'
        
    tooltipped = [] # things with tooltips
    tooltips = [('obj',        '@name'), # make it #name when we get to bokeh 0.13
                ('time',       '@abstime'), 
                ('altitude',   u"@alt\N{DEGREE SIGN}"), 
                ('airmass',    '@airmass')]

    for _, req in active.iterrows():
        req['ra'] = float(req['RA'])
        req['dec'] = float(req['DEC']) # iterrows doesn't preserve datatypes and turns ra, dec into decimals?
        color = alloc_color[req['allocation']]
        # vvv I got this formula from some website for the navy but forgot to copy the url
        alt = 180 / np.pi * np.arcsin(palo_cos_lat * \
              np.cos(np.pi/180 * (palo_long - req['ra'] + 15 * (18.697374558 + 24.06570982 * (delta_midnight.value/24. + approx_midnight - 2451545)))) * \
              np.cos(req['dec'] * np.pi/180) + palo_sin_lat * np.sin(req['dec'] * np.pi/180))
        airmass = 1./np.cos((90 - alt) * np.pi/180)
        source = ColumnDataSource(    dict(times=delta_midnight, 
                                             alt=alt,
                                         airmass=airmass,
                                         abstime=abstimes,
                                        priority=np.full(len(t), int(req['priority'])),
                                           alloc=np.full(len(t), req['allocation'][6:]),
                                            name=np.full(len(abstimes), req['object']))) # delete the name when we get to bokeh 0.13
        if len(active) == 1: # single object
            legend = req['object']
            line_width = 5
        else:
            legend = '{}'.format(programs[req['allocation']])
            #tooltips += [('priority',   '@priority'), ('allocation', '@alloc')]
            
            if req['status'] == 'COMPLETED': # plot that highlights observed part of the night
                # full path of the night
                dotted = p.line('times', 'alt', color=color, source=source, line_dash='2 2',
                                name=req['object'], line_width=1, legend=legend)
                # manually crop the source so only thick observed part has tooltips
                endtime = req['lastmodified']
                exptime = {req['obs_seq'][i]:req['exptime'][i] for i in range(len(req['obs_seq']))}['1ifu'] #TODO sometimes it's 2ifu or no ifu
                initime = endtime - exptime * u.second

                mask = np.logical_and(delta_midnight + midnight + utcoffset > initime,
                                      delta_midnight + midnight + utcoffset < endtime)
                source = ColumnDataSource(pd.DataFrame(source.data)[mask])
                line_width = int(req['priority'] + 3) # all it changes is the line width      
            else:
                line_width = int(req['priority'])
        
        path = p.line('times', 'alt', color=color, source=source, name=''.format(req['object']),
                      line_width=line_width, legend=legend)
        if not req['allocation'] == 'other':
            tooltipped.append(path)
        
    p.legend.click_policy = 'hide'
    p.legend.location = 'bottom_right'
    p.add_tools(HoverTool(renderers=tooltipped, tooltips=tooltips))
    
    curdoc().add_root(p)
    curdoc().title = 'Visibility plot'
    
    return components(p)

###############################################################################
# THIS SECTION IS THE WEATHER STATS SECTION.                                  #
# KEYWORD:WEATHER_STATS                                                       #
###############################################################################
def get_weather_stats(obsdate=None):
    message = ""
    if not obsdate:
        # get the weather stats
        statsfile, mydate = search_stats_file()
        stats_plot = plot_stats(statsfile, mydate)
        if stats_plot is None:
            message += " No statistics log found up to 100 days prior to today... Weather has been terrible lately!"
            script, div = None, None
        else:
            message += " Weather statistics for last opened day: %s" % (
                os.path.basename(os.path.dirname(os.path.dirname(statsfile))))
            script, div = components(stats_plot)
    else:
        mydate_in = obsdate.replace("-", "")

        # Just making sure that we have only allowed digits in the date
        mydate = re.findall(r"(2\d{3}[0-1]\d{1}[0-3]\d{1})", mydate_in)
        if len(mydate) == 0:
            message += "Incorrect format for the date! Your input is: %s. Shall be YYYYMMDD. \n" % mydate_in
            script, div = "", ""
        else:
            mydate = mydate[0]
            message = ""

            statsfile, mydate_out = search_stats_file(mydate)
            stats_plot = plot_stats(statsfile, mydate)
            if not statsfile:
                message = message + "No statistics log found for the date %s. Showing P18 data." % mydate
                script, div = components(stats_plot)

            else:
                stats_plot = plot_stats(statsfile, mydate)
                message = message + "Weather statistics for selected day: %s" % mydate
                script, div = components(stats_plot)

    return {'script': script, 'div': div, 'message': message}


def search_stats_file(mydate=None):
    '''
    Returns the last stats file that is present in the system according to the present date.
    It also returns a message stating what date that was.
    '''
    # If the date is specified, we will try to located the right file.
    # None will be returned if it does not exist.
    if mydate:
        s = os.path.join(phot_dir, mydate, "stats/stats.log")
        if os.path.isfile(s) and os.path.getsize(s) > 0:
            return s, mydate
        else:
            return None, None

    else:
        curdate = datetime.datetime.utcnow()
        # Try to find the stat files up to 100 days before today's date.
        i = 0
        while i < 100:
            newdate = curdate
            newdatedir = "%d%02d%02d" % (newdate.year, newdate.month, newdate.day)
            s = os.path.join(phot_dir, newdatedir, "stats/stats.log")

            if os.path.isfile(s) and os.path.getsize(s) > 0:
                return s, newdatedir
            i = i + 1
            curdate -= datetime.timedelta(days=1)
        return None, None


def load_p48seeing(obsdate):

    time, seeing = get_p18obsdata(obsdate)
    day_frac_diff = datetime.timedelta(
        np.ceil((datetime.datetime.now() - datetime.datetime.utcnow()).total_seconds()) / 3600 / 24)
    local_date = np.array(time) + day_frac_diff
    d = pd.DataFrame({'date': local_date, 'seeing': seeing})

    return d


def load_stats(statsfile='stats.log'):
    data = pd.read_csv(statsfile, header=None,
                       names=['path', 'obj', 'jd', 'ns', 'fwhm', 'ellipticity', 'bkg', 'airmass', 'in_temp', 'imtype',
                              'out_temp', 'in_hum'])

    jds = data['jd']
    t = Time(jds, format='jd', scale='utc')
    date = t.utc.datetime
    day_frac_diff = datetime.timedelta(
        np.ceil((datetime.datetime.now() - datetime.datetime.utcnow()).total_seconds()) / 3600 / 24)
    local_date = date + day_frac_diff

    data2 = data.assign(localdate=local_date)
    data2.set_index('localdate')

    return pd.DataFrame(
        {'date': data2['localdate'], 'ns': data2['ns'], 'fwhm': data2['fwhm'], 'ellipticity': data2['ellipticity'],
         'bkg': data2['bkg'], 'airmass': data2['airmass'], 'in_temp': data2['in_temp'], 'imtype': data2['imtype'],
         'out_temp': data2['out_temp'], 'in_hum': data2['in_hum']})


def plot_stats(statsfile, mydate):
    source = ColumnDataSource(
        data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[], in_temp=[], imtype=[], out_temp=[],
                  in_hum=[]))
    source_static = ColumnDataSource(
        data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[], in_temp=[], imtype=[], out_temp=[],
                  in_hum=[]))

    viewScience = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='SCIENCE')])
    viewAcquisition = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='ACQUISITION')])
    viewGuider = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='GUIDER')])
    viewFocus = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='FOCUS')])
    source_p48 = ColumnDataSource(data=dict(date=[], seeing=[]))

    def update(selected=None):

        if statsfile:
            data = load_stats(statsfile)
            source.data = source.from_df(data[['date', 'ns', 'fwhm', 'ellipticity', 'bkg', 'airmass', 'in_temp',
                                               'imtype', 'out_temp', 'in_hum']])
            source_static.data = source.data

        p48 = load_p48seeing(mydate)
        source_p48.data = source_p48.from_df(p48[['date', 'seeing']])
        source_static_p48.data = source_p48.data

    source_static_p48 = ColumnDataSource(data=dict(date=[], seeing=[]))
    tools = 'pan,box_zoom,reset'

    p48seeing = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
    p48seeing.circle('date', 'seeing', source=source_static_p48, color="black")
    p48seeing.title.text = "P18 seeing [arcsec]"

    if statsfile:
        ns = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        ns.line('date', 'ns', source=source_static)
        ns.circle('date', 'ns', size=1, source=source, color=None, selection_color="orange")
        ns.title.text = "Number of bright sources extracted"

        bkg = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        bkg.x_range = ns.x_range
        bkg.line('date', 'bkg', source=source_static)
        bkg.circle('date', 'bkg', size=1, source=source, color=None, selection_color="orange")
        bkg.title.text = "Background (counts)"

        temp = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        temp.x_range = ns.x_range
        temp.line('date', 'in_temp', source=source_static, color='blue', legend="Inside")
        temp.line('date', 'out_temp', source=source_static, color='green', legend="Outside")
        temp.circle('date', 'in_temp', size=1, source=source, color=None, selection_color="orange")
        temp.title.text = "Temperature [C]"

        fwhm = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        fwhm.x_range = ns.x_range
        fwhm.circle('date', 'fwhm', source=source_static, color="green", legend="Focus", view=viewFocus)
        fwhm.circle('date', 'fwhm', source=source_static, color="red", legend="Science", view=viewScience)
        fwhm.circle('date', 'fwhm', source=source_static, color="blue", legend="Acquisition", view=viewAcquisition)
        fwhm.circle('date', 'fwhm', source=source_static, color="black", legend="Guider", view=viewGuider)
        fwhm.circle('date', 'fwhm', size=1, source=source, color=None, selection_color="orange")
        fwhm.title.text = "P60 FWHM [arcsec]"

        airmass = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        airmass.x_range = ns.x_range
        airmass.line('date', 'airmass', source=source_static)
        airmass.circle('date', 'airmass', size=1, source=source, color=None, selection_color="orange")
        airmass.title.text = "Airmass"

        ellipticity = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime',
                             active_drag="box_zoom")
        ellipticity.x_range = ns.x_range
        ellipticity.line('date', 'ellipticity', source=source_static)
        ellipticity.circle('date', 'ellipticity', size=1, source=source, color=None, selection_color="orange")
        ellipticity.title.text = "Ellipticity"

        humidity = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        humidity.x_range = ns.x_range
        humidity.line('date', 'in_hum', source=source_static)
        humidity.circle('date', 'in_hum', size=1, source=source, color=None, selection_color="orange")
        humidity.title.text = "Inside Humidity [%]"

        p48seeing.x_range = ns.x_range

        left = column(fwhm, p48seeing, airmass)
        center = column(ellipticity, ns, bkg, )
        right = column(temp, humidity)
        layout = row(left, center, right)

    else:

        layout = row(column(p48seeing))

    # initialize
    update()

    curdoc().add_root(layout)
    curdoc().title = "Stats"

    return layout


def plot_not_found_message(day):
    not_found = figure(plot_width=900, plot_height=450, x_range=[0, 900], y_range=[0, 450])
    not_found.image(image=[np.zeros([900, 450]) + 0.1], x=0, y=0, dw=900, dh=450)
    citation = Label(x=50, y=225, x_units='screen', y_units='screen',
                     text='No statistics found for today \n (likely we were weathered out...)')
    not_found.add_layout(citation)
    not_found.title.text = "Statistics not found for day %s" % day

    layout = column(not_found)
    curdoc().add_root(layout)
    curdoc().title = "Stats not found"


###############################################################################
# THIS SECTION IS A COMMON UTILITIES SECTION                                  #
# KEYWORD:UTILITIES                                                           #
###############################################################################
def get_config_paths():
    return dict(path={
        'path_archive': redux_dir,
        'path_phot': phot_dir,
        'path_raw': raw_dir})


def make_dict_from_dbget(headers, data, decimal_to_float=True):
    """
    This function takes data from the returns of get_from_* returns and puts
    it in a dictionary form
    :param decimal_to_float: 
    :param headers: list of db header names
    :param data: tuples
    :return: 
    """

    if len(headers) != len(data):
        return {'error': 'headers and data are not of equal lengths'}

    return_dict = {}
    for i in range(len(headers)):
        if decimal_to_float and isinstance(data[i], Decimal):
            return_dict[headers[i]] = float(data[i])
        else:
            return_dict[headers[i]] = data[i]

    return return_dict


def get_filter_exptime(obsfilter, mag):
    """

    :param obsfilter: 
    :param mag: 
    :return: 
    """

    mag = float(mag)
    if mag > 18:
        ifu_exptime = 3600
        r_exptime = 180
        g_exptime = 180
        i_exptime = 180
        u_exptime = 300
    elif 15 > mag < 18:
        ifu_exptime = 2700
        r_exptime = 120
        g_exptime = 120
        i_exptime = 120
        u_exptime = 300
    elif 0 < mag < 5:
        ifu_exptime = 60
        r_exptime = 1
        g_exptime = 1
        i_exptime = 1
        u_exptime = 30
    elif 5 < mag < 10:
        ifu_exptime = 90
        r_exptime = 10
        g_exptime = 10
        i_exptime = 10
        u_exptime = 60
    elif 10 < mag < 12:
        ifu_exptime = 300
        r_exptime = 30
        g_exptime = 30
        i_exptime = 30
        u_exptime = 60
    elif 12 < mag < 13:
        ifu_exptime = 600
        r_exptime = 60
        g_exptime = 60
        i_exptime = 60
        u_exptime = 120
    elif 13 < mag < 15:
        ifu_exptime = 900
        r_exptime = 90
        g_exptime = 90
        i_exptime = 90
        u_exptime = 180

    else:
        ifu_exptime = 1800
        r_exptime = 90
        g_exptime = 90
        i_exptime = 90
        u_exptime = 90

    if obsfilter == 'ifu':
        return str(ifu_exptime)
    elif obsfilter == 'r':
        return str(r_exptime)
    elif obsfilter == 'g':
        return str(g_exptime)
    elif obsfilter == 'i':
        return str(i_exptime)
    elif obsfilter == 'u':
        return str(u_exptime)
    else:
        return str(0)


def get_p18obsdata(obsdate):
    """
    :param obsdate: Must be in "Year-Month-Day" or "YYYYMMDD" format
    :return: List of dates and average seeing
    """
    # 1. Create the URL to get the seeing for the requested night
    p18date = []
    p18seeing = []

    if "-" in obsdate:
        f = datetime.datetime.strptime(obsdate, "%Y-%m-%d") - datetime.timedelta(days=1)
    else:
        f = datetime.datetime.strptime(obsdate, "%Y%m%d") - datetime.timedelta(days=1)

    y, m, d = [f.strftime("%Y"), int(f.strftime("%m")), int(f.strftime("%d"))]
    p18obsdate = "%s-%s-%s" % (y, m, d)

    # 2. Get the data from the link
    page = requests.get('http://nera.palomar.caltech.edu/P18_seeing/seeing_log_%s.log' % p18obsdate)
    data = page.content.decode("ISO-8859-1")

    # 3. Split the page by newlines
    data = data.split('\n')

    # 4. Loop through the data and only use points that have 4 or more seeing values to average
    for i in data:
        i = i.split()

        if len(i) > 5 and int(i[5]) > 4:
            d = '%s %s' % (i[1], i[0])
            p18date.append(datetime.datetime.strptime(d, "%m/%d/%Y %H:%M:%S")
                           + datetime.timedelta(hours=8))
            p18seeing.append(float(i[4]))

    return p18date, p18seeing


if __name__ == "__main__":
    x = get_ifu_products('/scr7/rsw/sedm/redux/20180827/', 189)
    print(x)
