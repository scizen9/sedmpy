from werkzeug.security import check_password_hash
from db.SedmDb import SedmDB
import datetime
import os
import json
import re
import pandas as pd
import numpy as np
import requests
import glob
import time
from decimal import Decimal
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, Label, Span
from bokeh.models import CDSView, GroupFilter
from bokeh.models.tools import HoverTool
from bokeh.models.ranges import Range1d
from bokeh.models.axes import LinearAxis
from bokeh.models.annotations import BoxAnnotation
from bokeh.plotting import figure
from bokeh.embed import components
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun,\
    get_moon
from scheduler.scheduler import ScheduleNight

pd.options.mode.chained_assignment = None   # default='warn'

superuser_list = ['SEDm_admin', 2, 20180523190352189]

pd.set_option('display.max_colwidth', None)

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
                       'rc', 'do_r', 'do_g', 'do_i', 'do_u',
                       'r_exptime', 'g_exptime', 'i_exptime', 'u_exptime',
                       'r_repeats', 'g_repeats', 'i_repeats', 'u_repeats',
                       'ifu_exptime', 'priority', 'inidate', 'rc_use_mag',
                       'enddate', 'maxairmass', 'status', 'max_fwhm',
                       'min_moon_dist', 'max_moon_illum', 'max_cloud_cover',
                       'seq_repeats', 'seq_completed']

rc_filter_list = ['r', 'g', 'i', 'u']
schedule = ScheduleNight()

# this all needs to go in some sort of config file instead of changing the
# source code constantly
computer = os.uname()[1]  # a quick fix

port = 0
host = 'none'
if computer == 'pele':
    raw_dir = '/scr/rsw/sedm/raw/'
    phot_dir = '/scr/rsw/sedm/phot/'
    redux_dir = '/scr/rsw/sedm/data/redux/'
    new_phot_dir = '/scr/rsw/sedm/data/redux/phot/'
    status_dir = '/scr/rsw/'
    requests_dir = '/scr/rsw/'
    base_dir = '/scr/rsw/'
    host = 'minar.caltech.edu'
    port = 5432

elif computer == 'pharos':
    raw_dir = '/scr2/sedm/raw/'
    phot_dir = '/scr2/sedm/phot/'
    new_phot_dir = '/scr2/sedmdrp/redux/phot/'
    redux_dir = '/scr2/sedmdrp/redux/'
    status_dir = '/scr2/sedm/raw/telstatus/'
    requests_dir = '/scr2/sedm/logs/requests/'
    host = 'localhost'
    base_dir = '/scr2/sedmdrp/'
    port = 5432

elif computer == 'minar':
    raw_dir = '/data/sedmdrp/raw/'
    phot_dir = '/data/sedmdrp/redux/phot/'
    new_phot_dir = '/data/sedmdrp/redux/phot/'
    redux_dir = '/data/sedmdrp/redux/'
    status_dir = '/data/sedmdrp/raw/telstatus/'
    requests_dir = '/data/sedmdrp/logs/requests/'
    host = 'localhost'
    base_dir = '/data/sedmdrp/'
    port = 5432

elif computer == 'ether':
    raw_dir = '/home/rsw/sedm_data/raw/'
    phot_dir = '/home/rsw/sedm/phot/'
    redux_dir = '/home/rsw/sedm_data/redux/'
    new_phot_dir = '/home/rsw/sedm_data/redux/phot/'
    requests_dir = '/home/rsw/'
    base_dir = '/home/rsw/sedm_data/'
    host = 'localhost'
    port = 22222

print(computer, port, host, "inputs")
db = SedmDB(host=host, dbname='sedmdb', port=port)


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
        return {'message': df.to_html(escape=False, classes='table',
                                      index=False)}
    else:
        return obj_list


def fancy_request_table(df):
    """
    df: pandas dataframe
        intended for tables of requests from get_requests_for_user,
        ie with columns:
            ['allocation', 'object', 'RA', 'DEC', 'start date', 'end date',
             'priority','status', 'lastmodified', 'UPDATE']

    returns: IPython HTML object
        the html for df but with the following changes:
            -if 'RA' and 'Dec' mean it won't rise tonight, both fields are red
            -'name' column now has links to the growth marshal
            -table width is 100%,
             which is important for the fancy tables to display right
            -priority is a float to allow finer tuning of the scheduler
    """

    def highlight_set(hrow, color='#ff9999'):
        """
        makes 'RA' and 'DEC' fields highlighted if it won't get high when
        it's dark out meant for tables with both 'RA' and 'DEC' columns
        """
        red = 'background-color: {}'.format(color)
        try:
            # peak at longitude-based midnight, approx
            best_ra = ((Time.now().mjd - 58382.) * 360/365.)
            # more than 8h from midnight
            if 180 - abs(180 - (best_ra - float(hrow['RA'])) % 360) > 120:
                return [red if i == 'RA' else '' for i in hrow.index.values]
            # red if it'll never go above ~40deg
            if hrow['DEC'] < -15.:
                return [red if i == 'DEC' else '' for i in hrow.index.values]
            else:
                return ['' for _ in hrow.index.values]
        except KeyError:
            return ['' for _ in hrow.index.values]

    def improve_obs_seq(li):
        """
        takes a list like ['1ifu'], [180, 180, 180], etc and makes it take up
        less space and also be more human-readable
        """
        try:
            if type(li[0]) == str:
                # ie it's an obs_seq
                for i, val in enumerate(li):
                    if val[0] == '1':
                        li[i] = val[1:]
                    else:
                        li[i] = val[1:] + ' x' + val[0]
            else:  # ie exptime
                for i, val in enumerate(li):
                    if all([j == val for j in li[i:]]) and i < len(li) - 1:
                        # all the rest match, of which there's >1
                        li = li[:i] + ['{}ea'.format(val)]
                        break
                    else:
                        li[i] = str(val)
            return ', '.join(li)
        except:
            return "ERROR,PARSING_THE_FILTER_STRING"

    df['status'] = [i.lower() for i in df['status']]
    df['allocation'] = [i.replace('2018A-', '').replace('2018B-', '')
                        .replace('2019B-', '').replace('2021B-', '')
                        .replace('2022A-', '')
                        for i in df['allocation']]
    for col in ('obs_seq', 'exptime'):
        df[col] = [improve_obs_seq(i) for i in df[col]]

    if df.shape[0] > 0:
        styled = df.style\
               .apply(highlight_set, axis=1)\
               .format(
                {'object': '<a href="https://fritz.science/source/{0}">{0}</a>',
                 'RA': '{:.3f}', 'DEC': '{:.3f}', 'priority': '{:.1f}',
                 'start date': '{:%b %d}', 'end date': '{:%b %d}',
                 'lastmodified': '{:%b %d %H:%M}',
                 'UPDATE': '<a href="request?request_id={}">+</a>'})\
               .set_table_attributes('style="width:100%" '
                                     'class="dataframe_fancy table '
                                     'table-striped nowrap"')\
               .set_table_styles(
                    [{'selector': '.row_heading',
                      'props': [('display', 'none')]},
                     {'selector': '.blank.level0',
                      'props': [('display', 'none')]}])
    else:
        styled = df.style \
            .format(
             {'object': '<a href="https://fritz.science/source/{0}">{0}</a>',
              'RA': '{:.3f}', 'DEC': '{:.3f}', 'priority': '{:.1f}',
              'start date': '{:%b %d}', 'end date': '{:%b %d}',
              'lastmodified': '{:%b %d %H:%M}',
              'UPDATE': '<a href="request?request_id={}">+</a>'}) \
            .set_table_attributes('style="width:100%" '
                                  'class="dataframe_fancy table '
                                  'table-striped nowrap"') \
            .set_table_styles(
                 [{'selector': '.row_heading',
                   'props': [('display', 'none')]},
                  {'selector': '.blank.level0',
                   'props': [('display', 'none')]}])
    # .set_table_styles([{'text-align': 'left'}])\
    # this .replace() thing is super bad form but it's faster for now than
    # finding the right way
    return styled.to_html()\
                 .replace('RA</th>', '<a href="#" data-toggle="tooltip" '
                          'title="red if peaks >8h from midnight">RA</a></th>')\
                 .replace('DEC</th>', '<a href="#" data-toggle="tooltip" '
                          'title="red if peaks below 40deg">dec</a></th>')


def get_homepage(userid, username):
    sedm_dict = {'enddate':
                 datetime.datetime.utcnow() + datetime.timedelta(days=30),
                 'inidate':
                 datetime.datetime.utcnow() - datetime.timedelta(days=7,
                                                                 hours=8)}

    # 1. Get a dataframe of all requests for the current user
    reqs = get_requests_for_user(userid, sedm_dict['inidate'],
                                 sedm_dict['enddate'])

    # organize requests into dataframes by whether they are completed or not
    complete = reqs[(reqs['status'] == 'COMPLETED') |
                    (reqs['status'] == 'OBSERVED') |
                    (reqs['status'] == 'OBSERVED')]
    active = reqs[(reqs['status'] == 'ACTIVE')]
    pending = reqs[(reqs['status'] == 'PENDING')]
    expired = reqs[(reqs['status'] == 'EXPIRED')]
    failed = reqs[(reqs['status'] == 'FAILED')]

    # retrieve information about the user's allocations
    ac = get_allocations_user(userid)

    # Create html tables
    sedm_dict['active'] = {'table': fancy_request_table(active),
                           'title': 'Active Request'}

    sedm_dict['pending'] = {'table': fancy_request_table(pending),
                            'title': 'Pending Requests'}

    sedm_dict['complete'] = {'table': fancy_request_table(complete),
                             'title': 'Completed Requests in the last 7 days'}

    sedm_dict['expired'] = {'table': fancy_request_table(expired),
                            'title': 'Expired Requests in the last 7 days'}

    sedm_dict['failed'] = {'table': fancy_request_table(failed),
                           'title': 'Failed Exposures in the last 7 days'}

    sedm_dict['allocations'] = {
        'table': ac.to_html(escape=False, classes='table table-striped',
                            index=False, col_space=10),
        'title': 'Your Active Allocations'}

    sedm_dict['visibility'] = {'title': 'Visibilities for pending requests',
                               'url':   '/visibility'}

    # Make a greeting statement
    sedm_dict['greeting'] = 'Hello %s!' % username
    return sedm_dict


###############################################################################
# THIS SECTION HANDLES EXPIRED TARGETS.                                       #
# KEYWORD:EXPIRED                                                             #
###############################################################################
def show_expired(days=7):
    """

    :return:
    """
    print(days)


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE SCHEDULER PAGE.              #
# KEYWORD:SCHEDULER                                                           #
###############################################################################
def get_schedule():
    """
    :return:
    """

    with open('static/scheduler/scheduler.html', 'r') as myfile:
        data = myfile.read().replace('\n', '')

    return {'scheduler': data}


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
        choices = [(0, "You have no active allocations!")]
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

        ret_dict.update(parse_db_target_filters(ret_dict['obs_seq'],
                                                ret_dict['exptime']))
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
    # if 'last_obs_jd' in ret_dict:
    #    form.last_obs_jd.data = ret_dict['last_obs_jd']
    if 'do_ifu' in ret_dict:
        form.ifu.data = ret_dict['do_ifu']
    if 'ifu_exptime' in ret_dict:
        form.ifu_exptime.data = ret_dict['ifu_exptime']
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
    rc_filters = ['r', 'g', 'i', 'u']

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
        if flt in rc_filters:
            if 0 <= flt_exptime <= 1000:
                if 1 <= flt_repeat <= 100:
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

    return_dict = {'message': ''}

    # Add the object to the database
    if content['obj_ra'] and content['obj_dec']:
        ra = content['obj_ra']
        dec = content['obj_dec']
        if ":" in ra or ":" in dec:
            c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
            ra = c.ra.degree
            dec = c.dec.degree

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

    alloc = get_allocations_user(int(userid))

    if alloc is None or len(alloc) == 0:
        choices = [(0, "You have no active allocations!")]
    else:
        choices = [z for z in zip(alloc['id'], alloc['allocation'])]
        choices.insert(0, (0, '------'))

    form.allocation.choices = choices

    # 1. Let's start by making sure we have all the information needed for the
    #    object id

    if 'object_id' in content and content['object_id']:
        objid = content['object_id']
    else:
        message, objid = add_object_to_db(content)
        if not objid:
            return {**process_dict, **message}, form
        else:
            if 'message' in message:
                process_dict['message'] += message['message']

    request_dict['object_id'] = int(objid)

    # 2. Now let's put together the request
    # by getting all the values into a dictionary
    obs_seq_key_list = ['ifu', 'rc', 'do_r', 'do_i', 'do_u', 'do_g',
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
            # of the database format we must handle this data outside
            # the request dictionary.
            elif key in obs_seq_key_list:
                if key in content:
                    obs_seq_dict[key] = content[key]
                else:
                    obs_seq_dict[key] = False
            else:
                request_dict[key] = content[key]
        except Exception as e:
            print(str(e), key, 't')
            pass
    # print("I made it here")
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

        if 'request_id' in content and content['request_id']:

            request_dict['id'] = int(content['request_id'])
            request_dict.pop('request_id')
            for k, v in request_dict.items():
                if not v:
                    request_dict[k] = -1
            # ret = db.update_request(request_dict)
            db.update_request(request_dict)
        else:
            # print("I AM HERE NOW")
            if 'request_id' in request_dict:
                request_dict.pop('request_id')
            request_dict['user_id'] = int(request_dict['user_id'])
            # print(request_dict)
            if 'external_id' in content:
                request_dict['external_id'] = content['external_id']
            # ret = db.add_request(request_dict)
            db.add_request(request_dict)
            # print(ret)
    return process_dict, form


def get_add_csv(user_id, form, content):
    """

    :param user_id: 
    :param form: 
    :param content: 
    :return: 
    """

    return {'test': 'test'}, form, user_id, content


def process_add_csv(content, form, user_id):
    """

    :param content: 
    :param form: 
    :param user_id: 
    :return: 
    """

    return {'test': 'test'}, form, content, user_id


def make_obs_seq(obs_seq_dict):
    """

    :param obs_seq_dict: 
    :return: 
    """
    filters_list = []
    exptime_list = []
    ret_dict = {"proc_message": ""}

    if isinstance(obs_seq_dict['ifu'], bool):
        if obs_seq_dict['ifu']:
            obs_seq_dict['ifu'] = 'y'
        else:
            obs_seq_dict['ifu'] = 'n'

    if isinstance(obs_seq_dict['rc'], bool):
        if obs_seq_dict['rc']:
            obs_seq_dict['rc'] = 'y'
        else:
            obs_seq_dict['rc'] = 'n'

    if obs_seq_dict['ifu'].lower() in ['y', 'yes', 'true']:
        # There may be case in the future where people want more than one IFU
        # at a time.  In which case this code will need to be changed.
        if obs_seq_dict['ifu_use_mag']:
            if obs_seq_dict['ifu_exptime'] and \
                    int(obs_seq_dict['ifu_exptime']) > 0:
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
                                             "you. Otherwise there is something"
                                             " wrong with this '%s' value.  For"
                                             " the record here is the error "
                                             "message %s--" %
                                             (obs_seq_dict['obj_mag'], str(e)))

                ifu_exptime = False
        else:
            try:
                ifu_exptime = int(obs_seq_dict['ifu_exptime'])
                if 0 <= ifu_exptime <= 7200:
                    pass
                else:

                    ret_dict['proc_message'] += ("I don't know what you are "
                                                 "trying to do but %s is not an"
                                                 " acceptable IFU exposure "
                                                 "time. It's either less than "
                                                 " 0 or more than two hours.--"
                                                 % str(ifu_exptime))
                    ifu_exptime = False
            except Exception as e:
                ret_dict['proc_message'] += ("There is something wrong with "
                                             "your exposure time value.  '%s' "
                                             "is not a proper value.  Here is "
                                             "the error message return: %s--" %
                                             (obs_seq_dict['ifu_exptime'],
                                              str(e)))

                ifu_exptime = False

        if ifu_exptime:
            filters_list.append("1ifu")
            exptime_list.append(str(ifu_exptime))

    # print(obs_seq_dict)
    if obs_seq_dict['rc'].lower() in ['y', 'yes', 'true']:
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
                    if 0 <= exptime <= 1000:
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
        ret_dict["ERROR"] = "NO FILTERS COULD BE DETERMINED"
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
    res = db.execute_sql(""" SELECT a.id, a.designator, p.designator,
                        g.designator, a.time_allocated, a.time_spent
                        FROM allocation a, program p, groups g, usergroups ug
                        WHERE a.program_id = p.id AND p.group_id = g.id 
                        AND g.id = ug.group_id AND a.active is True AND
                        ug.user_id = %d""" % user_id)

    # create the dataframe and set the allocation names to be linked
    if return_type == 'list':
        data = []
        for i in res:
            data.append(i[0])
    else:
        data = pd.DataFrame(res, columns=['id', 'allocation', 'program',
                                          'group', 'time allocated',
                                          'time spent'])

    return data


def get_requests_for_user(user_id, inidate=None, enddate=None):
    """

    :param user_id: 
    :param inidate: 
    :param enddate: 
    :return: 
    """
    if not inidate:
        inidate = datetime.datetime.utcnow() - datetime.timedelta(days=7,
                                                                  hours=8)
    if not enddate:
        enddate = datetime.datetime.utcnow() + datetime.timedelta(days=1)

    request_query = ("""SELECT a.designator, o.name, o.ra, o.dec, r.inidate,
     r.enddate, r.priority, r.status, r.lastmodified, r.obs_seq, r.exptime, r.id 
    FROM request r, object o, allocation a 
    WHERE o.id = r.object_id AND a.id = r.allocation_id  
        AND ( r.enddate > DATE('%s') AND r.inidate <= DATE('%s') )
        AND r.allocation_id IN
        (SELECT a.id
            FROM allocation a, groups g, usergroups ug, users u, program p
            WHERE ug.user_id = u.id AND ug.group_id = g.id AND u.id = %d AND 
            p.group_id = g.id AND a.program_id = p.id
        ) ORDER BY r.lastmodified DESC;""" % (inidate, enddate, user_id))

    data = db.execute_sql(request_query)
    data = pd.DataFrame(data,
                        columns=['allocation', 'object', 'RA', 'DEC',
                                 'start date', 'end date', 'priority', 'status',
                                 'lastmodified', 'obs_seq', 'exptime',
                                 'UPDATE'])

    if user_id in superuser_list:
        pass
        # data['UPDATE'] = data['UPDATE'].apply(convert_to_link)
    else:
        data.drop(columns=['RA', 'DEC'])

    return data


def convert_to_link(reqid):
    return """http://minar.caltech.edu/request?request_id=%s""" % reqid


###############################################################################
# THIS SECTION HANDLES ALL THINGS RELATED TO THE OBJECT PAGE.                 #
# KEYWORD:OBJECT                                                              #
###############################################################################
def get_object(object_name, user_id):
    """
    
    :param object_name: 
    :param user_id:
    :return: 
    """

    if user_id:
        pass
    # 1. Start by getting the requested object
    objects = get_object_info(object_name, out_type='html')

    # 2. Check if there were and objects if not then go on
    if not objects:
        return {'message': 'Could not find any targets with that name under '
                           'your allocation'}
    else:
        return {'message': objects['message']}


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

    user_pass = db.get_from_users(['username', 'password', 'id'],
                                  {'username': username})

    # print(user_pass)
    if not user_pass:
        return False, 'Incorrect username or password!'

    if check_password_hash(user_pass[0][1], password=password):
        return True, user_pass[0][2]
    else:
        return False, 'Incorrect username or password!!'


def password_change(form, user_id):
    """
    :param form: 
    :param user_id: 
    :return: 
    """
    # check for correct password and change if true
    password = form.password.data
    new_password = form.pass_new.data
    new_password_conf = form.pass_conf.data
    user_pass = db.get_from_users(['username', 'password', 'id'],
                                  {'id': user_id})

    if not user_pass:
        return {'message': "User not found"}

    elif user_pass[0] == -1:
        message = user_pass[1]
        return {'message': message}

    elif check_password_hash(user_pass[0][1], password):
        if new_password == new_password_conf:
            db.update_user({'id': user_pass[0][2],
                            'password': new_password})
            return {'message': 'Password Changed!'}
    else:
        message = "Incorrect username or password!"
        return {'message': message}


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
    if user_id is None:
        res = db.get_from_allocation(["designator", "time_allocated",
                                      "time_spent"], {"active": True})
        df = pd.DataFrame(res, columns=["designator", "time_allocated",
                                        "time_spent"])

        alloc_hours = np.array([ta.total_seconds() / 3600.
                                for ta in df["time_allocated"]])
        spent_hours = np.array([ts.total_seconds() / 3600.
                                for ts in df["time_spent"]])
        free_hours = alloc_hours - spent_hours

        free_hours[np.where(free_hours < 0)] = 0.

        df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours,
                       free_hours=free_hours)

    else:
        if inidate is None or enddate is None:
            res = db.execute_sql(""" SELECT a.designator, a.time_allocated,
            a.time_spent
            FROM allocation a, program p, groups g, usergroups ug
            WHERE a.program_id = p.id AND p.group_id = g.id 
            AND g.id = ug.group_id AND a.active is True 
            AND ug.user_id = %d""" % user_id)

            df = pd.DataFrame(res, columns=["designator", "time_allocated",
                                            "time_spent"])

            alloc_hours = np.array([ta.total_seconds() / 3600.
                                    for ta in df["time_allocated"]])
            spent_hours = np.array([ts.total_seconds() / 3600.
                                    for ts in df["time_spent"]])
            free_hours = alloc_hours - spent_hours

            free_hours[np.where(free_hours < 0)] = 0.

            df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours,
                           free_hours=free_hours)

        else:
            res = db.execute_sql(""" SELECT DISTINCT a.id, a.designator,
            a.time_allocated
            FROM allocation a, program p, groups g, usergroups ug
            WHERE a.program_id = p.id AND p.group_id = g.id 
            AND g.id = ug.group_id AND a.active is True 
            AND ug.user_id = %d;""" % user_id)
            allocdes = []
            spent_hours = []
            alloc = []
            for ais in res:
                spent = db.get_allocation_spent_time(ais[0], inidate, enddate)
                allocdes.append(ais[1])
                spent_hours.append(int(spent) / 3600.)
                alloc.append(ais[2])
            res = np.array([allocdes, alloc, spent_hours])

            df = pd.DataFrame(res.T, columns=["designator", "time_allocated",
                                              "time_spent"])
            alloc_hours = np.array([ta.total_seconds() / 3600.
                                    for ta in df["time_allocated"]])
            free_hours = alloc_hours - spent_hours

            free_hours[np.where(free_hours < 0)] = 0.
            df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours,
                           free_hours=free_hours)

    df = df.sort_values(by=["alloc_hours"], ascending=False)

    alloc_names = df["designator"].values
    category = ["alloc_hours", "spent_hours", "free_hours"]

    data = {'allocations': alloc_names}

    for cat in category:
        data[cat] = df[cat]

    return data


def plot_stats_allocation(data):
    """
    Plots in the shape of bars the time available and spent for each active
    allocation.
    """

    data = {key: np.nan_to_num(data[key]) for key in data}

    # Create the first plot with the allocation hours
    alloc_names = data['allocations']
    categories = ["spent_hours", "free_hours"]
    colors = ["#e84d60", "darkgreen"]  # "#c9d9d3"

    n_names = len(alloc_names)

    source = ColumnDataSource(data=data)
    p = figure(x_range=alloc_names, plot_height=420, plot_width=80 * 8,
               title="Time spent/available for SEDM allocations this term",
               toolbar_location=None, tools="")

    p.vbar_stack(categories, x='allocations', width=0.9, color=colors,
                 source=source, legend_label=["Spent", "Available"])
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

    colors = n_names * ['#084594']
    '''for i, p in enumerate(percentage):
        if p<50: colors[i] = '#22A784'
        elif p>50 and p<75: colors[i] = '#FD9F6C'
        else: colors[i] = '#DD4968'''

    source = ColumnDataSource(data=dict(alloc_names=alloc_names,
                                        percentage=percentage, color=colors))

    p2 = figure(x_range=alloc_names, y_range=(0, 100), plot_height=420,
                plot_width=80 * 8,
                title="Percentage of time spent",
                toolbar_location=None, tools="")

    p2.vbar(x='alloc_names', top='percentage', width=0.9, color='color',
            source=source)

    p2.xgrid.grid_line_color = None
    # p2.legend.orientation = "horizontal"  # no legend for this plot
    # p2.legend.location = "top_center"
    p2.yaxis.axis_label = '% time spent'
    p2.xaxis.major_label_orientation = 0.3

    # Create the pie charts
    pie_colors = 10 * ["red", "green", "blue", "orange", "yellow", 'lime',
                       'brown', 'cyan', 'magenta', 'olive', 'black', 'teal',
                       'gold', 'crimson', 'moccasin', 'greenyellow', 'navy',
                       'ivory', 'lightpink']

    # First one with the time spent

    # define starts/ends for wedges from percentages of a circle
    percents_only = np.round(np.array(list(data["spent_hours"] /
                                           np.sum(data["spent_hours"])))
                             * 100, 1)
    percents = np.cumsum([0] + list(data["spent_hours"] /
                                    np.sum(data["spent_hours"])))
    starts = [per * 2 * np.pi for per in percents[:-1]]
    ends = [per * 2 * np.pi for per in percents[1:]]

    p3 = figure(x_range=(-1, 2.5), y_range=(-1.1, 1.1), plot_height=420,
                plot_width=600, title="% spent")

    # Add individual wedges:
    for i in range(n_names):
        p3.wedge(x=0, y=0, radius=.9, start_angle=starts[i], end_angle=ends[i],
                 color=pie_colors[i],
                 legend_label="[{0}%] {1}".format(percents_only[i], alloc_names[i]))

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
    percents_only = np.round(np.array(list(data["alloc_hours"] /
                                           np.sum(data["alloc_hours"])))
                             * 100, 1)
    percents = np.cumsum([0] + list(data["alloc_hours"] /
                                    np.sum(data["alloc_hours"])))
    starts = [per * 2 * np.pi for per in percents[:-1]]
    ends = [per * 2 * np.pi for per in percents[1:]]

    p4 = figure(x_range=(-1, 2.5), y_range=(-1.1, 1.1), plot_height=420,
                plot_width=600,
                title="% time allocated to each program")
    # Add individual wedges:
    for i in range(n_names):
        p4.wedge(x=0, y=0, radius=.9, start_angle=starts[i], end_angle=ends[i],
                 color=pie_colors[i],
                 legend_label="[{0}%] {1}".format(percents_only[i], alloc_names[i]))

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
def get_ab_what(obsdir):
    """get a pseudo what list for A/B cubes"""
    ablist = []
    cubes = glob.glob(os.path.join(obsdir, "e3d_crr_b_ifu*.fits"))
    for e3df in cubes:
        # get root filename
        rute = '_'.join(e3df.split('/')[-1].split('_')[1:7])
        # is this a standard single cube?
        crrf = glob.glob(os.path.join(obsdir, rute + '.fit*'))
        if len(crrf) > 0:
            continue
        fname = '_'.join(e3df.split('/')[-1].split('_')[3:7]) + '.fits'
        targ = e3df.split('/')[-1].split('_')[7].split('.fit')[0]
        ablist.append("   "+fname+" (1.000/0.1/1.0 s): " + targ + " [A]")
    return ablist


def get_ifu_products(obsdir=None, user_id=None, obsdate="", show_finder=True,
                     product_type='all', camera_type='ifu'):
    """

    :param obsdir: 
    :param user_id: 
    :param obsdate: 
    :param product_type:
    :param show_finder:
    :param camera_type:
    :return: 
    """
    # ifu_dict = {}
    if product_type:
        pass
    if camera_type:
        pass

    if not obsdate:
        obsdate = datetime.datetime.utcnow().strftime("%Y%m%d")
    else:
        obsdate = obsdate[0]

    if not obsdir:
        obsdir = '%s%s/' % (redux_dir, obsdate)
    else:
        obsdir = obsdir[0]

    # Look first to make sure there is a data directory.
    if not os.path.exists(obsdir):
        return {'message': 'No data directory could be located for %s UT' %
                           os.path.basename(os.path.normpath(obsdir)),
                'obsdate': obsdate}

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
    # print(calib_dict, 'calib products')

    if user_id == 2:    # SEDM_admin
        if os.path.exists(os.path.join(obsdir, 'report.txt')):
            ext_report = """<a href="http://minar.caltech.edu/data_r/redux/{0}/report.txt">Extraction</a>""".format(obsdate)
        else:
            ext_report = ""
        if os.path.exists(os.path.join(obsdir, 'report_ztf_fritz.txt')):
            frz_report = """<a href="http://minar.caltech.edu/data_r/redux/{0}/report_ztf_fritz.txt">Fritz</a>""".format(obsdate)
        else:
            frz_report = ""
        if os.path.exists(os.path.join(obsdir, 'report_ztf_growth.txt')):
            grw_report = """<a href="http://minar.caltech.edu/data_r/redux/{0}/report_ztf_growth.txt">Growth</a>""".format(obsdate)
        else:
            grw_report = ""
        if os.path.exists(os.path.join(obsdir, 'what.txt')):
            wha_report = """<a href="http://minar.caltech.edu/data_r/redux/{0}/what.txt" type="plain/text">What</a>""".format(obsdate)
        else:
            wha_report = ""
        div_str += """<div class="row">"""
        div_str += """<h4>Reports</h4>"""
        div_str += """{0} {1} {2} {3}""".format(ext_report, frz_report,
                                                grw_report, wha_report)
        div_str += "</div>"

    div_str += """<div class="row">"""
    div_str += """<h4>Calibrations</h4>"""
    for k, v in calib_dict.items():
        impath = "/data/%s/%s" % (obsdate, os.path.basename(v))
        impathlink = "/data/%s/%s" % (obsdate,
                                      os.path.basename(v.replace('.png',
                                                                 '.pdf')))
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

    # Go through the what list and return all non-calibration entries
    with open(os.path.join(obsdir, 'what.list')) as f:
        what_list = f.read().splitlines()

    if os.path.exists(os.path.join(obsdir, 'abpairs.tab')):
        what_list.extend(get_ab_what(obsdir))

    what_list.sort()

    science_list = []
    standard_list = []
    for targ in what_list:
        if 'Calib' in targ:
            pass
        elif '[A]' in targ or '[B]' in targ or 'STD' in targ:
            science_list.append(targ)
        elif 'STD' in targ:
            pass
            # standard_list.append(targ)
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
            object_id = False
            target_requests = False
            # Start by pulling up all request that match the science target
            targ_name = sci_targ.split(':')[1].split()[0]
            if 'STD' not in targ_name:
                # 1. Get the object id

                object_ids = db.get_object_id_from_name(targ_name)

                if len(object_ids) == 1:
                    object_id = object_ids[0][0]
                elif len(object_ids) > 1:
                    # TODO what really needs to happen here is that we need to
                    # TODO cont: find the id that is closest to the obsdate.
                    # TODO cont: For now I am just going to use last added
                    # print(object_ids)
                    object_id = object_ids[-1][0]
                elif not object_ids and ('at' in targ_name.lower()
                                         or 'sn' in targ_name.lower()):
                    # sometimes it's at 2018abc not at2018abc in the db
                    targ_name = targ_name[:2] + ' ' + targ_name[2:]
                    object_ids = db.get_object_id_from_name(targ_name)
                    try:
                        object_id = object_ids[-1][0]
                    except IndexError:
                        object_id = False
                        # print("There was an error. You can't see this")

                # If we are not the admin then we need to check
                # if the user can see the object

                if user_id not in [2, 20200227202025683]:

                    if object_id:
                        target_requests = db.get_from_request(
                            values=['allocation_id'],
                            where_dict={'object_id': object_id,
                                        'status': 'COMPLETED'})

                    if not target_requests:
                        target_requests = db.get_from_request(
                            values=['allocation_id'],
                            where_dict={'object_id': object_id,
                                        'status': 'OBSERVED'})

                    # print("Object id", object_id)
                    # Right now I am only seeing if there exists a match between
                    # allocations of all request.  It's possible the request
                    # could have been made by another group as another follow-up
                    # and thus the user shouldn't be able to see it.  This
                    # should be able to be fixed once all request are listed in
                    # the headers of the science images.
                    for req in target_requests:
                        # print(sci_targ, targ_name)
                        # print(allocation_id_list,
                        # "List of allocations this person can see")
                        if req[0] in allocation_id_list:
                            show_list.append((sci_targ, targ_name))
                        else:
                            print("You can't see this at allocation id list")
                else:
                    show_list.append((sci_targ, targ_name))
            else:
                targ_name = sci_targ.split(':')[1].split()[0].replace('STD-',
                                                                      '')
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
            # print(targ)
            targ_params = targ[0].split()
            fits_file = targ_params[0].replace('.fits', '')
            name = targ[1]

            image_list = (glob.glob('%sifu_spaxels_*%s*.png' % (obsdir,
                                                                fits_file)) +
                          glob.glob('%simage_%s*.png' % (obsdir, name)))

            spec_list = (glob.glob('%s%s_SEDM.png' % (obsdir, name)) +
                         glob.glob('%sspec_forcepsf*%s*.png' %
                                   (obsdir, fits_file)) +
                         glob.glob('%sspec_auto*%s*.png' % (obsdir, fits_file)))

            e3d_list = (glob.glob('%se3d*%s*.fits' % (obsdir, fits_file)))

            spec_ascii_list = (glob.glob('%sspec_forcepsf*%s*.txt' %
                                         (obsdir, fits_file)) +
                               glob.glob('%sspec_auto*%s*.txt' % (obsdir,
                                                                  fits_file)))

            fluxcals = (glob.glob('%sfluxcal_*%s*.fits' % (obsdir, fits_file)))

            if name not in science_dict:
                science_dict[name] = {'image_list': image_list,
                                      'spec_list': spec_list,
                                      'e3d_list': e3d_list,
                                      'spec_ascii_list': spec_ascii_list,
                                      'fluxcals': fluxcals}
            else:
                # We do this to handle cases where there are two or more of
                # the same object name
                science_dict[name+'_xRx_%s' % str(count)] = {
                    'image_list': image_list, 'spec_list': spec_list,
                    'e3d_list': e3d_list, 'spec_ascii_list': spec_ascii_list,
                    'fluxcals': fluxcals}
            count += 1
        # Alright now we build the table that will show the spectra, image file
        # and classification.

        # count = 0

        for obj, obj_data in science_dict.items():
            if '_xRx_' in obj:
                obj = obj.split('_xRx_')[0]

            if 'ZTF' in obj:
                obj_link = ('<a href="https://fritz.science/source/'
                            '%s">%s</a>' %
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
                if obj_data['spec_ascii_list']:
                    for j in obj_data['spec_ascii_list']:
                        impath = "/data/%s/%s" % (obsdate, os.path.basename(j))
                        div_str += ('<div class="col-md-{2}">'
                                    '<a href="%s">ASCII Spec File</a>'
                                    '</div>' % impath)
                if obj_data['fluxcals']:
                    for j in obj_data['fluxcals']:
                        impath = "/data/%s/%s" % (obsdate, os.path.basename(j))
                        div_str += ('<div class="col-md-{2}">'
                                    '<a href="%s">Flux calibration file</a>'
                                    '</div>' % impath)
            # ToDO: Grab data from somewhere to put in the meta data column
            if obj_data['image_list']:
                for i in obj_data['image_list']:

                    impath = "/data/%s/%s" % (obsdate, os.path.basename(i))
                    impathlink = "/data/%s/%s" % (
                        obsdate, os.path.basename(i.replace('.png', '.pdf')))
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
                # Check if finders exists in redux directory and if not then
                # log at the old phot directory location
                path1 = os.path.join(redux_dir, obsdate, 'finders')
                path2 = os.path.join(phot_dir, obsdate, 'finders')

                if os.path.exists(path1):
                    finder_path = path1
                else:
                    finder_path = path2

                if os.path.exists(finder_path):
                    finder_img = glob.glob(finder_path + '/*%s*.png' % obj)
                    if finder_img:
                        impathlink = "/data/%s/%s" % (
                            obsdate, os.path.basename(finder_img[-1]))
                        div_str += """<div class="col-md-{0}">
                                              <div class="thumbnail">
                                                <a href="{1}">
                                                  <img src="{2}" width="{3}px" height="{4}px">
                                                </a>
                                              </div>
                                      </div>""".format(4, impathlink,
                                                       impathlink, 250, 250)
            if obj_data['spec_list']:
                for i in obj_data['spec_list']:
                    impath = "/data/%s/%s" % (obsdate, os.path.basename(i))
                    impathlink = "/data/%s/%s" % (
                        obsdate, os.path.basename(i.replace('.png', '.pdf')))
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


def get_rc_products(obsdate=None, product=None, user_id=None, camera_type='rc'):
    """
    :param obsdate:
    :param product:
    :param user_id:
    :param camera_type:
    :return:
    """
    if user_id:
        pass
    if camera_type:
        pass
    # print(product, 'product')
    raw_png_dir = ['acquisition', 'bias', 'dome', 'focus',
                   'guider_images', 'guider_movies', 'twilight',
                   'science_raw']

    sedm_dict = {}
    if not obsdate:
        obsdate = datetime.datetime.utcnow().strftime("%Y%m%d")
        sedm_dict['obsdate'] = obsdate
    elif isinstance(obsdate, list):
        obsdate = obsdate[0]
        sedm_dict['obsdate'] = obsdate

    if not product:
        product = 'science'

    display_dict = {}
    ext = '*.png'

    sci_path = None
    if product.lower() == 'science':
        # print(new_phot_dir, obsdate)
        sci_path = os.path.join(new_phot_dir, obsdate, 'reduced', 'png')
        if not os.path.exists(sci_path):
            # print("Path doesn't exist", sci_path)
            sedm_dict['data'] = "No %s images found" % product
    elif product.lower() == 'acquisition':
        sci_path = os.path.join(new_phot_dir, obsdate, 'reduced', 'png')
        if not os.path.exists(sci_path):
            # print("Path doesn't exist", sci_path)
            sedm_dict['data'] = "No %s images found" % product
    elif product.lower() in raw_png_dir:
        # print(new_phot_dir, obsdate)
        if 'guider' in product.lower():
            p_split = product.split("_")
            if p_split[-1] == 'movies':
                ext = '*.gif'
            product = 'guider'

        sci_path = os.path.join(new_phot_dir, obsdate,
                                'pngraw', product.lower().replace('_raw', ''))
        # print(sci_path, "Science path in alt")
        if not os.path.exists(sci_path):
            print("Path doesn't exist")

    # print("Looking in directory:", sci_path)
    find_path = os.path.join(sci_path, ext)
    # print(find_path, 'find_path')
    files = glob.glob(find_path)

    # print("Files found", files)

    for file in files:
        base_name = os.path.basename(file).replace(".png", "")
        if product.lower() == 'science' and 'ACQ' in base_name:
            continue
        elif product.lower() == 'science':
            filters = base_name.split("_")
            if filters[-1] == "0":
                objfilt = filters[-3]
                imgfilt = filters[-2]
            else:
                objfilt = filters[-2]
                imgfilt = filters[-1]

            if objfilt == imgfilt:
                if 'data' in display_dict:
                    display_dict['data'].append(file)
                else:
                    display_dict['data'] = [file]
        elif product.lower() == 'acquisition':
            if 'ACQ' not in base_name:
                continue
            elif "_r_r" not in base_name and "_NA_r" not in base_name:
                continue
            else:
                if 'data' in display_dict:
                    display_dict['data'].append(file)
                else:
                    display_dict['data'] = [file]
        else:
            if 'data' in display_dict:
                display_dict['data'].append(file)
            else:
                display_dict['data'] = [file]

    div_str = ''

    if user_id == 2:    # SEDM_admin
        obsdir = os.path.join(new_phot_dir, obsdate)
        if os.path.exists(os.path.join(obsdir, 'rcwhat.txt')):
            wha_report = """<a href="http://minar.caltech.edu/data_r/redux/phot/{0}/rcwhat.txt" type="plain/text">RCWhat</a>""".format(obsdate)
            div_str += """<div class="row">"""
            div_str += """<h4>{0}</h4>""".format(wha_report)
            div_str += "</div>"

    if 'data' in display_dict:
        count = 100
        for fil in sorted(display_dict['data']):
            # fil = fil.replace(base_dir, '')
            impath = "/data_r/%s" % fil.replace(base_dir, '')

            if 'reduced' in fil:
                fits_suffix = '.fits'
                if os.path.exists(fil.replace('/png', '').replace('.png',
                                                                  '.fits.gz')):
                    fits_suffix = '.fits.gz'
                fil = fil.replace(base_dir, '')
                impathlink = "/data_r/%s" % fil.replace('/png/', '/').replace(
                    '.png', fits_suffix)
            elif 'pngraw' in fil and '.gif' not in fil:
                base_link = fil.replace(base_dir, '').split('/pngraw/')[0]
                fits_suffix = '.fits'
                png_suffix = '_all.png'
                if 'Bias' in fil or 'Flat' in fil:
                    png_suffix = '.png'
                if os.path.exists(
                        os.path.join(
                            new_phot_dir, obsdate,
                            os.path.basename(fil).replace(png_suffix,
                                                          '.fits.gz'))):
                    fits_suffix = '.fits.gz'
                fil = fil.replace(base_dir, '')
                impathlink = "/data_r/%s" % \
                             os.path.join(base_link,
                                          os.path.basename(fil).replace(
                                              png_suffix, fits_suffix))
            else:
                impathlink = "/data_r/%s" % fil.replace(base_dir, '')

            div_str += """<div class="col-sm-4"><div class="card">
                            <a href="{1}?image={4}" data-toggle="lightbox" data-gallery="example-gallery">
                                <img style="width:300px" class="card-img-top" src="{1}?image{4}" alt="Card image">
                            </a>
                            <div class="cardbody">
                                <h6 class="card-title">{2}</h6>
                                <a href="http://minar.caltech.edu{0}" class="btn btn-primary">
                                    Download
                                </a> 
                            </div>

                        </div></div>""".format(impathlink, impath,
                                               os.path.basename(fil), impath,
                                               count)

            count += 1

    else:
        div_str += "<p>No %s images found" % product

    sedm_dict['data'] = div_str

    # print(sedm_dict)
    return sedm_dict


###############################################################################
# THIS SECTION HANDLES THE ACTIVE_VISIBILITIES PAGE.                          #
# KEYWORD:VISIBILITIES #???                                                   #
###############################################################################
def get_pending_visibility(userid):
    sedm_dict = {'enddate': datetime.datetime.utcnow() +
                 datetime.timedelta(days=1),
                 'inidate': datetime.datetime.utcnow() -
                 datetime.timedelta(days=3, hours=8)}

    # 1. Get a dataframe of all requests for the current user
    reqs = get_requests_for_user(userid, sedm_dict['inidate'],
                                 sedm_dict['enddate'])

    # organize requests into dataframes by whether they are pending or not
    pending = reqs[(reqs['status'] == 'PENDING')]

    # retrieve information about the user's allocations
    # ac = get_allocations_user(userid)

    # Create html tables
    sedm_dict['pending'] = {'table': fancy_request_table(pending),
                            'title': 'Pending Requests'}

    sedm_dict['script'], sedm_dict['div'] = plot_visibility(userid, sedm_dict)
    return sedm_dict


def plot_visibility(userid, sedm_dict, obsdate=None):
    """
    plots visibilities for pending/active requests at the current date.
    Will be adapted to plot previous observations and arbitrary objects userid:
    user whose allocations will be shown in color with details. Others will be
    greyed out

    userid: <int>
    sedm_dict: <dict>
        should have ['active']['table'] and ['enddate'] and ['inidate']
    obsdate: <str> YYYYMMDD
        if "None", will use current date

    returns: components of a bokeh figure with the appropriate plot
    """

    allocpalette = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a',
                    '#b15928', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f',
                    '#cab2d6', '#ffff99']
    reqs = get_requests_for_user(2,
                                 sedm_dict['inidate'],
                                 sedm_dict['enddate'])  # admin
    active = reqs[(reqs['status'] == 'PENDING') | (reqs['status'] == 'ACTIVE')]
    # ['allocation', 'object', 'RA', 'DEC', 'start date', 'end date',
    # 'priority', 'status', 'lastmodified', 'obs_seq', 'exptime', 'UPDATE']

    allowed_allocs = get_allocations_user(userid)
    active['allocation'].mask(~np.in1d(active['allocation'],
                                       allowed_allocs['allocation']),
                              other='other', inplace=True)

    programs = {i['allocation']: i['program']
                for _, i in allowed_allocs.iterrows()}
    programs['other'] = 'other'
    # this needs to be alphabetical for the legend to look correct
    active.sort_values('allocation')

    p = figure(plot_width=700, plot_height=500, toolbar_location='above',
               y_range=(0, 90), y_axis_location="right")

    # setup with axes, sun/moon, frames, background
    # TODO Dima says to never ever use SkyCoord in production code
    palomar_mountain = EarthLocation(lon=243.1361*u.deg, lat=33.3558*u.deg,
                                     height=1712*u.m)
    utcoffset = -7 * u.hour  # Pacific Daylight Time

    # plotting a single object, or the pending objects in future
    if obsdate is None:
        time = (Time.now() - utcoffset).datetime  # date is based on local time
        time = Time(datetime.datetime(time.year, time.month, time.day))
    else:  # past observations on a particular night
        time = Time(datetime.datetime(int(obsdate[:4]), int(obsdate[4:6]),
                                      int(obsdate[6:8])))
        # all_requests = reqs[reqs['status'] == 'COMPLETED']
        # all_requests = all_requests[time - 12 * u.hour
        #                            <= all_requests['startdate']
        #                            < time + 12 * u.hour]
    midnight = time - utcoffset  # 7am local time of correct date, midnight UTC

    delta_midnight = np.linspace(-8, 8, 500) * u.hour
    t = midnight + delta_midnight
    abstimes = np.asarray([i.datetime.strftime('%I:%M %p')
                           for i in t + utcoffset])

    frame = AltAz(obstime=t, location=palomar_mountain)
    sun_alt = get_sun(t).transform_to(frame).alt
    moon_alt = get_moon(t).transform_to(frame).alt

    # shading for nighttime and twilight
    dark_times = delta_midnight[sun_alt < 0].value
    twilit_times = delta_midnight[sun_alt < -18 * u.deg].value
    plotted_times = delta_midnight[sun_alt < 5 * u.deg].value

    twilight = BoxAnnotation(left=min(twilit_times), right=max(twilit_times),
                             bottom=0, fill_alpha=0.15, fill_color='black',
                             level='underlay')
    night = BoxAnnotation(left=min(dark_times),    right=max(dark_times),
                          bottom=0, fill_alpha=0.25, fill_color='black',
                          level='underlay')
    earth = BoxAnnotation(top=0, fill_alpha=0.8, fill_color='sienna')

    p.add_layout(night)
    p.add_layout(twilight)
    p.add_layout(earth)

    # sun
    # p.line(delta_midnight, sun_alt,  line_color='red', name="Sun",
    #        legend='Sun', line_dash='dashed')
    # moon
    p.line(delta_midnight, moon_alt, line_color='yellow', line_dash='dashed',
           name="Moon", legend_label='Moon')
    # labels and axes
    p.title.text = "Visibility for %s UTC" % midnight
    p.xaxis.axis_label = "Hours from PDT Midnight"
    p.x_range.start = min(plotted_times)
    p.x_range.end = max(plotted_times)
    p.yaxis.axis_label = "Airmass"

    # primary airmass label on right
    airmasses = (1.01, 1.1, 1.25, 1.5, 2., 3., 6.)
    ticker = [90 - np.arccos(1./i) * 180/np.pi for i in airmasses]
    p.yaxis.ticker = ticker
    p.yaxis.major_label_overrides = {tick: str(airmasses[i])
                                     for i, tick in enumerate(ticker)}

    # add supplementary alt label on left
    p.extra_y_ranges = {"altitude": Range1d(0, 90)}
    p.add_layout(LinearAxis(y_range_name="altitude",
                            axis_label='Altitude [deg]'), 'left')

    ##########################################################################
    # adding data from the actual objects
    # objs = SkyCoord(np.array(ras,  dtype=np.float),
    #                np.array(decs, dtype=np.float), unit="deg")

    approx_midnight = int(Time.now().jd - .5) + .5 - utcoffset.value/24.
    palo_sin_lat = 0.549836545
    palo_cos_lat = 0.835272275
    palo_long = 243.1362

    alloc_color = {}
    for i, val in allowed_allocs.iterrows():
        alloc_color[val['allocation']] = allocpalette[i % len(allocpalette)]
    alloc_color['other'] = 'lightgray'

    tooltipped = []  # things with tooltips
    # make it #name when we get to bokeh 0.13
    tooltips = [('obj',        '@name'),
                ('time',       '@abstime'),
                ('altitude',   u"@alt\N{DEGREE SIGN}"),
                ('airmass',    '@airmass')]

    for _, req in active.iterrows():
        # iterrows doesn't preserve datatypes and turns ra, dec into decimals?
        req['ra'] = float(req['RA'])
        req['dec'] = float(req['DEC'])
        color = alloc_color[req['allocation']]
        # vvv I got this formula from some website for the navy
        # but forgot to copy the url
        alt = 180 / np.pi * np.arcsin(palo_cos_lat *
                                      np.cos(np.pi/180 *
                                             (palo_long - req['ra'] + 15 *
                                              (18.697374558 + 24.06570982 *
                                               (delta_midnight.value/24. +
                                                approx_midnight - 2451545)))) *
                                      np.cos(req['dec'] * np.pi/180) +
                                      palo_sin_lat * np.sin(req['dec'] *
                                                            np.pi/180))
        airmass = 1./np.cos((90 - alt) * np.pi/180)
        source = ColumnDataSource(dict(times=delta_midnight, alt=alt,
                                       airmass=airmass, abstime=abstimes,
                                       priority=np.full(len(t),
                                                        int(req['priority'])),
                                       alloc=np.full(len(t),
                                                     req['allocation'][6:]),
                                       name=np.full(len(abstimes),
                                                    req['object'])))
        # delete the name when we get to bokeh 0.13
        if len(active) == 1:  # single object
            legend = req['object']
            line_width = 5
        else:
            legend = '{}'.format(programs[req['allocation']])
            # tooltips += [('priority',   '@priority'),
            # ('allocation', '@alloc')]

            # plot that highlights observed part of the night
            if req['status'] == 'COMPLETED':
                # full path of the night
                dotted = p.line('times', 'alt', color=color, source=source,
                                line_dash='2 2', name=req['object'],
                                line_width=1, legend_label=legend)
                # manually crop the source so only thick observed
                # part has tooltips
                endtime = req['lastmodified']
                # TODO sometimes it's 2ifu or no ifu
                exptime = {req['obs_seq'][i]: req['exptime'][i]
                           for i in range(len(req['obs_seq']))}['1ifu']
                initime = endtime - exptime * u.second

                mask = np.logical_and(
                    delta_midnight + midnight + utcoffset > initime,
                    delta_midnight + midnight + utcoffset < endtime)
                source = ColumnDataSource(pd.DataFrame(source.data)[mask])
                # all it changes is the line width
                line_width = int(req['priority'] + 3)
            else:
                line_width = int(req['priority'])

        path = p.line('times', 'alt', color=color, source=source,
                      name=''.format(req['object']),
                      line_width=line_width, legend_label=legend)
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
        if statsfile is not None and mydate is not None:
            stats_plot = plot_stats(statsfile, mydate)
            if stats_plot is None:
                message += " No statistics log found up to 100 days prior to" \
                           " today... Weather has been terrible lately!"
                script, div = None, None
            else:
                message += " Weather statistics for last opened day: %s" % (
                    os.path.basename(os.path.dirname(os.path.dirname(
                        statsfile))))
                script, div = components(stats_plot)
        else:
            script, div = None, None
    else:
        mydate_in = obsdate.replace("-", "")

        # Just making sure that we have only allowed digits in the date
        mydate = re.findall(r"(2\d{3}[0-1]\d[0-3]\d)", mydate_in)
        if len(mydate) == 0:
            message += "Incorrect format for the date! Your input is: %s." \
                       " Shall be YYYYMMDD. \n" % mydate_in
            script, div = "", ""
        else:
            mydate = mydate[0]
            message = ""

            statsfile, mydate_out = search_stats_file(mydate)
            stats_plot = plot_stats(statsfile, mydate)
            if not statsfile:
                message = message + "No statistics log found for the date %s." \
                                    " Showing P18 data." % mydate
                script, div = components(stats_plot)

            else:
                stats_plot = plot_stats(statsfile, mydate)
                message = message + "Weather statistics for selected day: %s"\
                    % mydate
                script, div = components(stats_plot)

    return {'script': script, 'div': div, 'message': message}


def search_stats_file(mydate=None):
    """
    Returns the last stats file that is present in the system according to
    the present date. It also returns a message stating what date that was.
    """
    # If the date is specified, we will try to locate the right file.
    # None will be returned if it does not exist.
    if mydate:
        s = os.path.join(phot_dir, mydate, "stats/stats.log")
        if os.path.exists(s):
            if os.path.getsize(s) > 0:
                return s, mydate
            else:
                return None, None
        else:
            s = os.path.join(new_phot_dir, mydate, "stats/stats.log")
            if os.path.exists(s):
                if os.path.getsize(s) > 0:
                    return s, mydate
                else:
                    return None, None
            return None, None

    else:
        curdate = datetime.datetime.utcnow()
        # Try to find the stat files up to 100 days before today's date.
        i = 0
        while i < 100:
            newdate = curdate
            newdatedir = "%d%02d%02d" % (newdate.year, newdate.month,
                                         newdate.day)
            s = os.path.join(phot_dir, newdatedir, "stats/stats.log")
            s_new = os.path.join(new_phot_dir, newdatedir, "stats/stats.log")
            if os.path.exists(s):
                if os.path.getsize(s) > 0:
                    return s, newdatedir
                # else:
                #    return None, None
            elif os.path.exists(s_new):
                if os.path.getsize(s_new) > 0:
                    return s_new, newdatedir
                # else:
                #    return None, None
            i = i + 1
            curdate -= datetime.timedelta(days=1)
        return None, None


def load_p18seeing(obsdate):

    obtime, seeing = get_p18obsdata(obsdate)

    local_date = np.array(obtime)

    d = pd.DataFrame({'date': local_date, 'seeing': seeing})

    return d


def load_stats(statsfile='stats.log'):
    data = pd.read_csv(statsfile, header=None,
                       names=['path', 'obj', 'jd', 'ns', 'fwhm', 'ellipticity',
                              'bkg', 'airmass', 'in_temp', 'imtype', 'out_temp',
                              'in_hum'])

    jds = data['jd']
    t = Time(jds, format='jd', scale='utc')
    date = t.utc.datetime
    day_frac_diff = datetime.timedelta(
        np.ceil((datetime.datetime.now() -
                 datetime.datetime.utcnow()).total_seconds()) / 3600 / 24)
    local_date = date + day_frac_diff

    data2 = data.assign(localdate=local_date)
    data2.set_index('localdate')

    return pd.DataFrame(
        {'date': data2['localdate'], 'ns': data2['ns'], 'fwhm': data2['fwhm'],
         'ellipticity': data2['ellipticity'], 'bkg': data2['bkg'],
         'airmass': data2['airmass'], 'in_temp': data2['in_temp'],
         'imtype': data2['imtype'], 'out_temp': data2['out_temp'],
         'in_hum': data2['in_hum']})


def plot_stats(statsfile, mydate):

    source = ColumnDataSource(
        data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[],
                  in_temp=[], imtype=[], out_temp=[], in_hum=[]))
    source_static = ColumnDataSource(
        data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[],
                  in_temp=[], imtype=[], out_temp=[], in_hum=[]))

    view_science = CDSView(source=source,
                           filters=[GroupFilter(column_name='imtype',
                                                group='SCIENCE')])
    view_acquisition = CDSView(source=source,
                               filters=[GroupFilter(column_name='imtype',
                                                    group='ACQUISITION')])
    view_guider = CDSView(source=source,
                          filters=[GroupFilter(column_name='imtype',
                                               group='GUIDER')])
    view_focus = CDSView(source=source,
                         filters=[GroupFilter(column_name='imtype',
                                              group='FOCUS')])
    source_p18 = ColumnDataSource(data=dict(date=[], seeing=[]))

    def update(selected=None):

        if selected:
            pass

        if statsfile:
            data = load_stats(statsfile)
            source.data = source.from_df(data[['date', 'ns', 'fwhm',
                                               'ellipticity', 'bkg', 'airmass',
                                               'in_temp', 'imtype', 'out_temp',
                                               'in_hum']])
            source_static.data = dict(source.data)
            excl_focus = data[data['imtype'] != 'FOCUS']
            print("\nIFU min, max, median FWHM: ", excl_focus['fwhm'].min(),
                  excl_focus['fwhm'].max(), excl_focus['fwhm'].median())

        p18 = load_p18seeing(mydate)
        source_p18.data = source_p18.from_df(p18[['date', 'seeing']])
        source_p18.data['seeing_min'] = p18['seeing'].min()
        source_p18.data['seeing_max'] = p18['seeing'].max()
        source_p18.data['seeing_median'] = p18['seeing'].median()
        source_static_p18.data = dict(source_p18.data)
        print("\nP18 min, max, median FWHM: ", p18['seeing'].min(),
              p18['seeing'].max(), p18['seeing'].median())

    source_static_p18 = ColumnDataSource(data=dict(date=[], seeing=[]))
    tools = 'pan,box_zoom,reset'

    p18seeing = figure(plot_width=425, plot_height=250, tools=tools,
                       x_axis_type='datetime', active_drag="box_zoom")
    p18seeing.circle('date', 'seeing', source=source_static_p18, color="black")
    p18seeing.title.text = "P18 seeing [arcsec]"
    # seeing_median = Span(location='seeing_median', dimension='width',
    #                     source=source_static_p18, line_dash='dotted')
    # p18seeing.add_layout(seeing_median)

    if statsfile:
        ns = figure(plot_width=425, plot_height=250, tools=tools,
                    x_axis_type='datetime', active_drag="box_zoom")
        ns.line('date', 'ns', source=source_static)
        ns.circle('date', 'ns', size=1, source=source, color=None,
                  selection_color="orange")
        ns.title.text = "Number of bright sources extracted"

        bkg = figure(plot_width=425, plot_height=250, tools=tools,
                     x_axis_type='datetime', active_drag="box_zoom")
        bkg.x_range = ns.x_range
        bkg.line('date', 'bkg', source=source_static)
        bkg.circle('date', 'bkg', size=1, source=source, color=None,
                   selection_color="orange")
        bkg.title.text = "Background (counts)"

        temp = figure(plot_width=425, plot_height=250, tools=tools,
                      x_axis_type='datetime', active_drag="box_zoom")
        temp.x_range = ns.x_range
        temp.line('date', 'in_temp', source=source_static, color='blue',
                  legend_label="Inside")
        temp.line('date', 'out_temp', source=source_static, color='green',
                  legend_label="Outside")
        temp.circle('date', 'in_temp', size=1, source=source, color=None,
                    selection_color="orange")
        temp.title.text = "Temperature [C]"

        fwhm = figure(plot_width=425, plot_height=250, tools=tools,
                      x_axis_type='datetime', active_drag="box_zoom")
        fwhm.x_range = ns.x_range
        fwhm.circle('date', 'fwhm', source=source, color="green",
                    legend_label="Focus", view=view_focus)
        fwhm.circle('date', 'fwhm', source=source, color="red",
                    legend_label="Science", view=view_science)
        fwhm.circle('date', 'fwhm', source=source, color="blue",
                    legend_label="Acquisition", view=view_acquisition)
        fwhm.circle('date', 'fwhm', source=source, color="black",
                    legend_label="Guider", view=view_guider)
        fwhm.circle('date', 'fwhm', size=1, source=source, color=None,
                    selection_color="orange")
        fwhm.title.text = "P60 FWHM [arcsec]"

        airmass = figure(plot_width=425, plot_height=250, tools=tools,
                         x_axis_type='datetime', active_drag="box_zoom")
        airmass.x_range = ns.x_range
        airmass.line('date', 'airmass', source=source_static)
        airmass.circle('date', 'airmass', size=1, source=source, color=None,
                       selection_color="orange")
        airmass.title.text = "Airmass"

        ellipticity = figure(plot_width=425, plot_height=250, tools=tools,
                             x_axis_type='datetime',
                             active_drag="box_zoom")
        ellipticity.x_range = ns.x_range
        ellipticity.line('date', 'ellipticity', source=source_static)
        ellipticity.circle('date', 'ellipticity', size=1, source=source,
                           color=None, selection_color="orange")
        ellipticity.title.text = "Ellipticity"

        humidity = figure(plot_width=425, plot_height=250, tools=tools,
                          x_axis_type='datetime', active_drag="box_zoom")
        humidity.x_range = ns.x_range
        humidity.line('date', 'in_hum', source=source_static)
        humidity.circle('date', 'in_hum', size=1, source=source, color=None,
                        selection_color="orange")
        humidity.title.text = "Inside Humidity [%]"

        p18seeing.x_range = ns.x_range

        left = column(fwhm, p18seeing, airmass)
        center = column(ellipticity, ns, bkg, )
        right = column(temp, humidity)
        layout = row(left, center, right)

    else:

        layout = row(column(p18seeing))

    # initialize
    update()

    curdoc().add_root(layout)
    curdoc().title = "Stats"

    return layout


def plot_not_found_message(day):
    not_found = figure(plot_width=900, plot_height=450, x_range=[0, 900],
                       y_range=[0, 450])
    not_found.image(image=[np.zeros([900, 450]) + 0.1], x=0, y=0, dw=900,
                    dh=450)
    citation = Label(x=50, y=225, x_units='screen', y_units='screen',
                     text='No statistics found for today \n '
                          '(likely we were weathered out...)')
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
        'path_redux_phot': new_phot_dir,
        'path_raw': raw_dir,
        'path_requests': requests_dir})


def get_marshal_id(marshal='growth', request_id=None):
    """
    
    :param marshal: 
    :param request_id: 
    :return: 
    """

    try:
        request_id = int(request_id)
    except Exception as e:
        return {'error': str(e)}

    ret = db.get_from_request(values=['marshal_id', 'external_id'],
                              where_dict={'id': request_id})
    if not ret:
        return {'error': "No object found with that id number"}

    if marshal == 'growth':

        ret = make_dict_from_dbget(['marshal_id', 'external_id'], ret[0])
        if isinstance(ret['marshal_id'], int) and ret['marshal_id'] <= 100:
            return {'error': "Request is not a valid growth marshal request"}
        elif isinstance(ret['marshal_id'], str):
            return {'error': ret['marshal_id']}
        elif not ret['marshal_id']:
            return {'error': ret['marshal_id']}
        elif ret['external_id'] == 2:
            return {'error': "Not a growth request"}
        else:
            return ret


def get_user_observations(username, password, obsdate):
    """

    :param username:
    :param password:
    :param obsdate:
    :return:
    """
    # print(username, type(username))

    ret = check_login(username, password)

    # print(ret)
    if not ret[0]:
        return {'message': "User name and password do not match"}

    user_id = ret[1]

    obsdir = os.path.join(redux_dir, obsdate)
    obsdir += '/'

    calib_files = ['Xe.fits', 'Hg.fits', 'Cd.fits', 'dome.fits',
                   'bkgd_dome.fits', 'e3d_dome.fits', '%s_Flat.fits' % obsdate]

    pkl_list = (glob.glob('%s*.pkl' % obsdir))

    master_calib_list = []

    for file in calib_files:
        if os.path.exists(os.path.join(obsdir, file)):
            master_calib_list.append(os.path.join(obsdir, file))

    master_calib_list += pkl_list

    # print(master_calib_list, 'master')

    # Look first to make sure there is a data directory.
    if not obsdate:
        return {'message': 'No obsdate given in json request'}

    if not os.path.exists(obsdir):
        return {'message': 'No data directory could be located for %s UT' %
                           os.path.basename(os.path.normpath(obsdir)),
                'obsdate': obsdate}

    # sedm_dict = {'obsdate': obsdate, 'sci_data': ''}

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
    data_list = []

    for k, v in calib_dict.items():
        if not os.path.exists(v):
            remove_list.append(k)

    if remove_list:
        for i in remove_list:
            calib_dict.pop(i)
    # print(calib_dict, 'calib products')

    for v in master_calib_list:
        impath = "/data/%s/%s" % (obsdate, os.path.basename(v))
        data_list.append(impath)

    for k, v in calib_dict.items():
        impath = "/data/%s/%s" % (obsdate, os.path.basename(v))
        impathlink = "/data/%s/%s" % (obsdate,
                                      os.path.basename(v.replace('.png',
                                                                 '.pdf')))
        if not os.path.exists(impathlink):
            impathlink = impath

        data_list.append(impathlink)

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
            standard_list.append(targ)
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
            if user_id == 2:
                show_list.append((sci_targ, targ_name))
                continue
            if 'STD' not in targ_name:
                # 1. Get the object id

                object_ids = db.get_object_id_from_name(targ_name)

                object_id = None
                if len(object_ids) == 1:
                    object_id = object_ids[0][0]

                elif len(object_ids) > 1:
                    # TODO what really needs to happen here is that we need to
                    # TODO find the id that is closest to the obsdate.
                    # TODO For now I am just going to use last added
                    # print(object_ids)
                    object_id = object_ids[-1][0]
                elif not object_ids and ('at' in targ_name.lower()
                                         or 'sn' in targ_name.lower()):
                    # sometimes it's at 2018abc not at2018abc in the db
                    targ_name = targ_name[:2] + ' ' + targ_name[2:]
                    object_ids = db.get_object_id_from_name(targ_name)
                    try:
                        object_id = object_ids[-1][0]
                    except IndexError:
                        print("There was an error. You can't see this")

                target_requests = db.get_from_request(
                    values=['allocation_id'],
                    where_dict={'object_id': object_id, 'status': 'COMPLETED'})

                if not target_requests:
                    target_requests = db.get_from_request(
                        values=['allocation_id'],
                        where_dict={'object_id': object_id,
                                    'status': 'OBSERVED'})
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
                        print("You can't see this at target request")
            else:
                targ_name = sci_targ.split(':')[1].split()[0].replace('STD-',
                                                                      '')
                show_list.append((sci_targ, targ_name))

    if len(standard_list) >= 1:
        for std_targ in standard_list:
            targ_name = std_targ.split(':')[1].split()[0].replace('STD-', '')
            show_list.append((std_targ, targ_name))

    # We have our list of targets that we can be shown, now lets actually find
    # the files that we will show on the web page.  To make this backwards
    # compatible I have to look for two types of files
    # print(show_list, "Show list")
    if len(show_list) >= 1:
        science_dict = {}
        count = 0
        # div_str = ''

        for targ in show_list:
            # print(targ)
            targ_params = targ[0].split()
            fits_file = targ_params[0].replace('.fits', '')
            name = targ[1]
            # print(obsdir, fits_file)
            # print('%s%s_SEDM.png' % (obsdir, name))
            # print('%sspec_forcepsf*%s*.png' % (obsdir,fits_file))
            # print('%sspec_auto*%s*.png' % (obsdir, fits_file))

            image_list = (glob.glob('%sifu_spaxels_*%s*.png' % (obsdir,
                                                                fits_file)) +
                          glob.glob('%simage_%s*.png' % (obsdir, name)))

            spec_list = (glob.glob('%s%s_SEDM.png' % (obsdir, name)) +
                         glob.glob('%sspec_forcepsf*%s*.png' % (obsdir,
                                                                fits_file)) +
                         glob.glob('%sspec_auto*%s*.png' % (obsdir, fits_file)))

            spec_all_list = glob.glob("%sspec*%s*" % (obsdir, name))

            e3d_list = (glob.glob('%se3d*%s*.fits' % (obsdir, fits_file)))

            spec_ascii_list = (glob.glob('%sspec_forcepsf*%s*.txt'
                                         % (obsdir, fits_file)) +
                               glob.glob('%sspec_auto*%s*.txt' % (obsdir,
                                                                  fits_file)))

            fluxcals = (glob.glob('%sfluxcal_*%s*.fits' % (obsdir, fits_file)))

            background = (glob.glob('%sbkgd_crr_b_%s.fits' % (obsdir,
                                                              fits_file)))

            astrom_list = (glob.glob('%sguider_crr_b_%s_astrom.fits'
                                     % (obsdir, fits_file)))

            if name not in science_dict:

                science_dict[name] = {'image_list': image_list,
                                      'spec_list': spec_list,
                                      'e3d_list': e3d_list,
                                      'spec_ascii_list': spec_ascii_list,
                                      'fluxcals': fluxcals,
                                      'specall': spec_all_list,
                                      'background': background,
                                      'astrom': astrom_list}
            else:
                # We do this to handle cases where there are two or more of
                # the same object name
                science_dict[name+'_xRx_%s' % str(count)] = {
                    'image_list': image_list, 'spec_list': spec_list,
                    'e3d_list': e3d_list, 'spec_ascii_list': spec_ascii_list,
                    'fluxcals': fluxcals, 'specall': spec_all_list,
                    'background': background, 'astrom': astrom_list}
            count += 1
        # Alright now we build the table that will show the spectra, image file
        # and classification.

        # count = 0
        # print(science_dict)
        for obj, obj_data in science_dict.items():
            if '_xRx_' in obj:
                obj = obj.split('_xRx_')[0]

            if obj_data['e3d_list']:
                for j in obj_data['specall']:
                    if j.split('.')[-1] in ['fits', 'png', 'txt', 'pdf']:
                        data_list.append("/data/%s/%s" % (obsdate,
                                                          os.path.basename(j)))
                for j in obj_data['e3d_list']:
                    data_list.append("/data/%s/%s" % (obsdate,
                                                      os.path.basename(j)))

                if obj_data['spec_ascii_list']:
                    for j in obj_data['spec_ascii_list']:
                        data_list.append("/data/%s/%s" % (obsdate,
                                                          os.path.basename(j)))

                if obj_data['fluxcals']:
                    for j in obj_data['fluxcals']:
                        data_list.append("/data/%s/%s" % (obsdate,
                                                          os.path.basename(j)))
                if obj_data['background']:
                    for j in obj_data['background']:
                        data_list.append("/data/%s/%s" % (obsdate,
                                                          os.path.basename(j)))
                if obj_data['astrom']:
                    for j in obj_data['astrom']:
                        data_list.append("/data/%s/%s" % (obsdate,
                                                          os.path.basename(j)))

            # ToDO: Grab data from somewhere to put in the meta data column
            if obj_data['image_list']:
                for i in obj_data['image_list']:

                    impath = "/data/%s/%s" % (obsdate, os.path.basename(i))
                    impathlink = "/data/%s/%s" % \
                                 (obsdate, os.path.basename(i.replace('.png',
                                                                      '.pdf')))
                    if not os.path.exists(impathlink):
                        impathlink = impath
                    data_list.append(impathlink)

            # Check if finders exists in redux directory and if not then
            # log at the old phot directory location
            path1 = os.path.join(redux_dir, obsdate, 'finders')
            path2 = os.path.join(phot_dir, obsdate, 'finders')

            if os.path.exists(path1):
                finder_path = path1
            else:
                finder_path = path2

            if os.path.exists(finder_path):
                finder_img = glob.glob(finder_path + '/*%s*.png' % obj)
                if finder_img:
                    data_list.append("/data/%s/%s" %
                                     (obsdate,
                                      os.path.basename(finder_img[-1])))

            if obj_data['spec_list']:
                for i in obj_data['spec_list']:
                    impath = "/data/%s/%s" % (obsdate, os.path.basename(i))
                    impathlink = "/data/%s/%s" % \
                                 (obsdate, os.path.basename(i.replace('.png',
                                                                      '.pdf')))
                    if not os.path.exists(impathlink):
                        impathlink = impath

                    data_list.append(impathlink)
    return_dict = {'data': data_list}
    return return_dict


def get_status():
    """

    :return:
    """

    with open(status_dir+'telstatus.json') as json_file:
        try:
            data = json.load(json_file)
        except json.decoder.JSONDecodeError:
            print("JSON decode error, trying again")
            json_file.close()
            time.sleep(1)
            with open(status_dir + 'telstatus.json') as json_file2:
                try:
                    data = json.load(json_file2)
                except json.decoder.JSONDecodeError:
                    print("JSON decode error")
                    data = {}
    try:
        rc_start_time = datetime.datetime.strptime(data['rc_LastStartTime'],
                                                   '%Y-%m-%d %H:%M:%S.%f')
        rc_end_time = rc_start_time + datetime.timedelta(
            seconds=float(data['rc_ExposureTime']))
        data['rc_EndExposureTime'] = rc_end_time.strftime("%Y-%m-%d %H:%M:%S")
        data['rc_LastStartTime'] = rc_start_time.strftime("%Y-%m-%d %H:%M:%S")
    except:
        data['rc_EndExposureTime'] = "NA"
        data['rc_LastStartTime'] = "NA"

    try:
        ifu_start_time = datetime.datetime.strptime(data['ifu_LastStartTime'],
                                                    '%Y-%m-%d %H:%M:%S.%f')
        ifu_end_time = ifu_start_time + datetime.timedelta(
            seconds=float(data['ifu_ExposureTime']))
        data['ifu_EndExposureTime'] = ifu_end_time.strftime("%Y-%m-%d %H:%M:%S")
        data['ifu_LastStartTime'] = ifu_start_time.strftime("%Y-%m-%d %H:%M:%S")
    except:
        data['ifu_EndExposureTime'] = "NA"
        data['ifu_LastStartTime'] = "NA"

    # print("Last IFU exp start time: %s" % data['ifu_LastStartTime'])
    return data


def get_obstimes():
    times = schedule.get_observing_times(return_type='json')
    times['sciTime'] = '#'
    return times


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
        ifu_exptime = 2250
        r_exptime = 180
        g_exptime = 180
        i_exptime = 180
        u_exptime = 300
    elif 15 < mag < 18:
        ifu_exptime = 1800
        r_exptime = 120
        g_exptime = 120
        i_exptime = 120
        u_exptime = 300
    elif 13 < mag < 15:
        ifu_exptime = 1200
        r_exptime = 1
        g_exptime = 1
        i_exptime = 1
        u_exptime = 30
    elif 11 < mag < 13:
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
    if not obsdate:
        f = datetime.datetime.strptime(obsdate,
                                       "%Y%m%d") - datetime.timedelta(days=1)
        obsd = datetime.datetime.strptime(obsdate, "%Y%m%d")
    elif "-" in obsdate:
        f = datetime.datetime.strptime(obsdate,
                                       "%Y-%m-%d") - datetime.timedelta(days=1)
        obsd = datetime.datetime.strptime(obsdate, "%Y-%m-%d")
    else:
        f = datetime.datetime.strptime(obsdate,
                                       "%Y%m%d") - datetime.timedelta(days=1)
        obsd = datetime.datetime.strptime(obsdate, "%Y%m%d")

    y, m, d = [f.strftime("%Y"), int(f.strftime("%m")), int(f.strftime("%d"))]
    p18obsdate = "%s-%s-%s" % (y, m, d)

    # 2. Get the data from the link
    page = requests.get(
        'http://nera.palomar.caltech.edu/P18_seeing/seeing_log_%s.log'
        % p18obsdate)
    data = page.content.decode("ISO-8859-1")

    # 3. Split the page by newlines
    data = data.split('\n')

    # 4. Loop through the data and only use points that have
    # 4 or more seeing values to average
    for i in data:
        try:
            i = i.split()

            if len(i) > 5 and int(i[5]) > 4:
                d = '%s %s' % (i[1], i[0])
                p18date.append(datetime.datetime.strptime(d,
                                                          "%m/%d/%Y %H:%M:%S"))
                p18seeing.append(float(i[4]))
        except Exception as e:
            print(str(e))
            obsd = obsd.replace(hour=7)
            return [obsd], [0]
    return p18date, p18seeing


# if __name__ == "__main__":

#    x = get_ifu_products('/scr7/rsw/sedm/redux/20180827/', 189)
#    print(x)
