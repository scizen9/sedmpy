import json
# Required import to bypass the fact there is already a
# db module that gets imported from kpy
# import sys
# sys.path.append('/scr2/sedm/sedmpy/')
from string import Template
import datetime
import pandas as pd
import astroplan
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np
import os

computer = os.uname()[1]  # a quick fix
scheduler_path = None
if computer == 'pele':
    from db.SedmDb import SedmDB
    scheduler_path = \
        '/scr/rsw/sedm/projects/sedmpy/web/static/scheduler/scheduler.html'
    server = 'pharos.caltech.edu'
    port = 5432
elif computer == 'pharos':
    from db.SedmDb import SedmDB
    scheduler_path = '/scr2/sedm/sedmpy/web/static/scheduler/scheduler.html'
    server = 'localhost'
    port = 5432
elif computer == 'ether':
    from db.SedmDb import SedmDB
    scheduler_path = '/home/rsw/new/sedmpy/web/templates/scheduler_table.html'
    server = 'localhost'
    port = 22222
elif computer == 'modon':
    from db.SedmDb import SedmDB
    port = 5432

SITE_ROOT = os.path.realpath(os.path.dirname(__file__))
json_url = os.path.join(SITE_ROOT, 'config.json')


class ScheduleNight:
    def __init__(self, obsdatetime=None, config_file=json_url):
        # 1. Load config file
        with open(config_file, 'r') as fp:
            params = json.load(fp)

        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        self.site = params['site']['name']
        self.long = params['site']['longitude']
        self.lat = params['site']['latitude']
        self.elev = params['site']['elevation']
        self.obs_site = astroplan.Observer.at_site(site_name=self.site)
        self.obs_times = self.get_observing_times()
        self.running_obs_time = None
        self.conn = None
        self.cursor = None
        self.ph_db = SedmDB(host=server, port=port)
        self.target_frame = None
        self.tr_row = Template(
            """<tr id="${allocation}">
            <td>${obstime}</td>
            <td>${objname}</td>
            <td>${contact}</td>
            <td>${project}</td>
            <td>${ra}</td>
            <td>${dec}</td>
            <td>${ifu_exptime}</td>
            <td>Filters:${rc_seq}<br>Exptime:${rc_exptime}</td>
            <td>${total}</td>
            <td><a href='request?request_id=${request_id}'>+</a></td>
            </tr>""")

        self.req_query = Template(
            """SELECT r.id AS req_id, r.object_id AS obj_id,
            r.user_id, r.marshal_id, r.exptime, r.maxairmass,
            r.max_fwhm, r.min_moon_dist, r.max_moon_illum,
            r.max_cloud_cover, r.status,
            r.priority AS reqpriority, r.inidate, r.enddate,
            r.cadence, r.phasesamples, r.sampletolerance,
            r.filters, r.nexposures, r.obs_seq, r.seq_repeats,
            r.seq_completed, r.last_obs_jd, r.creationdate,
            r.lastmodified, r.allocation_id, r.marshal_id,
            o.id, o.name AS objname, o.iauname, o.ra, o."dec",
            o.typedesig, o.epoch, o.magnitude, o.creationdate,
            u.id, u.email, a.id AS allocation_id, a.inidate, a.enddate,
            a.time_spent, a.time_allocated, a.program_id, a.active,
            p.designator, p.name, p.group_id, p.pi,
            p.time_allocated, r.priority, p.inidate,
            p.enddate, pe.mjd0, pe.phasedays, pe.phi,
            r.phase, r.sampletolerance
            FROM "public".request r
            INNER JOIN "public"."object" o ON (r.object_id = o.id)
            INNER JOIN "public".users u ON (r.user_id = u.id)
            INNER JOIN "public".allocation a ON (r.allocation_id = a.id)
            INNER JOIN "public".program p ON (a.program_id = p.id)
            LEFT JOIN "public".periodic pe on (pe.object_id=o.id)
            ${where_statement}
            ${and_statement}
            ${group_statement}
            ${order_statement}"""
        )

    def _set_target_coordinates(self, target_df):
        """
        Add a column of SkyCoords to pandas dataframe
        :return: 
        """
        mask = (target_df['typedesig'] == 'f')
        target_df_valid = target_df[mask]

        target_df['SkyCoords'] = False
        target_df.loc[mask, 'SkyCoords'] = SkyCoord(ra=target_df_valid['ra'],
                                                    dec=target_df_valid['dec'],
                                                    unit="deg")

        target_df['FixedObject'] = target_df.apply(self._set_fixed_targets,
                                                   axis=1)
        # target_df['HA'] = target_df.apply(self._set_target_ha, axis=1)
        return target_df

    def _set_non_sidereal_target_coordinates(self, target_df):
        """

        :param target_df:
        :return:
        """
        # mask = target_df[target_df['objname'].str.contains("NONSID")]
        # mask = (target_df['typedesig'] == 'f')

    def _set_fixed_targets(self, row):
        """
        Add a column of SkyCoords to pandas dataframe
        :return: 
        """

        return astroplan.FixedTarget(name=row['objname'],
                                     coord=row['SkyCoords'])

    def _set_target_ha(self, row):
        """
        Add a column of SkyCoords to pandas dataframe
        :return: 
        """

        return self.obs_site.target_hour_angle(self.running_obs_time,
                                               row['FixedObject'])

    def _set_obs_seq(self, target):

        """
        Parse database target scheme

        :param target:
        :return:
        """
        # Prep the variables
        ifu = False
        rc = False
        rc_total = 0
        ifu_total = 0

        rc_filter_list = ['r', 'g', 'i', 'u']
        ifu_exptime = 0

        # 1. First we extract the filter sequence

        seq = list(target['obs_seq'])
        exptime = list(target['exptime'])
        print(seq, exptime)
        # 2. Remove ifu observations first if they exist
        index = [i for i, sq in enumerate(seq) if 'ifu' in sq]

        if index:
            for j in index:
                ifu = seq.pop(j)
                ifu_exptime = int(exptime.pop(j))

                if ifu_exptime == 0:
                    ifu = False
                elif ifu_exptime == 60:
                    ifu_exptime = 1800

        # 3. If the seq list is empty then there is no photmetry follow-up
        # and we should exit

        if not seq:
            ifu_total = ifu_exptime
            obs_seq_dict = {
                'ifu': ifu,
                'ifu_exptime': ifu_exptime,
                'ifu_total': ifu_total+47,
                'rc': rc,
                'rc_obs_dict': None,
                'rc_total': 0,
                'total': abs(ifu_total+47 + (rc_total * target['seq_repeats']))
            }
            return obs_seq_dict

        if ifu:
            ifu_total = ifu_exptime

        # 4. If we are still here then we need to get the photometry sequence
        obs_order_list = []
        obs_exptime_list = []
        obs_repeat_list = []

        for i in range(len(seq)):

            flt = seq[i][-1]
            flt_exptime = int(exptime[i])
            flt_repeat = int(seq[i][:-1])
            # 4a. After parsing the indivual elements we need to check that
            # they are
            # valid values
            if flt in rc_filter_list:

                if 0 < flt_exptime < 600:
                    if 0 < flt_repeat < 100:
                        obs_order_list.append(flt)
                        obs_exptime_list.append(str(flt_exptime))
                        obs_repeat_list.append(str(flt_repeat))
                        rc_total += ((flt_exptime+47) * flt_repeat)
            else:
                continue

        # 5. If everything went well then we should have three non empty list.

        if len(obs_order_list) >= 1:
            rc = True
            obs_dict = {
                'obs_order': ','.join(obs_order_list),
                'obs_exptime': ','.join(obs_exptime_list),
                'obs_repeat_filter': ','.join(obs_repeat_list),
                'obs_repeat_seq': target['seq_repeats']}
        else:
            rc = False
            obs_dict = None

        obs_seq_dict = {
            'ifu': ifu,
            'ifu_exptime': ifu_exptime,
            'ifu_total': ifu_total,
            'rc': rc,
            'rc_obs_dict': obs_dict,
            'rc_total': rc_total * target['seq_repeats'],
            'total': abs(ifu_total + (rc_total * target['seq_repeats']))
        }

        return obs_seq_dict

    def get_query(self, query, return_type=''):
        """
        Get a query return 
        :param return_type: str (will return pd.DataFrame if 'df')
        :param query: str (sql query)
        :return: dict or pd.DataFrame for results of query
        """

        if return_type == 'df':
            return pd.read_sql_query(query, self.ph_db.get_conn_sedmDB())
        else:
            return self.ph_db.execute_sql(query, return_type='dict')

    def load_targets(self, where_statement="",
                     and_statement="AND r.status = 'PENDING'",
                     group_statement="", order_statement="",
                     return_type=""):
        """
        retrieves all targets in the database that march self.req_query
        :param return_type: return_type kwarg of get_query()
        :param where_statement: 
        :param and_statement: 
        :param group_statement: 
        :param order_statement: 
        :return: pd.DataFrame if return_type is 'df', else dict
        """

        if not where_statement:
            enddate = (
                    datetime.datetime.utcnow() +
                    datetime.timedelta(days=1)).strftime("%Y-%m-%d")
            startdate = datetime.datetime.utcnow().strftime("%Y-%m-%d")

            where_statement = "WHERE r.enddate > '%s' AND r.object_id > 100 " \
                              "AND r.inidate <= '%s'" % (enddate, startdate)
            # where_statement = "WHERE r.enddate > '%s'
            # AND r.object_id > 100 " % (enddate)

        if not and_statement:
            and_statement = "AND r.status = 'PENDING'"

        query = self.req_query.substitute(where_statement=where_statement,
                                          and_statement=and_statement,
                                          group_statement=group_statement,
                                          order_statement=order_statement)
        print(query)
        results = self.get_query(query, return_type=return_type)
        print(results.name)
        return results

    def remove_setting_targets(self, target_list, start_time="", end_time="",
                               airmass=(1, 2.5), moon=(20, 180)):
        """
        
        :param target_list:
        :param start_time:
        :param end_time:
        :param airmass:
        :param moon:
        :return: 
        """
        constraint = [astroplan.AirmassConstraint(min=airmass[0],
                                                  max=airmass[1]),

                      astroplan.MoonSeparationConstraint(min=moon[0] * u.degree)
                      ]
        for row in target_list.itertuples():
            print(start_time.iso, end_time.iso)
            tx = astroplan.is_observable(constraint, self.obs_site,
                                         row.FixedObject,
                                         times=[start_time, end_time])
            ty = astroplan.is_event_observable(constraint, self.obs_site,
                                               row.FixedObject,
                                               times=[start_time, end_time])
            print(tx, ty)
            if not tx:
                print(row.objname)
                target_list = target_list[target_list.req_id != row.req_id]

        return target_list

    def simulate_night(self, start_time='', end_time='', do_focus=True,
                       do_standard=True, return_type='html'):
        """

        :param start_time: 
        :param end_time: 
        :param do_focus: 
        :param do_standard:
        :param return_type:
        :return: 
        """

        # 1. Get start and end times
        if not start_time:
            start_time = self.obs_times['evening_nautical']
            if datetime.datetime.utcnow() > start_time:
                start_time = Time(datetime.datetime.utcnow())
        if not end_time:
            end_time = self.obs_times['morning_astronomical']

        self.running_obs_time = start_time

        if return_type == 'html':
            html_str = """<table class='table'><tr><th>Expected Obs Time</th>
                              <th>Object Name</th>
                              <th>Start Date</th>
                              <th>End Date</th>
                              <th>Priority</th>
                              <th>Project ID</th>
                              <th>RA</th>
                              <th>DEC</th>
                              <th>IFU Exptime</th>
                              <th>RC Exptime</th>
                              <th>Total Exptime</th>
                              <th>Update Request</th>
                              </tr>"""
        else:
            html_str = ""

        # 2. Get all targets
        targets = self.load_targets(return_type='df')
        targets = self._set_target_coordinates(targets)

        # Remove non-observable targets
        print(len(targets), "before purge")
        # targets = self.remove_setting_targets(targets, start_time=start_time,
        #                                       end_time=end_time)
        print(len(targets), "after purge")

        targets['HA'] = targets.apply(self._set_target_ha, axis=1)
        targets['obs_dict'] = targets.apply(self._set_obs_seq, axis=1)

        targets = targets.sort_values(['priority', 'HA'], ascending=[False,
                                                                     False])
        print(targets[['priority', 'HA']])

        # 3. Go through all the targets until we fill up the night
        current_time = start_time

        while current_time <= end_time:
            print("Using input datetime of: ", current_time.iso)
            # Include focus time?
            if do_focus:
                current_time += TimeDelta(300, format='sec')
                do_focus = False
            if do_standard:
                current_time += TimeDelta(300, format='sec')
                do_standard = False

            time_remaining = end_time - current_time

            if time_remaining.sec <= 0:

                break

            self.running_obs_time = current_time
            targets['HA'] = targets.apply(self._set_target_ha, axis=1)

            targets = targets.sort_values(['priority', 'HA'], ascending=[False,
                                                                         False])
            # print(targets[['priority','HA']], current_time)
            z = self.get_next_observable_target(targets, obs_time=current_time,
                                                max_time=time_remaining,
                                                return_type=return_type)

            # targets = self.remove_setting_targets(
            # targets, start_time=current_time, end_time=end_time)

            idx, t = z

            if not idx:
                print(len(targets))
                if return_type == 'html':
                    html_str += self.tr_row.substitute(
                        {
                            'allocation': '',
                            'obstime': current_time.iso,
                            'objname': 'Standard [No other object available '
                                       'observation]',
                            'contact': 'SEDm-Calib',
                            'project': 'SEDm-Calib',
                            'ra': 'NA',
                            'dec': 'NA',
                            'ifu_exptime': 180,
                            'rc_seq': 'NA',
                            'rc_exptime': 'NA',
                            'total': 180,
                            'request_id': 'NA'
                        }
                    )
                current_time += TimeDelta(200, format='sec')
            else:
                if return_type == 'html':
                    html_str += t[1]
                    t = t[0]

                targets = targets[targets.req_id != idx]
                # Adding overhead
                current_time += TimeDelta(t['total']+60, format='sec')

        if return_type == 'html':
            html_str += "</table><br>Last Updated:%s UT" % \
                        datetime.datetime.utcnow()
            return html_str

    def get_next_observable_target(self, target_list=None, obs_time=None,
                                   max_time=-1, airmass=(1, 2.7),
                                   moon_sep=(0, 180), ignore_target=None,
                                   return_type=''):
        """

        :return:
            if no observable target:
                False, False
            elif return_type == 'html':
                req_id: float
                (obs_dict, html): tuple(dict, str)
            else:
                req_id: float
                obs_dict: dict
        """

        if max_time:
            pass
        if ignore_target:
            pass
        # If the target_list is empty then all we can do is return back no
        # target and do a standard for the time being

        if not isinstance(target_list, pd.DataFrame) and not target_list:
            print("Making a new target list")
            targets = self.load_targets(return_type='df')
            targets = self._set_target_coordinates(targets)
            targets['obs_dict'] = targets.apply(self._set_obs_seq, axis=1)
            target_list = targets.sort_values(['priority'], ascending=[False])

        if len(target_list) == 0:
            return False, False

        if not obs_time:
            obs_time = Time(datetime.datetime.utcnow())

        for row in target_list.itertuples():
            start = obs_time
            if start < self.obs_times['evening_astronomical']:
                if row.obs_dict['ifu']:
                    continue

            finish = start + TimeDelta(row.obs_dict['total'], format='sec')

            constraint = [
                astroplan.AirmassConstraint(min=airmass[0],
                                            max=airmass[1]),
                astroplan.MoonSeparationConstraint(min=moon_sep[0] * u.degree)
            ]
                          
            # we add an additional phase constraint if the object is periodic
            """if row.typedesig == 'v':
                if ~np.isnan(row['mjd0']):
                    epoch = Time(row['mjd0'], fmt='mjd')
                else:
                    # this is derived from a formula in Sesar et al 2017
                    epoch = Time(2400000-row['phi']*row['phasedays'], fmt='jd')
                    
                periodic_event = astroplan.PeriodicEvent(epoch=epoch,
                                        period=u.day * row['phasedays'])
                constraint.append(astroplan.PhaseConstraint(periodic_event,
                        min=(row['phase'] - row['sampletolerance']) % 1,
                        max=(row['phase'] + row['sampletolerance']) % 1))"""
                            
            if row.typedesig in 'vf':
                if astroplan.is_observable(constraint, self.obs_site,
                                           row.FixedObject,
                                           times=[obs_time, finish]):
                    print(row.objname, row.ra, row.dec)
                    if return_type == 'html':
                        if row.obs_dict['rc']:
                            rc_seq = row.obs_dict['rc_obs_dict']['obs_order'],
                            rc_exptime = row.obs_dict['rc_obs_dict'][
                                             'obs_exptime'],
                        else:
                            rc_seq = 'NA'
                            rc_exptime = 'NA'

                        html = self.tr_row.substitute(
                            {
                                'allocation': row.allocation_id,
                                'obstime': obs_time.iso,
                                'objname': row.objname,
                                'contact': row.priority,
                                'project': row.designator,
                                'ra': row.ra,
                                'dec': row.dec,
                                'ifu_exptime': row.obs_dict['ifu_exptime'],
                                'rc_seq': rc_seq,
                                'rc_exptime': rc_exptime,
                                'total': row.obs_dict['total'],
                                'startdate': row.inidate,
                                'enddate': row.enddate,
                                'request_id': row.req_id
                            }
                        )
                        return row.req_id, (row.obs_dict, html)
                    else:
                        return row.req_id, row.obs_dict

        # If we make it back out of the loop then there was no observable target
        return False, False

    def get_observing_times(self, obsdatetime=None, return_type=''):
        """

        :param return_type: 
        :param obsdatetime: 
        :return: 
        """
        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        obstime = Time(self.obsdatetime)
        sun_set = self.obs_site.sun_set_time(obstime, which="nearest")

        # if sun set is greater than input time we are
        if sun_set > obstime:
            print('Using next')
            which = 'next'
        else:
            print('Using nearest')
            which = 'nearest'

        ret = {
            'sun_set':
                self.obs_site.sun_set_time(obstime, which=which),
            'sun_rise':
                self.obs_site.sun_rise_time(obstime, which=which),
            'evening_civil':
                self.obs_site.twilight_evening_civil(obstime, which=which),
            'evening_nautical':
                self.obs_site.twilight_evening_nautical(obstime, which=which),
            'evening_astronomical':
                self.obs_site.twilight_evening_astronomical(obstime,
                                                            which=which),
            'morning_civil':
                self.obs_site.twilight_morning_civil(obstime, which=which),
            'morning_nautical':
                self.obs_site.twilight_morning_nautical(obstime, which=which),
            'morning_astronomical':
                self.obs_site.twilight_morning_astronomical(obstime,
                                                            which=which),
            'moon_set':
                self.obs_site.moon_set_time(obstime),
        }

        if return_type == 'json':
            json_dict = {k: v.iso.split()[-1] for k, v in ret.items()}
            json_dict['obsdate'] = ret['evening_astronomical'].iso.split()[0]
            json_dict['moon_illumination'] = self.obs_site.moon_illumination(obstime)
            json_dict['moon_rise'] = self.obs_site.moon_rise_time(obstime).iso.split()[-1]
            return json_dict
        else:
            return ret

    def get_lst(self, obsdatetime=None):
        if obsdatetime:
            self.obsdatetime = obsdatetime

        obstime = Time(self.obsdatetime)

        return self.obs_site.local_sidereal_time(obstime)

    def get_sun(self, obsdatetime=None):
        if obsdatetime:
            self.obsdatetime = obsdatetime

        obstime = Time(self.obsdatetime)

        return self.obs_site.sun_altaz(obstime)

    def get_twilight_coords(self, obsdatetime=None, dec=33.33):
        """

        :param obsdatetime: 
        :param dec: 
        :return: 
        """
        if obsdatetime:
            self.obsdatetime = obsdatetime

        # Get sidereal time
        lst = self.get_lst()

        ra = lst.degree - 45

        if ra < 0:
            ra = 360 + ra

        return {'ra': round(ra, 4),
                'dec': dec}

    def get_twilight_exptime(self, obsdatetime=None, camera='rc'):
        """

        :param obsdatetime: 
        :param camera:
        :return: 
        """
        if obsdatetime:
            self.obsdatetime = obsdatetime

        # Get sun angle
        sun_pos = self.get_sun()
        sun_angle = sun_pos.alt.degree

        if -10 >= sun_angle >= -12:
            exptime = 120
        elif -8 >= sun_angle >= -10:
            exptime = 60
        elif -6 >= sun_angle >= -8:
            exptime = 30
        elif -4 >= sun_angle >= -6:
            exptime = 1
        else:
            exptime = 1

        if camera == 'ifu':
            exptime *= 1.5

        return {'exptime': exptime}

    def get_focus_coords(self, obsdatetime=None, dec=33.33):
        """

        :param obsdatetime: 
        :param dec: 
        :return: 
        """
        if obsdatetime:
            self.obsdatetime = obsdatetime

        # Get sidereal time
        lst = self.get_lst()

        ra = lst.degree

        return {'ra': round(ra, 4),
                'dec': dec}

    def get_standard_request_id(self, name="", exptime=180):
        """

        :param name: 
        :param exptime: 
        :return: bool, id
        """
        object_id = self.ph_db.get_object_id_from_name(name)
        for obj in object_id:
            if obj[1].lower() == name.lower():
                object_id = obj[0]
                break

        print(object_id)

        start = datetime.datetime.utcnow()
        end = datetime.timedelta(days=1) + start
        request_dict = {'obs_seq': '{1ifu}',
                        'exptime': '{%s}' % int(exptime),
                        'object_id': object_id,
                        'marshal_id': '-1',
                        'user_id': 2,
                        'allocation_id': '20180131224646741',
                        'priority': '-1',
                        'inidate': start.strftime("%Y-%m-%d"),
                        'enddate': end.strftime("%Y-%m-%d"),
                        'maxairmass': '2.5',
                        'status': 'PENDING',
                        'max_fwhm': '10',
                        'min_moon_dist': '30',
                        'max_moon_illum': '1',
                        'max_cloud_cover': '1',
                        'seq_repeats': '1',
                        'seq_completed': '0'}

        return {'object_id': object_id,
                'request_id': self.ph_db.add_request(request_dict)[0]}

    def filter_offsets(self, ra, dec, offset_filter=None):

        """
        Get a filter offset dictionary
        :param ra:
        :param dec:
        :param offset_filter:
        :return:
        """
        try:
            filter_dict = {
                'r': {'ra': None, 'dec': None},
                'g': {'ra': None, 'dec': None},
                'i': {'ra': None, 'dec': None},
                'u': {'ra': None, 'dec': None},
                'ifu': {'ra': None, 'dec': None}
            }

            if not offset_filter:
                offset_filter = {
                    'r': {'ra': 130, 'dec': 110},
                    'g': {'ra': -325, 'dec': -328},
                    'i': {'ra': -320, 'dec': 93},
                    'u': {'ra': 111, 'dec': -328},
                    'ifu': {'ra': -97, 'dec': -98}
                }

            obj = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
            for k, v in offset_filter.items():
                offra = (Angle(offset_filter[k]['ra'], unit=u.arcsec) /
                         np.cos(obj.dec.to('radian')))
                offdec = Angle(offset_filter[k]['dec'], unit=u.arcsec)

                new_pos = SkyCoord(obj.ra + offra, obj.dec + offdec,
                                   frame='icrs')
                filter_dict[k]['ra'] = round(new_pos.ra.value, 6)
                filter_dict[k]['dec'] = round(new_pos.dec.value, 6)

            return True, filter_dict
        except Exception as e:
            print(str(e))
            return False, str(e)

    def get_calib_request_id(self, camera='ifu', nexp=1, object_id="",
                             exptime=0):
        """

        :param camera: 
        :param nexp:
        :param object_id: 
        :param exptime: 
        :return: 
        """

        if camera == 'ifu':
            pass
        elif camera == 'rc':
            camera = 'r'
        else:
            return {'request_id': ''}

        start = datetime.datetime.utcnow()
        end = datetime.timedelta(days=1) + start
        request_dict = {'obs_seq': '{%s%s}' % (nexp, camera),
                        'exptime': '{%s}' % int(exptime),
                        'object_id': object_id,
                        'marshal_id': '-1',
                        'user_id': 2,
                        'allocation_id': '20180131224646741',
                        'priority': '-1',
                        'inidate': start.strftime("%Y-%m-%d"),
                        'enddate': end.strftime("%Y-%m-%d"),
                        'maxairmass': '2.5',
                        'status': 'PENDING',
                        'max_fwhm': '10',
                        'min_moon_dist': '30',
                        'max_moon_illum': '1',
                        'max_cloud_cover': '1',
                        'seq_repeats': '1',
                        'seq_completed': '0'}

        return {'request_id': self.ph_db.add_request(request_dict)[0]}

    def reset_targets(self, start_date=None, end_date=None,
                      start_status='ACTIVE', end_status='PENDING'):
        """

        :param start_date:
        :param end_date:
        :param start_status:
        :param end_status:
        :return:
        """

        if end_date:
            pass
        if start_status:
            pass
        if end_status:
            pass

        if not start_date:
            start_date = datetime.datetime.utcnow() - datetime.timedelta(days=7)

        query = """UPDATE request SET status = 'PENDING' WHERE 
                   status = 'ACTIVE' AND inidate > '%s'""" % start_date

        ret = self.get_query(query, return_type='')
        print(ret)


if __name__ == "__main__":
    import time

    x = ScheduleNight()
    s = time.time()

    # print(x.get_observing_times(return_type='json'))
    # print(x.get_standard_request_id(name='HZ44', exptime=90))
    r = x.simulate_night()
    x.reset_targets()

    data = open(scheduler_path, 'w')
    data.write(r)
    data.close()

    print(time.time() - s)
