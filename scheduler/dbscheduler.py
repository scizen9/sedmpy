import json
from string import Template
import datetime
import pandas as pd
import astroplan
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
import astropy.units as u
import os
import psycopg2.extras
import psycopg2
import time

import obstimes
import sqlite3

import sedmpy_version

# Get scheduler configuration
# Find config file: default is sedmpy/config/schedulerconfig.json
try:
    configfile = os.environ["SEDMSCHEDULERCONFIG"]
except KeyError:
    configfile = os.path.join(sedmpy_version.CONFIG_DIR, "schedulerconfig.json")
with open(configfile) as config_file:
    db_cfg = json.load(config_file)

_dbconn = {"host": db_cfg["dbconn"]["host"],
           "port": db_cfg["dbconn"]["port"],
           "dbname": db_cfg["dbconn"]["dbname"],
           "user": db_cfg["dbconn"]["user"],
           "password": db_cfg["dbconn"]["password"]
           }

_standard_db = db_cfg["paths"]["standard_db"]
_target_dir = db_cfg["paths"]["target_dir"]
_schedule = db_cfg["paths"]["schedule"]


class Scheduler:
    def __init__(self, site_name='Palomar', obsdatetime=None,
                 save_as="targets.json"):

        self.scheduler_config_file = configfile
        # Open the file with the correct encoding

        self.standard_db = _standard_db
        self.target_dir = _target_dir
        self.standard_dict = {}
        self.standard_star_list = []

        self.site_name = site_name
        self.times = obstimes.ScheduleNight()
        self.obs_times = self.times.get_observing_times_by_date()
        self.running_obs_time = None
        self.site = EarthLocation.of_site(self.site_name)
        self.obs_site_plan = astroplan.Observer.at_site(
            site_name=self.site_name)
        self.obsdatetime = obsdatetime
        self.save_as = save_as
        self.dbconn = psycopg2.connect(**_dbconn)
        # self.ph_db = sedmpy_import.dbconnect()
        self.query = Template(
            "SELECT r.id AS req_id, r.object_id AS obj_id, \n"
            "r.user_id, r.marshal_id, r.exptime, r.maxairmass,\n"
            "r.max_fwhm, r.min_moon_dist, r.max_moon_illum, \n"
            "r.max_cloud_cover, r.status, \n"
            "r.priority AS reqpriority, r.inidate, r.enddate,\n"
            "r.cadence, r.phasesamples, r.sampletolerance, \n"
            "r.filters, r.nexposures, r.obs_seq, r.seq_repeats,\n"
            "r.seq_completed, r.last_obs_jd, r.creationdate,\n"
            "r.lastmodified, r.allocation_id, r.marshal_id, \n"
            "o.id AS obj_id, o.name AS objname, o.iauname, o.ra, o.\"dec\",\n"
            "o.typedesig, o.epoch, o.magnitude, o.creationdate, \n"
            "u.id, u.email, a.id AS allocation_id, \n"
            "a.inidate, a.enddate, a.time_spent, a.designator as p60prid, \n"
            "a.time_allocated, a.program_id, a.active, \n"
            "p.designator, p.name, p.group_id, p.pi,\n"
            "p.time_allocated, r.priority, p.inidate,\n"
            "p.enddate, pe.mjd0, pe.phasedays, pe.phi,\n"
            "r.phase, r.sampletolerance\n"
            "FROM \"public\".request r\n"
            "INNER JOIN \"public\".\"object\" o ON (r.object_id = o.id)\n"
            "INNER JOIN \"public\".users u ON (r.user_id = u.id)\n"
            "INNER JOIN \"public\".allocation a ON (r.allocation_id = a.id)\n"
            "INNER JOIN \"public\".program p ON (a.program_id = p.id)\n"
            "LEFT JOIN \"public\".periodic pe on (pe.object_id=o.id)\n"
            "${where_statement}\n"
            "${and_statement}\n"
            "${group_statement}\n"
            "${order_statement}"
        )

        self.tr_row = Template(
            """<tr id="${allocation}">
               <td>${obstime}</td>
               <td>${objname}</td> 
               <td>${priority}</td> 
               <td>${project}</td> 
               <td>${ra}</td> 
               <td>${dec}</td> 
               <td>${start_ha}</td> 
               <td>${end_ha}</td> 
               <td>${ifu_exptime}</td> 
               <td>Filters:${rc_seq}<br>Exptime:${rc_exptime}</td> 
               <td>${total}</td> 
               <td><a href='request?request_id=${request_id}'>+</a></td> 
               <td>${rejects}</td> 
               </tr>"""
        )

    def __load_standards_from_db(self):
        """
        Open the sqlite database of targets

        :return:
        """
        print(self.standard_db)
        conn = sqlite3.connect(self.standard_db)
        cur = conn.cursor()

        results = cur.execute("SELECT * FROM standards")
        standards = results.fetchall()
        for s in standards:
            name, ra, dec, exptime = s[0].rstrip().encode('utf8'), s[3], s[4],\
                                     s[5]

            if name.upper() == 'LB227':
                continue
            coords = SkyCoord(ra=ra, dec=dec, unit='deg')
            obj = astroplan.FixedTarget(name=name, coord=coords)
            self.standard_star_list.append(obj)
            self.standard_dict[name] = {
                'name': name,
                'ra': ra,
                'dec': dec,
                'exptime': exptime
            }

    def get_standard(self, name=None, obsdate=None):
        """
        If the name is not given find the closest standard star to zenith
        :param name:
        :param obsdate:
        :return:
        """
        print("Choosing a standard")
        start = time.time()
        self.__load_standards_from_db()
        if not obsdate:
            obsdate = datetime.datetime.utcnow()

        if not name:
            name = 'zenith'
        print("Finding standard: %s" % name)
        if name.lower() == 'zenith':
            sairmass = 100
            for standard in self.standard_star_list:
                print(standard)
                airmass = self.obs_site_plan.altaz(obsdate, standard).secz

                if sairmass > airmass > 0:
                    target = standard
                    sairmass = airmass
                    name = target.name
        std = self.standard_dict[name]
        std['name'] = std['name'].decode('utf-8')
        return {'elaptime': time.time() - start,
                'data': std}

    def _set_obs_seq(self, row):
        """
        Parse database target scheme

        :param row:
        :return:
        """

        obs_seq_list = row['obs_seq']
        exp_time_list = row['exptime']
        repeat = 1

        # Prep the variables
        ifu = False
        rc = False
        rc_total = 0
        ifu_total = 0

        rc_filter_list = ['r', 'g', 'i', 'u']
        ifu_exptime = 0

        # 1. First we extract the filter sequence

        seq = list(obs_seq_list)
        exptime = list(exp_time_list)

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
                'ifu_total': ifu_total + 47,
                'rc': rc,
                'rc_obs_dict': None,
                'rc_total': 0,
                'total': abs(ifu_total + 47 + (rc_total * repeat))
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
            # 4a. After parsing the individual elements we need to check that
            # they are
            # valid values
            if flt in rc_filter_list:

                if 0 < flt_exptime < 600:
                    if 0 < flt_repeat < 100:
                        obs_order_list.append(flt)
                        obs_exptime_list.append(str(flt_exptime))
                        obs_repeat_list.append(str(flt_repeat))
                        rc_total += ((flt_exptime + 47) * flt_repeat)
            else:
                continue

        # 5. If everything went well then we should have three non empty list.

        if len(obs_order_list) >= 1:
            rc = True
            obs_dict = {
                'obs_order': ','.join(obs_order_list),
                'obs_exptime': ','.join(obs_exptime_list),
                'obs_repeat_filter': ','.join(obs_repeat_list),
                'obs_repeat_seq': repeat}
        else:
            rc = False
            obs_dict = None

        obs_seq_dict = {
            'ifu': ifu,
            'ifu_exptime': ifu_exptime,
            'ifu_total': ifu_total,
            'rc': rc,
            'rc_obs_dict': obs_dict,
            'rc_total': rc_total * repeat,
            'total': abs(ifu_total + (rc_total * repeat))
        }

        return obs_seq_dict

    def _set_fixed_targets(self, row):
        """
        Add a column of SkyCoords to pandas dataframe
        :return:
        """

        return astroplan.FixedTarget(name=row['objname'],
                                     coord=row['SkyCoords'])

    def _set_end_time(self, row):
        return row['start_obs'] + TimeDelta(row['obs_seq']['total'],
                                            format='sec')

    def _set_start_altaz(self, row):
        return row['SkyCoords'].transform_to(AltAz(obstime=row['start_obs'],
                                                   location=self.site)).alt

    def _set_end_altaz(self, row):
        return row['SkyCoords'].transform_to(AltAz(obstime=row['end_obs'],
                                                   location=self.site)).alt

    def _set_start_ha(self, row):
        return self.obs_site_plan.target_hour_angle(row['start_obs'],
                                                    row['fixed_object'])

    def _set_end_ha(self, row):
        return self.obs_site_plan.target_hour_angle(row['end_obs'],
                                                    row['fixed_object'])

    def _set_rise_time(self, row):
        return self.obs_site_plan.target_rise_time(row['start_obs'],
                                                   row['fixed_object'],
                                                   horizon=15 * u.degree,
                                                   which="next")

    def _set_set_time(self, row):
        return self.obs_site_plan.target_set_time(row['start_obs'],
                                                  row['fixed_object'],
                                                  horizon=15 * u.degree,
                                                  which="next")

    def _convert_row_to_json(self, row):    # , fields=('name', 'ra', 'dec',
                                            #     'obj_id', 'req_id',
                                            #     'email', 'objname', 'pi')):
        """

        :param row:
        :param fields:
        :return:
        """

        return dict(name=row.objname, p60prnm=row.name,
                    p60prid=row.p60prid, p60prpi=row.pi,
                    ra=row.ra, dec=row.dec, equinox=row.epoch,
                    req_id=row.req_id, obj_id=row.obj_id,
                    obs_dict=row.obs_seq, marshal_id=row.marshal_id)

    def simulate_night(self, start_time='', end_time='', do_focus=True,
                       do_standard=True, target_list=None,
                       get_current_observation=True,
                       return_type='html',
                       sort_columns=('priority', 'start_alt'),
                       sort_order=(False, False), ):
        """

        :param start_time:
        :param end_time:
        :param do_focus:
        :param do_standard:
        :param target_list:
        :param get_current_observation:
        :param return_type:
        :param sort_columns:
        :param sort_order:
        :return:
        """

        start = time.time()

        # 1. Get start and end times
        if not start_time:
            start_time = self.obs_times['evening_nautical']
            if datetime.datetime.utcnow() > start_time:
                start_time = Time(datetime.datetime.utcnow())
        if not end_time:
            end_time = self.obs_times['morning_astronomical']

        self.running_obs_time = start_time

        if not isinstance(target_list, pd.DataFrame) and not target_list:
            print("Making a new target list")
            ret = self.get_active_targets()
            if 'data' in ret:
                targets = ret['data']
            else:
                return ret

            target_list = self.initialize_targets(targets)['data']

        if len(target_list) == 0:
            return {'data': False, 'elaptime': time.time() - start}

        obsdatetime = self.running_obs_time

        if return_type == 'html':
            html_str = """<table class='table'><tr><th>Expected Obs Time</th>
                              <th>Object Name</th>
                              <th>Priority</th>
                              <th>Project ID</th>
                              <th>RA</th>
                              <th>DEC</th>
                              <th>Start HA</th>
                              <th>End HA</th>
                              <th>IFU Exptime</th>
                              <th>RC Exptime</th>
                              <th>Total Exptime</th>
                              <th>Update Request</th>
                              <th>Priority 4+ reject<br>reasons</th>
                              </tr>"""
        else:
            html_str = ""

        # 2. Get all targets
        targets = target_list

        # 3. Go through all the targets until we fill up the night
        current_time = start_time

        while current_time <= end_time:
            targets = self.update_targets_coords(targets, current_time)['data']
            targets = targets.sort_values(list(sort_columns),
                                          ascending=list(sort_order))
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

            z = self.get_next_observable_target(targets,
                                                obsdatetime=current_time,
                                                update_coords=False,
                                                return_type=return_type,
                                                do_sort=False)

            # targets = self.remove_setting_targets(targets,
            #                                       start_time=current_time,
            #                                       end_time=end_time)

            idx, t = z

            if not idx:
                print(len(targets))
                if return_type == 'html':
                    html_str += self.tr_row.substitute(
                        {
                            'allocation': "",
                            'obstime': current_time.iso,
                            'objname': "Standard",
                            'priority': "",
                            'project': "Calib",
                            'ra': "",
                            'dec': "",
                            'start_ha': "",
                            'end_ha': "",
                            'ifu_exptime': 300,
                            'rc_seq': "",
                            'rc_exptime': "",
                            'total': 300,
                            'request_id': "NA",
                            'rejects': ""
                        }
                    )
                current_time += TimeDelta(300, format='sec')
            else:
                if return_type == 'html':
                    html_str += t[1]
                    t = t[0]

                targets = targets[targets.req_id != idx]
                # Adding overhead
                current_time += TimeDelta(t['total']*1.2 + 120, format='sec')

        if return_type == 'html':
            html_str += "</table><br>Last Updated:%s UT" \
                        % datetime.datetime.utcnow()
            return html_str

    def get_active_targets(self, startdate=None, enddate=None,
                           where_statement="", and_statement="",
                           group_statement="", order_statement="",
                           save_copy=True):

        start = time.time()

        if not startdate:
            if datetime.datetime.utcnow().hour >= 14:
                startdate = (datetime.datetime.utcnow() +
                             datetime.timedelta(days=1))
            else:
                startdate = datetime.datetime.utcnow()

            startdate = startdate.replace(hour=23, minute=59, second=59)

        if not enddate:
            enddate = (datetime.datetime.utcnow() +
                       datetime.timedelta(days=1)).strftime("%Y-%m-%d")

        if not where_statement:
            where_statement = ("WHERE r.enddate >= '%s' AND r.object_id > 100 "
                               "AND r.inidate <= '%s'" % (enddate, startdate))

        if not and_statement:
            and_statement = "AND r.status = 'PENDING'"

        q = self.query.substitute(where_statement=where_statement,
                                  and_statement=and_statement,
                                  group_statement=group_statement,
                                  order_statement=order_statement)

        df = pd.read_sql_query(q, self.dbconn)

        if save_copy:
            df.to_csv(self.save_as)

        return {"data": df, "elaptime": time.time() - start}

    def initialize_targets(self, target_df, obstime=''):
        start = time.time()
        mask = (target_df['typedesig'] == 'f')
        target_df_valid = target_df[mask]

        target_df['SkyCoords'] = False
        target_df.loc[mask, 'SkyCoords'] = SkyCoord(ra=target_df_valid['ra'],
                                                    dec=target_df_valid['dec'],
                                                    unit="deg")
        target_df['start_obs'] = False
        if not obstime:
            obstime = datetime.datetime.utcnow()

        target_df.loc[mask, 'start_obs'] = Time(obstime)

        target_df['obs_seq'] = target_df.apply(self._set_obs_seq, axis=1)
        target_df['end_obs'] = target_df.apply(self._set_end_time, axis=1)
        target_df['start_alt'] = target_df.apply(self._set_start_altaz,
                                                 axis=1)
        target_df['end_alt'] = target_df.apply(self._set_end_altaz, axis=1)
        target_df['fixed_object'] = target_df.apply(self._set_fixed_targets,
                                                    axis=1)
        target_df['start_ha'] = target_df.apply(self._set_start_ha, axis=1)
        target_df['end_ha'] = target_df.apply(self._set_end_ha, axis=1)
        target_df['rise_time'] = target_df.apply(self._set_rise_time, axis=1)
        target_df['set_time'] = target_df.apply(self._set_set_time, axis=1)

        return {'data': target_df, 'elaptime': time.time() - start}

    def update_targets_coords(self, df, obstime=None):
        start = time.time()
        df['start_obs'] = False

        if not obstime:
            obstime = datetime.datetime.utcnow()

        df['start_obs'] = Time(obstime)
        df['end_obs'] = df.apply(self._set_end_time, axis=1)
        df['start_alt'] = df.apply(self._set_start_altaz, axis=1)
        df['end_alt'] = df.apply(self._set_end_altaz, axis=1)
        df['start_ha'] = df.apply(self._set_start_ha, axis=1)
        df['end_ha'] = df.apply(self._set_end_ha, axis=1)
        return {'data': df, 'elaptime': time.time() - start}

    def look_for_new_targets(self, df, startdate=None, enddate=None,
                             where_statement="", and_statement="",
                             group_statement="", order_statement="",
                             field='req_id'):

        start = time.time()

        ret = self.get_active_targets(startdate=startdate,
                                      enddate=enddate,
                                      where_statement=where_statement,
                                      and_statement=and_statement,
                                      group_statement=group_statement,
                                      order_statement=order_statement)

        if 'data' in ret:
            new_df = ret['data']
        else:
            return ret

        new_targets = (list(set(new_df[field]) - set(df[field])))

        dropped_targets = (list(set(df[field]) - set(new_df[field])))

        if len(new_targets) >= 1:
            new_df = new_df[new_df['req_id'].isin(new_targets)]
            ret = self.initialize_targets(new_df)

            if 'data' in ret:
                df = df.append(ret['data'])

        if len(dropped_targets) >= 1:
            df = df[-df["req_id"].isin(dropped_targets)]

        return {'data': df, 'elaptime': time.time() - start}

    def get_next_observable_target(self, target_list=None, obsdatetime=None,
                                   airmass=(1, 3.0), moon_sep=(30, 180),
                                   altitude_min=15, ha=(18.75, 5.75),
                                   return_type='', do_airmass=True,
                                   do_sort=True, do_moon_sep=True,
                                   sort_columns=('priority', 'set_time'),
                                   sort_order=(False, False), save=False,
                                   save_as='',
                                   check_end_of_night=True, update_coords=True):
        """

        :return:
        """

        s = time.time()
        # If the target_list is empty then all we can do is return back no
        # target and do a standard for the time being
        next_target = False

        if not isinstance(target_list, pd.DataFrame) and not target_list:
            print("Making a new target list")
            ret = self.get_active_targets()
            if 'data' in ret:
                targets = ret['data']
            else:
                return ret

            target_list = self.initialize_targets(targets)['data']

        if len(target_list) == 0:
            return {'data': False, 'elaptime': time.time() - s}

        if not obsdatetime:
            obsdatetime = Time(datetime.datetime.utcnow())
        else:
            obsdatetime = Time(obsdatetime)

        print("Using time:", obsdatetime.iso)

        if update_coords:
            target_list = self.update_targets_coords(target_list,
                                                     obsdatetime)['data']

        # Remove targets outside of HA range
        if do_sort:
            target_list = target_list.sort_values(list(sort_columns),
                                                  ascending=list(sort_order))
        rej_html = ""
        target_reorder = False
        print(target_list['typedesig'])
        for row in target_list.itertuples():
            start = obsdatetime

            finish = start + TimeDelta(row.obs_seq['total'],
                                       format='sec')
            if row.priority <= 2 and not target_reorder:
                print(target_list.keys())
                target_list = target_list.sort_values('start_ha',
                                                      ascending=False)
                target_reorder = True
                continue

            # Force altitude constraint
            constraint = [astroplan.AltitudeConstraint(
                min=altitude_min * u.deg)]

            if do_airmass:
                constraint.append(astroplan.AirmassConstraint(min=airmass[0],
                                                              max=airmass[1]))

            if do_moon_sep:
                constraint.append(astroplan.MoonSeparationConstraint(
                    min=moon_sep[0] * u.degree))

            if row.typedesig == 'f':
                print(row.objname)
                if astroplan.is_observable(constraint, self.obs_site_plan,
                                           row.fixed_object,
                                           times=[start, finish],
                                           time_grid_resolution=0.1 * u.hour):
                    print(type(row.start_ha))
                    s_ha = float(row.start_ha.to_string(unit=u.hour,
                                                        decimal=True))
                    e_ha = float(row.end_ha.to_string(unit=u.hour,
                                                      decimal=True))

                    if 18.75 > s_ha > 5.75:
                        continue
                    if 18.75 > e_ha > 5.75:
                        continue

                    print(row.objname, row.priority, row.ra, row.dec,
                          row.start_ha, row.end_ha, row.start_obs)
                    if return_type == 'html':
                        if row.obs_seq['rc']:
                            rc_seq = row.obs_seq['rc_obs_dict']['obs_order'],
                            rc_exptime = row.obs_seq[
                                             'rc_obs_dict']['obs_exptime'],
                        else:
                            rc_seq = 'NA'
                            rc_exptime = 'NA'

                        html = self.tr_row.substitute(
                            {
                                'allocation': row.allocation_id,
                                'obstime': start.iso,
                                'objname': row.objname,
                                'priority': row.priority,
                                'project': row.designator,
                                'ra': row.ra,
                                'dec': row.dec,
                                'start_ha': row.start_ha,
                                'end_ha': row.end_ha,
                                'ifu_exptime': row.obs_seq['ifu_exptime'],
                                'rc_seq': rc_seq,
                                'rc_exptime': rc_exptime,
                                'total': row.obs_seq['total'],
                                'request_id': row.req_id,
                                'rejects': rej_html
                            }
                        )
                        return row.req_id, (row.obs_seq, html)
                    elif return_type == 'json':
                        print(rej_html)
                        print("No Regj")
                        targ = self._convert_row_to_json(row)

                        if save:
                            if not save_as:
                                save_as = os.path.join(
                                    self.target_dir,
                                    "next_target_%s.json" %
                                    datetime.datetime.utcnow().strftime(
                                        "%Y%m%d"))

                                with open(save_as, 'w') as outfile:
                                    outfile.write(json.dumps(targ))

                        return {"elaptime": time.time() - s, "data": targ}
                    else:
                        return row.req_id, row.obs_seq
                else:
                    print("I am here")
                    if row.priority >= 4:
                        count = 1
                        num = []
                        for i in constraint:
                            ret = astroplan.is_observable(
                                [i], self.obs_site_plan, row.fixed_object,
                                times=[start, finish],
                                time_grid_resolution=0.1 * u.hour)
                            print(ret, i)
                            if ret:
                                num.append(str(count))
                            count += 1
                        if return_type == 'html' and len(num) >= 1:
                            rej_html += """%s: %s<br>""" % (row.objname, ','.join(num))
                        elif return_type == 'json' and len(num) >= 1:
                            rej_html += ','.join(num)

        if return_type == 'json':
            return {"elaptime": time.time() - s, "error": "No targets found"}
        return False, False

    def get_lst(self, obsdatetime=None):
        start = time.time()
        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        obstime = Time(self.obsdatetime)
        lst = self.obs_site_plan.local_sidereal_time(obstime)
        return {"elaptime": time.time() - start,
                "data": lst}

    def get_sun(self, obsdatetime=None):
        start = time.time()
        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        obstime = Time(self.obsdatetime)

        sun = self.obs_site_plan.sun_altaz(obstime)

        return {"elaptime": time.time() - start,
                "data": sun}

    def get_twilight_coords(self, obsdatetime=None, dec=33.33):
        """

        :param obsdatetime:
        :param dec:
        :return:
        """
        start = time.time()
        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        # Get sidereal time
        lst = self.get_lst(self.obsdatetime)

        ra = lst['data'].degree

        return {'elaptime': time.time() - start,
                'data': {'ra': round(ra, 4),
                         'dec': dec}
                }

    def get_twilight_exptime(self, obsdatetime=None, camera='rc'):
        """

        :param obsdatetime:
        :param camera:
        :return:
        """
        start = time.time()

        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        # Get sun angle
        sun_pos = self.get_sun(self.obsdatetime)
        sun_angle = sun_pos['data'].alt.degree

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

        return {'elaptime': time.time() - start, 'data': {'exptime': exptime}}

    def get_focus_coords(self, obsdatetime=None, dec=33.33):
        """

        :param obsdatetime:
        :param dec:
        :return:
        """
        start = time.time()

        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        # Get sidereal time
        lst = self.get_lst(self.obsdatetime)

        ra = lst['data'].degree

        return {'elaptime': time.time() - start, 'data': {'ra': round(ra, 4),
                                                          'dec': dec}}

    def get_standard_request_id(self, name="", exptime=180):
        """

        :param name:
        :param exptime:
        :return: bool, id
        """
        start = time.time()

        object_id = self.ph_db.get_object_id_from_name(name)
        for obj in object_id:
            if obj[1].lower() == name.lower():
                object_id = obj[0]
                break

        start_date = datetime.datetime.utcnow()
        end_date = start_date + datetime.timedelta(days=1)
        request_dict = {'obs_seq': '{1ifu}',
                        'exptime': '{%s}' % int(exptime),
                        'object_id': object_id,
                        'marshal_id': '-1',
                        'user_id': 2,
                        'allocation_id': '20211011220000019',
                        'priority': '-1',
                        'inidate': start_date.strftime("%Y-%m-%d"),
                        'enddate': end_date.strftime("%Y-%m-%d"),
                        'maxairmass': '2.5',
                        'status': 'PENDING',
                        'max_fwhm': '10',
                        'min_moon_dist': '30',
                        'max_moon_illum': '1',
                        'max_cloud_cover': '1',
                        'seq_repeats': '1',
                        'seq_completed': '0'}
        request_id = self.ph_db.add_request(request_dict)[0]
        return {'elaptime': time.time() - start,
                'data': {'object_id': object_id,
                         'request_id': request_id}}

    def get_calib_request_id(self, camera='ifu', N=1, object_id="", exptime=0):
        """

        :param camera:
        :param N:
        :param object_id:
        :param exptime:
        :return:
        """
        start = time.time()
        print(camera, "right here at them moment")
        if camera == 'ifu':
            pass
        elif camera == 'rc':
            camera = 'r'
        else:
            return {'request_id': ''}

        start_date = datetime.datetime.utcnow()
        end_date = start_date + datetime.timedelta(days=1)
        request_dict = {'obs_seq': '{%s%s}' % (N, camera),
                        'exptime': '{%s}' % int(exptime),
                        'object_id': object_id,
                        'marshal_id': '-1',
                        'user_id': 2,
                        'allocation_id': '20211011220000019',
                        'priority': '-1',
                        'inidate': start_date.strftime("%Y-%m-%d"),
                        'enddate': end_date.strftime("%Y-%m-%d"),
                        'maxairmass': '2.5',
                        'status': 'PENDING',
                        'max_fwhm': '10',
                        'min_moon_dist': '30',
                        'max_moon_illum': '1',
                        'max_cloud_cover': '1',
                        'seq_repeats': '1',
                        'seq_completed': '0'}

        ret_id = self.ph_db.add_request(request_dict)[0]

        return {'elaptime': time.time() - start, 'data': ret_id}

    def update_request(self, request_id, status="PENDING"):
        """

        :param request_id:
        :param status:
        :return:
        """
        start = time.time()
        ret = self.ph_db.update_request({'id': request_id,
                                         'status': status})

        print(ret)

        return {'elaptime': time.time()-start, 'data': ret['data']}


if __name__ == "__main__":

    # s = time.time()
    planner = Scheduler()

    schedule = planner.simulate_night(do_focus=True, do_standard=True)
    data = open(_schedule, 'w')
    data.write(schedule)
    data.close()

    data = open("/scr2/sedm/scheduler/test_scheduler.%s.html"
                % datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S"), 'w')
    data.write(schedule)
    data.close()
