import json
from string import Template
import datetime
import psycopg2
import psycopg2.extras
import pandas as pd
import astroplan
import os
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord

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
        self.conn_str = ("dbname=%(db)s user=%(user)s host=%(host)s "
                         "password=%(password)s port=%(port)s" % (params['database']))

        self.target_frame = None
        self.req_query = Template("""SELECT r.id AS req_id, r.object_id AS obj_id, 
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
                                    u.id, u.email, a.id, a.inidate, a.enddate, a.time_spent, 
                                    a.time_allocated, a.program_id, a.active, 
                                    p.designator, p.name, p.group_id, p.pi,
                                    p.time_allocated, r.priority, p.inidate,
                                    p.enddate 
                                    FROM "public".request r
                                    INNER JOIN "public"."object" o ON ( r.object_id = o.id  )  
                                    INNER JOIN "public".users u ON ( r.user_id = u.id  )  
                                    INNER JOIN "public".allocation a ON ( r.allocation_id = a.id  )  
                                    INNER JOIN "public".program p ON ( a.program_id = p.id  )
                                    ${where_statement} 
                                    ${and_statement}
                                    ${group_statement}
                                    ${order_statement}
                                    """)

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

        target_df['FixedObject'] = target_df.apply(self._set_fixed_targets, axis=1)
        target_df['HA'] = target_df.apply(self._set_target_ha, axis=1)
        return target_df

    def _set_fixed_targets(self, row):
        """
        Add a column of SkyCoords to pandas dataframe
        :return: 
        """

        return astroplan.FixedTarget(name=row['objname'], coord=row['SkyCoords'])


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
        ifu_exptime = 1800

        # 1. First we extract the filter sequence
        seq = list(target['obs_seq'])
        exptime = list(target['exptime'])

        # 2. Remove ifu observations first if they exist
        index = [i for i, s in enumerate(seq) if 'ifu' in s]

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
                'ifu_total': ifu_total,
                'rc': rc,
                'rc_obs_dict': None,
                'rc_total': 0,
                'total': ifu_total + (rc_total * target['seq_repeats'])
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
            flt_exptime =int(exptime[i])
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
                        rc_total += (flt_exptime * flt_repeat)
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
            'rc_total': rc_total*target['seq_repeats'],
            'total': ifu_total + (rc_total*target['seq_repeats'])
        }

        return obs_seq_dict

    def db(self):
        return psycopg2.connect(self.conn_str)



    def get_query(self, query, return_type=''):
        """
        Get a query return 
        :param return_type: 
        :param query: 
        :return: 
        """

        if not self.conn:
            self.conn = self.db()

        if return_type == 'df':
            df = pd.read_sql_query(query, self.conn)
            return df
        else:
            self.cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            self.cursor.execute(query)
            return self.cursor.fetchall()

    def load_targets(self, where_statement="", and_statement="",
                     group_statement="", order_statement="",
                     return_type=""):
        """
        
        :param return_type: 
        :param where_statement: 
        :param and_statement: 
        :param group_statement: 
        :param order_statement: 
        :return: 
        """

        if not where_statement:
            enddate = datetime.datetime.utcnow()
            where_statement = "WHERE r.enddate > '%s' " % enddate

        if not and_statement:
            and_statement = "AND r.status = 'PENDING'"

        query = self.req_query.substitute(where_statement=where_statement,
                                          and_statement=and_statement,
                                          group_statement=group_statement,
                                          order_statement=order_statement)

        results = self.get_query(query, return_type=return_type)

        return results

    def simulate_night(self, start_time='', end_time='', do_focus=True, do_standard=True,
                       block_list='P,P,P,C', add_columns=True, save_as='', return_type=''):
        """
        
        :param start_time: 
        :param end_time: 
        :param do_focus: 
        :param do_standard: 
        :param block_list: 
        :param add_columns: 
        :param save_as: 
        :return: 
        """

        # 1. Get start and end times
        if not start_time:
            start_time = self.obs_times['evening_nautical']

        if not end_time:
            end_time = self.obs_times['morning_nautical']

        self.running_obs_time = start_time

        # 2. Get all targets
        targets = self.load_targets(return_type='df')
        targets = self._set_target_coordinates(targets)
        targets['obs_dict'] = targets.apply(self._set_obs_seq, axis=1)
        print(targets['allocation_id'])
        targets = targets.sort_values(['priority', 'ra'], ascending=[False, True])

        # 3. Go through all the targets until we fill up the night
        current_time = start_time
        count = 0
        obslist = []

        block = 0
        if isinstance(block_list, str):
            block_list = block_list.split(",")

        while current_time <= end_time and count < 1000:

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

            # Get the current observing block
            if block > len(block_list):
                block = 0

            observing_block = block_list[block]

            idx, t = self.get_next_observable_target(targets, obs_time=current_time,
                                                     max_time=time_remaining)

            if idx:
                targets = targets[targets.req_id != idx]
                obslist.append({**t, **{'obs_time': current_time.datetime}})
            elif time_remaining.sec >= 300:
                obslist.append({'name': "STANDARD OBS", "ifu_total": 300,
                                'priority': -1, "rc_total": 0,
                                "obs_time": current_time.datetime})
                t = {'total': 300}
            else:
                break
            targets = self._set_target_coordinates(targets)
            current_time += TimeDelta(t['total'], format='sec')

            count += 1

        if return_type == 'table':
            table_str = """<table class="table table-striped"><tr>
            <th>Projected Time</th><th>Priority</th><th>Name</th>
            <th>IFU Exptime</th><th>RC Exptime</th></tr>"""
            for i in obslist:
                table_str += ('<tr><td>%s</td><td>%s</td><td>%s</td>'
                              '<td>%s</td><td>%s</td></tr>' %
                              (i['obs_time'], i['priority'],
                               i['name'], i['ifu_total'], i['rc_total']))
            table_str += "</table>"

        return table_str


    def get_next_observable_target(self, target_list, obs_time, max_time=-1,
                                   airmass=(1, 2.2),  moon_sep=(30, 180),
                                   block_type=''):
        """
        
        :return: 
        """

        # If the target_list is empty then all we can do is return back no
        # target and do a standard for the time being
        next_target = False

        if target_list.empty:
            return next_target, ""

        for row in target_list.itertuples():

            constraint = astroplan.AirmassConstraint(min=airmass[0],
                                                     max=airmass[1])

            if row.typedesig == 'f':
                if astroplan.is_observable(constraint, self.obs_site,
                                           row.FixedObject,
                                           times=[obs_time]):

                    return_dict = {'name': row.objname,
                                   'priority': row.priority}

                    return row.req_id, {**row.obs_dict, **return_dict}

        return False, False


    def get_observing_times(self, obsdatetime=None):
        """
        
        :param obsdatetime: 
        :return: 
        """
        if obsdatetime:
            self.obsdatetime = obsdatetime

        obstime = Time(self.obsdatetime)
        sun_set = self.obs_site.sun_set_time(obstime, which="nearest")

        # if sun set is greater than input time we are
        if sun_set > obstime:
            print('Using next')
            which = 'next'
        else:
            print('Using nearest')
            which = 'nearest'

        return {'sun_set':
                     self.obs_site.sun_set_time(obstime,
                                                which=which),
                     'sun_rise':
                         self.obs_site.sun_rise_time(obstime,
                                                     which=which),
                     'evening_civil':
                         self.obs_site.twilight_evening_civil(obstime,
                                                              which=which),
                     'evening_nautical':
                         self.obs_site.twilight_evening_nautical(obstime,
                                                                 which=which),
                     'evening_astronomical':
                         self.obs_site.twilight_evening_astronomical(obstime,
                                                                     which=which),
                     'morning_civil':
                         self.obs_site.twilight_morning_civil(obstime,
                                                              which=which),
                     'morning_nautical':
                         self.obs_site.twilight_morning_nautical(obstime,
                                                                 which=which),
                     'morning_astronomical':
                         self.obs_site.twilight_morning_astronomical(obstime,
                                                                     which=which)}




if __name__ == "__main__":
    import time

    x = ScheduleNight()
    s = time.time()
    print(x.obs_times)
    print(x.simulate_night(return_type='table'))
    print(time.time()-s)