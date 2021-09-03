import datetime
import astroplan
from astropy.time import Time, TimeDelta


class ScheduleNight:
    def __init__(self, obsdatetime=None, config_file="json_url"):
        # 1. Load config file
        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        self.site = 'palomar'
        #self.long = params['site']['longitude']
        #self.lat = params['site']['latitude']
        #self.elev = params['site']['elevation']
        self.obs_site = astroplan.Observer.at_site(site_name=self.site)
        #self.obs_times = self.get_observing_times()

    def get_observing_times_by_date(self, obsdatetime="", return_type=""):
        """

        :param return_type:
        :param obsdatetime:
        :return:
        """

        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        if isinstance(self.obsdatetime, datetime.datetime):
            print(self.obsdatetime, 'time')

            if datetime.datetime.utcnow().hour > 14:
                self.obsdatetime += datetime.timedelta(days=1)
                self.obsdatetime = self.obsdatetime.replace(hour=7)
            else:
                self.obsdatetime = self.obsdatetime.replace(hour=7)

        obstime = Time(self.obsdatetime)

        print(obstime.iso)

        which = 'nearest'

        ret = {'sun_set':
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

        if return_type == 'json':
            json_dict = {k: v.iso for k, v in ret.items()}
            return json_dict
        else:
            return ret

    def get_observing_times(self, obsdatetime=None, return_type='',
                            use_ut_date=True):
        """

        :param return_type:
        :param obsdatetime:
        :return:
        """

        if obsdatetime:
            self.obsdatetime = obsdatetime
        else:
            self.obsdatetime = datetime.datetime.utcnow()

        if isinstance(self.obsdatetime, datetime.datetime):
            self.obsdatetime.strftime('')

        obstime = Time(self.obsdatetime)
        if use_ut_date:
            obstime = Time()
        sun_set = self.obs_site.sun_set_time(obstime, which="nearest")

        # if sun set is greater than input time we are
        if sun_set > obstime:
            print('Using next')
            which = 'next'
        else:
            print('Using nearest')



            obstime = self.obs_site.twilight_evening_astronomical(obstime,
                                                                  which="nearest")

            which = 'nearest'

        ret = {'sun_set':
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

        if return_type == 'json':
            json_dict = {k: v.iso for k, v in ret.items()}
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


if __name__ == '__main__':
    import time
    import json

    #obs_times = Times()
    obs = ScheduleNight()

    #start = time.time()
    #obs_times.get_observing_times()print(obs.get_observing_times(return_type='json'))
    #print(time.time()-start)
    obsdict = {}
    obsstart = datetime.datetime.strptime('2017-01-01', "%Y-%m-%d")
    start = time.time()
    datestr = "2017-01-01"
    x = obs.get_observing_times('2019-09-12')
    print(x['evening_nautical'].iso)
    print(x['morning_nautical'].iso)
    print((x['morning_nautical']-x['evening_nautical']).to_datetime().seconds)
    while datestr != '2021-10-31':
        datestr = obsstart.strftime("%Y-%m-%d")
        print(datestr)
        x = obs.get_observing_times(datestr)
        print(x['evening_nautical'].iso)
        print(x['morning_nautical'].iso)
        print((x['morning_nautical'] - x['evening_nautical']).to_datetime().seconds)
        obsdict[datestr] = {'morning': x['morning_nautical'].iso, 'evening': x['evening_nautical'].iso,
                            'total': (x['morning_nautical'] - x['evening_nautical']).to_datetime().seconds}
        obsstart += datetime.timedelta(days=1)
    with open('data.txt', 'w') as outfile:
        json.dump(obsdict, 'science_dates.json')
    print(time.time()-start)