# generate ephemerides for SEDM Solar System observations
# usage: geneph.py object_name YYYYMMDD
# (c) Quanzhi Ye

"""
SEDM's TCS requires the following metadata (from Richard Walters, 02/26/19)
Name: (object name)
Start RA: (decimal hours)
Start DEC: (decimal degrees)
equinox: (decimal years)
ra_rate: (RA rate in arcsec/hr)
dec_rate: (DEC rate in arcsec/hr)
epoch: (Epoch of coordinates in decimal hours)
"""

import numpy as np
import sys
from astroquery.jplhorizons import Horizons
from astropy.time import Time, TimeDelta
import datetime
import astroplan
from dateutil import parser

# start of config
time_intv = 2  # time interval of the output [min]
# end of config

obs_site = astroplan.Observer.at_site(site_name='Palomar')


def get_observing_times(obsdatetime=None, return_type=''):
    """

    :param return_type: 
    :param obsdatetime: 
    :return: 
    """
    
    if not obsdatetime:
        obsdatetime = datetime.datetime.utcnow()

    obstime = Time(obsdatetime)
    sun_set = obs_site.sun_set_time(obstime, which="nearest")

    # if sun set is greater than input time we are
    if sun_set > obstime:
        which = 'next'
    else:
        which = 'nearest'

    ret = {
        'sun_set':
            obs_site.sun_set_time(obstime, which=which),
        'sun_rise':
            obs_site.sun_rise_time(obstime, which=which),
        'evening_civil':
            obs_site.twilight_evening_civil(obstime, which=which),
        'evening_nautical':
            obs_site.twilight_evening_nautical(obstime, which=which),
        'evening_astronomical':
            obs_site.twilight_evening_astronomical(obstime, which=which),
        'morning_civil':
            obs_site.twilight_morning_civil(obstime, which=which),
        'morning_nautical':
            obs_site.twilight_morning_nautical(obstime, which=which),
        'morning_astronomical':
            obs_site.twilight_morning_astronomical(obstime, which=which)
    }

    if return_type == 'json':
        json_dict = {k: v.iso for k, v in ret.items()}
        return json_dict
    else:
        return ret


def get_ephem(target, utdate, use_visible_window=True, return_nearest=True,
              input_time=None, print_all=False):
    """

    :param target:
    :param utdate:
    :param use_visible_window:
    :param return_nearest:
    :param input_time:
    :param print_all:
    :return:
    """
    if use_visible_window:
        obstimes = get_observing_times()
        t = Time(obstimes['evening_nautical'])
        tt = [t.jd]

        if datetime.datetime.utcnow() < obstimes['evening_nautical']:
            start = obstimes['evening_nautical']
        else:
            start = Time(datetime.datetime.utcnow())

        count = 0
        print(obstimes['evening_nautical'].to_datetime())
        while start < obstimes['morning_nautical'] and count < 3600:
            t.jd += time_intv / 60 / 24
            tt.append(t.jd)
            start += TimeDelta(time_intv*60, format='sec')
            count += 1

    else:
        t = Time('%sT00:00:00' % utdate, format='isot', scale='utc')
        tt = [t.jd]
        for ti in np.arange(0, 1, time_intv / 60 / 24):
            tt.append(t.jd + ti)

    horizons = Horizons(target, location='675', epochs=tt).ephemerides()

    print('%10s %10s %10s %07s %10s %10s %10s %06s'
          % ('name', 'RA', 'Dec', 'equinox', 'RA rate', 'Dec rate', 'time',
             'V'))

    nearest = 1000
    index = -1
    for i, _ in enumerate(horizons['V']):
        if horizons['solar_presence'][i] is True or horizons['EL'][i] < 12 or \
                horizons['airmass'][i] > 2.3:
            continue

        if return_nearest:
            if not input_time:
                input_time = datetime.datetime.utcnow()

            ephem_date = parser.parse(horizons['datetime_str'][i])
            if input_time > ephem_date:
                dt = (input_time-ephem_date).seconds
            else:
                dt = (ephem_date-input_time).seconds

            if dt < nearest:
                nearest = dt
                index = i

        if print_all:
            print('%10s %10f %10f %07f %10f %10f %s %06f' % (
                target, horizons['RA'][i], horizons['DEC'][i], 2000.0,
                horizons['RA_rate'][i],
                horizons['DEC_rate'][i],
                horizons['datetime_str'][i], horizons['V'][i]))

    if nearest != 1000 and index != -1:
        print('%10s %10f %10f %07f %10f %10f %s %06f' % (
            target, horizons['RA'][index], horizons['DEC'][index], 2000.0,
            horizons['RA_rate'][index],
            horizons['DEC_rate'][index],
            horizons['datetime_str'][index], horizons['V'][index]))

        uthour = datetime_to_hour(horizons['datetime_str'][index])
        print(uthour)


def datetime_to_hour(time_str=""):
    """

    :param time_str:
    :return:
    """
    time_str = time_str.split()[-1]
    time_str = time_str.split(":")
    time_float = [float(i) for i in time_str]
    h, m, s = time_float
    return round(h+(m/60)+(s/3600), 3)


def get_ephem2(target, utdate):
    """

    :param target:
    :param utdate:
    :return:
    """
    t = Time('%sT00:00:00' % utdate, format='isot', scale='utc')

    tt = [t.jd]
    for ti in np.arange(0, 2, 1):
        tt.append(t.jd + ti)
    print(tt)
    horizons = Horizons(target, location='675', epochs=tt).ephemerides()

    print('%10s %10s %10s %07s %10s %10s %10s %06s' % (
        'name', 'RA', 'Dec', 'equinox', 'RA rate', 'Dec rate', 'time', 'V'))

    for i, _ in enumerate(horizons['V']):
        if horizons['solar_presence'][i] is True or horizons['EL'][i] < 12 or \
                horizons['airmass'][i] > 2.3:
            continue

        print('%10s %10f %10f %07f %10f %10f %10f %06f' % (
            target, horizons['RA'][i] / 24, horizons['DEC'][i], 2000.0,
            horizons['RA_rate'][i],
            horizons['DEC_rate'][i],
            (horizons['datetime_jd'][i] - np.floor(horizons['datetime_jd'][i]))
            * 24, horizons['V'][i]))


def command_line():
    t = Time('%s-%s-%sT00:00:00' % (sys.argv[2][0:4], sys.argv[2][4:6],
                                    sys.argv[2][6:8]), format='isot',
             scale='utc')

    tt = [t.jd]
    for ti in np.arange(0, 1, time_intv / 60 / 24):
        tt.append(t.jd + ti)
    print(tt)
    horizons = Horizons(sys.argv[1], location='675', epochs=tt).ephemerides()

    print('%10s %10s %10s %07s %10s %10s %10s %06s' % ('name', 'RA', 'Dec',
                                                       'equinox', 'RA rate',
                                                       'Dec rate', 'time', 'V'))

    for i, _ in enumerate(horizons['V']):
        if horizons['solar_presence'][i] is True or horizons['EL'][i] < 0:
            continue

        print('%10s %10f %10f %07f %10f %10f %10f %06f' % (
            sys.argv[1], horizons['RA'][i] / 24, horizons['DEC'][i], 2000.0,
            horizons['RA_rate'][i], horizons['DEC_rate'][i],
            (horizons['datetime_jd'][i] - np.floor(horizons['datetime_jd'][i]))
            * 24, horizons['V'][i]))


if __name__ == "__main__":
    input_time_str = '2019-05-03 05:00:00'

    get_ephem('2019 CH', '2019-05-03', input_time=parser.parse(input_time_str))
