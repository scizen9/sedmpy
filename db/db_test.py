import SedmDb
import SedmDb_tools as db_tools
import numpy as np


db = SedmDb.SedmDB()


def basictest():
    #add_user',
    user = db.add_user({'username': 'db_test_user', 'name': 'new_user', 'email': '', 'password': 'cvs'})
    assert user[0] >= 0
    # get_from_users'
    assert db.get_from_users(['username', 'name', 'password'], {'id': user[0]}) == [('db_test_user', 'new_user', 'cvs')]
    #update_user',
    db.update_user({'name': 'updated_user', 'id': user[0]})
    assert db.get_from_users(['name'], {'id': user[0]})
    #remove_user',
    #add_usergroup',
    #get_from_usergroups',
    #add_group',
    #remove_from_group',
    #add_program',
    #get_from_program',
    #update_program',
    #add_allocation'
    all = db.add_allocation()
    #get_from_allocation',
    #update_allocation',
    #add_atomicrequest',
    #get_from_atomicrequest',
    #update_atomicrequest',
    #add_phot_calib',
    #add_spec_calib',
    #get_from_phot_calib',
    #get_from_spec_calib',
    #update_phot_calib',
    #update_spec_calib',
    #add_classification',
    #get_from_classification',
    #update_classification',
    #add_flexure',
    #get_from_flexure',
    #add_elliptical_heliocentric',
    #add_hyperbolic_heliocentric',
    #add_parabolic_heliocentric',
    #get_from_elliptical_heliocentric',
    #get_from_hyperbolic_heliocentric',
    #get_from_parabolic_heliocentric',
    #get_object_id_from_name',
    #get_objects_near',
    #add_object',
    #get_from_object',
    #add_observation',
    #get_from_observation',
    #update_observation',
    #_add_planet_satellite_orbit',
    #add_metrics_phot',
    #add_phot',
    #get_from_metrics_phot',
    #get_from_phot',
    #add_request',
    #get_from_request',
    #update_request',
    #expire_requests',
    #add_earth_satellite',
    #get_from_earth_satellite',
    #get_conn_sedmDB',
    #add_spec',
    #get_from_spec',
    #execute_sql',
    #add_telescope_stats',
    #get_from_telescope_stats',



# TODO: make tests that cause IntegrityError (giving a value in of the wrong format e.g. string for a decimal parameter) and ProgrammingError (attempting to insert/update a column that doesn't exist)
# TODO: make the functions reject pardic keys that aren't columns

# there are groups 'default group' and 'second default group'
# (user_id = 1, username=default_user) can be tied to the requests

    
# TODO: more robust testing of the objects
fixed_object = {'marshal_id': 60, 'name': 'test_obj', 'ra': 20., 'dec': 40., 'typedesig': 'f', 'epoch': 2000.}

def test_object_creation():
    # TODO: test all types of objects
    db.add_object(fixed_object)


request_dict = {'object_id': 2, 'user_id': 1, 'program_id': 1, 'marshal_id': 90, 
                'exptime': '{2700, 600}', 'maxairmass': 3.5, 'priority': 3.5, 'inidate': '2016-11-10',
                'enddate': '2016-11-14', 'nexposures': '{1,2,1,2,0}', 'ordering': '{2g,2r,1ifu,2g,2r}'}
# TODO: add to request_dict to test empty and invalid keys

def test_request_manipulation():
    assert db.add_request(request_dict) == (-1, "ERROR: nexposures and ordering are inconsistent!")
    request_dict.pop('nexposures', None)
    db.add_request(request_dict)  # this re-adds nexposures by generating from ordering
    assert db.execute_sql("SELECT id FROM request WHERE nexposures = '{1,0,4,4,0}' AND status != 'EXPIRED'")
    request_dict.pop('ordering', None)
    request_dict.pop('nexposures', None)
    assert db.add_request(request_dict) == (-1, "ERROR: nexposures or ordering is required!")
    request_dict['ordering'] = '{2g,2r,1ifu,2g,2r}'  # return ordering to the dict

    request_dict.pop('priority', None)
    assert db.add_request(request_dict) == (-1, "ERROR: priority not in dictionary!")
    request_dict['priority'] = 3.5  # re-add it to the dict

    request_dict['object_id'] = 0
    assert db.add_request(request_dict) == (-1, "ERROR: object does not exist!")
    request_dict['object_id'] = 1


def test_request_update():
    obj, status, priority = db.execute_sql("SELECT object_id, status, priority FROM request WHERE id='1'")[0]
    db.update_request({'id': 1, 'object_id': int(obj+1), 'status': 'NEW', 'priority': float(priority+1)})
    new_obj, new_stat, new_priority = db.execute_sql("SELECT object_id, status, priority FROM request WHERE id='1'")[0]
    # object_id isn't allowed to be updated and the given status isn't valid
    assert new_obj == obj and new_stat == status and new_priority == priority+1
    assert db.update_request({'object_id': obj+1, 'status': 'ACTIVE'}) == (-1, "ERROR: no id provided!")
    assert db.update_request({'id': 0, 'maxairmass': 3}) == (-1, "ERROR: request does not exist!")
    assert db.update_request({'id': 1, 'creationdate': '2016-02-04'}) == (-1, "ERROR: no parameters given to update!")
    db.update_request({'id': 1, 'object_id': 1, 'status': 'PENDING', 'priority': 6.})  # reset to defaults


def test_request_expiration():
    test_request = db.execute_sql("SELECT id FROM request WHERE enddate = '2016-11-11';")
    if not test_request:
        request_dict['enddate'] = '2016-11-11'
        db.add_request(request_dict)
    else:
        db.update_request({'id': int(test_request[0][0]), 'status': 'PENDING'})
    assert 'PENDING' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11';")[0][0]
    db.expire_requests()
    assert 'EXPIRED' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11';")[0][0]


def test_cancel_request():
    test_request = db.execute_sql("SELECT * FROM request")
    if not test_request:
        db.add_request(request_dict)
    else:
        db.update_request({'id': 1, 'status': 'PENDING'})
    assert 'PENDING' == db.execute_sql("SELECT status FROM request WHERE id=1")[0][0]
    db.cancel_scheduled_request(1)
    assert 'CANCELED' == db.execute_sql("SELECT status FROM request WHERE id=1")[0][0]
    assert db.cancel_scheduled_request(0) == (-1, "ERROR: request does not exist!")

atomic_dict = {'exptime': 1800.0, 'filter': 'f', 'priority': 32., 'inidate': '2016-12-12',
               'enddate': '2017-1-12', 'object_id': 1, 'order_id': 1, 'when': '54'}
# 'when' should do nothing (if it does something sql will error)
# lack of 'request_id' should cause a failure


def test_atomicrequest_manipulation():
    # TODO: test for add
    assert db.add_atomicrequest(atomic_dict) == (-1, "ERROR: request_id not provided!")
    atomic_dict['request_id'] = 3.2
    print db.add_atomicrequest(atomic_dict)
    assert db.add_atomicrequest(atomic_dict) == (-1, "ERROR: request_id must be of type 'int'!")
    atomic_dict['request_id'] = db.get_from_atomicrequest(['request_id'], {'object_id': 1})
    assert db.add_atomicrequest(atomic_dict)[0] == 0


    # TODO: test for update
    pass


def test_request_atomic_requests():
    # create_request_atomic_requests is called by add_request already, so its tests above test valid cases
    assert db_tools.create_request_atomic_requests(0) == (-1, "ERROR: request does not exist!")
    # TODO: find a way to test other situations that doesn't run into (atomicrequests for that request already exist)?


def test_get_request_atomic_requests():
    pass


def test_fits_header_parse():
    fitsfile = '/scr2/sedm/phot/20161012/rc20161012_12_36_34.fits'
    db_tools.add_observation_fitsfile(fitsfile)


def baseline_tests():
    # db.add_user
    # db.get_from_users
    # db.remove_user
    # db.add_group
    # db.add_usergroup
    # db.get_from_usergroups
    # db.remove_from_group

#    print db.add_object({'name': 'test_obj', 'ra': 40., 'dec': 30., 'epoch': 2000., 'typedesig': 'f'})
#    print db.get_from_object(['ra', 'dec'], {'name': 'test_obj'}) , (40, 30)
#    print db.get_object_id_from_name('defa') , 2
#    print db.add_elliptical_heliocentric({'object_id': 1, 'inclination': 0.1, 'perihelion_o': 0.5, 'a': 3.2, 'n': 20,
#                                          'e': 0.01, 'M': 1., 'mjdepoch': 2000, 'D': 2009, 'M1': 1., 'M2': 2.,
#                                          'longascnode_O': 3.})
#    print db.get_from_elliptical_heliocentric(['e', 'mjdepoch'], {'D': 2009}), (0, 2000)
#    print db.add_hyperbolic_heliocentric({'object_id': 1, 'T': '2000-01-01', 'inclination': 0.1, 'longascnode_O': 0.5,
#                                         'perihelion_o': 2., 'e': 3.2, 'q': 20., 'D': 2009, 'M1': 1., 'M2': 2.})
#    print db.get_from_hyperbolic_heliocentric(['e', 'q'], {'D': 2009}), (3.2, 20)
#    print db.add_parabolic_heliocentric({'object_id': 1, 'T': '2000-01-01', 'inclination': 0.1, 'longascnode_O': 0.5,
#                                         'perihelion_o': 2., 'q': 20., 'D': 2009, 'M1': 1., 'M2': 2.})
#    print db.get_from_parabolic_heliocentric(['object_id', 'inclination'], {'M1': 1.}), (1, 0.1)
#    print db.add_earth_satellite({'object_id': 1, 'T': '2000-01-01', 'inclination': 0.1, 'ra': 0.5, 'e': 2., 'pedigree': 2.,
#                                  'M': 3.2, 'n': 20., 'decay': 1., 'reforbit': 10})
#    print db.get_from_earth_satellite(['object_id', 'M'], {'n': 20.}), (1, 3.2)

    print db.add_request({'object_id': 1, 'user_id': 1, 'program_id': 1, 'exptime': '{2500, 300}', 'priority': 3.,
                          'inidate': '2017-01-01', 'enddate': '2017-01-02', 'ordering': '{1ifu, 1g, 1r, 1i}'})
    ids = db.get_from_request(['id'], {'enddate': '2017-01-02'})
    print ids
    print db.update_request({'id': int(ids[0][0]), 'priority': 5.})
    print db.get_from_request(['priority'], {'id': int(ids[0][0])})
    print db.expire_requests()
    print db.get_from_request(['status'], {'enddate': '2017-01-02'})

    print db.add_atomicrequest({'request_id': 1, 'exptime': 300., 'filter': 'g', 'priority': 5., 'inidate':'2017-01-01',
                                'enddate': '2017-04-03'})
    ids = db.get_from_atomicrequest(['id'], {'enddate': '2017-04-03'})
    print ids
    print db.update_atomicrequest({'id': ids[0][0], 'status': 'CANCELED'})  # why doesn't have 'id' long fail?
    print db.get_from_atomicrequest(['status'], {'id': int(ids[0][0])})

    obs_dict = {'object_id': 1, 'request_id': 1, 'atomicrequest_id': 1, 'mjd': 23452.34, 'airmass': 2.2,
                'exptime': 299., 'fitsfile': 'made_up_file', 'lst': 'lst', 'ra': 20., 'dec': 10., 'tel_ra': '42',
                'tel_dec': '99', 'tel_az': 12., 'tel_el': 32., 'tel_pa': 42., 'ra_off': .1, 'dec_off': .2}
    print db.add_observation(obs_dict)
    # TODO: find out why it returned None
    print db.update_observation
    print db.get_from_observation(['fitsfile'], {'dec_off': .2})
    print db.add_telescope_stats
    print db.get_from_telescope_stats

    print db.add_phot({'observation_id': 1, 'astrometry': 'false', 'filter': 'g', 'reducedfile': 'nofile', 'sexfile': 'c',
                       'biasfile': 'd', 'maskfile': 'e', 'flatfile': 'f', 'pipeline': 'phot', 'marshal_phot_id': 1})
    # test add_phot as update
    print db.get_from_phot(['biasfile'], {'filger': 'g'})
    print db.add_spec({'observation_id': 1, 'reducedfile': 'a', 'sexfile': 'b', 'biasfile': 'c', 'flatfile': 'd',
                       'imgset': 'new', 'quality': 3, 'cubefile': 'e', 'standardfile': 'f', 'skysub': 'true'})
    # tested add_spec as update as well
    print db.get_from_spec(['imgset', 'skysub'], {'sexfile': 'b'})

    print db.add_metrics_phot({'phot_id': 1, 'fwhm': 4., 'background': 102., 'zp': 2., 'zperr': 8.,
                               'ellipticity': 0., 'nsources': 40})
    # test add_metrics_phot as update
    print db.get_from_metrics_phot(['nsources', 'fwhm'], {'phot_id': 1})
    print db.add_metrics_spec({'spec_id': 1, 'fwhm': 2., 'background': 1.})
    # test add_metrics_spec as update
    print db.get_from_metrics_spec(['life_fwhm', 'background'], {'spec_id': 1})
    obs_dict['atomicrequest_id']= 2
    db.add_observation(obs_dict)
    db.add_spec({'observation_id': 2, 'reducedfile': 'a', 'sexfile': 'b', 'biasfile': 'c', 'flatfile': 'd',
                 'imgset': 'new', 'quality': 3, 'cubefile': 'e', 'standardfile': 'f', 'skysub': 'true'})
    print db.add_flexure({'rms': 1., 'spec_id_1': 1, 'spec_id_2': 2, 'timestamp1': '2017-01-12 10:20:20',
                          'timestamp2': '2017-01-12 10:20:40'})
    print db.get_from_flexure(['rms'], {'spec_id_1': 1})
    print db.add_classification({'spec_id': 1, 'object_id': 1, 'classification': 'none', 'redshift': 0.1,
                                 'redshift_err': 1., 'classifier': 'none', 'score': 0.})
    print db.update_classification
    print db.get_from_classification(['classification', 'score'], {'classifier': 'none'})


if __name__ == '__main__':
    baseline_tests()

