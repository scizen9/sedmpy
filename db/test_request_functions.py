import SedmDb
import datetime

# default entries (shouldn't be changed by functions)
# table: users; 'id': 1, 'username': 'default_user', 'name': 'apo', 'email': 'kiwi'
# table: gropus; 'id': 1, 'designator': 'default group'
# table: groups; 'id': 2, 'designator': 'second default group'
# table: object; 'id': 1, 'name': 'vega', 'ra': 279.23, 'dec': 38.783, 'typedesig': 'f', 'epoch': 2000.
# table: request; 'id': 1, 'object_id': 1, 'user_id': 1, 'exptime': '{ 2400, 360}', 'priority': 6.,
#                 'inidate': '2016-12-30', 'enddate': '2017-02-12', 'nexposures': '{1, 2, 2, 2, 2}'
# table: atomicrequest; 'id': 1 through 'id': 9
db = SedmDb.SedmDB()


def test_request_manipulation():
    # test successful add_request
    req = {'object_id': 1, 'user_id': 1, 'program_id': 4, 'exptime': '{2400, 360}', 'priority': 5.6,
           'inidate': '2017-02-03', 'enddate': '2017-04-25', 'ordering': '{2g,2r,1ifu,2g,2r}', 'maxairmass': 3.2}
    assert db.add_request(req) == (0, 'Request added')
    # test get_from and update
    assert (db.get_from_request(['user_id', 'inidate', 'nexposures'], {'enddate': '2017-04-25'})[0] ==
            (1, datetime.date(2017, 02, 03), [1, 0, 4, 4, 0]))
    id = db.get_from_request(['id'], {'enddate': '2017-04-25'})
    upd = {'id': id[0][0], 'status': 'ACTIVE', 'maxairmass': 2.2, 'user_id': 2} # user_id isn't allowed
    assert db.update_request(upd) == (0, "Requests updated")
    assert (db.get_from_request(['status', 'maxairmass', 'user_id'], {'enddate': '2017-04-25'})[0] ==
            ('ACTIVE', 2.2, 1))
    # test get_from failures
    assert db.get_from_request([], {})[0] == (-1, "ERROR: no valid values requested!")
    assert db.get_from_request(['greetings'], {})[0] == (-1, "ERROR: no valid values requested!")
    assert db.get_from_request(['status'], {'maxairmass': 'r'}) == (-1, "ERROR: maxairmass must be of type 'float'!")
    # test add_request failure
    req['object_id'] = 0
    assert db.add_request(req) == (-1, "ERROR: object does not exist!")
    req['object_id'] = 1
    req['user_id'] = 0
    assert db.add_request(req) == (-1, "ERROR: user does not exist!")
    req.pop('enddate')
    req['user_id'] = 1
    assert db.add_request(req) == (-1, "ERROR: enddate not in dictionary!")
    # test update_request failure
    upd = {'id': id[0][0], 'user_id': 2, 'object_id': 3}
    assert db.update_request(upd) == (-1, "ERROR: no parameters given to update!")
    upd = {'id': id[0][0], 'maxairmass': 'r'}
    assert db.update_request(upd) == (-1, "ERROR: maxairmass must be of type 'float'!")
    # test expire_requests
    req['inidate'] = '2016-02-03'
    req['enddate'] = '2016-04-25'
    db.add_request(req)
    assert db.expire_requests() == (0, "Requests expired")
    assert db.get_from_request(['status'], {'enddate': '2016-04-25'})[0][0] == 'EXPIRED'


def test_atomicrequest_manipulation():
    # test adding an atomicrequest
    areq = {'request_id':1, 'exptime': 240, 'filter': 'g', 'priority': 3, 'inidate': '2017-01-01',
            'enddate': '2017-02-23', 'object_id': 1}
    assert db.add_atomicrequest(areq) == (0, "Request added")
    # test mismatch object_id
    areq['object_id'] = 2
    assert db.add_atomicrequest(areq) == (-1, "ERROR: object_id given doesn't match request_id!")
    # test invalid filter
    areq['object_id'] = 1
    areq['filter'] = 'n'
    assert db.add_atomicrequest(areq) == (-1, "ERROR: invalid filter given!")
    # test other failure
    areq.pop('filter')
    assert db.add_atomicrequest(areq) == (-1, "ERROR: filter not provided!")
    # test successful get_from and update
    id = db.get_from_atomicrequest(['id'], {'enddate': '2017-02-23'})[0][0]
    assert db.get_from_atomicrequest(['exptime', 'filter', 'object_id'], {'id': id})[0] == (240, 'g', 1)
    upd = {'id': id, 'exptime': 250, 'object_id': 2}  # object_id doesn't update
    assert db.update_atomicrequest(upd) == (0, "Atomicrequest updated")
    assert db.get_from_atomicrequest(['exptime', 'filter', 'object_id'], {'id': id})[0] == (250, 'g', 1)
    # test update failure
    assert db.update_atomicrequest({'id': 0, 'priority': 3.3}) == (-1, "ERROR: atomicrequest does not exist!")
    assert db.update_atomicrequest({'id': 1, 'priority': 'l'}) == (-1, "ERROR: priority must be of type 'float!")
    # test get_from failure
    assert db.get_from_atomicrequest(['status', 'priority'], {'inidate': 43}) == (-1, "ERROR: inidate must be of the "
                                                                                      "format 'year-month-day'!")
# TODO: test create_atomicrequests_from_request
import SedmDb_tools

if __name__ == '__main__':
    test_request_manipulation()
    test_atomicrequest_manipulation()
