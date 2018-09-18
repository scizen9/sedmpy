import SedmDb
import SedmDb_tools


db = SedmDb.SedmDB()

# add default user
db.add_user({'username': 'default_user', 'name': 'apo', 'email': 'kiwi'})

# add default groups
db.add_group({'designator': 'default group'})
db.add_group({'designator': 'second default group'})

#  add default objects
vega = {'ra': 279.23, 'dec': 38.783, 'name': 'vega', 'typedesig': 'f', 'epoch': 2000.}
db.add_object(vega)
default_object = {'ra': 20., 'dec': 40., 'name': 'default_object', 'typedesig': 'f', 'epoch': 2000.,
                  'marshal_id': 90}
db.add_object(default_object)

# add default requests
request = {'object_id': 1, 'user_id': 1, 'program_id': 1, 'exptime': '{ 2400, 360}', 'priority': 6.,
           'inidate': '2016-12-30', 'enddate': '2017-02-12', 'nexposures': '{1, 2, 2, 2, 2}'}
db.add_request(request)

# add default atomic requests
SedmDb_tools.create_request_atomic_requests(1)


# add default observations?
