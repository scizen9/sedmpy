
import numpy as np
from astropy.time import Time
from astropy.io import fits

class DbTools(object):
    def __init__(self, db_class):
        self.db = db_class
    # TODO: refactor and move functions into the class

    def submit_asteroids_edb(self, file = '/home/sedm/kpy/ephem/asteroids.edb', delimiter=','):
        f = open(file)
        for line in f:
            if line[0] != '#':
                par = line.split(',')
                if par[-1][-1] == 'n':  # remove /n from the end of the line
                    par[-1] = par[-1][:-1]
                obj_dict = {'name': par[0], 'typedesig': par[1]}
                obj_id = self.db.add_object(obj_dict)[0]
                if obj_id == -1:
                    print("object %s failed to be created" % (par[0],))
                    continue
                if par[1] == 'e':
                    time=par[9].split('/')
                    pardic = {'object_id': obj_id, 'inclination': float(par[2]), 'longascnode_O': float(par[3]), 'perihelion_o': float(par[4]),
                    'a': float(par[5]), 'n': float(par[6]), 'e': float(par[7]), 'M': float(par[8]), 'mjdepoch': Time(time[2]+'-'+time[0]+'-'+time[1]).mjd, 'D': float(par[10]), 'M1': float(par[11]),
                    'M2': float(par[12])}
                    if len(par) == 14:
                        pardic['s'] = float(par[13])
                    add = self.db.add_elliptical_heliocentric(pardic)
                    print(add)
                """
                if par[1] == 'h':
                    pardic = {'object_id': obj_id, 'T': par[2], 'e': par[6], 'inclination': par[3], 'longascnode_O': par[4],
                       'perihelion_o': par[5], 'q': par[7], 'D': float(par[8]), 'M1': par[9], 'M2': par[10]}
                    if len(par) == 12:
                        pardic['s'] =par[11]
                    add=self.db.add_hyperbolic_heliocentric(pardic)
                    print add
                if par[1] == 'p':
                    pardic = {'object_id': obj_id, 'T': par[2], 'inclination': par[3], 'longascnode_O': par[5],
                       'perihelion_o': par[4], 'q': par[5], 'D': float(par[7]), 'M1': par[8], 'M2': par[9]}
                    if len(par) == 11:
                        pardic['s'] =par[10]
                    add=self.db.add_parabolic_heliocentric(pardic)
                    print add
                if par[1] == 'E':
                    pardic = {'object_id': obj_id, 'T': par[2], 'e': par[5], 'inclination': par[3], 'ra': par[4],
                       'pedigree': par[6], 'M': par[7], 'n': par[8], 'decay': par[9], 'reforbit': par[10], 'drag': par[11]}
                    add=self.db.add_earth_satellite(pardic)
                    print add
                """

    def get_object_parameters(self, object_id):
        """
        Takes an object_id and returns the object's information and orbit parameters for SSO
        Args:
            object_id (int): id of the object

        Returns:
            (-1, "ERROR...") if an error occurred

            [('marshal_id', 'name', 'iauname', 'ra', 'dec', typedesig', 'epoch'), None] for fixed objects, or if there are
            no orbit parameters in the database

            [('marshal_id', 'name', 'iauname', 'ra', 'dec', typedesig', 'epoch'), (orbit params arranged in ephem ordering)]
        """  # TODO: test if the orbit_params can be input directly into an ephem call
        object = self.db.get_from_object(['marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch'],
                                    {'id': object_id})
        if not object:
            return (-1, "ERROR: no object found matching the given object_id")
        elif object[0] == -1:  # object[0] should be positive or empty if successful
            return object  # this should carry the error message
        # format the response depending on the typedesig
        if object[0][5] == 'f':
            return [object[0], None]

        elif object[0][5] == 'e':
            orbit_params = self.db.get_elliptical_orbit(['inclination', 'longascnode_0', 'perihelion_o', 'a', 'n', 'e', 'M',
                                                    'mjdepoch', 'D', 'M1', 'M2', 's'], {'object_id': object_id})
            if not orbit_params:
                return [object[0], None]
            elif orbit_params[0] == -1:  # without an error orbit_params[0] should be a tuple
                return orbit_params
            else:
                return [object[0], orbit_params[0]]

        elif object[0][5] == 'h':
            orbit_params = self.db.get_hyperbolic_orbit(['T', 'inclination', 'longascnode_0', 'perihelion_o', 'e', 'q', 'D',
                                                    'M1', 'M2', 's'], {'object_id': object_id})
            if not orbit_params:
                return [object[0], None]
            elif orbit_params[0] == -1:  # without an error orbit_params[0] should be a tuple
                return orbit_params
            else:
                return [object[0], orbit_params[0]]

        elif object[0][5] == 'p':
            orbit_params = self.db.get_parabolic_orbit(['T', 'inclination', 'perihelion_o', 'q', 'longascnode_0', 'D',
                                                   'M1', 'M2', 's'], {'object_id': object_id})
            if not orbit_params:
                return [object[0], None]
            elif orbit_params[0] == -1:  # without an error orbit_params[0] should be a tuple
                return orbit_params
            else:
                return [object[0], orbit_params[0]]

        elif object[0][5] == 'E':
            orbit_params = self.db.get_earth_satellite_orbit(['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n', 'decay',
                                                         'reforbit', 'drag'], {'object_id': object_id})
            if not orbit_params:
                return [object[0], None]
            elif orbit_params[0] == -1:  # without an error orbit_params[0] should be a tuple
                return orbit_params
            else:
                return [object[0], orbit_params[0]]

        elif object[0][5] == 'P':
            raise NotImplementedError  # TODO: add this

        else:
            return (-1, "ERROR: object has an invalid typedesig '%s'" % (object[0][5],))

    def program_time_used(self, start_date, end_date, program_id):

        # TODO: return to SedmDb.py because of how much sql "understanding" it requires?
        # TODO: evaluate assumption that the most recent change to the request will be setting 'status' to 'COMPLETED'
        where_dict = {'program_id': program_id, 'lastmodified': '>' + start_date, 'lastmodified': '<' + end_date}
        sql = ("SELECT exptime, nexposures, phasesamples FROM request WHERE program_id='%s', status='COMPLETED', "
               "lastmodified>'%s', lastmodified<'%s';" % (program_id, start_date, end_date))
        # TODO: fix this, it's really bad that it doesn't go through the checks of SedmDb functions
        # TODO: allow where_dict values to be lists (if it is a list iterate through)?
        program_requests = self.db.execute_sql(sql)
        if not program_requests:
            return 0.0
        else:
            raise NotImplementedError("exptime isn't guarenteed a particular format, and calculations aren't implemented")
        # TODO: solidify what "exptime" is, then implement calculations

    def get_user_requests(self, user_id=None, username=None, values=['id', 'object_id', 'program_id', 'status']):
        """
        Get the parameters for all requests associated with a user

        Args:
            user_id (int): id of the user
            username (str): username of the user
                only needed if user_id is not provided
            values (list): values to be returned
                defaults to ['id', 'object_id', 'program_id', 'status']

        Returns:
            list of tuples conaining the ``values`` for each request of the user
            (-1, "ERROR ... ") if there was an issue
            [] if no requests are associated with the user
        """

        if user_id:
            where_dict = {'user_id': user_id}
        elif username:
            user_id = self.db.get_from_users(['id'], {'username': username})
            if not user_id:
                return (-1, 'ERROR: not user with that username!')
            where_dict = {'user_id': user_id[0][0]}
        else:
            return (-1, "ERROR: neither username nor user_id were provided!")
        user_requests = self.db.get_from_request(values, where_dict)
        return user_requests

    def get_active_requests(self, values=['id', 'object_id', 'user_id', 'program_id', 'marshal_id', 'exptime',
                            'maxairmass', 'priority', 'cadence', 'phasesamples', 'sampletolerance',
                            'filters', 'nexposures', 'ordering']):
        """
        Get the parameters for all 'ACTIVE' requests
        Args:
            values (list): the values to get for each request
                defaults to ['id', 'object_id', 'user_id', 'program_id', 'marshal_id', 'exptime',
                             'maxairmass', 'priority', 'cadence', 'phasesamples', 'sampletolerance',
                             'filters', 'nexposures', 'ordering']
        Returns:
            list of tuples conaining the ``values`` for each active request
            (-1, "ERROR ... ") if there was an issue
            [] if no requests are 'ACTIVE'
        """
        # TODO: test?

        active_requests = self.db.get_from_request(values, {'status': 'ACTIVE'})
        return active_requests

    def cancel_request(self, requestid):
        """
        Changes the status of the request and any related atomicrequests to "CANCELED"
        """
        # TODO: return to SedmDb.py because of how much sql "understanding" it requires?
        self.db.update_request({'id': requestid, 'status': 'CANCELED'})
        # cancel the associated atomicrequests
        # TODO: allow more nuanced update function inputs (e.g. add a where_dict)?
        self.db.execute_sql("UPDATE atomicrequest SET status='CANCELED' WHERE request_id='%s'" % (requestid,))
        return (0, "Request canceled")

    def create_atomic_requests(self):
        """
        Finds requests without atomicrequets and creates their atomicrequests

        Returns:
            (0, "...") containing a list of the requests and which, if any, failed
        """
        # get the ids of the requests with no atomicrequest
        request_ids = self.db.execute_sql('SELECT id FROM request WHERE NOT EXISTS '
                                     '(SELECT id FROM atomicrequest WHERE request.id = atomicrequest.request_id);')
        # TODO: make this only create for requests with status 'PENDING'/'ACTIVE'?
        failed = []
        for request in request_ids:
            success = self.create_request_atomic_requests(request)
            if success[0] == -1:
                failed.append(request)
        for request in failed:
            request_ids.remove(request)
        if failed:
            return (0, "Added atomicrequests for requests %s, failed to add for request  %s" % (request_ids, failed))
        else:
            return (0, "Added atomicrequests for requests %s" % (request_ids,))

    def create_request_atomic_requests(self, request_id):
        """
        create atomicrequest entries for a given request

        Args:
            request_id: int
                id of the request that needs atomicrequests

        Returns:
            (-1, "ERROR: ...") if there was an issue
            (0, "Added (#) atomic requests for request (#)") if it was successful
        """
        status = self.db.execute_sql("SELECT status FROM request WHERE id=%s" % (request_id,))
        if not status:
            return (-1, "ERROR: request does not exist!")
        elif status == 'CANCELED':
            return (-1, "ERROR: request has been canceled!")
        elif status == 'EXPIRED':
            return (-1, "ERROR: request has expired!")

        if self.db.execute_sql("SELECT id FROM atomicrequest WHERE request_id='%s';" % (request_id,)):
            return (-1, "ERROR: atomicrequests already exist for that request!")
        request = self.db.get_from_request(['object_id', 'exptime', 'priority', 'inidate', 'enddate', 'cadence', 'phasesamples',
                                        'sampletolerance', 'filters', 'nexposures', 'ordering'], {'id': request_id})[0]
        # TODO: implement cadence/phasesamples/sampletolerance (I have no idea how they interact with nexposures)
        pardic = {'object_id': int(request[0]), 'priority': float(request[2]), 'inidate': str(request[3]),
                  'enddate': str(request[4]), 'request_id': int(request_id)}
        obs_order = []
        if request[10]:
            for num_fil in request[10]:
                for n in range(int(num_fil[0])):  # the number should be single digit
                    obs_order.append(num_fil[1:])
        elif request[9]:
            for filter_idx in range(len(request[8])):
                for n in range(request[9][filter_idx]):
                    obs_order.append(request[8][filter_idx])
        else:
            return (-1, "ERROR: request contains neither nexposures nor ordering!")  # add_request should prevevnt this
        ifus = np.where(np.array(obs_order) == 'ifu')[0]
        if len(ifus) == 2:  # TODO: rewrite, have another way of indicating a/b
            obs_order[ifus[0]] = 'ifu_a'
            obs_order[ifus[1]] = 'ifu_b'
        if any([(filt not in ['u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b']) for filt in obs_order]):
            return (-1, "ERROR: either filters or ordering has an invalid entry!")
        for n, filter_des in enumerate(obs_order):
            pardic['filter'] = filter_des
            pardic['order_id'] = n
            # TODO: do exptime modifications here per filter and ab vs single exposure
            if 'ifu' in filter_des:
                pardic['exptime'] = float(request[1][0])  # TODO: make sure the sql returns the proper format
            else:
                pardic['exptime'] = float(request[1][1])
            add_return = self.db.add_atomicrequest(pardic)
            if add_return[0] == -1:
                return (-1, "ERROR: adding atomicrequest #%s, filter:%s failed with message, '%s'!"
                        % (n + 1, filter_des, add_return))
        return (0, "Added %s atomic requests for request %s" % (len(obs_order), request_id))

    def get_request_atomicrequests(self, request_id, values):
        """
        return the atomicreqests associated with a single request

        Args:
            request_id (int):
                The request_id of the desired request
            values (list):
                list of values to be retreived about each atomicrequest

        Returns:
            list of tuples:
                (id, object_id, order_id, exptime, filter, status, priority)
                for each atomicreqest associated with the desired request
            [] if there were no atomicrequests matching the given request_id
            (-1, "ERROR...") if there was an issue

        """
        if not isinstance(request_id, int):
            return []
        # TODO: test
        atomic_requests = self.db.get_from_atomicrequest(values, {'request_id': request_id})
        return atomic_requests

    def add_observation_fitsfile(self, fitsfile, atomicrequest_id):
        """
        Adds an observation from a fitsfile, sets the atomicobservation to 'OBSERVED'. Checks if all the atomicobservations
        associated with its request are observed, if so it sets the request's status to 'COMPLETED'.

        Args:
            fitsfile: str (fits file path)
            atomicrequest_id: int (temporary until it is included in the header)
        """
        hdulist = fits.open(fitsfile)
        header = hdulist[0].header
        # TODO: make sure all of the parameters are of the right type (so SedmDb.py doesn't kill it)
        header_dict = {'mjd': header['JD'] - 2400000.5, 'airmass': header['AIRMASS'], 'exptime': header['EXPTIME'],
                       'fitsfile': fitsfile, 'lst': header['LST'], 'ra': ra_to_decimal(header['RA']),
                       'dec': dec_to_decimal(header['DEC']), 'tel_ra': header['TEL_RA'], 'tel_dec': header['TEL_DEC'],
                       'tel_az': header['TEL_AZ'], 'tel_el': header['TEL_EL'], 'tel_pa': header['TEL_PA'],
                       'ra_off': header['RA_OFF'], 'dec_off': header['DEC_OFF'], 'camera': header['CAM_NAME'],
                       'atomicrequest_id': int(atomicrequest_id)}  # TODO: remove atomicrequest_id arg from function
        # header_dict['atomicrequest_id'] = int(header['ATOM_ID'])
        ids = self.db.get_from_atomicrequest(['request_id', 'object_id'], {'id': header_dict['atomicrequest_id']})[0]
        if ids:
            header_dict['request_id'] = int(ids[0][0])
            header_dict['object_id'] = int(ids[0][1])
        else:
            return (-1, "ERROR: no atomicrequests found with an id matching the header's ATOM_ID!")

        obs_added = self.db.add_observation(header_dict)
        # TODO: add ATOM_ID to the header (above)
        # TODO: add 'imtype', generate from fitsfile name?

        tel_stats = {'date': header['OBSDATE'], 'dome_status': header['DOMEST'], 'in_temp': header['IN_AIR'],
                     'in_humidity': header['IN_HUM'], 'in_dew': header['IN_DEW'], 'out_temp': header['OUT_AIR'],
                     'out_humidity': header['OUT_HUM'], 'out_dew': header['OUT_HUM'], 'wind_dir': header['WIND_DIR'],
                     'wsp_cur': header['WSP_CUR'], 'wsp_avg': header['WSP_AVG'], 'mir_temp': header['MIR_TEMP'],
                     'top_air': header['TOP_AIR'], 'pri_temp': header['PRI_TEMP'], 'sec_temp': header['SEC_TEMP'],
                     'flo_temp': header['FLO_TEMP'], 'bot_temp': header['BOT_TEMP'], 'mid_temp': header['MID_TEMP'],
                     'top_temp': header['TOP_TEMP']}
        observation_id = self.db.get_from_observations(['id'], {'atomicrequest_id': header_dict['atomicrequest_id']})
        if observation_id:
            tel_stats['observation_id'] = int(observation_id[0][0])
        else:
            return obs_added

        stats_added = self.db.add_tel_stats(tel_stats)

        # check if all of the request's atomicrequests are 'OBSERVED', if so set the request's status to 'COMPLETED'
        # TODO: must test with connected request/atomicrequest/fitsfile
        if obs_added[0] == 0:
            print(self.db.update_atomicrequest(
                {'id': header_dict['atomicrequest_id'], 'status': 'OBSERVED'}))  # better way than print?
            request_id = self.db.execute_sql("SELECT request_id FROM atomicrequest WHERE id='%s'"
                                        % (header_dict['atomicrequest_id'],))
            atomic_requests = self.db.get_request_atomicrequests(request_id)
            if all(request[5] == 'OBSERVED' for request in atomic_requests):
                self.db.update_request({'id': request_id, 'status': 'COMPLETED'})
        return obs_added


def ra_to_decimal(ra):
    hms = ra.split(':')
    return float(hms[0])*15+float(hms[1])/4+float(hms[2])/240


def dec_to_decimal(dec):
    dms = dec.split(':')
    if dms[0][0] == '-':
        return -(abs(float(dms[0]))+float(dms[1])/60+float(dms[2])/3600)
    else:
        return float(dms[0])+float(dms[1])/60+float(dms[2])/3600
