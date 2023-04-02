from marshals import interface
import json
import db.SedmDb
from datetime import datetime


sedmdb = db.SedmDb.SedmDB()


def process_new_request(request, isfile=True, add2db=True,
                        check_rejection=False, create_request=True):
    """

    :param request: request string path
    :param isfile:
    :param add2db:
    :param check_rejection:
    :param create_request:

    :return:
    """
    # 1. Open the request
    req_dict = interface.read_request(request, isfile=isfile)

    # 2. Check if there were any issues reading in the request
    if 'iserror' in req_dict:
        print("Error reading in the request file")
        return False, req_dict['msg']
    else:
        print("Request is in proper format")

    # 3. Check any rejection criteria
    # TODO: Check database to determine if group has allocation time
    if check_rejection:
        ret = interface.checker(req_dict, check_source=False)
        if 'iserror' in ret:
            return False, ret['msg']
        else:
            print("Target passed checker")

    # 4. Determine if a new request or if we need to update an old request
    if req_dict['status'].lower() == 'delete':
        delete_request_entry(req_dict)
        create_request = False
    elif req_dict['status'].lower() == 'edit':
        delete_request_entry(req_dict)

    # 4. Now check to see if we are creating a request entry
    if create_request:
        if create_request_entry(req_dict):
            print("new request succeeded")
        else:
            print("new request failed")

    # 5. Add the raw request to database
    if add2db:
        _ = add_request_to_db(req_dict)


def delete_request_entry(request_dict):
    # a. start by trying to get the request id from the marshal id
    sedm_requestid = sedmdb.get_from_request(
        values=['id'], where_dict={'marshal_id': request_dict['requestid'],
                                   'external_id': 2})
    print("Canceling SedmDB request ids:\n", sedm_requestid)

    if not sedm_requestid:
        print("No request matching that id")
    else:
        for i in sedm_requestid:
            try:
                print("Canceling: ", i[0])
                print(sedmdb.update_request({'id': i[0], 'status': 'CANCELED'}))
            except Exception as e:
                print(str(e))


def add_request_to_db(request_dict):
    """

    :param request_dict:
    :return:
    """

    # Convert previous photmetry to json string
    try:
        request_dict['prior_photometry'] = json.dumps(
            request_dict['prior_photometry'])
    except Exception as e:
        print(str(e))
        request_dict['prior_photometry'] = '{"error": "conversion error"}'

    columns = ', '.join(str(x).replace('/', '_').lower()
                        for x in request_dict.keys())
    values = ', '.join("'" + str(x).replace('/', '_') + "'"
                       for x in request_dict.values())

    # This table is not currently used, AFAIK: neill 2021-Oct-04
    insert_statement = 'INSERT INTO marshal_requests (%s) VALUES (%s)' % (
        columns, values)
    try:
        ret = sedmdb.execute_sql(insert_statement)
        return ret
    except Exception as e:
        print(str(e), "Unable to add request to database")


def create_request_entry(request, custom_dict=None,
                         fill_by_custom_dict=False,
                         add_raw_to_db=True):
    """

    :param request:
    :param custom_dict:
    :param fill_by_custom_dict:
    :param add_raw_to_db:
    :return:
    """

    if add_raw_to_db:
        pass

    # 0. If delete request, cancel it from the database if it has an entry
    if request['status'].lower() == 'delete':
        # a. start by trying to get the request id from the marshal id
        sedm_requestid = sedmdb.get_from_request(
            values=['id'], where_dict={'marshal_id': request['requestid']})

        # b. if no request id was found then
        # it means the target was not in there
        if not sedm_requestid:
            print("No request matching that id")
            return False

        # c. even though there should only be one request id we have seen
        # cases where more than one request have been made with the same
        # marshal id.  Therefore we delete all request matching the id.
        else:
            for i in sedm_requestid:
                try:
                    print(sedmdb.update_request(
                        {'id': i[0], 'status': 'CANCELED'}))
                    return True
                except Exception as e:
                    print(str(e))
                    return False

    # 1. If we are filling the request by custom dict then we don't need to
    # look up the info.
    shareid = 0
    if fill_by_custom_dict:
        user_id = custom_dict['user_id']
        object_id = custom_dict['object_id']
        priority = custom_dict['priority']
        allocation_id = custom_dict['allocation_id']
        objdict = {'magnitude': 18.5}
    else:
        # 1a. Get the username
        try:
            ret = sedmdb.get_from_users(
                values=['id'],
                where_dict={'username': request['username']})
            if len(ret) > 0:
                user_id = ret[0][0]
                print("Found user: %s" % request['username'])
            else:
                udict = {'username': request['username'],
                         'email': request['email'],
                         'password': datetime.now().strftime(
                             "%Y-%m-%dT%H:%M:%S")}
                ret = sedmdb.add_user(udict)
                user_id = ret[0]
                print("Added user: %s" % request['username'])
        except Exception as e:
            print(str(e))
            print("Using default user!")
            user_id = 254
            # 1b. Create an entry in the database
        objdict = {
            'name': request['sourcename'],
            'ra': request['ra'],
            'dec': request['dec'],
            'typedesig': 'f',
            'epoch': 2000,
            'magnitude': get_prior_mag(request['prior_photometry'])
        }

        object_id, msg = sedmdb.add_object(objdict)

        if object_id == -1:
            print(msg)
            object_id = msg.split()[-1]
        priority = request['priority']

        # TODO: Replace this and only use the programname to add request
        # request['programname'] is the 'name' in the FRITZ group record
        # name is the 'name' field in the SEDM db table 'program'
        # shareid is 2 for ZTF programs and 3 for Caltech programs
        # ZTF programs
        if request['programname'] == 'Asteroids':
            name = 'Astrometric Follow-up of 10 m scale asteroids'
            shareid = 2
        elif request['programname'] == 'Infant Supernovae':
            name = 'Infant Supernovae 2021'
            shareid = 2
        elif request['programname'] == 'Nuclear Transients':
            name = 'A complete spectroscopic sample of bright TDEs'
            shareid = 2
        elif request['programname'] == 'Sollerman Research Group':
            name = 'An unbiased stripped-envelope supernova sample from ZTF'
            shareid = 2
        elif request['programname'] == 'SEDM Team':
            name = 'SEDM Czar Discretionary Time'
            shareid = 2
        elif request['programname'] == 'Type Ia Supernovae':
            name = 'SEDm follow-up of Type Ia supernovae with early spectra'
            shareid = 2
        elif request['programname'] == 'Neutrino follow-up':
            name = 'Multi-Messenger Astronomy'
            shareid = 2
        elif request['programname'] == 'Lensed SNe':
            name = 'Identifying gravitationally lensed supernovae'
            shareid = 2
        elif request['programname'] == 'Cataclysmic Variables':
            name = 'Identification Spectra of Eruptive variables and transients in our galaxy'
            shareid = 2
        # Caltech programs
        elif request['programname'] == 'Redshift Completeness Factor':
            name = 'The ZTF-II Bright Transient Survey Caltech'
            shareid = 3
        elif request['programname'] == 'Red Transients':
            name = 'The ZTF-II Bright Transient Survey Caltech'
            shareid = 3
        elif request['programname'] == 'Census of the Local Universe Caltech':
            name = 'The ZTF-II Bright Transient Survey Caltech'
            shareid = 3
        elif request['programname'] == 'Palomar Gattini-IR':
            name = 'Exploring the dynamic infrared sky with Palomar Gattini-IR'
            shareid = 3
        elif request['programname'] == 'AMCVn ':
            name = 'Finding AM CVn systems in alert-streams'
            shareid = 3
        elif request['programname'] == 'cyclotron WDs':
            name = 'Director Discretionary Time'
            shareid = 3
        # elif request['programname'] == 'Electromagnetic Counterparts to Gravitational Waves':
        #    name = 'Gravitational Wave Follow-Up'
        #    shareid = 3     # inactive
        else:
            name = 'SEDm program'
            shareid = 0

        program_id = sedmdb.get_from_program(values=['id'],
                                             where_dict={'name': name})[0][0]

        try:
            allocation_id = sedmdb.get_from_allocation(
                values=['id'],
                where_dict={'program_id': program_id, 'active': True})
            print(allocation_id)
            allocation_id = allocation_id[0][0]
        except Exception as e:
            print(str(e))
            # TODO add an email alert indicating an unknown
            # program has been added
            allocation_id = 1

    start_date = request['startdate']
    end_date = request['enddate']

    follow_up = "%s|%s" % (request['Followup'], request['Filters'])
    obs_seq, exptime = get_observing_sequence(follow_up, objdict['magnitude'],float(request['exptime']))

    if 'origins_url' in request:
        print(request['origins_url'])
        if 'private' in request['origins_url']:
            external_id = 2
        elif 'skipper' in request['origins_url']:
            external_id = 3
        else:
            external_id = 0
    else:
        external_id = 0

    request_dict = {
        'user_id': int(user_id),
        'object_id': int(object_id),
        'priority': float(priority),
        'allocation_id': int(allocation_id),
        'exptime': exptime,
        'obs_seq': obs_seq,
        'inidate': start_date,
        'enddate': end_date,
        'external_id': 2,
        'shareid': shareid,
        'marshal_id': int(request['requestid']),
        'maxairmass': float(request['maxairmass']),
        'max_fwhm': float(request['max_fwhm'])
    }
    print(external_id)
    print(request_dict)
    print(sedmdb.add_request(request_dict))

    return True


def get_prior_mag(mag_dict):
    """
    Get the last prior magnitude
    :param mag_dict:
    :return:
    """
    print("GETTING MAGNITUDE")
    print(mag_dict)
    print(type(mag_dict))
    mag = 17
    if not isinstance(mag_dict, dict):
        print("Not a dictionary so using mag=17")
        return mag

    for k, v in mag_dict.items():
        mag = v['mag']

    try:
        mag = float(mag)
    except Exception as e:
        print(str(e), "Error getting magnitude")
        mag = 17
    print(mag)
    return mag


def get_observing_sequence(obs_str, mag=17, exptime=0):
    """
    Get the observing sequence from the predefined followups and or filters

    :param obs_str:
    :param mag:
    :param exptime:
    :return:
    """
    sequence, filters = obs_str.split('|', 1)

    sequence_list = []
    exptime_list = []

    if sequence:
        print("Found sequence")
        if sequence == 'IFU':
            sequence_list.append('1ifu')
            if exptime <= 0:
                exptime_list.append(get_exptime(obsfilter='ifu', mag=mag))
            else:
                exptime_list.append(str(exptime))
        elif sequence == 'Four Shot (r,g,i,u)':
            obs_list = ['1r', '1g', '1i', '1u']
            for o in obs_list:
                sequence_list.append(o)
                # if exptime <= 0:
                exptime_list.append(get_exptime(obsfilter=o[1], mag=mag))
                # else:
                #    exptime_list.append(str(exptime))
        elif sequence == 'Three Shot (r,g,i)':
            obs_list = ['1r', '1g', '1i']
            for o in obs_list:
                sequence_list.append(o)
                # if exptime <= 0:
                exptime_list.append(get_exptime(obsfilter=o[1], mag=mag))
                # else:
                #    exptime_list.append(str(exptime))
        elif sequence == 'Fourshot + IFU':
            obs_list = ['1r', '1g', '1i', '1u']
            sequence_list.append('1ifu')
            if exptime <= 0:
                exptime_list.append(get_exptime(obsfilter='ifu', mag=mag))
            else:
                exptime_list.append(str(exptime))
            for o in obs_list:
                sequence_list.append(o)
                # if exptime <= 0:
                exptime_list.append(get_exptime(obsfilter=o[1], mag=mag))
                # else:
                #    exptime_list.append(str(exptime))

    if filters:
        filters = filters.split(',')
        for f in filters:
            sequence_list.append('1%s' % f)
            # if exptime <= 0:
            exptime_list.append(get_exptime(obsfilter=f.lower(), mag=mag))
            # else:
            #    exptime_list.append(str(exptime))

    print(sequence_list)
    print(exptime_list)

    sequence_str = '{%s}' % ','.join(sequence_list)
    exptime_str = '{%s}' % ','.join(exptime_list)

    return sequence_str, exptime_str


def get_exptime(obsfilter, mag):
    """
    Get exposure times for filters
    :param obsfilter:
    :param mag:
    :return:
    """

    mag = float(mag)
    print("Finding exptime for filter:%s at magnitude:%.2f" % (obsfilter, mag))
    if mag > 18:
        ifu_exptime = 2250
        r_exptime = 180
        g_exptime = 180
        i_exptime = 180
        u_exptime = 300
    elif 15 > mag < 18:
        ifu_exptime = 2000
        r_exptime = 120
        g_exptime = 120
        i_exptime = 120
        u_exptime = 300
    elif 0 < mag < 5:
        ifu_exptime = 60
        r_exptime = 1
        g_exptime = 1
        i_exptime = 1
        u_exptime = 30
    elif 5 < mag < 10:
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
        print("No valid filter")
        return str(0)


if __name__ == "__main__":
    request_test = "/home/rsw/PycharmProjects/sedmpy/growth/archived/test.json"
    request_test_dict = interface.read_request(request_test)
    # ret = process_new_request(request, isfile=True,
    # check_rejection=True, add2db=True)
    print(delete_request_entry(request_test_dict))

    """
    # 4. Check if we are adding it to the SEDm DB if not continue
    if add2db:
        print("Adding request to database")
        #print(add_request_to_db(req_dict))
        print("Done adding request to database")
        pass

    # 4. Assuming it passed the above conditions send back the response
    ret = growth.update_request(request_id=req_dict['requestid'], status=status)
    print(ret.text)

    # 5. If the update was successful add it to the archive
    if archive:
        current_archive = archive_dir + request_date
        if not os.path.exists(current_archive):
            os.system('mkdir %s' % current_archive)
        os.system('mv %s %s' % (request, current_archive))
    return ret"""
