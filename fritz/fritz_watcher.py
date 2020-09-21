import glob
import os
import time
import datetime
import db.SedmDb
import fritz.fritz as fritz

sedmdb = db.SedmDb.SedmDB()


def check_for_new_request(request_path="requests/", request_date="",
                          archived_dir="archived/"):
    """
    
    :return: 
    """
    # 1. Check for a directory path, if not set then use the default path
    if not request_path:
        request_path = 'requests/'

    # 2. Determine if we are only looking for request on a certain date.  If
    # not then get all the requests
    if request_date:
        requests = sorted(glob.glob(request_path + '*%s*' % request_date))
    else:
        requests = glob.glob(request_path + '*')

    # 3. If there are no request then we are done and can exit
    if len(requests) == 0:
        print("No new requests have been sent")
        return False

    # 4. If there are request then we need to check and see if they are new
    # or are parts of an older request set by looking through the archive or DB

    # TODO: Look through the SEDm DB
    request_list = []
    prev_req_list = [os.path.basename(x)
                     for x in glob.glob(archived_dir + "%s/*" % request_date)]

    for req in requests:
        
        if os.path.basename(req) not in prev_req_list:
            request_list.append(req)

    # 5. Finally send back a list of all new requests
    if len(request_list) != 0:
        return request_list
    else:
        return False


def pull_request_from_remote(site='rsw@nera.palomar.caltech.edu',
                             remote_dir='/tmp/', remote_files='request',
                             local_dir='requests/'):
    """
    
    :param site: 
    :param remote_dir: 
    :param remote_files: 
    :param local_dir: 
    :return: 
    """
    remote_cmd = "scp -q %s:%s%s* %s" % (site, remote_dir, remote_files,
                                         local_dir)

    os.system(remote_cmd)


def process_new_request(request, status='ACCEPTED', add2db=False,
                        check_rejection=False, archive=True,
                        archive_dir='archived/', request_date=''):
    """
    
    :param add2db: 
    :param check_rejection: 
    :param archive: 
    :param archive_dir: 
    :param request_date: 
    :param request: request string path
    :param status: 
    :return: 
    """
    # 1. Open the request
    req = fritz.read_request(request)

    # 2. Check any rejection criteria
    # TODO: Check database to determine if group has allocation time
    if check_rejection:
        pass

    # 3. Check if we are adding it to the SEDm DB if not continue
    if add2db:

        print("Adding request to database")
        print(add_request_to_db(req))
        print("Done adding request to database")
        pass

    # 4. Assuming it passed the above conditions send back the response
    ret = fritz.update_request(request_id=req['requestid'], status=status)
    print(ret.text)

    # 5. If the update was successful add it to the archive
    if archive:
        current_archive = archive_dir + request_date
        if not os.path.exists(current_archive):
            os.system('mkdir %s' % current_archive)
        os.system('mv %s %s' % (request, current_archive))
    return ret


def add_request_to_db(request, custom_dict=None, fill_by_custom_dict=False):
    """
    
    :param request:
    :param custom_dict:
    :param fill_by_custom_dict:
    :return: 
    """

    # 0. If delete request, cancel it from the database if it has an entry
    if request['status'].lower() == 'delete':
        # a. start by trying to get the request id from the marshal id
        sedm_requestid = sedmdb.get_from_request(
            values=['id'], where_dict={'marshal_id': request['requestid']})
         
        if not sedm_requestid:
            print("No request matching that id")
            return False
        else:
            
            for i in sedm_requestid:
                try:
                    print(sedmdb.update_request(
                        {'id': i[0], 'status': 'CANCELED'}))
                    return True
                except Exception as e:
                    print(str(e))
                    return False
    # TODO ADD UPDATE CHECK

    # 1. If we are filling the request by custom dict then we don't need to
    # look up the info.
    if fill_by_custom_dict:
        user_id = custom_dict['user_id']
        object_id = custom_dict['object_id']
        priority = custom_dict['priority']
        allocation_id = custom_dict['allocation_id']

    else:
        # 1a. Get the username

        try:        
            user_id = sedmdb.get_from_users(
                values=['id'],
                where_dict={'username': request['username']})[0][0]
        except Exception:
            user_id = 189                                            
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
            print("TARGET IN DB ALREADY")
            print(msg)
            object_id = sedmdb.get_from_object(
                values=['id'],
                where_dict={'name': request['sourcename'].lower()})[0][0]
        priority = request['priority']
        print(priority)
        # TODO: Replace this and only use the programname to add request
        if request['programname'] == 'ZTF Science Validation':
            name = 'GROWTH classification program'
        else:
            name = 'ZTF commissioning'

        program_id = sedmdb.get_from_program(
            values=['id'], where_dict={'name': name})[0][0]

        allocation_id = sedmdb.get_from_allocation(
            values=['id'],
            where_dict={'program_id': program_id, 'active': True})[0][0]

        print(allocation_id)
    start_date = request['startdate']
    end_date = request['enddate']

    follow_up = "%s|%s" % (request['Followup'], request['Filters'])
    obs_seq, exptime = get_observing_sequence(follow_up, objdict['magnitude'])

    request_dict = {
        'user_id': int(user_id),
        'object_id': int(object_id),
        'priority': float(priority),
        'allocation_id': int(allocation_id),
        'exptime': exptime,
        'obs_seq': obs_seq,
        'inidate': start_date,
        'enddate': end_date,
        'marshal_id': int(request['requestid'])
    }
    print(request_dict)
    print(sedmdb.add_request(request_dict))


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


def get_observing_sequence(obs_str, mag=17):
    """
    Get the observing sequence from the predefined followups and or filters
    
    :param obs_str:
    :param mag:
    :return: 
    """
    print(obs_str)
    sequence, filters = obs_str.split('|', 1)
    print(sequence)
    print(filters)
    sequence_list = []
    exptime_list = []

    if sequence:
        print("Found sequence")
        if sequence == 'IFU':
            sequence_list.append('1ifu')
            exptime_list.append(get_exptime(obsfilter='ifu', mag=mag))
        elif sequence == 'Four Shot (r,g,i,u)':
            obs_list = ['1r', '1g', '1i', '1u']
            for o in obs_list:
                sequence_list.append(o)
                exptime_list.append(get_exptime(obsfilter=o[1], mag=mag))
        elif sequence == 'Three Shot (r,g,i)':
            obs_list = ['1r', '1g', '1i']
            for o in obs_list:
                print(o)
                sequence_list.append(o)
                exptime_list.append(get_exptime(obsfilter=o[1], mag=mag))
        elif sequence == 'Fourshot + IFU':
            obs_list = ['1r', '1g', '1i', '1u']
            sequence_list.append('1ifu')
            exptime_list.append(get_exptime(obsfilter='ifu', mag=mag))
            for o in obs_list:
                sequence_list.append(o)
                exptime_list.append(get_exptime(obsfilter=o[1], mag=mag))

    if filters:
        filters = filters.split(',')
        for f in filters:
            sequence_list.append('1%s' % f)
            exptime_list.append(get_exptime(obsfilter=f.lower(), mag=mag))

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
    print("Finding exptime for filter:%s at magnitude:%s" % (obsfilter, mag))
    if mag > 18:
        ifu_exptime = 3600
        r_exptime = 180
        g_exptime = 180
        i_exptime = 180
        u_exptime = 300
    elif 15 > mag < 18:
        ifu_exptime = 2700
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


def start_watcher():
    """
    Check for new targets and update database
    :return: 
    """
    while True:
        request_date = datetime.datetime.utcnow().strftime("%Y%m%d")
        pull_request_from_remote(remote_files="*%s*" % request_date)
        new_requests = check_for_new_request(request_date=request_date)
        if not new_requests:
            time.sleep(5)
            continue

        # noinspection PyTypeChecker
        for r in new_requests:
            print("Processing %s" % r)
            try:
                ret = process_new_request(r, request_date=request_date,
                                          add2db=True)
                print(ret)
            except:
                os.system('cp -r %s /home/sedm/fritz_marshal/archived/failed/'
                          % r)
                os.system('cp -r %s /home/sedm/fritz_marshal/archived/%s/' %
                          (r, request_date))

        print("Waiting %ss before checking for new request" % 5)
        time.sleep(5)


if __name__ == '__main__':
    start_watcher()
