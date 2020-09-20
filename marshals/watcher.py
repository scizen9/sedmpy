from marshals import interface


def process_new_request(request, isfile=True, status='ACCEPTED', add2db=False,
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

if __name__ == "__main__":
    request = "/home/rsw/PycharmProjects/sedmpy/growth/archived/test.json"
    ret = process_new_request(request, isfile=True, check_rejection=True)
    print(ret)

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
