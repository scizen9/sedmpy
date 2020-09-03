import os
import json

SITE_ROOT = os.path.abspath(os.path.dirname(__file__)+'/../..')

with open(os.path.join(SITE_ROOT, 'config', 'marshals.json')) as data_file:
    params = json.load(data_file)

print(params)


def read_request(request, isfile=True):
    """
    Read an incoming request from external marshal.  All request are expected
    to be in a json format that can be converted to a python dictionary

    :param isfile: Is the request a file
    :param request: Expected type is pathfile,
    :return: dictionary of request
    """

    # 1. If isfile is true then check for the file and assume it is
    # in json format
    if isfile and os.path.exists(request):
        try:
            with open(request) as f:
                output = json.load(f)
        except Exception as e:
            output = {'iserror': True,
                      'msg': str(e)}
        return output

    # 2. If it is not a file but a string then assume it is json format
    # and read it in
    if not isfile and isinstance(request, str):
        try:
            output = json.loads(request)
        except Exception as e:
            output = {'iserror': True,
                      'msg': str(e)}
        return output

    # 3. If it is not a file or string but instead a dictionary then just go
    # ahead and return it back
    if not isfile and isinstance(request, dict):
        return request

    # 4. Finally if none of those work then just send back an error message
    return {'iserror': True,
            'msg': "Request is in an unknown format:%s" % type(request)}


def check_field(key, input_dict):
    if key not in input_dict:
        return "Missing key: %s in request" % key
    if not input_dict[key]:
        return "Key: %s is empty or null" % key
    return False


def checker(request, check_source=True, check_followup=True,
            check_program=False, check_token=False, check_user=True,
            check_email=True, check_status=True, check_id=True,
            check_dates=True):

    """
    The following fields are considered to be vital for ingest into the
    SEDm database.  Any fields that are missing or does not pass the
    constraints should result in an invalid request.

    :param request:
    :param check_id:
    :param check_dates:
    :param check_source: Check if there is a source in the request
    :param check_followup: Check if there is a followup type or filters listed
    :param check_program: Check that there is a valid follow-up program
    :param check_token: Check that token matches url sender
    :param check_user: Check if a valid user is in the request
    :param check_email: Check if a valid email is in the request
    :param check_status: Check that we have a status of new, edit, delete
    :return:
    """
    # Make sure all dictionary keys are lowercase
    if not isinstance(request, dict):
        return {"iserror": True, "msg": "Can't check request when not in py dict form"}

    request = dict((k.lower(), v) for k, v in request.items())
    msg_list = []

    # 1. Check that this request comes from a valid source
    if check_source:
        if "origins_url" in request:
            url = request["orgins_url"]
            if not any(x in url for x in params['sources']):
                msg_list.append("Source of request '%s' not regonized" % url)
        else:
            msg_list.append("origins_url missing from request")

    if check_followup:
        followup = False
        filters = False

        if 'followup' in request and request['followup']:
            followup = True
        if 'filters' in request:
            filters = True

        if not followup and not filters:
            msg_list.append("No followup or filters found in request")

    if check_program:
        ret = check_field('programname', request)
        if ret:
            msg_list.append(ret)
        elif request['programname'] not in params['programs']:
            msg_list.append("Progam: %s is not a valid program")

    if check_token:
        ret = check_field('token', request)
        if ret:
            msg_list.append(ret)

    if check_user:
        ret = check_field('username', request)
        if ret:
            msg_list.append(ret)

    if check_email:
        ret = check_field('email', request)
        if ret:
            msg_list.append(ret)

    if check_status:
        ret = check_field('status', request)
        if ret:
            msg_list.append(ret)
        elif request['status'].lower() not in params['status']:
            msg_list.append("Status: %s is not recognized" % request['status'])

    if check_id:
        ret = check_field('requestid', request)
        if ret:
            msg_list.append(ret)

    if check_dates:
        ret = check_field('startdate', request)
        if ret:
            msg_list.append(ret)
        ret = check_field('enddate', request)
        if ret:
            msg_list.append(ret)

    if msg_list:
        return {'iserror': True, 'msg': msg_list}
    else:
        return True

