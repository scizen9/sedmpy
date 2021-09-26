import requests
import argparse

from marshals.interface import api, update_status_request
import db.SedmDb

# URL constants
fritz_base_url = 'https://fritz.science/'
fritz_view_source_url = fritz_base_url + 'source'
fritz_req_url = fritz_base_url + 'api/followup_request/'

# Use the database
sedmdb = db.SedmDb.SedmDB()


def get_fritz_req_status(marshal_id):
    """Query fritz for followup request status"""
    try:
        res = api('GET', fritz_req_url+str(marshal_id)).json()
    except (requests.exceptions.ConnectionError,
            requests.exceptions.RetryError):
        res = {'status': 'Error', 'message': 'ConnectionError', 'data': None}

    if 'success' in res['status']:
        data = res['data']
    else:
        data = None

    return data


def ifu_present(sequence):
    ret = False
    for item in sequence:
        if 'ifu' in item:
            ret = True
    return ret


def update_fritz_status(request_id=None, status='COMPLETED',
                        ifu=False, testing=False):

    if ifu:
        chan = 'ifu'
    else:
        chan = 'rc'
    # Return values
    image_ret = None
    status_ret = None
    return_link = None
    if request_id:
        print("Searching SedmDB")

        # Search for target in the database
        try:
            res = sedmdb.get_from_request(["marshal_id",
                                           "object_id",
                                           "user_id",
                                           "external_id",
                                           "obs_seq"],
                                          {"id": request_id})[0]
        except IndexError:
            print("Unable to retrieve ids from database")
            return return_link, image_ret, status_ret
        marshal_id = res[0]
        object_id = res[1]
        user_id = res[2]
        external_id = res[3]
        obs_seq = res[4]
        # does ifu status override?
        if ifu_present(obs_seq) and not ifu:
            print("IFU status overrides RC status!")
            return return_link, image_ret, status_ret
        if not ifu_present(obs_seq) and ifu:
            print("No IFU obs in this request!")
            return return_link, image_ret, status_ret
        # is this a Fritz object?
        if external_id != 2 and external_id != 4:
            print("Not a Fritz object!")
            return return_link, image_ret, status_ret
        else:
            if external_id == 4:
                print("AMPEL trigger")
            else:
                print("Fritz trigger")
        # get source name
        try:
            res = sedmdb.get_from_object(["name"], {"id": object_id})[0]
        except IndexError:
            print("Unable to retrieve object_name from database")
            return return_link, image_ret, status_ret
        object_name = res[0]
        # get user name and email
        try:
            res = sedmdb.get_from_users(["name", "email"],
                                        {"id": user_id})[0]
        except IndexError:
            print("Unable to retrieve username, email from database")
            return return_link, image_ret, status_ret
        username = res[0]
        email = res[1]
    else:
        print("no request_id given!")
        return return_link, image_ret, status_ret

    # Did we get a marshal ID?
    if marshal_id is None:
        print("Unable to find marshal id for target %s" % object_name)
    else:
        print("Updating %s request status for %s using marshal id %d" %
              (chan, object_name, marshal_id))
        try:
            status_ret = update_status_request(status, marshal_id, 'fritz',
                                               testing=testing)
        except requests.exceptions.ConnectionError:
            status_ret = None

        return_link = fritz_view_source_url + "/%s" % object_name

        print("Send to %s at %s\nRequest status = %s\n%s" %
              (username, email, status, return_link))

    return return_link, image_ret, status_ret


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
                         
Updates IFU/RC results to the fritz marshal.
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to marshal (for testing)')
    parser.add_argument('-r', '--request_id', type=int, default=None,
                        help="Request ID to update (int)")
    parser.add_argument('-s', '--status', type=str, default=None,
                        help="New status for request (str)")
    parser.add_argument('-i', '--ifu', action="store_true", default=False,
                        help='Update IFU status on fritz (else update RC)')
    args = parser.parse_args()

    # resu = sedmdb.get_from_request(["id", "obs_seq"],
    #                               {"status": "COMPLETED", "external_id": 2})

    # print("Found %d rows" % len(resu))

    # for r in resu:
    #    if ifu_present(r[1]):
    #        continue
    #    update_rc_status(r[0], testing=args.testing)

    if args.request_id is not None:
        if args.status is not None:
            update_fritz_status(args.request_id, status=args.status,
                                ifu=args.ifu)
        else:
            update_fritz_status(args.request_id, ifu=args.ifu)
    else:
        print("Must specify request id!")
