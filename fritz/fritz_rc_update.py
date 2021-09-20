import requests
import argparse

from marshals.interface import update_status_request
import db.SedmDb

# Path constants
pharos_spec_dir = '/scr2/sedmdrp/redux/'
pharos_phot_dir = '/scr2/sedmrp/redux/phot/'
# URL constants
fritz_base_url = 'https://fritz.science/'
fritz_spec_url = fritz_base_url + 'api/spectrum/ascii'
fritz_phot_url = fritz_base_url + 'api/photometry'
fritz_view_source_url = fritz_base_url + 'source'
fritz_alloc_url = fritz_base_url + 'api/allocation'

default_id = 37
instrument_id = 2
telescope_id = 37

# Use the database

sedmdb = db.SedmDb.SedmDB()


def ifu_present(sequence):
    ret = False
    for item in sequence:
        if 'ifu' in item:
            ret = True
    return ret


def update_rc_status(request_id=None, status='COMPLETED', testing=False):

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
        if ifu_present(obs_seq):
            print("IFU status overrides RC status!")
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
        print("Updating rc request status for %s using marshal id %d" %
              (object_name, marshal_id))

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
                         
Updates RC results to the fritz marshal.
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to marshal (for testing)')
    args = parser.parse_args()

    rows = sedmdb.get_from_request(["request_id", "obs_seq"],
                                   {"status": "COMPLETED", "external_id": 2})[0]

    rids = rows[0]
    seq_list = rows[1]
    print("Found %d rows" % len(rids))

    for i, seq in enumerate(seq_list):
        if ifu_present(seq):
            continue
        update_rc_status(rids[i], testing=args.testing)
