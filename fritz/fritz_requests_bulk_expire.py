import requests
import argparse

from marshals.interface import api, update_status_request
import db.SedmDb

import datetime

# URL constants
fritz_base_url = 'https://fritz.science/'
fritz_view_source_url = fritz_base_url + 'source'
fritz_req_url = fritz_base_url + 'api/followup_request/'

# Use the database
sedmdb = db.SedmDb.SedmDB()


def get_all_fritz_reqests():
    """Query fritz for followup requests"""
    try:
        res = api('GET', fritz_req_url).json()
    except requests.exceptions.ConnectionError:
        res = {'status': 'Error', 'message': 'ConnectionError', 'data': None}

    if 'success' in res['status']:
        data = res['data']
    else:
        data = None

    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
                         
Expires outdated requests on the fritz marshal.
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to marshal (for testing)')
    args = parser.parse_args()

    print("Getting all the requests from Fritz...")
    reqs = get_all_fritz_reqests()
    print("Done")

    for r in reqs:
        if 'submitted' in r['status']:
            pay = r['payload']
            end_date = pay['end_date']
            if datetime.date(*[int(du) for du in
                               end_date.split('-')]) < datetime.date.today():
                print("Expiring outdated request for %s" % r['obj_id'])
                marshal_req_id = r['id']
                res = update_status_request('Expired', marshal_req_id, 'fritz',
                                            testing=True)
                break
