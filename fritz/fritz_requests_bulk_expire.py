import requests
import argparse

from marshals.interface import api, update_status_request
import json

import datetime


def get_all_fritz_reqests(fnam):
    f = open(fnam)

    data = json.load(f)

    if 'success' in data['status']:
        ret = data['data']
    else:
        print(data['message'])
        ret = None

    return ret


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
                         
Expires outdated requests on the fritz marshal.
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to marshal (for testing)')
    parser.add_argument('--filename', type=str, default='followup_request.json',
                        help='Followup request JSON file')

    args = parser.parse_args()

    print("Getting all the requests from Fritz...")
    reqs = get_all_fritz_reqests(args.filename)
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
                                            testing=args.testing)
                break
