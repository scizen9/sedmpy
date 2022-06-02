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

Download the file from the web url:

https://fritz.science/api/followup_request

and save the output to the file followup_request.json

in /data/sedmdrp/fritz

cd into that directory and run:

~/spy ~/sedmpy/fritz/fritz_requests_bulk_expire.py
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to marshal (for testing)')
    parser.add_argument('--filename', type=str, default='followup_request.json',
                        help='Followup request JSON file')

    args = parser.parse_args()

    print("Getting all the Fritz requests from %s" % args.filename)
    reqs = get_all_fritz_reqests(args.filename)

    today = datetime.date.today()

    n_expired = 0
    for r in reqs:
        if 'submitted' in r['status']:
            pay = r['payload']
            end_date = pay['end_date']
            if datetime.date(*[int(du) for du in end_date.split('-')]) < today:
                print("Expiring outdated (%s) request for %s" %
                      (end_date, r['obj_id']))
                marshal_req_id = r['id']
                res = update_status_request('Expired', marshal_req_id, 'fritz',
                                            testing=args.testing)
                if 'success' in res['status']:
                    n_expired += 1
                else:
                    print(res)

    print("Expired %d outdated requests" % n_expired)
