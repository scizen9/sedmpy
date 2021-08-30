import json
import os
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

DEFAULT_TIMEOUT = 5  # seconds

# URL constants
fritz_base_url = 'https://fritz.science/'
fritz_alloc_url = fritz_base_url + 'api/allocation'

default_id = 37
instrument_id = 2
telescope_id = 37

SITE_ROOT = os.path.abspath(os.path.dirname(__file__)+'/../..')

with open(os.path.join(SITE_ROOT, 'sedmpy', 'marshals', 'config',
                       'marshals.json')) as data_file:
    params = json.load(data_file)

token = params['marshals']['fritz']['alloc_token']


class TimeoutHTTPAdapter(HTTPAdapter):
    def __init__(self, *args, **kwargs):
        self.timeout = DEFAULT_TIMEOUT
        if "timeout" in kwargs:
            self.timeout = kwargs["timeout"]
            del kwargs["timeout"]
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        timeout = kwargs.get("timeout")
        if timeout is None:
            kwargs["timeout"] = self.timeout
        return super().send(request, **kwargs)


session = requests.Session()
session_headers = {'Authorization': 'token {}'.format(token)}
retries = Retry(
    total=5,
    backoff_factor=2,
    status_forcelist=[405, 429, 500, 502, 503, 504],
    method_whitelist=["HEAD", "GET", "PUT", "POST", "PATCH"]
)
adapter = TimeoutHTTPAdapter(timeout=5, max_retries=retries)
session.mount("https://", adapter)
session.mount("http://", adapter)


def api(method, endpoint, data=None, verbose=False):
    headers = {'Authorization': 'token {}'.format(token)}
    response = session.request(method, endpoint, json=data, headers=headers)
    print('HTTP code: {}, {}'.format(response.status_code, response.reason))
    if response.status_code in (200, 400) and verbose:
        print(response.text)
        # print('JSON response: {}'.format(response.json()))

    return response


def create_allocation(pi=None, proposal_id=None,
                      start_date=None, end_date=None, hours_allocated=None,
                      group_id=None, inst_id=2, testing=False):
    """
    Create an allocation on the fritz marshal

    :param pi:
    :param proposal_id:
    :param start_date:
    :param end_date:
    :param hours_allocated:
    :param group_id:
    :param inst_id:
    :param testing:

    :return:
    """
    if not hours_allocated or not group_id or not inst_id:
        print("Must have hours allocated, group id, and instrument id")
        return None

    submission_dict = {}

    # create payload
    submission_dict.update(
        {'hours_allocated': hours_allocated,
         'group_id': group_id,
         'instrument_id': inst_id}
    )
    if pi is not None:
        submission_dict['pi'] = pi
    if proposal_id is not None:
        submission_dict['proposal_id'] = proposal_id
    if start_date is not None:
        submission_dict['start_date'] = start_date
    if end_date is not None:
        submission_dict['end_date'] = end_date
    # Are we just testing?
    if testing:
        print(fritz_alloc_url)
        print(submission_dict)
        ret = {"message": "string", "status": "success"}
        return ret
    else:
        # post the new allocation
        try:
            ret = api("POST", fritz_alloc_url, data=submission_dict)
        except requests.exceptions.ConnectionError:
            ret = {
                'status': 'Error', 'message': 'ConnectionError',
                'data': None
            }

    return ret.json()


def update_allocation(alloc_id, pi=None, proposal_id=None,
                      start_date=None, end_date=None, hours_allocated=None,
                      group_id=None, inst_id=2, testing=False):
    """
    Update an allocation on the fritz marshal

    :param alloc_id:
    :param pi:
    :param proposal_id:
    :param start_date:
    :param end_date:
    :param hours_allocated:
    :param group_id:
    :param inst_id:
    :param testing:

    :return:
    """
    if not alloc_id or not hours_allocated or not group_id or not inst_id:
        print("Must have allocation id, hours allocated, group id,"
              " and instrument id")
        return None

    alloc_url = fritz_alloc_url + "/%d" % alloc_id

    submission_dict = {}

    # create payload
    submission_dict.update(
        {'hours_allocated': hours_allocated,
         'group_id': group_id,
         'instrument_id': inst_id}
    )
    if pi is not None:
        submission_dict['pi'] = pi
    if proposal_id is not None:
        submission_dict['proposal_id'] = proposal_id
    if start_date is not None:
        submission_dict['start_date'] = start_date
    if end_date is not None:
        submission_dict['end_date'] = end_date
    # Are we just testing?
    if testing:
        print(alloc_url)
        print(submission_dict)
        ret = {"message": "string", "status": "success"}
        return ret
    else:
        # put the update
        try:
            ret = api("PUT", alloc_url, data=submission_dict)
        except requests.exceptions.ConnectionError:
            ret = {
                'status': 'Error', 'message': 'ConnectionError',
                'data': None
            }

    return ret.json()
