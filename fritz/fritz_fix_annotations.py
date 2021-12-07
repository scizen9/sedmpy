import os
import json
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

DEFAULT_TIMEOUT = 5  # seconds

# URL constants
fritz_base_url = 'https://fritz.science/api/'

SITE_ROOT = os.path.abspath(os.path.dirname(__file__)+'/../..')

with open(os.path.join(SITE_ROOT, 'sedmpy', 'marshals', 'config',
                       'marshals.json')) as data_file:
    params = json.load(data_file)

token = params['marshals']['fritz']['token']


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
    method_whitelist=["HEAD", "GET", "PUT", "POST", "PATCH", "DELETE"]
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
    ret = response.json()

    return ret


def get_source_autoannot(obj_id, testing=False):
    """
    gets autoannotations

    obj_id: 'ZTF18aaaaaa', for example
    testing: <bool> are we testing only?

    return: True if success, False if not
    """

    aids = []
    adat = []
    aorg = []
    spid = []
    if testing:
        print("TESTING get_source_autoannot(): no data sent to marshal")
    else:
        fritz_annotation_url = fritz_base_url + \
                               'sources/%s/annotations' % obj_id
        r = api("GET", fritz_annotation_url)
        if 'success' in r['status']:
            for ann in r['data']:
                if 'SNIascore:' in ann['origin'] or 'sedm:' in ann['origin']:
                    aids.append(ann['id'])
                    adat.append(ann['data'])
                    spid.append(int(ann['origin'].split(':')[-1].split(
                        'spc')[-1]))
                    if 'SNIascore:' in ann['origin']:
                        aorg.append('sedm:SNIascore')
                    else:
                        aorg.append('sedm:SNID')
        else:
            print('error getting annotations for %s' % obj_id)
            print(r['status'])
            print(r['message'])

    return aids, spid, aorg, adat


def delete_annotation(annot_id, res_type='sources', res_id=None,
                      testing=False):
    if 'sources' in res_type:
        fritz_delete_ann_url = fritz_base_url + '%s/%s/annotations/%d' % (
            res_type, res_id, annot_id)
    else:
        fritz_delete_ann_url = fritz_base_url + '%s/%d/annotations/%d' % (
            res_type, res_id, annot_id)

    if testing:
        print(fritz_delete_ann_url)
    else:
        r = api("DELETE", fritz_delete_ann_url)
        if 'success' in r['status']:
            print("Annotation at %s successfully deleted!" %
                  fritz_delete_ann_url)
        else:
            print("error deleting annotation at %s" % fritz_delete_ann_url)
            print(r['status'])
            print(r['message'])


def add_spec_autoannot(obj_id, andic, spec_id=None, origin=None, testing=False):
    """
    adds an autoannotation without attachment to a particular SEDM spectrum
    on the view_spec page, which will also appear elsewhere.

    obj_id: 'ZTF18aaaaaa', for example
    andic: <dict> of annotations
    spec_id: <int>
    origin: <str> giving origin of annotations
    testing: <bool> are we testing only?

    return: True if success, False if not
    """

    ddict = {'obj_id': obj_id,  # This goes
             'origin': origin,
             'data': andic}

    new_anid = -1
    if testing:
        print("TESTING add_spec_autoannot(): no data sent to marshal")
        print(ddict)
    else:
        fritz_annotation_url = fritz_base_url + \
                               'spectra/%d/annotations' % spec_id
        r = api("POST", fritz_annotation_url, data=ddict)
        if 'success' in r['status']:
            r_data = r['data']
            if 'annotation_id' in r_data:
                new_anid = int(r_data['annotation_id'])
                print("annotation id = %d" % new_anid)
            print('{}: {} posted'.format(obj_id, origin))
        else:
            print('error submitting comment')
            print(r['status'])
            print(r['message'])

    return new_anid


def verify_spec_annot(spec_id, andic, new_anid, verbose=False):
    fritz_annotation_url = fritz_base_url + \
                           'spectra/%d/annotations/%d' % (spec_id, new_anid)
    r = api("GET", fritz_annotation_url)
    verified = False
    if 'success' in r['status']:
        new_dict = r['data']['data']
        num_to_verify = len(andic)
        for k, v in andic.items():
            if k in new_dict:
                if v == new_dict[k]:
                    num_to_verify -= 1
                    if verbose:
                        print("%s verified as " % k, v)
        verified = (num_to_verify <= 0)
    else:
        print("error verifying annotation %d" % new_anid)
        print(r['status'])
        print(r['message'])

    return verified


def fix_annotations(source_id, testing=False):

    anids, spids, anorig, andicts = get_source_autoannot(source_id,
                                                         testing=testing)
    for i, ani in enumerate(anids):
        new_anid = add_spec_autoannot(source_id, andicts[i], spec_id=spids[i],
                                      origin=anorig[i], testing=testing)
        if not testing:
            if verify_spec_annot(spids[i], andicts[i], new_anid):
                delete_annotation(ani, res_id=source_id)
            else:
                print("Could not verify annotation %d" % new_anid)


if __name__ == "__main__":
    print("here")
    anids, anspec, anorig, andata = get_source_autoannot('ZTF21abtsoky')
    print(anids)
    print(anspec)
    print(anorig)
    print(andata)
    # print("%d annotations found")
