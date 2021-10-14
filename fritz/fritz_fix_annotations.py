
from marshals.interface import api

fritz_base_url = 'https://fritz.science/api/'


def get_source_autoannot(obj_id, testing=False):
    """
    gets autoannotations

    obj_id: 'ZTF18aaaaaa', for example
    testing: <bool> are we testing only?

    return: True if success, False if not
    """

    aids = []
    if testing:
        print("TESTING get_source_autoannot(): no data sent to marshal")
    else:
        fritz_annotation_url = fritz_base_url + \
                               'sources/%s/annotations' % obj_id
        r = api("GET", fritz_annotation_url)
        if 'success' in r['status']:
            annotations = r['data']['annotations']
            for ann in annotations:
                if 'SNIascore:' in ann['origin'] or 'sedm:' in ann['origin']:
                    aids.append(ann['id'])
        else:
            print('error getting annotations for %s' % obj_id)
            print(r['status'])
            print(r['message'])

    return aids


if __name__ == "__main__":
    anids = get_source_autoannot('ZTF21abtsoky')
    print("%d annotations found")
