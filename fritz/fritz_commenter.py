import re
import base64
from glob import glob
from getpass import getpass
from pprint import pprint

from marshals.interface import api

fritz_base_url = 'https://fritz.science/api/'
fritz_comment_url = fritz_base_url + 'comment'
fritz_annotation_url = fritz_base_url + 'annotation'


def add_spec_attachment(obj_id, comment, fname, spec_id=None, testing=False):
    """
    adds a comment (not autoannotation) with attachment to a particular SEDM
    spectrum on the view_spec page, which will also appear elsewhere.

    obj_id: <str>
    comment: <str>
    fname: <str>
    spec_id: <int>
    testing: <bool> are we just testing?

    return: True if success, False if not
    """

    if obj_id is None or spec_id is None:
        print("ERROR - Unable to get info required to post comment")
        return False

    # read in file
    with open(fname, 'rb') as image_file:
        encoded = base64.b64encode(image_file.read()).decode('ascii')
    # create payload
    ddict = {'obj_id': obj_id,   # 'commentable_id': spec_id,
             'text': comment,
             'attachment': {'body': encoded,
                            'name': fname.split('/')[-1]}}
    if testing:
        print("TESTING add_spec_attachment(): no data sent to marshal")
        print("%s: %s encoded with length %d" % (obj_id, fname.split('/')[-1],
                                                 len(encoded)))
        return True
    else:
        r = api("POST", fritz_comment_url, data=ddict).json()
        if 'success' in r['status']:
            print('{} uploaded'.format(fname.split('/')[-1]))
            return True
        else:
            print('error submitting comment with attachment')
            print(r['status'])
            return False


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

    ddict = {'obj_id': obj_id,  # 'commentable_id': spec_id,
             'origin': origin,
             'data': andic}

    if testing:
        print("TESTING add_spec_autoannot(): no data sent to marshal")
        print(ddict)
        return True
    else:
        r = api("POST", fritz_annotation_url, data=ddict).json()

        if 'success' in r['status']:
            print('{}: {} posted'.format(obj_id, origin))
            return True
        else:
            print('error submitting comment')
            print(r['status'])
            return False


def add_SNID_pysedm_autoannot(fname, object_id=None, spec_id=None,
                              testing=False):
    """
    if z < 0.3 and rlap > 5.0
        adds autoannotations with SNID rlap, z, type, etc
        adds a comment with the SNID plot attached

    fname: '*ZTF18aaaaaaa.txt' that has a bunch of
            "# SNIDMATCH[something]: [val]" in the header
    cred: ('username', 'password')
    reducedby: (str)
    testing: (bool)

    returns: True if all four comments/attachments works, False
            (and it'll exit early) otherwise
    """

    file_ext = fname.split('.')[-1]
    assert file_ext == 'txt' or file_ext == 'ascii'

    with open(fname) as f:
        header = {line.split(':', 1)[0][1:].strip().lower():
                  line.split(':', 1)[-1].strip()
                  for line in f if line[0] == '#'}

    # Upload pysedm_report
    try:
        pysedm_report = glob(fname.replace('spec',
                                           'pysedm_report').replace('.txt',
                                                                    '.png'))[0]
        pr_posted = add_spec_attachment(object_id,
                                        'pysedm_report:spc%d' % spec_id,
                                        pysedm_report, spec_id=spec_id,
                                        testing=testing)
    except IndexError:
        print('no pysedm_report for {}?'.format(header['name']))
        pr_posted = False

    # SNID RESULTS
    if 'snidmatchtype' not in header:
        print(fname, "never run through snid?")
        return False
        
    if header['snidmatchtype'].lower() == 'none':
        print('no match')
        return False

    elif float(header['snidmatchrlap']) < 5:
        print('bad rlap, only {}'.format(header['snidmatchrlap']))
        return False

    elif (header['snidmatchtype'][0] == 'I') \
            and not (0.01 <= float(header['snidmatchredshift']) <= 0.3):
        print('bad redshift, {snidmatchredshift} '
              'for {snidmatchtype}'.format(**header))
        return False

    if header['snidmatchsubtype'] == '-':
        header['snidmatchmatch'] = header['snidmatchtype']
    else:
        header['snidmatchmatch'] = '-'.join([header['snidmatchtype'],
                                             header['snidmatchsubtype']])

    # construct annotations dictionary
    andic = {'match': 'None', 'rlap': 0., 'redshift': 0., 'age': 0.}
    for key in andic:
        andic[key] = header['snidmatch' + key]
    # construct origin
    origin = 'sedm:spc%d' % spec_id

    if not add_spec_autoannot(object_id, andic, origin=origin, testing=testing):
        return False

    if pr_posted:
        return True  # we already have an attachment comment so don't overwrite

    # SNID PLOT
    # NOTE: this makes a major assumption about the naming scheme of snid plots
    image_filename = fname.replace('.txt',
                                   '_{}.png'.format(header['snidmatchtype']))
    if not glob(image_filename):
        return False
    ret = add_spec_attachment(object_id, 'AUTO_SNID_plot', image_filename,
                              spec_id=spec_id, testing=testing)
    return ret


if __name__ == "__main__":
    auth = (input('GROWTH marshal username:'), getpass())

    successes = []
    for filename in glob('*ZTF?????????.txt'):
        if add_SNID_pysedm_autoannot(filename, auth):
            print('success!')
            successes.append(re.findall(r'ZTF\d{2}[a-z]{7}', filename)[-1])
            break
    pprint(successes)
