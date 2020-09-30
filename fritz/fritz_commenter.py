import re
from glob import glob
from getpass import getpass
from pprint import pprint

from marshals.interface import api

fritz_base_url = 'http://private.caltech.edu/'
fritz_comment_url = fritz_base_url + 'comment'


def add_spec_attachment(obj_id, comment, fname, spec_id=None, testing=False):
    """
    adds a comment (not autoannotation) with attachment to a particular SEDM
    spectrum on the view_spec page, which will also appear elsewhere.

    ztfname: 'ZTF18aaaaaa', for example
    comment: <str>
    fname: <str> spec file name
    cred: (<str>, <str>) credentials
    sourceid: <int> around ~4000
    specid: <int> around ~1700. faster if you provide it
    obsdate: 'YYYYMMDD', necessary for finding the correct spectrum. Assumes
            exactly 1 SEDM spectrum that night
    reducedby: <str> who did the reduction, defaults to 'auto'
    testing: <bool> are we just testing?

    return: True if success, False if not
    """

    if obj_id is None or spec_id is None:
        print("ERROR - Unable to get info required to post comment")
        return False

    # TODO: put in base64-encoded contents from file as attachement
    ddict = {'obj_id': obj_id,
             'text': comment,
             'attachment': fname}

    if testing:
        print("TESTING add_spec_attachment(): no data sent to marshal")
        print(fname)
        print(ddict)
        return True
    else:
        r = api("POST", fritz_comment_url, data=ddict)
        if 'success' in r['status']:
            print('{} uploaded'.format(fname.split('/')[-1]))
            return True
        else:
            print('error submitting comment')
            print(r.status)
            raise Exception
            # return False


def add_spec_autoannot(obj_id, text, spec_id=None, testing=False):
    """
    adds an autoannotation without attachment to a particular SEDM spectrum
    on the view_spec page, which will also appear elsewhere.

    obj_id: 'ZTF18aaaaaa', for example
    text: <str>
    spec_id: <int>
    testing: <bool> are we testing only?

    return: True if success, False if not
    """

    ddict = {'obj_id': obj_id,
             'text': text}

    if testing:
        print("TESTING add_spec_autoannot(): no data sent to marshal")
        print(ddict)
        return True
    else:
        r = api("POST", fritz_comment_url, data=ddict)

        if 'success' in r.status:
            print('{}: {} posted'.format(obj_id, text))
            return True
        else:
            print('error submitting comment')
            print(r.status)
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
    # TODO if we look at the _snid.output we can get more info, eg phase
    file_ext = fname.split('.')[-1]
    assert file_ext == 'txt' or file_ext == 'ascii'

    with open(fname) as f:
        header = {line.split(':', 1)[0][1:].strip().lower():
                  line.split(':', 1)[-1].strip()
                  for line in f if line[0] == '#'}

    # PYSEDM_REPORT
    # must be posted after the SNID plot or else it'll be overwritten
    try:
        pysedm_report = glob(fname.replace('spec',
                                           'pysedm_report').replace('.txt',
                                                                    '.png'))[0]
        pr_posted = add_spec_attachment(object_id, 'pysedm_report',
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

    dtypes = {'match': 'STRING', 'rlap': 'FLOAT', 'redshift': 'FLOAT'}
    for key in dtypes:
        if 'STRING' in dtypes[key]:
            r = add_spec_autoannot(object_id, '[AUTO_SNID_' + key + '] ' +
                                   header['snidmatch' + key], testing=testing)
        elif 'FLOAT' in dtypes[key]:
            r = add_spec_autoannot(object_id, '[AUTO_SNID_' + key + ']  %.3f' %
                                   header['snidmatch' + key], testing=testing)
        else:
            r = False
        return r

    if pr_posted:
        return True  # we already have an attachment comment so don't overwrite

    # SNID PLOT
    # NOTE: this makes a major assumption about the naming scheme of snid plots
    image_filename = fname.replace('.txt',
                                   '_{}.png'.format(header['snidmatchtype']))
    if not glob(image_filename):
        return False
    add_spec_attachment(object_id, 'AUTO_SNID_plot', image_filename,
                        spec_id=spec_id, testing=testing)

    return True


if __name__ == "__main__":
    auth = (input('GROWTH marshal username:'), getpass())

    successes = []
    for filename in glob('*ZTF?????????.txt'):
        if add_SNID_pysedm_autoannot(filename, auth):
            print('success!')
            successes.append(re.findall(r'ZTF\d{2}[a-z]{7}', filename)[-1])
            break
    pprint(successes)
