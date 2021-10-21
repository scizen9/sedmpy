import re
import base64
import math
from glob import glob
from getpass import getpass
from pprint import pprint

from marshals.interface import api
import tns.sedm_auto_tns as tns

fritz_base_url = 'https://fritz.science/api/'
fritz_classification_url = fritz_base_url + 'classification'
fritz_redshift_update_url = fritz_base_url + 'sources/'


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
    ddict = {'obj_id': obj_id, 'spectrum_id': spec_id,  # these go
             'text': comment,
             'attachment': {'body': encoded,
                            'name': fname.split('/')[-1]}}
    if testing:
        print("TESTING add_spec_attachment(): no data sent to marshal")
        print("%s: %s encoded with length %d" % (obj_id, fname.split('/')[-1],
                                                 len(encoded)))
        return True
    else:
        fritz_comment_url = fritz_base_url + 'spectra/%d/comments' % spec_id
        r = api("POST", fritz_comment_url, data=ddict)
        if 'success' in r['status']:
            r_data = r['data']
            if 'comment_id' in r_data:
                print("Comment id = %d" % int(r_data['comment_id']))
            print('{} uploaded'.format(fname.split('/')[-1]))
            return True
        else:
            print('error submitting comment with attachment')
            print(r['status'])
            print(r['message'])
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

    ddict = {'origin': origin, 'data': andic}

    if testing:
        print("TESTING add_spec_autoannot(): no data sent to marshal")
        print(ddict)
        return True
    else:
        # fritz_annotation_url = fritz_base_url + \
        #                       'spectra/%d/annotations' % spec_id
        fritz_annotation_url = fritz_base_url + \
                               'sources/%s/annotations' % obj_id
        r = api("POST", fritz_annotation_url, data=ddict)
        if 'success' in r['status']:
            r_data = r['data']
            if 'annotation_id' in r_data:
                print("annotation id = %d" % int(r_data['annotation_id']))
            print('{}: {} posted'.format(obj_id, origin))
            return True
        else:
            print('error submitting annotation')
            print(r['status'])
            print(r['message'])
            return False


def add_SNIascore_pysedm_autoannot(fname, object_id=None, spec_id=None,
                                   testing=False, upload_tns=True):
    """
    adds autoannotations with SNIASCORE and error
    if SNIASCORE > 0.9, also adds SNIASCORE redshift and error

    fname: '*ZTF18aaaaaaa.txt' that has a bunch of
            "# SNIASCORE[something]: [val]" in the header
    object_id: (str)
    spec_id: (int)
    testing: (bool)
    upload_tns: (bool)

    returns: True if autoannotation works, False
            (and it'll exit early) otherwise
    """

    tns_upl = False
    file_ext = fname.split('.')[-1]
    assert file_ext == 'txt' or file_ext == 'ascii'

    with open(fname) as f:
        header = {line.split(':', 1)[0][1:].strip():
                  line.split(':', 1)[-1].strip()
                  for line in f if line[0] == '#'}

    # SNIascore RESULTS
    if 'SNIASCORE' not in header:
        print(fname, "never run through SNIascore?")
        return False

    if float(header['SNIASCORE']) < 0:
        print('no score')
        return False

    # construct annotations dictionary
    if float(header['SNIASCORE']) >= 0.9:
        # Post classification
        if add_SNIascore_classification(fname, object_id=object_id,
                                        testing=testing):
            print("POSTed Ia classification to fritz")
            # Attempt to post to TNS
            if upload_tns:
                try:
                    if tns.sedm_tns_classify(fname, ztfname=object_id,
                                             specid=spec_id,
                                             testing=testing):
                        print("Uploaded SNIa classification to TNS")
                        tns_upl = True
                    else:
                        print("Unable to upload SNIa classification to TNS")
                except:
                    print("Problems connecting")
        else:
            print("Unable to post Ia classification to fritz")
        # Generate annotation dictionary
        andic = {
            'SNIascore': header['SNIASCORE'],
            'SNIascore_err': header['SNIASCORE_ERR'],
            'SNIa_z': header['SNIASCORE_Z'],
            'SNIa_z_err': header['SNIASCORE_ZERR']
        }
    else:
        andic = {
            'SNIascore': header['SNIASCORE'],
            'SNIascore_err': header['SNIASCORE_ERR']}
    # construct origin
    origin = 'SNIascore:spc%d' % spec_id  # origin = 'SNIascore'

    print(andic)

    return add_spec_autoannot(object_id, andic, spec_id=spec_id,
                              origin=origin, testing=testing), tns_upl


def add_SNIascore_classification(fname, object_id=None, testing=False):
    """
    adds SNIASCORE "Ia" classification if SNIASCORE > 0.9

    fname: '*ZTF18aaaaaaa.txt' that has a bunch of
            "# SNIASCORE[something]: [val]" in the header
    object_id: (str)
    testing: (bool)

    returns: True if classification works, False
            (and it'll exit early) otherwise
    """

    file_ext = fname.split('.')[-1]
    assert file_ext == 'txt' or file_ext == 'ascii'

    with open(fname) as f:
        header = {line.split(':', 1)[0][1:].strip():
                  line.split(':', 1)[-1].strip()
                  for line in f if line[0] == '#'}

    # SNIascore RESULTS
    if 'SNIASCORE' not in header:
        print(fname, "never run through SNIascore?")
        return False

    if float(header['SNIASCORE']) < 0:
        print('no score')
        return False

    # construct annotations dictionary
    if float(header['SNIASCORE']) >= 0.9:
        cldict = {
            "obj_id": object_id,
            "classification": "Ia",
            "taxonomy_id": 3,
            "probability": float(header['SNIASCORE'])
        }
        if testing:
            print("TESTING add_SNIascore_classification():"
                  " no data sent to marshal")
            print(cldict)
            return True
        else:
            r = api("POST", fritz_classification_url, data=cldict)
            if 'success' in r['status']:
                r_data = r['data']
                if 'classification_id' in r_data:
                    print("classification id = %d" %
                          int(r_data['classification_id']))
                print('{}: Ia classification posted'.format(object_id))
                # now add redshift
                if 'SNIASCORE_Z' in header and 'SNIASCORE_ZERR' in header:
                    # What is the current redshift set to?
                    rc = api("GET", fritz_redshift_update_url + object_id)
                    if 'success' in rc['status']:
                        rc_data = rc['data']
                        current_redshift = None
                        if 'redshift' in rc_data:
                            current_redshift = rc_data['redshift']
                        # Only set redshift if it is not already set
                        if current_redshift is None:
                            new_redshift = float(header['SNIASCORE_Z'])
                            new_redshift_error = float(header['SNIASCORE_ZERR'])
                            try:
                                new_z_round = math.ceil(abs(
                                    math.log10(new_redshift_error)))
                            # Handle negative, NaN, Inf, None and <str> values
                            except (ValueError, OverflowError, TypeError):
                                new_z_round = 1
                            new_z = round(new_redshift,
                                          1 if new_z_round < 1 else new_z_round)
                            new_error = round(new_redshift_error,
                                              1 if new_z_round < 1 else
                                              new_z_round)
                            rsdict = {"redshift": new_z,
                                      "redshift_error": new_error}
                            rr = api("PATCH",
                                     fritz_redshift_update_url + object_id,
                                     data=rsdict)
                            if 'success' in rr['status']:
                                print("redshift for %s updated to %.4f +- %.4f"
                                      % (object_id, rsdict['redshift'],
                                         rsdict['redshift_error']))
                            else:
                                print('error updating %s redshift' % object_id)
                                print(rr['status'])
                                print(rr['message'])
                        else:
                            print('Redshift for %s already set to %.4f' %
                                  (object_id, float(rc_data['redshift'])))
                    else:
                        print('error getting current redshift for %s' %
                              object_id)
                else:
                    print('No SNIascore redshift records found for %s ' %
                          object_id)
                return True
            else:
                print('error submitting classification')
                print(r['status'])
                print(r['message'])
                return False

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
    origin = 'sedm:spc%d' % spec_id  # origin = 'sedm'

    if not add_spec_autoannot(object_id, andic, spec_id=spec_id,
                              origin=origin, testing=testing):
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
