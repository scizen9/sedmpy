import numpy as np
from astropy.time import Time
import os
import sys
import datetime
import requests
import json
import argparse
from time import sleep

from marshals.interface import api

global TOKEN, BASEURL
GETTOKEN = ''      # Fritz API Key, input your TOKEN from Fritz
BASEURL = 'https://fritz.science/'                     # Fritz base url

API_KEY = "54916f1700966b3bd325fc1189763d86512bda1d"     # TNS API Key

# TNS URLs for real uploads
TNS_BASE_URL = "https://www.wis-tns.org/api/"
upload_url = "https://www.wis-tns.org/api/file-upload"
report_url = "https://www.wis-tns.org/api/bulk-report"
reply_url = "https://www.wis-tns.org/api/bulk-report-reply"

# SANDBOX URLs for TNS upload trials
SAND_TNS_BASE_URL = "https://sandbox-tns.org/api/"
SAND_upload_url = "https://sandbox-tns.org/api/"
SAND_report_url = "https://sandbox-tns.org/api/bulk-report"
SAND_reply_url = "https://sandbox-tns.org/api/bulk-report-reply"


def get_source_api(ztf_name):
    """ Info : Query a single source, takes input ZTF name
        Returns : all basic data of that source (excludes photometry and
         spectra, includes redshift, classification, comments, etc.)
    """
    url = BASEURL+'api/sources/'+ztf_name+'?includeComments=true'
    resp = api('GET', url).json()
    return resp['data']


def get_group_ids(groupnames=None):
    """ Info : Query group ids of groups specified
        Input : Name or names of groups in an array []
        Returns : List of group  names and their group ids
    """

    if groupnames is None:
        groupnames = ['Redshift Completeness Factor',
                      'Census of the Local Universe Caltech']

    url = BASEURL+'api/groups'
    headers = {'Authorization': 'token {}'.format(GETTOKEN)}
    groupnames = np.atleast_1d(groupnames)
    grpids = []
    for grpname in groupnames:
        resp = requests.request('GET', url, params={'name': grpname},
                                headers=headers).json()
        answer = str(grpname)+' = '+str(resp['data'][0]['id'])
        grpids.append(answer)

    return grpids


def get_number_of_sources(group_id, date):
    """ Info : Query number of sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : Number of sources saved after a given date to
         the specified group
    """

    url = BASEURL+'api/sources?saveSummary=true&group_ids=' + group_id + \
        '&savedAfter='+date+'T00:00:00.000001'
    resp = api('GET', url).json()
    return len(resp['data']['sources'])


def get_group_sources(group_id, date):
    """ Info : Query all sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : List of jsons of all sources in group(s)
        Comment : Takes a little time based on the date
    """

    srces = []

    for si in range(get_number_of_sources(group_id, date)):

        url = BASEURL+'api/sources?saveSummary=true&group_ids=' + group_id + \
              '&savedAfter='+date+'T00:00:00.000001'
        resp = api('GET', url).json()
        ztf_name = resp['data']['sources'][si]['obj_id']
        srces.append(ztf_name)

    return srces


def get_total_number_of_sources(group_id):
    """ Info : Query total number of sources saved in a group
        Input : group id
        Returns : Total number of sources saved in a group
    """

    url = BASEURL + 'api/sources?saveSummary=true&group_ids=' + group_id
    resp = api('GET', url).json()
    return len(resp['data']['sources'])


def get_all_group_sources(group_id):
    """ Info : Query all sources saved in a group
        Input : group id
        Returns : List of jsons of all sources in group(s)
        Comment : Takes a long time
    """

    srces = []

    for si in range(get_number_of_sources(group_id)):

        url = BASEURL + 'api/sources?saveSummary=true&group_ids=' + group_id
        resp = api('GET', url).json()
        ztf_name = resp['data']['sources'][si]['obj_id']
        srces.append(ztf_name)

    return srces


def get_iau_name(ztf_name):

    """ Info : Query the TNS name for any source
        Input : ZTFname
        Returns : ATname
    """

    url = BASEURL + 'api/alerts/ztf/' + ztf_name + '/aux'
    resp = api('GET', url).json()
    return resp["data"]["cross_matches"]["TNS"]


def get_classification(ztf_name):

    """ Info : Query the classification and classification date for any source
        Input : ZTFname
        Returns : Classification and Classification date
        Comment : You need to choose the classification
                if there are multiple classifications
    """

    url = BASEURL + 'api/sources/' + ztf_name + '/classifications'
    resp = api('GET', url).json()
    output = resp['data']

    classification = None
    classification_date = None

    if len(output) < 1:
        classification = "No Classification found"
        classification_date = "None"

    if len(output) == 1:

        classification = resp['data'][0]['classification']
        classification_date = resp['data'][0]['created_at'].split('T')[0]

    if len(output) > 1:

        classification = []
        classification_date = []

        for si in range(len(output)):

            clss = resp['data'][si]['classification']
            classify_date = resp['data'][si]['created_at']

            classification.append(clss)
            classification_date.append(classify_date)

        for si in range(len(classification)):

            print((si+1), ")", "Classification: ", classification[si],
                  "\t Classification date:",
                  classification_date[si].split('T')[0])

        class_pick = input("Choose classification: ")

        classification = classification[int(class_pick)-1]
        classification_date = classification_date[
            int(class_pick)-1].split('T')[0]

    return classification, classification_date


def get_redshift(ztf_name):

    """ Info : Query the redshift for any source
        Input : ZTFname
        Returns : redshift
    """

    url = BASEURL + 'api/sources/' + ztf_name
    resp = api('GET', url).json()

    redshift = resp['data']['redshift']

    if redshift is None:
        redshift = "No redshift found"

    return redshift


def get_tns_information(ztf_name):

    # url = BASEURL+'api/sources/'+ztf_name
    # resp = api('GET', url).json()

    iau = get_iau_name(ztf_name)

    if not iau:
        iau = "Not reported to TNS"

    else:
        iau = iau[0]['name']

    clss = get_classification(ztf_name)

    if clss[1] == 'None':
        clss = "Not classified yet"

    else:
        clss = ('Classification: ' + str(clss[0]) + ',' +
                ' Classification date: ' + str(clss[1]))

    redshift = get_redshift(ztf_name)

    if redshift is None:
        redshift = "No redshift found"

    else:
        redshift = ('redshift:'+str(redshift))

    return ztf_name, iau, clss, redshift


def convert_to_jd(date):

    dd = Time(date, format='fits')
    return dd.jd


def get_spectrum_api(spectrum_id):
    """ Info : Query all spectra corresponding to a source, takes input ZTF name
        Returns : list of spectrum jsons
    """
    url = BASEURL + 'api/spectrum/' + str(spectrum_id)
    resp = api('GET', url).json()
    return resp


def get_all_spectra_len(ztf_name):

    url = BASEURL+'api/sources/'+ztf_name+'/spectra'
    resp = api('GET', url).json()
    return len(resp['data']['spectra'])


def get_all_spectra_id(ztf_name):
    """ Info : Query all spectra corresponding to a source, takes input ZTF name
        Returns : list of spectrum jsons
    """

    spec_id = []

    for si in range(get_all_spectra_len(ztf_name)):

        url = BASEURL + 'api/sources/' + ztf_name + '/spectra'
        resp = api('GET', url).json()

        sid = resp['data']['spectra'][si]['id']
        spec_id.append(sid)

    return spec_id


def get_required_spectrum_id(ztf_name, spec_file):

    spec_id = None

    specfn = os.path.basename(spec_file)

    nspec = (get_all_spectra_len(ztf_name))

    if nspec > 0:

        spec_ids = get_all_spectra_id(ztf_name)

        for sid in spec_ids:
            spec = get_spectrum_api(sid)
            if specfn in spec['data']['original_file_filename']:
                spec_id = sid

    return spec_id


def write_ascii_file(ztf_name, specid, path=None):

    specfn = None

    if specid is None:
        print("ERROR: no spec_id")

    else:
        a = get_spectrum_api(specid)

        inst = (a['data']['instrument_name'])

        if inst == 'SEDM':

            header = (a['data']['altdata'])

            specfn = (ztf_name + '_' + str(header['OBSDATE']) +
                      '_' + str(inst) + '.ascii')

            with open(os.path.join(path, specfn), 'w') as f:
                f.write(a['data']['original_file_string'])

        else:
            print('ERROR: not an SEDM spectrum!')
            specfn = None

    return specfn


def pprint(*args, **kwargs):
    """
    slightly more convenient function instead of print(get_pprint)

    params:
        *args (arguments to pass to get_pprint)
        **kwargs (keyword arguments to pass to get_pprint)
    """
    print(get_pprint(*args, **kwargs))


def post_comment(ztf_name, text):

    data = {"obj_id": ztf_name,
            "text": text,
            }

    url = BASEURL + 'api/comment'

    resp = api('POST', url, data=data).json()

    return resp


def pprint(*args, **kwargs):
    """
    slightly more convenient function instead of print(get_pprint)

    params:
        *args (arguments to pass to get_pprint)
        **kwargs (keyword arguments to pass to get_pprint)
    """
    print(get_pprint(*args, **kwargs))


def get_pprint(item, indent=0, tab=' '*4, maxwidth=float('inf')):
    """
    it's just like 'from pprint import pprint', except instead of
    having dictionaries use hanging indents dependent on the length
    of their key, if the value is a list or dict it prints it indented
    by the current indent plus tab

    params:
        item <di or li> (the thing to be printed)
        indent <int> (the number of times it's been indented so far)
        tab <str> (how an indent is represented)
        maxwidth <int> (maximum characters per line in ouptut)

    returns:
        result <str>
    """
    def get_pprint_di(di, indent, tab=' '*4):
        """
        pprints a dictionary

        params:
            di <dict>
            indent <int> (the number of indents so far)

        returns:
            di_str <str>
        """
        di_str = ''
        for i, (key, item) in enumerate(di.items()):
            di_str += tab*indent
            di_str += repr(key) + ': ' + get_pprint(item, indent, tab)
            if i+1 < len(di):
                # everything until the last item has a trailing comma
                di_str += ',\n'
            else:
                di_str += '\n'
        return di_str

    def get_pprint_li(li, indent, tab=' '*4):
        """
        pprints a list

        params:
            li <list>
            indent <int> (the number of indents so far)

        returns:
            current_result <str>
        """
        li_str = ''
        for i, item in enumerate(li):
            li_str += tab*indent
            pprint(item, indent, tab)
            if i+1 < len(li):
                li_str += ',\n'
            else:
                li_str += '\n'
        return li_str

    result = ''
    if isinstance(item, dict):
        result += '{\n'
        result += get_pprint_di(item, indent+1, tab)
        result += tab*indent + '}'
    elif isinstance(item, list):
        result += '[\n'
        result += get_pprint_li(item, indent+1, tab)
        result += tab*indent + ']'
    else:
        result += repr(item)

    # this gets rid of too-long lines, but only supports space tabs
    lines = result.split('\n')
    for i, line in enumerate(lines):
        while max([len(li) for li in line.split('\n')]) > maxwidth:
            tabs = line[:-len(line.lstrip())]
            if len(tabs) > maxwidth - 8:
                break  # giving up
            line = line[:78] + '\\\n' + tabs + 2*tab + line[78:]
            lines[i] = line
    result = '\n'.join(lines)

    return result


def get_number(group_id, date):
    """ Info : Query number of sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : Number of sources saved after a given date to
                the specified group
    """

    url = BASEURL + 'api/sources?saveSummary=true&group_ids=' + group_id + \
        '&savedAfter='+date+'T00:00:00.000001'
    response = api('GET', url).json()
    return len(response['data']['sources'])


def get_sources(group_id, date):
    """ Info : Query all sources saved in a group after a certain date
        Input : group id, date [yyyy-mm-dd]
        Returns : List of jsons of all sources in group(s)
        Comment : Takes a little time based on the date
    """

    sources = []

    for i in range(get_number(group_id, date)):

        url = BASEURL + 'api/sources?saveSummary=true&group_ids=' + group_id + \
              '&savedAfter='+date+'T00:00:00.000001'
        response = api('GET', url).json()
        ztf_name = response['data']['sources'][i]['obj_id']
        sources.append(ztf_name)

    return sources


def get_tns_classification_id(classification):

    class_ids = {'Afterglow': 23, 'AGN': 29, 'CV': 27, 'Galaxy': 30, 'Gap': 60,
                 'Gap I': 61, 'Gap II': 62, 'ILRT': 25, 'Kilonova': 70,
                 'LBV': 24, 'M dwarf': 210, 'Nova': 26, 'Novae': 26, 'QSO': 31,
                 'SLSN-I': 18, 'SLSN-II': 19, 'SLSN-R': 20, 'SN': 1, 'I': 2,
                 'Type I': 2, 'I-faint': 15, 'I-rapid': 16, 'Ia': 3,
                 'Ia-norm': 3, 'Ia-91bg': 103, 'Ia-91T': 104, 'Ia-CSM': 106,
                 'Ia-pec': 100, 'Ia-SC': 102, 'Ia-02cx': 105, 'Ib': 4,
                 'Ib-norm': 4, 'Ib-Ca-rich': 8, 'Ib-pec': 107, 'Ib/c': 6,
                 'SN Ibn': 9, 'Ic': 5, 'Ic-norm': 5, 'Ic-BL': 7, 'Ic-pec': 108,
                 'II': 10, 'Type II': 10, 'II-norm': 10, 'II-pec': 110,
                 'IIb': 14, 'IIL': 12, 'IIn': 13, 'IIn-pec': 112, 'IIP': 11,
                 'SN impostor': 99, 'Std-spec': 50, 'TDE': 120, 'Varstar': 28,
                 'WR': 200, 'WR-WC': 202, 'WR-WN': 201, 'WR-WO': 203,
                 'Other': 0}

    # keys = np.array(class_ids.keys())
    for keys in class_ids:
        if keys == classification:
            classkey = class_ids[keys]
            return classkey


def get_tns_instrument_id(inst):

    inst_ids = {'DBSP': 1, 'ALFOSC': 41, 'LRIS': 3, 'DIS': 70, 'SEDM': 149,
                'SPRAT': 156, 'GMOS': 6, 'Lick-3m': 10, 'LFC': 2, 'TSPEC': 109}

    if inst in inst_ids:
        return inst_ids[inst]
    else:
        return None


class TNSClassificationReport:
    def __init__(self):
        self.name = ''
        self.fitsName = ''
        self.asciiName = ''
        self.classifierName = ''
        self.classificationID = ''
        self.redshift = ''
        self.classificationComments = ''
        self.obsDate = ''
        self.instrumentID = ''
        self.expTime = ''
        self.observers = ''
        self.reducers = ''
        self.specTypeID = ''
        self.spectrumComments = ''
        self.groupID = ''
        self.spec_proprietary_period_value = ''
        self.spec_proprietary_period_units = ''

    def fill(self):
        spectrumdict = {
            'obsdate': self.obsDate,
            'instrumentid': self.instrumentID,
            'exptime': self.expTime,
            'observer': self.observers,
            'reducer': self.reducers,
            'spectypeid': self.specTypeID,
            'ascii_file': self.asciiName,
            'fits_file': self.fitsName,
            'remarks': self.spectrumComments,
            'spec_proprietary_period': self.spec_proprietary_period_value}

        classification_dict = {
            'classification_report': {
                '0': {
                    'name': self.name,
                    'classifier': self.classifierName,
                    'objtypeid': self.classificationID,
                    'redshift': self.redshift,
                    'groupid': self.groupID,
                    'remarks': self.classificationComments,
                    'spectra': {
                        'spectra-group': {
                            '0': spectrumdict
                        }
                    }
                }
            }
        }

        return classification_dict

    def classification_json(self):
        return json.dumps(self.fill())


def upload_to_tns(filename, base_url=upload_url, api_key=API_KEY,
                  filetype='ascii'):
    """
    uploads a file to TNS and returns the response json
    """
    url = base_url
    data = {'api_key': api_key}

    if filetype is 'ascii':
        files = [('files[]', (filename, open(filename), 'text/plain'))]

    elif filetype is 'fits':
        files = [('files[0]', (filename, open(filename, 'rb'),
                               'application/fits'))]
    else:
        files = None

    if files is not None:
        response = requests.post(url, data=data, files=files)
        try:
            return response.json()
        except:
            print(url, data, files, response.content, sep='\n')
            return False
    else:
        return False


def tns_classify(classification_report, base_url=report_url, api_key=API_KEY):
    """
    submits classification report to TNS and returns the response json
    """
    url = base_url
    data = {'api_key': api_key,
            'data': classification_report.classification_json()}
    resp = requests.post(url, data=data).json()
    if not resp:
        return False

    res_code = resp['id_code']
    reprt_id = resp['data']['report_id']
    print("ID:", reprt_id)
    print(res_code, resp['id_message'], "reporting finished")
    if res_code == 200:
        return reprt_id
    else:
        print("Result reporting didn't work")
        pprint(resp)
        print("re-submit classification, but don't re-upload files")
        return False


def tns_feedback(reprt_id):
    data = {'api_key': API_KEY, 'report_id': reprt_id}
    resp = requests.post(TNS_BASE_URL + 'bulk-report-reply',
                         data=data).json()
    feedback_code = resp['id_code']
    print(feedback_code, resp['id_message'], "feedback finished")
    if feedback_code == 200:
        return True
    elif feedback_code == 404:
        print("Waiting and retrying...")
        sleep(2)
        try:
            return tns_feedback(reprt_id)
        except KeyboardInterrupt:
            return False
    elif feedback_code == 400:
        print(resp)
        return False
    else:
        # error receiving the feedback from TNS about the upload
        print("Something went wrong with the feedback, but the report may",
              "still have been fine?")
        return False


def sedm_tns_classify(spec_file, ztfname=None, specid=None, testing=False):
    """Verify the input source, prepare a classification report and
    upload to TNS"""

    if not os.path.exists(spec_file):
        print("ERROR: File not found!: %s" % spec_file)
        return False

    with open(spec_file) as f:
        header = {line.split(':', 1)[0][1:].strip():
                  line.split(':', 1)[-1].strip()
                  for line in f if line[0] == '#'}

    if ztfname is None:
        ztfname = header['NAME']

    comments = [c['text'] for c in get_source_api(ztfname)['comments']]
    if 'Uploaded to TNS' in comments:
        print("Already uploaded to TNS")
        return False

    info = get_tns_information(ztfname)

    if info[2] == 'Not classified yet':         # Check if classified
        print(info[2])
        return False

    class_date = info[2].split(',')[-1].split(':')[-1].strip()
    classify = info[2].split(',')[0].split(':')[-1].strip()

    # print(info)

    path = os.path.dirname(spec_file)

    if specid is None:
        specid = get_required_spectrum_id(ztfname, spec_file)

    spectrum_name = write_ascii_file(ztfname, specid, path=path)

    if spectrum_name is None:
        print("No spectrum found")
        return False

    specfile = os.path.join(path, spectrum_name)

    classifiers = 'A. Dahiwale, C. Fremling(Caltech) on behalf of the ' \
                  'Zwicky Transient Facility (ZTF)'
    source_group = 48  # Require source group id from drop down list, 0 is for
                       # None

    proprietary_period = '0'
    proprietary_units = "years"
    spec_comments = ''
    classification_comments = ''
    spectype = 'object'
    spectype_id = ['object', 'host', 'sky', 'arcs',
                   'synthetic'].index(spectype) + 1

    a = get_spectrum_api(specid)
    header = (a['data']['altdata'])
    obsdate = str((header['UTC']).split('T')[0]) + \
        ' ' + str((header['UTC']).split('T')[1])

    classification_report = TNSClassificationReport()
    classification_report.name = get_iau_name(ztfname)[0]['name'][3:]
    classification_report.fitsName = ''
    classification_report.asciiName = spectrum_name
    classification_report.classifierName = classifiers
    classification_report.classificationID = get_tns_classification_id(classify)
    classification_report.redshift = get_redshift(ztfname)
    classification_report.classificationComments = classification_comments
    classification_report.obsDate = obsdate
    classification_report.instrumentID = get_tns_instrument_id('SEDM')
    classification_report.expTime = (header['EXPTIME'])
    classification_report.observers = 'SEDmRobot'
    classification_report.reducers = (header['REDUCER'])
    classification_report.specTypeID = spectype_id
    classification_report.spectrumComments = spec_comments
    classification_report.groupID = source_group
    classification_report.spec_proprietary_period_value = proprietary_period
    classification_report.spec_proprietary_period_units = proprietary_units

    pprint(classification_report.fill(), tab='  ')

    # ASCII FILE UPLOAD
    if not testing:
        print("\n")
        response = upload_to_tns(specfile)
        print(response)

        if not response:
            print("File upload didn't work")
            print(response)
            return False

        print(response['id_code'], response['id_message'],
              "\nSuccessfully uploaded ascii spectrum")
        # classification_report.asciiName = response['data'][-1]

        report_id = tns_classify(classification_report)
        post_comment(ztfname, 'Uploaded to TNS')
        tns_feedback(report_id)
    else:
        print(classification_report)

    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""

Uploads classification report to the TNS website.

""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--data_file', type=str, default=None,
                        help='Data file to upload: '
                             '/scr2/sedmdrp/redux/YYYYMMDD/spec_*.txt')
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to TNS (for testing)')
    args = parser.parse_args()

    infile = args.data_file

    # Check input
    if not os.path.exists(infile):
        print("File not found: %s" % infile)
    else:
        print("Uploading from %s" % infile)
        sedm_tns_classify(infile, testing=args.testing)
