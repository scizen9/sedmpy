import re
import requests
from glob import glob
from getpass import getpass
from pprint import pprint
import json

growth_base_url = 'http://skipper.caltech.edu:8080/cgi-bin/growth/'
growth_source_summary_url = growth_base_url + 'source_summary.cgi?sourcename='

try:
    user, pwd = open('/home/sedm/.growth_creds.txt', 'r').readlines()[0].split()
    auth = (user, pwd)
    # config = json.load(open('config.json'))
    # auth = tuple(config['growth_auth']['adugas'])
except FileNotFoundError:
    print("ERROR - could not find credentials file!")
    auth = None


def get_missing_info(ztfname, obsdate, sourceid, specid, reducedby=None):
    """
    #TODO this is currently being called 5 times per object, which may
        include a lot of downloading the same thing repeatedly.
        Should probably make a ztf_object class, then save specid etc to that,
        then call add_spec_attachment, etc methods to it. Another day '\_0_/'
    """

    try:
        source_summary = requests.get(growth_source_summary_url
                                      + ztfname, auth=auth).json()
    except json.decoder.JSONDecodeError:
        print("ERROR - json could not decode")
        return None, None

    if not sourceid:
        try:
            sourceid = source_summary['id']
        except IndexError as e:
            print(e)
            print("ERROR - count not get id from source summary")
        except KeyError:
            print("ERROR - could not obtain source summary")

    if not specid:
        if not reducedby:
            reducedby = 'auto'
        obsdate = obsdate.replace('-', '')  # YYYYMMDD, YYYY-MM-DD, Y-YY-YM-MDD
        try:
            specid = [spec['specid']
                      for spec in source_summary['uploaded_spectra']
                      if spec['obsdate'].replace('-', '') == obsdate
                      and spec['instrumentid'] == 65
                      and spec['reducedby'].strip() == reducedby][-1]
        except KeyError:
            print("ERROR - could not obtain source summary")
        except IndexError as e:
            print(e)
            print("Target OBSDATE: %s" % obsdate)
            print("Target reducedby: %s" % reducedby)
            pprint([(spec['reducedby'], spec['obsdate'])
                    for spec in source_summary['uploaded_spectra']
                    if spec['instrumentid'] == 65])
    return sourceid, specid


def add_spec_attachment(ztfname, comment, fname, cred, sourceid=None,
                        specid=None, obsdate=None, reducedby=None,
                        testing=False):
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

    sourceid, specid = get_missing_info(ztfname, obsdate, sourceid, specid,
                                        reducedby=reducedby)

    if sourceid is None or specid is None:
        print("ERROR - Unable to get info required to post comment")
        return False

    ddict = {'comment':    comment,
             'tablename': 'spec',
             'tableid':    specid,
             'camefrom':  'view_spec',
             'name':       ztfname,
             'sourceid':   sourceid,
             'specid':     specid,
             'commit':    'yes'}

    if testing:
        print("TESTING add_spec_attachment(): no data sent to marshal")
        print(fname)
        print(ddict)
        return True
    else:
        with open(fname, 'rb') as att:
            r = requests.post(growth_base_url + "add_spec.cgi", auth=cred,
                              data=ddict, files={"attachment": att})

        if 'Updating Database' in r.text or 'copied file successfully' in r.text:
            print('{} uploaded'.format(fname.split('/')[-1]))
            return True
        else:
            print('error submitting comment')
            print(r.text)
            raise Exception
            # return False


def add_spec_autoannot(ztfname, value, annot_type, datatype, cred,
                       sourceid=None, specid=None, obsdate=None,
                       reducedby=None, testing=False):
    """
    adds an autoannotation without attachment to a particular SEDM spectrum
    on the view_spec page, which will also appear elsewhere.

    ztfname: 'ZTF18aaaaaa', for example
    comment: <str>
    datatype: 'STRING' or 'FLOAT' or 'BOOL'
    cred: (growth username, growth password) as (<str>, <str>)
    sourceid: <int> around ~4000
    specid: <int> around ~1700. faster if you provide it
    obsdate: 'YYYYMMDD' or 'YYYY-MM-DD', necessary for finding the correct
            spectrum if no specid. Assumes exactly 1 SEDM spectrum that night
    testing: <bool> are we testing only?

    return: True if success, False if not
    """

    sourceid, specid = get_missing_info(ztfname, obsdate, sourceid, specid,
                                        reducedby=reducedby)

    ddict = {'comment':    value,
             'datatype':   datatype,
             'tablename': 'spec',
             'type':       annot_type,
             'tableid':    specid,
             'table':     'autoannotations',
             'camefrom':  'view_spec',
             'name':       ztfname,
             'sourceid':   sourceid,
             'specid':     specid,
             'commit':    'yes'}

    if testing:
        print("TESTING add_spec_autoannot(): no data sent to marshal")
        print(ddict)
        return True
    else:
        r = requests.post(growth_base_url + "add_spec.cgi", auth=cred,
                          data=ddict)

        if 'Updating Database' in r.text or 'copied file successfully' in r.text:
            print('{}: {} posted'.format(annot_type, value))
            return True
        else:
            print('error submitting comment')
            print(r.text)
            return False


def add_SNID_pysedm_autoannot(fname, cred, reducedby=None, testing=False):
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

    # get the specid by comparing all the header info
    sourcename = header['name']
    try:
        source_summary = requests.get(growth_base_url +
                                  'source_summary.cgi?'
                                  'sourcename=%s' % sourcename,
                                  auth=cred).json()
    except json.decoder.JSONDecodeError:
        print("ERROR: could not decode results")
        return False

    if 'uploaded_spectra' in source_summary:
        if source_summary['uploaded_spectra']:
            for spec in source_summary['uploaded_spectra']:
                f = requests.get('http://skipper.caltech.edu:8080/growth-data/'
                                 + spec['datapath'],
                                 auth=cred).text
                fheader = {line.split(':', 1)[0][1:].strip().lower():
                           line.split(':', 1)[-1].strip()
                           for line in f.split('\n') if '#' in line}
                if header == fheader:
                    # specid = spec['specid']
                    break
        else:
            print("ERROR: No uploaded spectra for %s" % sourcename)
            return False

    # Update obsdate
    if ':' not in header['obsdate']:
        obsdate = header['obsdate'].strip() + " " + \
                  header['obstime'].strip().split('.')[0]

    # PYSEDM_REPORT
    # must be posted after the SNID plot or else it'll be overwritten
    try:
        # TODO use os.path.dir or something
        pysedm_report = glob(fname.replace('spec',
                                           'pysedm_report').replace('.txt',
                                                                    '.png'))[0]
        pr_posted = add_spec_attachment(header['name'], 'pysedm_report',
                                        pysedm_report, cred,
                                        obsdate=obsdate,
                                        reducedby=reducedby,
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
        if not add_spec_autoannot(header['name'], header['snidmatch' + key],
                                  'AUTO_SNID_' + key, dtypes[key], cred,
                                  obsdate=obsdate,
                                  reducedby=reducedby,
                                  testing=testing):
            return False

    if pr_posted:
        return True  # we already have an attachment comment so don't overwrite

    # SNID PLOT
    # NOTE: this makes a major assumption about the naming scheme of snid plots
    image_filename = fname.replace('.txt',
                                      '_{}.png'.format(header['snidmatchtype']))
    if not glob(image_filename):
        return False
    add_spec_attachment(header['name'], 'AUTO_SNID_plot', image_filename,
                        cred,  obsdate=obsdate,
                        reducedby=reducedby, testing=testing)

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
