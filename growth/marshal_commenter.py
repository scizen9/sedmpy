import re, requests
from glob import glob
from getpass import getpass
from pprint import pprint

auth = ('adugas', getpass())
sourcelist = requests.get("http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi", data={'programidx':-1}, auth=auth).json()

def get_missing_info(ztfname, obsdate, sourceid, sourcelist, specid):
    '''#TODO this is currently being called four times per object, which may include a lot of downloading the same thing repeatedly.
    have it save everything to a dictionary or something that gets passed around? posssssibly as a global?
    '''
    if not sourceid:
        if not sourcelist:
            sourcelist = requests.get("http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi", data={'programidx':-1}, auth=auth).json()
        try:
            sourceid = [source['id'] for source in sourcelist if source['name'] == ztfname][-1]
        except IndexError:
            print("you probably don't have permissions to view this source")
            raise

    if not specid:
        source_summary = requests.get("http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi?sourceid={}".format(sourceid), auth=auth).json() # why can't I pass the sourceid as data?
        try:
            specid = [spec['specid'] for spec in source_summary['uploaded_spectra'] if spec['obsdate'].replace('-', '') == obsdate.replace('-', '') and spec['instrumentid'] == 65 and spec['reducedby'].strip() == 'auto'][-1]
        except IndexError as e:
            print(e)
            pprint([(spec['reducedby'], spec['obsdate']) for spec in source_summary['uploaded_spectra'] if spec['instrumentid'] == 65])
            raise
    return sourceid, sourcelist, specid
            
    

def add_spec_attachment(ztfname, comment, filename, auth, sourceid=None, specid=None, obsdate=None, sourcelist=None):
    '''
    adds a comment (not autoannotation) with attachment to a particular SEDM spectrum on the view_spec page, which will also appear elsewhere.
    
    ztfname: 'ZTF18aaaaaa', for example
    comment: <str>
    sourceid: <int> around ~4000
    specid: <int> around ~1700. faster if you provide it
    obsdate: 'YYYYMMDD', necessary for finding the correct spectrum. Assumes exactly 1 SEDM spectrum that night
    sourcelist: list of sources in your program, from list_program_sources.cgi. If you don't know sourceid it's faster if you provide it
    
    return: True if success, False if not
    '''

    sourceid, sourcelist, specid = get_missing_info(ztfname, obsdate, sourceid, sourcelist, specid)

    with open(filename, 'rb') as att:
        r = requests.post("http://skipper.caltech.edu:8080/cgi-bin/growth/add_spec.cgi", auth=auth,
                      data={'comment':   comment,
                            'tablename': 'spec',
                            'tableid':    specid,
                            'camefrom':  'view_spec',
                            'name':       ztfname,
                            'sourceid':   sourceid,
                            'specid':     specid,
                            'commit':    'yes'
                           }, files={"attachment":att})

    if 'Updating Database' in r.text or 'copied file successfully' in r.text:
        return True
    else:
        print('error submitting comment')
        print(r.text)
        raise Exception
        return False
        
def add_spec_autoannot(ztfname, value, annot_type, datatype, auth, sourceid=None, specid=None, obsdate=None, sourcelist=None):
    '''
    adds an autoannotation without attachment to a particular SEDM spectrum on the view_spec page, which will also appear elsewhere.
    
    ztfname: 'ZTF18aaaaaa', for example
    comment: <str>
    datatype: 'STRING' or 'FLOAT' or 'BOOL'
    auth: (growth username, growth password) as (<str>, <str>)
    sourceid: <int> around ~4000
    specid: <int> around ~1700. faster if you provide it
    obsdate: 'YYYYMMDD' or 'YYYY-MM-DD', necessary for finding the correct spectrum if no specid. Assumes exactly 1 SEDM spectrum that night
    sourcelist: list of sources in your program, from list_program_sources.cgi. If you don't know sourceid it's faster if you provide it
    
    return: True if success, False if not
    '''

    sourceid, sourcelist, specid = get_missing_info(ztfname, obsdate, sourceid, sourcelist, specid)

    r = requests.post("http://skipper.caltech.edu:8080/cgi-bin/growth/add_spec.cgi", auth=auth,
                  data={'comment':    value,
                        'datatype':   datatype,
                        'tablename': 'spec',
                        'type':       annot_type,
                        'tableid':    specid,
                        'table':     'autoannotations',
                        'camefrom':  'view_spec',
                        'name':       ztfname,
                        'sourceid':   sourceid,
                        'specid':     specid,
                        'commit':    'yes'
                       })

    if 'Updating Database' in r.text or 'copied file successfully' in r.text:
        return True
    else:
        print('error submitting comment')
        print(r.text)
        return False

def add_SNID_autoannot(filename, auth, sourcelist=None):
    '''
    if z < 0.3 and rlap > 5.0
        adds autoannotations with SNID rlap, z, type, etc
        adds a comment with the SNID plot attached
    
    filename: '*ZTF18aaaaaaa.txt' that has a bunch of"# SNIDMATCH[something]: [val]" in the header
    auth: ('username', 'password')
    sourcelist: dictionary of sources in your program, from list_program_sources.cgi
    
    returns: True if all four comments/attachments works, False (and it'll exit early) otherwise
    '''
    # TODO if we look at the _snid.output we can get more info, eg phase which would be pretty helpful to me personally
    file_ext = filename.split('.')[-1]
    assert file_ext == 'txt' or file_ext == 'ascii'

    with open(filename) as f:
        header = {line.split(':')[0][1:].strip().lower(): line.split(':', 1)[-1].strip() for line in f if line[0] == '#'}
    header['snidmatchmatch'] = '{snidmatchtype}-{snidmatchsubtype}'.format(**header)


    if header['snidmatchtype'].lower() == 'none':
        print('no match')
        return False
        
    elif float(header['snidmatchrlap']) < 5:
        print('bad rlap, only {}'.format(header['rlap']))
        return False
        
    elif (header['snidmatchtype'][0] == 'I') and not (0.01 <= float(header['snidmatchredshift']) <= 0.3):
        print('bad redshift, {snidmatchredshift} for {snidmatchtype}'.format(**header))
        return False

    dtypes = {'match': 'STRING', 'rlap': 'FLOAT', 'redshift': 'FLOAT'}
    for key in dtypes:
        if not add_spec_autoannot(header['name'], header['snidmatch' + key], 'AUTO_SNID_' + key, 
                           dtypes[key], auth, obsdate=header['obsdate'], sourcelist=sourcelist):
            return False
        
    # NOTE: this makes a major assumption about the naming scheme of snid plots
    image_filename = filename.replace('.txt', '_{}.png'.format(header['snidmatchtype']))
    assert glob(image_filename)
    if not add_spec_attachment(header['name'], 'AUTO_SNID_plot', image_filename, auth, 
                        obsdate=header['obsdate'], sourcelist=sourcelist):
        return False
    
    return True

if __name__ == "__main__":
    auth = (raw_input('GROWTH marshal username:'), getpass())
    sourcelist = requests.get("http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi", data={'programidx':-1}, auth=auth).json()
    
    successes = []
    for filename in glob('*ZTF?????????.txt'):
        if add_SNID_autoannot(filename, auth, sourcelist=sourcelist):
            print('success!')
            successes.append(re.findall(r'ZTF\d{2}[a-z]{7}', filename)[-1])
            break
    pprint(successes)

