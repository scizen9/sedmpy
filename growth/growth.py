import json
import glob
import requests
import subprocess
import argparse
import os
import datetime
import sys
try:
    from marshal_commenter import add_SNID_pysedm_autoannot as add_annots
    from marshal_commenter import sourcelist, auth
except ImportError:
    from growth.marshal_commenter import add_SNID_pysedm_autoannot as add_annots
    from growth.marshal_commenter import sourcelist, auth


# Path constants
add_target_url = 'http://nera.palomar.caltech.edu/cgi-bin/' \
                 'telescopes/p60/sedm/exFollowup.cgi'

pharos_spec_dir = '/scr2/sedmdrp/redux/'
pharos_phot_dir = '/scr2/sedm/phot/'

growth_base_url = 'http://skipper.caltech.edu:8080/cgi-bin/growth/'
growth_inst_url = growth_base_url + 'update_followup_config.cgi'
growth_stat_url = growth_base_url + 'update_followup_status.cgi'
growth_spec_url = growth_base_url + 'add_spec_auto.cgi'
growth_phot_url = growth_base_url + 'edit_phot_auto.cgi'
growth_view_source_url = growth_base_url + 'view_source.cgi?'

default_id = 65
if auth:
    user, pwd = auth
else:
    user = None
    pwd = None


def write_json_file(pydict, output_file):
    """
    Write the python dictionary to a json file
    :param pydict: 
    :param output_file: 
    :return: json file path
    """

    json_file = open(output_file, 'w')
    json_file.write(json.dumps(pydict))
    json_file.close()

    return output_file


def timestamp():
    """
    UTC timestamp.  Use this when saving files
    :return: 
    """
    return datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S")


def send_instrument_configuration(instrument_id="",
                                  options_programs=None,
                                  post_url="", save=False):
    """
    Update the instrument configuration.  Updating this file will change the
    available options shown on the growth marshal page

    :param instrument_id:
    :param options_programs:
    :param post_url:
    :param save:
    :return:
    """

    if not instrument_id:
        instrument_id = default_id

    if not post_url:
        post_url = add_target_url

    # 1. Create the main json dictionary file
    inst_config = {
        'instrument_id': instrument_id,
        'post_url': post_url,
        'options': []
    }

    # 2. Check the options dictionary for values if none use default for SEDm
    if options_programs:
        inst_config['options'] = options_programs
    else:
        fourshot = {
            'name': 'Followup',
            'type': 'select',
            'value': 'Four Shot (r,g,i,u)'
        }

        three_shot = {
            'name': 'Followup',
            'type': 'select',
            'value': 'Three Shot (r,g,i)'
        }

        ifu = {
            'name': 'Followup',
            'type': 'select',
            'value': 'IFU'
        }

        ifu_fourshot = {
            'name': 'Followup',
            'type': 'select',
            'value': 'Fourshot + IFU'
        }

        check1 = {"name": "Filters", "type": "check", "value": "r"}

        check2 = {"name": "Filters", "type": "check", "value": "g"}

        check3 = {"name": "Filters", "type": "check", "value": "i"}

        check4 = {"name": "Filters", "type": "check", "value": "ifu"}

        inst_config['options'] = [fourshot, three_shot, ifu, ifu_fourshot,
                                  check1, check2, check3, check4]

    # 3. Create a json file with the request and read it in to memory
    json_file = open(write_json_file(inst_config, 'config.txt'), 'r')

    # 4. Send json file request to growth marshal
    ret = requests.post(growth_inst_url,
                        auth=(user, pwd),
                        files={'jsonfile': json_file})

    # 5. Close the file and save it if needed
    json_file.close()

    if not save:
        os.remove('config.txt')

    # 6. Return the request response for user to determine if it was a success
    return ret


def update_request(status, request_id, instrument_id='',
                   output_dir='targets', filename='', save=True):
    """
    Function to update the status of any request as long as it has
    not been deleted. The new status will show up on the status section
    of the request on the growth marshal.
    
    :param filename: 
    :param status: 
    :param instrument_id: 
    :param request_id: 
    :param output_dir: 
    :param save: 
    :return: 
    """

    # 1. Make sure that we have the required fields
    if not request_id:
        print("Can't update a request without the request id")
        return False

    if not instrument_id:
        instrument_id = default_id

    if not filename:
        filename = "%s_%s.txt" % (str(request_id), timestamp())

    output_file = os.path.join(output_dir, filename)

    # 2. Create the new status dictionary
    status_config = {'instrument_id': instrument_id,
                     'request_id': request_id,
                     'new_status': status}

    # 3. Write and read in the json file to memory
    json_file = open(write_json_file(status_config, output_file), 'r')

    # 4. Send the request, close the file, and save if needed
    ret = requests.post(growth_stat_url, auth=(user, pwd),
                        files={'jsonfile': json_file})
    json_file.close()

    if not save:
        os.remove(output_file)

    # 5. Print the request response for the user
    print(ret)

    return True


def get_keywords_from_file(inputfile, keywords, sep=':'):
    """
    Get keywords from file.  It is dependent on files having a specific format
    where the keyword is on the left and the value on the right by some common
    seperator

    :param inputfile: input spectrum text file
    :param keywords: dictionary of keywords to get from file
    :param sep: separator character
    :return:
    """
    return_dict = {}

    for k, v in keywords.items():
        try:
            out = subprocess.check_output('grep %s %s' % (v, inputfile),
                                          shell=True, universal_newlines=True)
            if k.upper() == 'EXPTIME':
                outstr = out.split(sep, 1)[-1]
                return_dict[k] = int(outstr.split('.')[0])
            elif v.upper() == 'OBSDATE':
                date_str = out.split(sep, 1)[-1]
                out = subprocess.check_output('grep OBSTIME %s' % inputfile,
                                              shell=True,
                                              universal_newlines=True)
                date_str += " " + out.split(sep, 1)[-1]
                return_dict[k] = date_str
            else:
                return_dict[k] = out.split(sep, 1)[-1]
        except subprocess.CalledProcessError:
            print("Not found: %s" % k)

    return return_dict


def upload_phot(phot_file, instrument_id=65, request_id=''):
    """

    :param phot_file:
    :param instrument_id:
    :param request_id:
    :return:
    """

    with open(phot_file, 'r') as photometryFile:
        photometry = photometryFile.read()

    photometry = photometry.split('\n')
    column_names = photometry[0].split(',')
    photometry = photometry[1:-1]

    photometry_list = []
    for entry in photometry:
        new_dict = {}
        photometry_point = entry.split(',')
        for index, column in enumerate(column_names):
            data = photometry_point[index]
            if '"' in data:
                data = data.replace('"', '')
            elif data == 't':
                data = True
            elif data == 'f':
                data = False
            else:
                data = float(data)
            if data == 'None':
                data = None
            new_dict[column] = data
        photometry_list.append(new_dict)

        submission_dict = {
            'photometry_list': photometry_list, 'instrument_id': instrument_id,
            'request_id': request_id
        }

        json_file = open('photometryExample.txt', 'w')
        json_file.write(json.dumps(submission_dict))
        json_file.close()

        json_file = open('photometryExample.txt', 'r')
        ret = requests.post(growth_phot_url, auth=(user, pwd),
                            files={'jsonfile': json_file})
        json_file.close()
        return ret


def upload_spectra(spec_file, fill_by_file=False, instrument_id=65,
                   request_id='', exptime=3600, observer='SEDmRobot',
                   reducedby="auto", obsdate="", output_dir='targets/',
                   format_type='ascii',
                   check_quality=True, quality=1, min_quality=2):
    """
    Add spectra to the growth marshal.  If the fill_by_file is selected then
    most of the keywords will be filled from the spectra file itself.  If
    the request has been canceled for some reason then it will not be possible
    to update the request
    
    :param spec_file:
    :param fill_by_file:
    :param instrument_id:
    :param request_id: 
    :param exptime:
    :param observer:
    :param reducedby: 
    :param obsdate: 
    :param output_dir: 
    :param format_type:
    :param check_quality:
    :param quality:
    :param min_quality:

    :return: 
    """
    if not request_id:
        print("Can't update without a request id")
        return False

    basefile = os.path.basename(spec_file)
    output = os.path.join(output_dir, basefile)
    output_file = output.replace('.txt', '.json')

    # 1. Create the mandatory keyword dictionary
    if not fill_by_file:
        submission_dict = {'exptime': exptime,
                           'obsdate': obsdate,
                           'reducedby': reducedby}
    else:
        keywords_dict = {'reducedby': 'REDUCER',
                         'obsdate': 'OBSDATE',
                         'exptime': 'EXPTIME',
                         'quality': 'QUALITY'}
        if '_SEDM' in spec_file:
            keywords_dict.update({'obsdate': 'OBSUTC'})

        submission_dict = get_keywords_from_file(spec_file, keywords_dict)

        quality = int(submission_dict['quality'])
        del submission_dict['quality']

    # print(type(format_type), type(request_id), type(instrument_id))
    submission_dict.update({'format': format_type.rstrip().lstrip(),
                           'instrument_id': instrument_id,
                            'request_id': request_id,
                            'observer': observer.rstrip().lstrip(),
                            'proprietary': 1})

    # 1a. [Optional] Check the quality if it is smaller or equal to min quality
    # then upload spectra
    
    if check_quality and quality > min_quality:
        print("Spectra quality does not pass")
        return False
    
    # 2. Open the configuration and spec file for transmission
    json_file = open(write_json_file(submission_dict, output_file), 'r')
    upfile = open(spec_file, 'r')

    # 3. Send the request
    ret = requests.post(growth_spec_url, auth=(user, pwd),
                        files={'jsonfile': json_file, 'upfile': upfile})
    print(ret)

    # 4. Close files and send request response
    upfile.close()
    json_file.close()

    return True


def read_request(request_file):
    """
    Read in request file to python dictionary
    :param request_file: 
    :return: 
    """
    return json.load(open(request_file, 'r'))


def update_target_by_object(objname, add_spectra=False, spectra_file='',
                            add_status=False, status='Completed',
                            pull_requests=False, request_id=None,
                            add_phot=False, phot_file='', search_db=None,
                            target_dir='requests/', target_base_name='request'):
    """
    Go through the request and find the one that matches the objname
    :param objname:
    :param add_spectra:
    :param spectra_file:
    :param add_status:
    :param status:
    :param pull_requests:
    :param request_id:
    :param add_phot:
    :param phot_file:
    :param search_db:
    :param target_dir:
    :param target_base_name:
    :return: 
    """

    marshal_id = None
    object_name = None
    username = None
    email = None
    # Look in the SEDM Db
    if search_db:
        if request_id:
            print("Searching SedmDB")

            # 1a. Search for target in the database
            res = search_db.get_from_request(["marshal_id",
                                              "object_id",
                                              "user_id"],
                                             {"id": request_id})[0]
            marshal_id = res[0]
            object_id = res[1]
            user_id = res[2]
            res = search_db.get_from_object(["name"], {"id": object_id})[0]
            object_name = res[0]
            res = search_db.get_from_users(["name", "email"],
                                           {"id": user_id})[0]
            username = res[0]
            email = res[1]
            out_dir = ''
        else:
            print("No request id provided")
    else:
        # 1. Start by looking at all files in the target directory
        # this is currently the directory that holds all incoming request
        # dictionary from growth.
        if pull_requests:
            print("Gathering all request files")
            import growth.growth_watcher as growth_watcher
            growth_watcher.pull_request_from_remote()
        match_list = []
        files = glob.glob('%s%s*' % (target_dir, target_base_name))

        for i in files:
            targ = read_request(i)
            # Add any request file that contains the objname.  There may be more
            # than one if an update to the original request has been sent.
            # Ignore any request with the status delete
            # as we can not update those
            if (targ['sourcename'].lower() == objname.lower() and
                    targ['status'] != 'delete'):
                match_list.append(targ)
         
        if len(match_list) == 1:
            print(match_list)
            target = match_list[0]
            marshal_id = target['requestid']
            object_name = target['sourcename']
            username = target['username']
            print("Uploading target %s files" % objname)
        elif len(match_list) == 0:
            print("Could not match name with any request file")
        else:
            print('Multiple matches have been made for target: %s' % objname)
            request_id_list = []
            for j in match_list:
                request_id_list.append(j['requestid'])

            # If all the request id matches then there is no problem and we can
            # send data
            if all(x == request_id_list[0] for x in request_id_list):
                target = match_list[0]

            else:
                # TODO: Handle case in which an update has been sent
                # For now we just use highest value
                request_id = sorted(request_id_list)[-1]
                print(request_id_list)
                for j in match_list:
                    if j['requestid'] == request_id:
                        target = j
                        break
            marshal_id = target['requestid']
            object_name = target['sourcename']
            username = target['username']
            out_dir = 'targets/'

    # Return values
    spec_ret = None
    phot_ret = None
    status_ret = None
    return_link = None
    spec_stat = ''
    phot_stat = ''
    # Did we get a marshal ID?
    if marshal_id is None:
        print("Unable to find marshal id for target %s" % objname)
    elif marshal_id <= 0:
        print("Not an object from the marshal")
    else:
        print("Updating target %s using id %d" % (objname, marshal_id))

        if add_spectra:

            spec_ret = upload_spectra(spectra_file, fill_by_file=True,
                                      request_id=marshal_id,
                                      output_dir=out_dir)
            if not spec_ret:
                spec_stat = 'IFU: Failed'
            else:
                spec_stat = 'IFU: Complete'
                annots_posted = add_annots(spectra_file, auth, 
                                           sourcelist=sourcelist)
                if annots_posted:
                    print("Annotations successfully posted")
                else:
                    print("Warning: Annotations encountered a problem")

        if add_phot:
            phot_ret = upload_phot(phot_file, request_id=marshal_id)
            if not phot_ret:
                phot_stat = 'RC: Failed'
            else:
                phot_stat = 'RC: Complete'

        if add_status:
            if add_spectra and add_phot:
                status = spec_stat + ', ' + phot_stat
            elif add_spectra:
                status = spec_stat
            elif add_phot:
                status = phot_stat
            status_ret = update_request(status, request_id=marshal_id,
                                        output_dir=out_dir)

        return_link = growth_view_source_url + "name=%s" % object_name

        print("Send to %s at %s\nRequest status = %s\n%s" %
              (username, email, status, return_link))

    return return_link, spec_ret, phot_ret, status_ret

          
def parse_ztf_by_dir(target_dir, upfil=None, dbase=None):
    """Given a target directory get all files that have ztf or ZTF as base 
       name

       :param target_dir:
       :param upfil:
       :param dbase:
       """

    if target_dir[-1] != '/':
        target_dir += '/'

    # files = glob.glob('%sZTF*.txt' % target_dir)
    # files += glob.glob('%sztf*.txt' % target_dir)
    # files += glob.glob('%sspec_*ZTF*.txt' % target_dir)
    files = glob.glob('%sspec_*.txt' % target_dir)

    started = os.path.exists(os.path.join(target_dir, "report_ztf.txt"))
    out = open(target_dir + "report_ztf.txt", "a")
    if not started:
        out.write("\nZTF Upload report for %s started on %s\n\n" %
                  (target_dir.split('/')[-2],
                   datetime.datetime.now().strftime("%c")))
    pr = True
    for fi in files:
        # Has it already been uploaded?
        if os.path.exists(fi.split('.')[0] + ".upl"):
            print("Already uploaded: %s" % fi)
            continue

        # Is it flux calibrated?
        if "notfluxcal" in fi:
            print("Not flux calibrated: %s" % fi)
            continue

        # Extract request ID
        req_id = subprocess.check_output(('grep', 'REQ_ID', fi),
                                         universal_newlines=True)
        req_id = req_id.split(':', 1)[-1].strip()
        if not req_id:
            print("No REQ_ID found: %s" % fi)
            continue
        # Extract object name
        if "spec" in fi:
            objname = os.path.basename(fi).split('_')[-1].split('.')[0]
        else:
            objname = os.path.basename(fi).split('_')[0]
        # Are we uploading only one file?
        if upfil is not None:
            # if this is not the file, skip
            if upfil not in fi:
                continue
        # Upload
        r, spec, phot, stat = update_target_by_object(objname,
                                                      add_status=True,
                                                      status='Completed',
                                                      add_spectra=True,
                                                      spectra_file=fi,
                                                      request_id=req_id,
                                                      search_db=dbase,
                                                      pull_requests=pr)
        # Mark as uploaded
        os.system("touch " + fi.split('.')[0] + ".upl")
        # Only need to pull requests the first time
        pr = False
        # log upload
        out.write(" %s: " % objname)
        # Was a spectrum uploaded?
        if spec:
            out.write("OK ")
        else:
            out.write("NO ")
        # Was status updated?
        if stat:
            out.write("OK ")
        else:
            out.write("NO ")
        print("URL: " + r)
        out.write("%s\n" % r)

    # Close log file
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
                         
Uploads results to the growth marshal.
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('indate', type=str, default=None,
                        help='input directory date (UT date as YYYYMMDD)')
    parser.add_argument('--data_file', type=str, default=None,
                        help='Data file to upload.')
    parser.add_argument('--no_usedb', action="store_true", default=False,
                        help='Do not use SEDM database')
    args = parser.parse_args()

    # Check environment
    try:
        reddir = os.environ["SEDMREDUXPATH"]
    except KeyError:
        print("please set environment variable SEDMREDUXPATH")
        sys.exit(1)

    # Get source dir
    if args.indate:
        utc = args.indate
    else:
        utc = datetime.datetime.utcnow().strftime("%Y%m%d")
    srcdir = reddir + '/' + utc + '/'

    # Check source dir
    if not os.path.exists(srcdir):
        print("Dir not found: %s" % srcdir)
    else:
        print("Uploading from %s" % srcdir)

        # Will we use the database?
        if args.no_usedb:
            # Check requests, targets dirs
            reqdir = srcdir + 'requests'
            trgdir = srcdir + 'targets'
            if not os.path.exists(reqdir):
                os.mkdir(reqdir)
            if not os.path.exists(trgdir):
                os.mkdir(trgdir)
            parse_ztf_by_dir(srcdir, upfil=args.data_file)
        else:
            import db.SedmDb
            sedmdb = db.SedmDb.SedmDB()
            parse_ztf_by_dir(srcdir, upfil=args.data_file, dbase=sedmdb)
