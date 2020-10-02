import json
import glob
import requests
import subprocess
import argparse
import os
import datetime
import sys
from marshals.interface import api, update_status_request
try:
    from fritz_commenter import add_SNID_pysedm_autoannot as add_annots
except ImportError:
    from fritz.fritz_commenter import add_SNID_pysedm_autoannot as add_annots

# Path constants
add_target_url = 'http://private.caltech.edu/api/'

pharos_spec_dir = '/scr2/sedmdrp/redux/'
pharos_phot_dir = '/scr2/sedmrp/redux/phot/'

fritz_base_url = 'http://private.caltech.edu/'
fritz_inst_url = fritz_base_url + 'UNKNOWN'     # change request params?
fritz_stat_url = fritz_base_url + 'followup_request/request_id'
fritz_spec_url = fritz_base_url + 'api/spectrum/ascii'
fritz_phot_url = fritz_base_url + 'api/photometry'
fritz_view_source_url = fritz_base_url + 'source'
fritz_token = 'cf127f93-19ef-4ba2-a692-b754c16412b8'

default_id = 37
instrument_id = 2
telescope_id = 37


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
                return_dict[k] = float(outstr)
            elif v.upper() == 'OBSDATE':
                date_str = out.split(sep, 1)[-1].split()[-1]
                out = subprocess.check_output('grep OBSTIME %s' % inputfile,
                                              shell=True,
                                              universal_newlines=True)
                date_str += "T" + out.split(sep, 1)[-1].split()[-1]
                date_str = date_str.split('.')[0] + "Z"
                return_dict[k] = date_str
            else:
                return_dict[k] = out.split(sep, 1)[-1].split()[-1].strip()
        except subprocess.CalledProcessError:
            print("Not found: %s" % k)

    return return_dict


def upload_phot(phot_file, instrument_id=65, request_id='', testing=False):
    """

    :param phot_file:
    :param instrument_id:
    :param request_id:
    :param testing:
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

        if testing:
            ret = "TESTING upload_phot(): no data sent to marshal"
        else:
            json_file = open('photometryExample.txt', 'r')
            ret = requests.post(fritz_phot_url, files={'jsonfile': json_file})
            json_file.close()
        return ret


def upload_spectra(spec_file, request_id=None, sourceid=None, inst_id=2,
                   observer='SEDmRobot', check_quality=True, min_quality=2,
                   testing=False):
    """
    Add spectra to the fritz marshal.  If the fill_by_file is selected then
    most of the keywords will be filled from the spectra file itself.  If
    the request has been canceled for some reason then it will not be possible
    to update the request
    
    :param spec_file:
    :param inst_id:
    :param request_id:
    :param observer:
    :param sourceid:
    :param check_quality:
    :param min_quality:
    :param testing:

    :return: 
    """
    if not request_id or not sourceid:
        print("Can't update without a request id and a source id")
        return None

    # Create the mandatory keyword dictionary payload
    keywords_dict = {   # 'reduced_by': 'REDUCER',
                     'observed_at': 'OBSDATE',
                     'quality': 'QUALITY'}
    if '_SEDM' in spec_file:
        keywords_dict.update({'observed_at': 'OBSUTC'})

    submission_dict = get_keywords_from_file(spec_file, keywords_dict)

    # Retrieve quality
    quality = int(submission_dict['quality'])
    del submission_dict['quality']
    # Check the quality if it is smaller or equal to min quality
    if check_quality and quality > min_quality:
        print("Spectra quality does not pass")
        return None

    # create payload
    with open(spec_file) as sfh:
        contents = sfh.read()
    submission_dict.update({'filename': spec_file,
                            'obj_id': sourceid,
                            'instrument_id': inst_id,
                            # 'followup_request_id': request_id,
                            # 'observer': observer.rstrip().lstrip(),
                            'ascii': contents,
                            'group_ids': []
                            })
    # Are we just testing?
    if testing:
        print(submission_dict)
        ret = {"message": "string", "status": "success", "data": {"id": -1}}
    else:
        # post the spectrum
        ret = api("POST", fritz_spec_url, data=submission_dict)
    print(ret)

    return ret


def read_request(request_file):
    """
    Read in request file to python dictionary
    :param request_file: 
    :return: 
    """
    return json.load(open(request_file, 'r'))


def update_target_by_request_id(request_id, add_spectra=False, spectra_file='',
                                add_status=False, status='Completed',
                                search_db=None, reducedby=None, testing=False):
    """
    Go through the request and find the one that matches the objname
    :param request_id:
    :param add_spectra:
    :param spectra_file:
    :param add_status:
    :param status:
    :param search_db:
    :param reducedby:
    :param testing:
    :return: 
    """

    marshal_id = None
    object_name = None
    username = None
    email = None
    # Return values
    spec_ret = None
    status_ret = None
    return_link = None
    spec_stat = ''
    # Look in the SEDM Db
    if search_db:
        print("Searching SedmDB")

        # Search for target in the database
        try:
            res = search_db.get_from_request(["marshal_id",
                                              "object_id",
                                              "user_id",
                                              "external_id"],
                                             {"id": request_id})[0]
        except IndexError:
            print("Unable to retrieve ids from database")
            return return_link, spec_ret, status_ret
        marshal_id = res[0]
        object_id = res[1]
        user_id = res[2]
        external_id = res[3]
        # is this a Fritz object?
        if external_id != 2:
            print("Not a Fritz object!")
            return return_link, spec_ret, status_ret
        try:
            res = search_db.get_from_object(["name"], {"id": object_id})[0]
        except IndexError:
            print("Unable to retrieve object_name from database")
            return return_link, spec_ret, status_ret
        object_name = res[0]
        try:
            res = search_db.get_from_users(["name", "email"],
                                           {"id": user_id})[0]
        except IndexError:
            print("Unable to retrieve username, email from database")
            return return_link, spec_ret, status_ret
        username = res[0]
        email = res[1]
    else:
        print("no dbase given!")

    # Did we get a marshal ID?
    if marshal_id is None:
        print("Unable to find marshal id for target %s" % object_name)
    else:
        print("Updating target %s using id %d" % (object_name, marshal_id))

        now = datetime.datetime.now()
        ts_str = "%4d%02d%02d %02d_%02d_%02d" % (now.year, now.month,
                                                 now.day, now.hour,
                                                 now.minute,
                                                 now.second)
        if add_spectra:

            spec_ret = upload_spectra(spectra_file, request_id=marshal_id,
                                      sourceid=object_name, testing=testing)
            if not spec_ret:
                spec_stat = 'IFU: Failed ' + ts_str
            else:
                spec_stat = 'IFU: Complete ' + ts_str
                ret_data = spec_ret['data']
                spec_id = ret_data['id']
                annots_posted = add_annots(spectra_file, object_id=object_name,
                                           spec_id=spec_id, testing=testing)
                if annots_posted:
                    print("Annotations successfully posted")
                else:
                    print("Warning: Annotations encountered a problem")

        if add_status:
            status_ret = update_status_request(spec_stat, marshal_id, 'fritz',
                                               testing=testing)

        return_link = fritz_view_source_url + "/%s" % object_name

        print("Send to %s at %s\nRequest status = %s\n%s" %
              (username, email, status, return_link))

    return return_link, spec_ret, status_ret

          
def parse_ztf_by_dir(target_dir, upfil=None, dbase=None, reducedby=None,
                     testing=False):
    """Given a target directory get all files that have ztf or ZTF as base 
       name

       :param target_dir:
       :param upfil:
       :param dbase:
       :param reducedby:
       :param testing:
       """

    if target_dir[-1] != '/':
        target_dir += '/'

    # files = glob.glob('%sZTF*.txt' % target_dir)
    # files += glob.glob('%sztf*.txt' % target_dir)
    # files += glob.glob('%sspec_*ZTF*.txt' % target_dir)

    # list of all spectra in directory
    fls = glob.glob('%sspec_*.txt' % target_dir)
    # scrape out unneeded files or find upfil in list
    files = []
    for fi in fls:
        # are we uploading a specific file?
        if upfil is not None:
            # is this our file?
            if upfil in fi:
                files.append(fi)
            else:
                continue
        # uploading all files in directory
        else:
            # skip uncalibrated spectra
            if "notfluxcal" in fi:
                print("Not flux calibrated: %s" % fi)
                continue
            # skip contsep extractions
            if "contsep" in fi:
                print("Not uploading contsep extraction yet: %s" % fi)
                continue
            # add all others
            files.append(fi)

    started = os.path.exists(os.path.join(target_dir, "report_fritz.txt"))
    out = open(target_dir + "report_fritz.txt", "a")
    if not started:
        out.write("\nZTF Upload report for %s started on %s\n\n" %
                  (target_dir.split('/')[-2],
                   datetime.datetime.now().strftime("%c")))
    for fi in files:
        # Has it already been uploaded?
        if os.path.exists(fi.split('.')[0] + ".upl"):
            print("Already uploaded: %s" % fi)
            continue

        # Extract request ID
        req_id = subprocess.check_output(('grep', 'REQ_ID', fi),
                                         universal_newlines=True)
        req_id = req_id.split(':', 1)[-1].strip()
        if not req_id:
            print("No REQ_ID found: %s" % fi)
            continue
        # Extract object name
        tname = fi.split('_ifu')[-1].split('_')[4:]
        if len(tname) > 1:
            objname = '_'.join(tname).split('.txt')[0]
        else:
            objname = tname[0].split('.txt')[0]
        # Extract observation id
        fname = os.path.basename(fi)
        if 'ifu' in fname:
            obs_id = ":".join(fname.split('_ifu')[-1].split('_')[1:4])
        elif 'rc' in fname:
            obs_id = ":".join(fname.split('_rc')[-1].split('_')[1:4])
        else:
            obs_id = "..:..:.."
        # Are we uploading only one file?
        if upfil is not None:
            # if this is not the file, skip
            if upfil not in fi:
                continue
        # Upload
        r, spec, stat = update_target_by_request_id(req_id, add_status=True,
                                                    status='Completed',
                                                    add_spectra=True,
                                                    spectra_file=fi,
                                                    search_db=dbase,
                                                    reducedby=reducedby,
                                                    testing=testing)
        # Mark as uploaded
        os.system("touch " + fi.split('.')[0].replace(" ", "\ ") + ".upl")
        # log upload
        out.write("%s %s: " % (obs_id, objname))
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
        if r:
            print("URL: " + r)
            out.write("%s\n" % r)
        else:
            print("URL: None")
            out.write("None\n")

    # Close log file
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
                         
Uploads results to the fritz marshal.
                         
""",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('indate', type=str, default=None,
                        help='input directory date (UT date as YYYYMMDD)')
    parser.add_argument('--data_file', type=str, default=None,
                        help='Data file to upload.')
    parser.add_argument('--reducedby', type=str, default=None,
                        help='reducer (defaults to auto)')
    parser.add_argument('--testing', action="store_true", default=False,
                        help='Do not actually post to marshal (for testing)')
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
        # Use the database
        import db.SedmDb
        sedmdb = db.SedmDb.SedmDB()
        parse_ztf_by_dir(srcdir, upfil=args.data_file, dbase=sedmdb,
                         reducedby=args.reducedby, testing=args.testing)
