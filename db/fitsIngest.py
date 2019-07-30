# Script to add all fits files into the database
from astropy.io import fits
import glob
import os
import json
import datetime
import pprint
from db.SedmDb import SedmDB
import astropy.units as u
from astropy.coordinates import SkyCoord
import re
from string import Template

class ingestNight:
    """
    Ingest a directory of files or given a specific night ingest all raw and
    processed images for that night.
    """
    def __init__(self, raw_dir='/scr2/sedm/raw/',
                 ifu_proc_dir='/scr/rsw/sedm/data/redux/',
                 rc_proc_dir='/scr2/sedm/phot/',
                 ifu_prexfix='ifu', rc_prefix='rc',
                 ingest_file='ingest.json',
                 problem_file='/scr/rsw/problem_files.json'):
        """

        :param raw_dir:
        :param ifu_proc_dir:
        :param rc_proc_dir:
        """
        self.raw_dir = raw_dir
        self.ifu_proc_dir = ifu_proc_dir
        self.rc_proc_dir = rc_proc_dir
        self.ifu_prefix = ifu_prexfix
        self.rc_prefix = rc_prefix
        self.ingest_file = ingest_file
        self.problem_file = problem_file
        self.required_obs_header = ['obj_id', 'req_id', 'mjd_obs',
                                    'airmass', 'exptime',
                                    'lst', 'ra', 'dec', 'tel_az',
                                    'tel_el', 'tel_pa', 'ra_off',
                                    'dec_off']
        """'airmass_end'(float),
                'parang'(float),
                'parang_end'(float),
                'ra_off'(float),
                'dec_off'(float),
                'imtype'(str),
                'time_elapsed'(float)
                'filter'(str), options - 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b', 'NA'
                'camera'(str)"""

        self.optional_obs_header = {'endair': 'airmass_end', 'elaptime':
                                    'time_elapsed', 'imgtype': 'imtype',
                                    'filter': 'filter', 'tel_pa': 'parang',
                                    'end_pa': 'parang_end', 'ser_no': 'camera'}


        self.telescope_stats_dict = dict(BITPIX="int", NAXIS="int", NAXIS1="int", NAXIS2="int", EXTEND="bool",
                                         BZERO="int", EXPTIME="float", ADCSPEED="float", TEMP="float", GAIN_SET="float",
                                         ADC="float", MODEL="str", INTERFC="str", SNSR_NM="str", SER_NO="str",
                                         TELESCOP="str", LST="str", MJD_OBS="float", JD="float", APPEQX="float",
                                         EQUINOX="float", TEL_HA="str", RA="str", TEL_RA="str", DEC="str",
                                         TEL_DEC="str", TEL_AZ="float", TEL_EL="float", AIRMASS="float", TEL_PA="float",
                                         RA_RATE="float", DEC_RATE="float", RA_OFF="float", DEC_OFF="float",
                                         TELHASP="float", TELDECSP="float", RA_REFR="float", DEC_REFR="float",
                                         FOCPOS="float", IFUFOCUS="float", IFUFOC2="float", DOMEST="str", DOMEMO="str",
                                         DOME_GAP="float", DOMEAZ="float", WSCRMO="str", TELCONT="str", LAMPSTAT="str",
                                         LAMPCUR="float", HG_LAMP="str", XE_LAMP="str", CD_LAMP="str", TELPOWST="str",
                                         OILSTAT="str", WEASTAT="str", SUNSTAT="str", REMOTST="str", TELRDST="str",
                                         HAAX_ST="str", FOCSTAT="str", DEC_AX="str", OBJECT="str", OBJTYPE="str",
                                         IMGTYPE="str", OBJNAME="str", OBJEQX="str", OBJRA="str", OBJDEC="str",
                                         ORA_RAT="float", ODEC_RAT="float", SUNRISE="str", SUNSET="str", TEL_MO="str",
                                         WSCR_EL="str", SOL_RA="time", SOL_DEC="str", WIND_DIR="float", WSP_CUR="float",
                                         WSP_AVG="float", OUT_AIR="float", OUT_HUM="float", OUT_DEW="float",
                                         IN_AIR="float", IN_HUM="float", IN_DEW="float", MIR_TEMP="float",
                                         TOP_AIR="float", PRI_TEMP="float", SEC_TEMP="float", FLO_TEMP="float",
                                         BOT_TEMP="float", MID_TEMP="float", TOP_TEMP="float", WETNESS="float",
                                         FILTER="str", ABPAIR="str", IMGSET="str", NAME="str", UTC="datetime",
                                         OBSDATE="date", OBSTIME="time", P60PRID="str", P60PRNM="str", P60PRPI="str",
                                         EMAIL="str", INSTRUME="str", REQ_ID="int", OBJ_ID="int", GAIN="float",
                                         CAM_NAME="str", CRPIX1="float", CRPIX2="float", CDELT1="float", CDELT2="float",
                                         CTYPE1="str", CTYPE2="str", CRVAL1="float", CRVAL2="float", CCDTEMP="float",
                                         END_SHUT="datetime", ENDAIR="float", ENDDOME="str", END_RA="str",
                                         END_DEC="str", END_PA="str", ELAPTIME="float", OBSERVATION_ID="int",
                                         UTCDAY="int", OBJ_RA_DEG="float", OBJ_DEC_DEG="float",
                                         TEL_RA_DEG="float", TEL_DEC_DEG="float", HA_HOUR_DEG="float")

        self.conversion_dict = {'TEMP': "CCDTEMP", "OBJECT": "OBJNAME"}

        self.standard_file_dict = {
            'adr_fit': 'adr_fit',
            'calibch': 'calibcheck',
            'fluxcal': 'fluxcal',
            'ifu_spa': 'ifu_spaxels',
            'psfprof': 'psfprofile',
            'pysedm_': 'pysedm_report',
            'spaxels': 'spaxels',
            'spec_ap': 'spec_aperture',
            'spec_au': 'spec_auto',
            'verify_': 'verify'
        }

        self.science_file_dict = {
            'flex_sodi': 'flex_sodium',
            'flexuretr': 'flexure',
            'adr_fit_a': 'adr_fit',
            'forcepsf_': 'forcepsf_fitted',
            'forcepsfm': 'forcepsf_model',
            'guider_cr': 'guider_stack',
            'ifu_spaxe': 'ifu_spaxels',
            'maskcrr_b': 'mask_crr_b',
            'psfprofil': 'psfprofile',
            'pysedm_re': 'pysedm_report',
            'spaxels_s': 'spaxels',
            'spec_auto': 'spec_auto',
            'spec_aper': 'spec_aperture',
            'verify_au': 'verify'
        }

        self.calibration_file_dict = {
            'what.list': {'product_type': 'what.list', 'ext_type': 'txt'},
            'Makefile': {'product_type': 'makefile', 'ext_type': 'txt'},
            'bias0.1.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'bias0.1.lst': {'product_type': 'calibration', 'ext_type': 'lst'},
            'bias2.0.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'bias2.0.lst': {'product_type': 'calibration', 'ext_type': 'lst'},
            'dome.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'dome.lst': {'product_type': 'calibration', 'ext_type': 'lst'},
            'Hg.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'Hg.lst': {'product_type': 'calibration', 'ext_type': 'lst'},
            'Cd.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'Cd.lst': {'product_type': 'calibration', 'ext_type': 'lst'},
            'Xe.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'Xe.lst': {'product_type': 'calibration', 'ext_type': 'lst'},
            'bkgd_dome.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'e3d_dome.fits': {'product_type': 'calibration', 'ext_type': 'fits'},
            'snid.param': {'product_type': 'snid_param', 'ext_type': 'txt'},
            'report.txt': {'product_type': 'report', 'ext_type': 'txt'},
            'report_ztf.txt': {'product_type': 'report_ztf', 'ext_type': 'txt'},
        }


        self.dataproduct = Template("""INSERT INTO dataproducts (observation_id,
                                  obsdate, product_name, product_type, camera, basedir, ext_type) 
                                  VALUE 
                                  (${observation_id}, ${obsdate}, ${product_name}, ${product_type}, 
                                  ${camera}, ${basedir}, ${product_type})""")

        self.db = SedmDB(host="pharos", dbname='sedmdb', port=5432)

    def _insert_dataproducts_dict(self, obsdate, observation_id=0, product_name="",
                                  product_type=):
        """

        :param productdict:
        :param observation_id:
        :return:
        """

        query = self.dataproduct.substitute(p)


    def _insert_telescope_stats_dict(self, headerdict, observation_id):
        new_dict = {}

        # Check if their
        print(observation_id, "<--This is the observation id")
        pprint.pprint(headerdict)
        try:
            sql = """INSERT INTO telescope_stats (observation_id) VALUES(%s)""" % observation_id

            self.db.execute_sql(sql)
        except Exception as e:
            print(e)
            print(str(e))
            print("HERE")
            if "already exists" in str(e):
                pass

        # Convert needed entreies
        headerdict['UTCDAY'] = headerdict['UTC'].split(":")[1]
        headerdict['UTC'] = "%s %s" % (headerdict['OBSDATE'], headerdict['UTC'].split(":", 2)[-1])

        for k in ['OBJ', 'TEL_']:
            if headerdict[k+'RA'] and headerdict[k+'DEC']:
                ra = headerdict[k+'RA']
                dec = headerdict[k+'DEC']
                if ":" in ra or ":" in dec:
                    k = k.replace('_',"")
                    c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
                    headerdict[k+'_RA_DEG'] = c.ra.degree
                    headerdict[k+'_DEC_DEG'] = c.dec.degree


        for k, i in headerdict.items():
            print(k, i, type(i))
            if k.upper() in self.conversion_dict:
                k = self.conversion_dict[k]

            if k.upper() not in self.telescope_stats_dict:
                print(k, "Is not in the dictionary")
                continue
            elif k.upper() in ['EXPTIME', 'SIMPILE']:
                print("EXCLUDING:", k.upper())
                continue
            elif self.telescope_stats_dict[k.upper()] == 'str':
                try:
                    new_dict[k.lower()] = str(i)
                except:
                    print("String Exception Class", i)
                    new_dict[k.lower()] = 'StrTypeFailure'
            elif self.telescope_stats_dict[k.upper()] == 'float':
                try:
                    new_dict[k.lower()] = float(i)
                except:
                    print("Float Exception Class", k, i)
                    new_dict[k.lower()] = -9999.999
            elif self.telescope_stats_dict[k.upper()] == 'int':
                try:
                    new_dict[k.lower()] = int(i)
                except:
                    print("Int Exception Class", i)
                    new_dict[k.lower()] = -9999
            elif self.telescope_stats_dict[k.upper()] == 'time':
                try:
                    new_dict[k.lower()] = str(i)
                except:
                    print("Time Exception Class", i)
                    new_dict[k.lower()] = '00:00:00'
            elif self.telescope_stats_dict[k.upper()] == 'date':
                try:
                    new_dict[k.lower()] = str(i)
                except:
                    print("Date Exception Class", i)
                    new_dict[k.lower()] = '2000-00-00'
            elif self.telescope_stats_dict[k.upper()] == 'datetime':
                try:
                    new_dict[k.lower()] = str(i)
                except:
                    print("Datetime Exception Class", i)
                    new_dict[k.lower()] = "2000-01-01 00:00:00"
            elif self.telescope_stats_dict[k.upper()] == 'bool':
                try:
                    print("BOOL Character = ", i, type(i))
                    if isinstance(i, bool):
                        new_dict[k.lower()] = i
                    elif i.lower() in ['n', 'no', 'false']:
                        print("False element")
                        new_dict[k.lower()] = False
                    else:
                        print("True Element")
                        new_dict[k.lower()] = True
                except:
                    print("Bool Exception Class", i)
                    new_dict[k.lower()] = False
            else:
                print(k, "Has no class")
                print(self.telescope_stats_dict[k.upper()])
                continue
            if isinstance(new_dict[k.lower()], str):
                sql = """UPDATE telescope_stats SET %s = '%s' WHERE observation_id=%s""" % (k.lower(), new_dict[k.lower()], observation_id)
            else:
                sql = """UPDATE telescope_stats SET %s = %s WHERE observation_id=%s""" % (k.lower(), new_dict[k.lower()], observation_id)
            self.db.execute_sql(sql)

        return new_dict

    def gather_files(self, directory, prefix='all', suffix='all',
                     recursive=False):
        """
        Given a directory gather the selected files in a directory
        matching the prefix and suffix given.  If recursive is true
        then all the corresponding files in subdirectories will be
        added as well.

        :param recursive:
        :param suffix:
        :param directory:
        :param prefix:
        :return:
        """

        # Generate a search string to be used with the glob command
        if prefix == 'all':
            prefix = ''
        if suffix == 'all':
            suffix = ''

        search_str = os.path.join(directory, prefix)

        search_str += '*'
        search_str += suffix
        print(search_str)

        # Run the search
        files = glob.glob(search_str)

        return files

    def __read_in_json_file(self, jsonFile):
        """

        :return:
        """
        if not os.path.exists(jsonFile):
            return False

        with open(jsonFile) as json_data_file:
            return json.load(json_data_file)

    def add_files_to_db(self):
        """

        :return:
        """

        pass

    def remove_ingested_files(self):
        """

        :return:
        """

        pass

    def ingest_nightly_data_products(self, night, camera="", directory="",
                                     include_subdirectories=True, ignore_raw_files=True):
        """

        :param night:
        :param directory:
        :param include_subdirectories:
        :return:
        """

        if not directory and camera.lower() == 'ifu':
            directory = os.path.join(self.ifu_proc_dir, night)

        # Start by getting all files in a directory
        if camera == 'ifu':
            print("Searching %s for data products" % directory)
            # Start by getting all the observation id for a night
            obs_list = self.db.execute_sql("SELECT id, fitsfile FROM observation "
                                           "WHERE fitsfile LIKE '%%ifu%s%%'"  % night)

            # Get the list of all files in directory tree at given path
            listOfFiles = list()
            for (dirpath, dirnames, filenames) in os.walk(directory):
                listOfFiles += [os.path.join(dirpath, file) for file in filenames]

            # Now look to see if the dataproduct is a fits file and already in the obs_list
            if listOfFiles == 0:
                return "No products to add"

            non_fits_id_list = []
            for f in listOfFiles:
                # Start by checking if it has a fits id associated with it
                print("Parsing file: %s" % f)
                base_file = os.path.basename(f)
                match = re.search(r'\d{8}_\d{2}_\d{2}_\d{2}', base_file)
                if not match:
                    if base_file in self.calibration_file_dict:
                        product_type = self.calibration_file_dict[base_file]['product_type']
                        ext_type = self.calibration_file_dict[base_file]['ext_type']
                    elif base_file[-3:] == 'pkl':
                        product_type = 'pickle_file'
                        ext_type = 'pickle'
                    elif "wavesolution_dispersionmap" in base_file:
                        product_type = 'dispersionmap'
                        ext_type = base_file.split('.')[-1]
                    elif "wavesolution_stats" in base_file:
                        product_type = 'wavesolution_stats'
                        ext_type = base_file.split('.')[-1]
                    elif "_flat3d" in base_file:
                        product_type = 'flat3d'
                        ext_type = base_file.split('.')[-1]
                    elif "%s_Flat.fits" % night == base_file:
                        product_type = "calibration"
                        ext_type = 'fits'
                    elif "%s_dome_stats.txt" % night == base_file:
                        product_type = "dome_stats"
                        ext_type = 'txt'
                    else:
                        print("No match for %s" % f)
                        non_fits_id_list.append(f)
                        continue
#                    self._insert_dataproducts_dict()
                    continue
                # Check if it is a raw fits file
                print("Match '%s' found in file %s" % (match.group(), base_file))
                match_str = match.group()
                raw_filename = "ifu%s.fits" % match_str
                #print(filename, os.path.basename(f))
                if raw_filename == os.path.basename(f):
                    print("MATCH!!!")
                    continue
                elif re.search(r'\d{4}_\d{8}_\d{2}_\d{2}_\d{2}.txt', base_file):
                    product_type = 'growth_file'
                    ext_type = 'json'
                    #self._insert_dataproducts_dict()

                # If we are here then we have a file that is associated with a particular observation
                filename = os.path.basename(f)

                # Try and find the observation_id associated with the match string
                observation_id = ""
                for obs in obs_list:
                    if match_str in obs[1]:
                        print("Found observation! Using id number:%s" % obs[0])
                        observation_id = obs[0]
                        break

                # Calibration Group
                if len(filename) == 27 and filename[:2] == 'b_':
                    print("Bias corrected file")
                    product_type = 'calibration'
                    ext_type = filename.split('.')[-1]
                elif len(filename) == 31 and filename[:5] == 'crr_b':
                    print("Corrected file")
                    product_type = 'calibration'
                    ext_type = filename.split('.')[-1]
                elif len(filename) == 35 and filename[:9] == 'maskcrr_b':
                    print("Mask Corrected file")
                    product_type = 'calibration'
                    ext_type = filename.split('.')[-1]
                elif len(filename) == 36 and filename[:4] == 'bkgd':
                    print("Background Corrected file")
                    product_type = 'calibration'
                    ext_type = filename.split('.')[-1]
                elif len(filename) > 30 and filename[:4] == 'guid':
                    print("Guider Corrected file")
                    product_type = 'calibration'
                    ext_type = filename.split('.')[-1]


                # Standard group
                elif 'STD' in filename:
                    if "e3" in filename:
                        print("E3D file")
                    elif filename[:7] in self.standard_file_dict:
                        typ = self.standard_file_dict[filename[:7]]
                        print("File type: %s" % typ)
                    else:
                        print("UNKOWN FILETYPE Standard", f)

                    base, ext = filename.split('STD')
                    print("Base:%s\nExt:%s" % (base, ext))

                else:
                    if "e3" in filename:
                        print("E3D file")
                    elif "snid" in filename:
                        print("SNID output")
                    elif filename[:9] in self.science_file_dict:
                        typ = self.science_file_dict[filename[:9]]
                        print("File type: %s" % typ)
                    elif filename == 'rc%s.fits' % match_str:
                        print("Raw filts files")
                    elif 'finder' in filename:
                        print("Finder file")
                    elif 'rc%s' % night in filename:
                        print("Astrometry file")
                    else:
                        print("UNKOWN FILETYPE Science", f)

                    #base, ext = filename.split('STD')
                    #print("Base:%s\nExt:%s" % (base, ext))
                    pass
            print(non_fits_id_list)

    def ingest_nightly_raw_files(self, night, directory='', prefix='rc', updateIngest=False,
                                 add_to_inges_list=True, check_ingest_list=True,
                                 check_by='night', ingest_file="", add_telescope_stats=True):
        """

        :return:
        """
        if not directory:
            directory = os.path.join(self.raw_dir, night)
        else:
            directory = directory

        # 1. Start by gathering all the requested files in the raw directory
        raw_files = self.gather_files(directory, prefix=prefix, suffix='fits')

        # 2. Check if the file has already been added to the ingest list
        if check_ingest_list:
            if not ingest_file:
                ingest_file = self.ingest_file

            ingest_list = self.__read_in_json_file(ingest_file)

            if ingest_list:
                nightly_data = ingest_list[night]['rc']['raw']

                for f in raw_files:
                    if f in nightly_data:
                        raw_files.pop(f)

        # 3. Now check the database for all matches
        needs_to_be_removed_list = []
        if check_by == 'night':
            query = "SELECT fitsfile,id FROM observation WHERE fitsfile LIKE '%%%s%%'" % night
            print(query)
            ret = self.db.execute_sql(query, return_type='')

            for idx, val in enumerate(raw_files):
                base_file = os.path.basename(val)
                if "a_" in base_file:
                    continue
                for idx2, sublist in enumerate(ret):
                    #print(sublist, base_file)
                    if base_file in sublist[0]:
                        print("Found it!", sublist, base_file)
                        needs_to_be_removed_list.append(idx)
                        if add_telescope_stats:
                            hdu = fits.getheader(val)
                            tstats_dict = dict(hdu)
                            self._insert_telescope_stats_dict(tstats_dict, sublist[1])

                        break

            for index in sorted(needs_to_be_removed_list, reverse=True):
                del raw_files[index]
        # 4. At this point we should only have files that have not been added to the database

        import time
        #time.sleep(1000)
        # 5. Add files to the database
        for f in sorted(raw_files):

            # Open the header file
            hdu = fits.getheader(f)
            print(f, "here")
            print(hdu.keys())
            if "a_" in f:
                continue
            c = SkyCoord(ra=hdu['ra'], dec=hdu['dec'], unit=(u.hourangle, u.deg))

            if all(key in hdu for key in self.required_obs_header):
                # Create dictionary of data
                header_dict = {'fitsfile': os.path.basename(f)}
                all_fields_accounted_for = True
                missing_value_list = []

                # This is a special case section.  These keys do not
                # match the keys in the header and need to be pulled
                # out one x one.  They are all required
                if hdu['obj_id']:
                    header_dict['object_id'] = hdu['obj_id']
                else:
                    print("Missing value in header")
                    missing_value_list.append(["object_id",  hdu['obj_id']])
                    all_fields_accounted_for = False

                if hdu['req_id']:
                    header_dict['request_id'] = hdu['req_id']
                else:
                    print("Missing value in header")
                    missing_value_list.append(["request_id", hdu['req_id']])
                    all_fields_accounted_for = False

                if hdu['mjd_obs']:
                    header_dict['mjd'] = hdu['mjd_obs']
                else:
                    print("Missing value in header")
                    missing_value_list.append(["mjd", hdu['mjd_obs']])
                    all_fields_accounted_for = False

                # These keywords are optional and not required for ingest into
                # the db
                for hkey, dkey in self.optional_obs_header.items():
                    if hkey in hdu:
                        header_dict[dkey] = hdu[hkey]
                    else:
                        print("Missing %s" % hkey)


                # The following keywords match the database keywords and header
                header_list = ['airmass', 'exptime', 'lst', 'ra', 'dec',
                               'tel_az', 'tel_el', 'tel_pa', 'ra_off',
                               'dec_off']

                for h in header_list:
                    if h == 'lst':
                        header_dict[h] = str(hdu[h])
                    elif h == 'ra':
                        header_dict[h] = c.ra.degree
                    elif h == 'dec':
                        header_dict[h] = c.dec.degree
                    else:
                        header_dict[h] = hdu[h]

            else:
                print("IN MISSING HEADER SECTION")
                all_fields_accounted_for = False
                for k in self.required_obs_header:
                    if k not in hdu:
                        print(k)
                        print("Missing header value: %s" % k)
                        missing_file = open('missing_header_list.txt', 'w')
                        missing_file.write('%s: %s\n' % (f,k))
                continue

            # Now we check for optional values

            # At this point we check to see if we have everything we need for
            # adding the target to the database
            ret = None

            if all_fields_accounted_for:
                print("Adding %s to the database" % f)
                #ret = self.db.add_observation(header_dict)
                print(ret)

            if add_telescope_stats:
                print("Adding header information for %s "
                      "into the telescope_stats database" % f)

                if not ret:
                    ret = self.db.execute_sql("SELECT id FROM observation WHERE fitsfile = '%s'" % os.path.basename(f))
                tstats_dict = dict(hdu)
                self._insert_telescope_stats_dict(tstats_dict, ret[0])



                # Start by looking if their is a ret value for the observation id

        # Now add the rest of the header information into the telescope stats
        # database
        """
        Adds an observation

        Args:
            header_dict (dict):
                required:
                    'object_id' (int/long),
                    'request_id' (int/long),
                    'mjd' (float),
                    'airmass' (float),
                    'exptime' (float),
                    'fitsfile' (abspath str),
                    'lst' (str),
                    'ra' (float),
                    'dec' (float),
                    'tel_az' (float),
                    'tel_el' (float),
                    'tel_pa' (float),
                    'ra_off' (float),
                    'dec_off' (float),
                optional:
                    'airmass_end' (float),
                    'parang' (float),
                    'parang_end' (float),
                    'ra_off' (float),
                    'dec_off' (float),
                    'imtype' (str),
                    'time_elapsed' (datetime.timedelta object or
                                    float/int seconds),
                    'filter' (str),
                        options - 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b',
                                    'NA'
                    'camera' (str)

        Returns:
            (-1, "ERROR...") if there was an issue

            (id (long), "Observation added") if it completed successfully
        """

if __name__ == '__main__':
    x = ingestNight(raw_dir='/scr/rsw/sedm/data/raw/')

    #x.ingest_nightly_raw_files('20181204', prefix='all')
    x.ingest_nightly_data_products('20181204', 'ifu')
    #hdu = fits.getheader('/scr/rsw/sedm/data/raw/20181204/rc20181204_05_02_24.fits')
    #tstats_dict = dict(hdu)
    #print(x._insert_telescope_stats_dict(tstats_dict, 20190717231726230))
    #print(x.gather_files('/scr/rsw/sedm/redux/20190629', suffix=''))




