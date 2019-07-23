# Script to add all fits files into the database
from astropy.io import fits
import glob
import os
import json
from db.SedmDb import SedmDB
import astropy.units as u
from astropy.coordinates import SkyCoord


class ingestNight:
    """
    Ingest a directory of files or given a specific night ingest all raw and
    processed images for that night.
    """
    def __init__(self, raw_dir='/scr2/sedm/raw/',
                 ifu_proc_dir='/scr2/sedmdrp/redux/',
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
                                         CTYPE1="str", CTYPE2="str", CRVAL1="float", CRVAL2="float",
                                         END_SHUT="datetime", ENDAIR="float", ENDDOME="str", END_RA="str",
                                         END_DEC="str", END_PA="str", ELAPTIME="float", observation_id="int")

        self.db = SedmDB(host="pharos", dbname='sedmdb', port=5432)

    def _format_telescope_stats_dict(self, headerdict):

        for k, i in headerdict.items():
            print(k, i)

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
        print(len(raw_files), raw_files)
        if check_by == 'night':
            query = "SELECT fitsfile,id FROM observation WHERE fitsfile LIKE '%%%s%%'" % night
            print(query)
            ret = self.db.execute_sql(query, return_type='')

            for idx, val in enumerate(raw_files):
                base_file = os.path.basename(val)
                for idx2, sublist in enumerate(ret):
                    #print(sublist, base_file)
                    if base_file in sublist[0]:
                        print("Found it!", sublist, base_file)
                        needs_to_be_removed_list.append(idx)
                        break

            print(needs_to_be_removed_list)
            for index in sorted(needs_to_be_removed_list, reverse=True):
                del raw_files[index]
        # 4. At this point we should only have files that have not been added to the database
        print(len(raw_files), raw_files)


        # 5. Add files to the database
        for f in raw_files:

            # Open the header file
            hdu = fits.getheader(f)
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


            # Now we check for optional values

            # At this point we check to see if we have everything we need for
            # adding the target to the database
            ret = None

            if all_fields_accounted_for:
                print("Adding %s to the database" % f)
                ret = self.db.add_observation(header_dict)
                print(ret)

            if add_telescope_stats:
                print("Adding header information for %s "
                      "into the telescope_stats database" % f)
                tstats_dict = dict(hdu)

                

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

    x.ingest_nightly_raw_files('20181204', prefix='all')

    #print(x.gather_files('/scr/rsw/sedm/redux/20190629', suffix=''))




