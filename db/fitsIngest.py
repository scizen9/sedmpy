# Script to add all fits files into the database
from astropy.io import fits
import glob
import os
import json
from db.SedmDb import SedmDB


class ingestNight:
    """
    Ingest a directory of files or given a specific night ingest all raw and
    processed images for that night.
    """
    def __init__(self, raw_dir='/scr2/sedm/raw/',
                 ifu_proc_dir='/scr2/sedmdrp/redux/',
                 rc_proc_dir='/scr2/sedm/phot/',
                 ifu_prexfix='ifu', rc_prefix='rc',
                 ingest_file='ingest.json'):
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
        self.db = SedmDB(host="pharos", dbname='sedmdb', port=5432)

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
                                 check_by='night', ingest_file=""):
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
        needs_to_added_list = []

        if check_by == 'night':
            query = "SELECT fitsfile,id FROM observation WHERE fitsfile LIKE '%%%s%%'" % night
            print(query)
            ret = self.db.execute_sql(query, return_type='')

            for i in raw_files:
                for sublist in ret:
                    print(sublist, i)
                    if i in sublist[0]:
                        print("Found it!", sublist)
                        raw_files.pop(i)
                        break

        print(raw_files)

if __name__ == '__main__':
    x = ingestNight(raw_dir='/scr/rsw/sedm/data/raw/')

    x.ingest_nightly_raw_files('20181204')

    #print(x.gather_files('/scr/rsw/sedm/redux/20190629', suffix=''))




