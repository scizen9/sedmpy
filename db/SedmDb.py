import sqlalchemy.pool as pool
from sqlalchemy import exc
import psycopg2
import numpy as np
import subprocess
from astropy.time import Time
from datetime import timedelta
from werkzeug.security import generate_password_hash
import os
import sys
import psycopg2.extras
import psycopg2.errors

import smtplib

from marshals.interface import update_status_request

from email.message import EmailMessage
# from email.utils import make_msgid

if sys.version_info > (3,):
    long = int


SITE_ROOT = os.path.realpath(os.path.dirname(__file__))
json_url = os.path.join(SITE_ROOT, 'config.json')
fritz_base_url = 'https://fritz.science/'
fritz_view_source_url = fritz_base_url + 'source/'
growth_base_url = 'http://skipper.caltech.edu:8080/cgi-bin/growth/'
growth_view_source_url = growth_base_url + 'view_source.cgi?name='


# Singleton/SingletonPattern.py

class SedmDB:
    class __SedmDB:
        def __init__(self, dbname, host, port=5432, supply_pass=False,
                     passwd=None):
            """
            Creates the instance of db connections.
            Needs the username as a parameter.
            The password for SedmDB must be stored in ~/.pgpass
            """
            cmd = "cat ~/.pgpass | grep sedmdbtest | grep -v '#' | " \
                  "awk -F ':' '{print $4}' | head -n 1"

            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            self.user_sedmdb = p.stdout.read().decode("utf-8").replace('\n', '')
            self.dbname = dbname
            self.host = host
            self.port = port
            self.supply_pass = supply_pass
            self.passwd = passwd
            self.pool_sedmdb = pool.QueuePool(self.__getSedmDBConn__,
                                              max_overflow=10, pool_size=2,
                                              recycle=True)

        def __str__(self):
            return repr(self)

        def __getSedmDBConn__(self):
            """
            Creates the connection to SedmDB.
            """
            if not self.supply_pass:
                sedmdbcon = psycopg2.connect(host=self.host, port=self.port,
                                             dbname=self.dbname,
                                             user=self.user_sedmdb)
            else:
                sedmdbcon = psycopg2.connect(host=self.host, port=self.port,
                                             dbname=self.dbname,
                                             user=self.user_sedmdb,
                                             password=self.passwd)

            return sedmdbcon

    instance = None

    def __init__(self, dbname='sedmdb', host='localhost', port=5432,
                 supply_pass=False, passwd=None):
        """
        Makes sure only one instance is created.
        """
        if not SedmDB.instance:
            SedmDB.instance = SedmDB.__SedmDB(dbname, host, port, supply_pass,
                                              passwd)

        self.sso_objects = None
        # Email templates

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def execute_sql(self, sql, return_type='list'):
        """
        Runs the SedmDB sql query in a safe way through the DBManager.

        Returns the object with the results.
        """
        conn = self.pool_sedmdb.connect()
        if return_type == 'list':
            cursor = conn.cursor()
        else:
            cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        try:
            cursor.execute(sql)
        except exc.DBAPIError as e:
            # an exception is raised, Connection is invalidated.
            if e.connection_invalidated:
                print("Connection was invalidated!")

        if 'SELECT' in sql[:8]:
            obj = cursor.fetchall()
            # print(obj)
            return obj
        else:
            cursor.execute("commit;")
            return []

    def get_conn_sedmDB(self):
        """
        Runs the WSDB sql query in a safe way through the DBManager.

        Returns the object with the results.
        """
        conn = self.pool_sedmdb.connect()

        return conn

    def add_user(self, pardic):
        """
        adds a new user

        Args:
            pardic (dict):
                required:
                    'username' (str),
                    'name' (str),
                    'email' (str),
                    'password' (str) (will be hashed+salted)

        Returns:
            (-1, "ERROR...") if there is an issue

            (id (long), "User added") if the user was added
        """
        # no need to check parameter value types as they are all strings
        uid = _id_from_time()
        pardic['id'] = uid
        keys = list(pardic.keys())
        if 'username' not in keys:
            return -1, "ERROR: no username provided!"
        # check for duplicate username
        usernames = [user[0] for user in
                     self.execute_sql('SELECT username FROM users')]
        if pardic['username'] in usernames:
            return -1, "ERROR: user with that username exists!"
        if 'password' in keys:
            pardic['password'] = generate_password_hash(pardic['password'])
        for key in reversed(keys):  # remove group keys and any other bad keys
            if key not in ['id', 'username', 'name', 'email', 'password']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        sql = _generate_insert_sql(pardic, keys, 'users')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_user sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_user sql command failed with " \
                       "a ProgrammingError!"
        return uid, "User added"

    def update_user(self, pardic):
        """
        updates a user

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'name' (str),
                    'email' (str),
                    'password' (str) (will be hashed+salted)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "User updated, columns 'column_names'")
             if the user is updated successfully
        """
        param_types = {'id': int, 'name': str, 'email': str, 'password': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"

        elif pardic['id'] not in [ux[0] for ux in
                                  self.execute_sql('SELECT id FROM users;')]:
            return -1, "ERROR: no user with the id!"
        keys.remove('id')
        if 'password' in keys:
            pardic['password'] = generate_password_hash(pardic['password'])

        # remove any keys that are invalid or not allowed to be updated
        for key in reversed(keys):
            if key not in ['name', 'email', 'password']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_update_sql(pardic, keys, 'users')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_user sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_user sql command failed with " \
                       "a ProgrammingError!"
        return pardic['id'], "User updated, columns " + str(keys)[1:-1]

    def remove_user(self, pardic):
        """
        Removes an existing user

        Args:
            pardic (dict):
                required:
                    'username' (str)
                    OR
                    'id': (str)

        Returns:
            (-1, "ERROR...") if there was an issue
            (user doesn't exist, not enough information in pardic)

            (0, "User removed") if the removal was successful
        """
        if 'username' in pardic.keys():
            user_id = self.get_from_users(['id'],
                                          {'username': pardic['username']})
            if user_id:
                if user_id[0] == -1:  # if get_from_users failed
                    return user_id
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'"
                                 % (user_id[0][0],))
                self.execute_sql("DELETE FROM users WHERE id='%s';"
                                 % (user_id[0][0],))
                return 0, "User removed"
            else:
                return -1, "ERROR: no user with that username!"
        elif 'id' in pardic.keys():
            if not (isinstance(pardic['id'], int) or isinstance(pardic['id'],
                                                                long)):
                return -1, "ERROR: id must be of type 'int'"
            if pardic['id'] in [ux[0] for ux in
                                self.execute_sql('SELECT id FROM users;')]:
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'"
                                 % (pardic['id'],))
                self.execute_sql("DELETE FROM users WHERE id='%s';"
                                 % (pardic['id'],))
                return 0, "User removed"
            else:
                return -1, "ERROR: no user with that id!"
        else:
            return -1, "ERROR: username or id required!"

    def get_from_users(self, values, where_dict=None, compare_dict=None):
        """
        select values from `users`

        Args:
            values: list of str
                values to be returned
            where_dict (dict):
                ``'param':'value'`` to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'username' (str),
                'name' (str),
                'email' (str),
                'password' (str)

        Returns:
            list of tuples with the values for each user matching the criteria
            
            empty list if no users match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'username': str, 'name': str, 'email': str,
                          'password': str}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'users')
        # print(sql)
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql
        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in designator

        Args:
            pardic (dict):
                required:
                    'designator' (str)

        Returns:
            (-1, "ERROR...") if no designator was provided or
             there is already a group with it

            (id (long), "Group added") if the adding was successful
        """
        gid = _id_from_time()
        pardic['id'] = gid
        if 'designator' not in pardic.keys():
            return -1, 'ERROR: no group designator provided!'
        groups = [des[0] for des in
                  self.execute_sql('SELECT designator FROM groups;')]
        if pardic['designator'] not in groups:
            sql = ("INSERT INTO groups (id, designator) VALUES ('%s', '%s')"
                   % (pardic['id'], pardic['designator']))
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_group sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_group sql command failed with " \
                           "a ProgrammingError!"
            return gid, "Group added"
        else:
            return -1, "ERROR: group exists!"

    def add_usergroup(self, user, group):
        """
        Adds the user as member of the group. Checks for duplicates in name.

        Args:
            user (int/long):
                id of the user in the 'users' Table
            group (int/long):
                id of the group in the 'groups' Table

        Returns:
            (-1, "ERROR...") if there was areason for failure

            (id of user, "User added to group") if the adding was successful
        """
        if user not in [user_id[0] for user_id in
                        self.execute_sql('SELECT id FROM users')]:
            return -1, "ERROR: user does not exist!"
        if group not in [group_id[0] for group_id in
                         self.execute_sql('SELECT id FROM groups')]:
            return -1, "ERROR: group does not exist!"
        usergroups = self.execute_sql(
            'SELECT user_id, group_id FROM usergroups')
        if (user, group) in usergroups:
            return -1, "ERROR: user already in group!"
        else:
            sql = "INSERT INTO usergroups (user_id, group_id) VALUES " \
                  "('%s','%s')" % (user, group)
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_to_group sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_to_group sql command failed with " \
                           "a ProgrammingError!"
            return user, "User added to group"

    def remove_from_group(self, user, group):
        """
        removes the user from the group.

        Args:
            user (int/long):
                id of the user in the 'users' Table
            group: int
                id of the group in the 'groups' Table

        Returns:
            (-1, "ERROR...") if there was areason for failure

            (0, "User removed from group") if the removal was successful
        """
        if user not in [user_id[0] for user_id in
                        self.execute_sql('SELECT id FROM users')]:
            return -1, "ERROR: user does not exist!"
        if group not in [group_id[0] for group_id in
                         self.execute_sql('SELECT id FROM groups')]:
            return -1, "ERROR: group does not exist!"
        usergroups = self.execute_sql(
            'SELECT user_id, group_id FROM usergroups')
        if (user, group) not in usergroups:
            return -1, "ERROR: user not in group!"
        else:
            sql = "DELETE FROM usergroups WHERE " \
                  "user_id='%s' AND group_id='%s'" % (user, group)
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return -1, "ERROR: remove_from_group sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: remove_from_group sql command failed with " \
                           "a ProgrammingError!"
            return 0, "User removed from group"

    def get_from_usergroups(self, values, where_dict=None, compare_dict=None):
        """
        select values from `usergroups`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options: ['user_id', 'group_id']

        Returns:
            list of tuples containing the values for each usergroup
            matching the criteria

            empty list if no usergroups match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'user_id': int, 'group_id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'usergroups')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: get_from_usergroups sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: get_from_usergroups sql command failed with " \
                       "a ProgrammingError!"
        return results

    def add_program(self, pardic):
        """
        creates a new program

        Args:
            pardic (dict):
                required:
                    'designator' (str) (unique)
                    'name' (str)
                    'group_id' (int/long)
                    'pi' (str)
                optional:
                    'time_allocated' (datetime.timedelta object or
                                      float/int seconds)
                    'priority' (float)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)

        Returns:
            (-1, "ERROR...") if it failed to add

            (id (long), "Program added") if the program is added successfully
        """
        param_types = {'id': int, 'name': str, 'designator': str,
                       'group_id': int, 'pi': str, 'color': str,
                       'time_allocated': 'timedelta', 'priority': float,
                       'inidate': 'datetime', 'enddate': 'datetime'}
        pid = _id_from_time()
        pardic['id'] = pid
        keys = list(pardic.keys())

        if 'designator' in keys:
            if pardic['designator'] in [obj[0] for obj in
                                        self.execute_sql(
                                            'SELECT designator FROM program')]:
                return -1, "ERROR: a program with that designator " \
                           "already exists!"
        # check for required keys
        for key in ['designator', 'name', 'group_id', 'pi']:
            if key not in keys:
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(keys):  # remove any extraneous keys
            if key not in ['id', 'name', 'designator', 'group_id', 'pi',
                           'color', 'time_allocated', 'priority', 'inidate',
                           'enddate']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        program_sql = _generate_insert_sql(pardic, keys, 'program')
        try:
            self.execute_sql(program_sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_program sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_program sql command failed with " \
                       "a ProgrammingError!"
        return pid, "Program added"

    def update_program(self, pardic):
        """
        updates a pragram

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'name' (str)
                    'PI' (str)
                    'time_allocated' (datetime.timedelta object or
                                      float/int seconds)
                    'priority' (float)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Program updated, columns 'column_names'")
            if the program is updated successfully
        """
        param_types = {'id': int, 'time_allocated': 'timedelta', 'name': str,
                       'PI': str, 'priority': float, 'inidate': 'datetime',
                       'enddate': 'datetime', 'color': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"

        elif pardic['id'] not in [px[0] for px in
                                  self.execute_sql('SELECT id FROM program;')]:
            return -1, "ERROR: no program with the id!"
        keys.remove('id')
        # remove any keys that are invalid or not allowed to be updated
        for key in reversed(keys):
            if key not in ['time_allocated', 'name', 'PI', 'priority',
                           'inidate', 'enddate', 'color']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_update_sql(pardic, keys, 'program')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_program sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_program sql command failed with " \
                       "a ProgrammingError!"
        return pardic['id'], "Program updated, columns " + str(keys)[1:-1]

    def update_object(self, pardic):
        """
        updates a pragram

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'marshal_id' (int/long)
                    'name' (str)
                    'iauname' (str)
                    'ra' (float)
                    'dec' (float)
                    'epoch' (str)
                    'magnitude' (float)
                    'creationdate' ('year-month-day hour:minute:second')

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Object updated, columns 'column_names'")
            if the program is updated successfully
        """
        param_types = {'id': int, 'marshal_id': int, 'name': str,
                       'iauname': str, 'ra': float, 'dec': float, 'epoch': str,
                       'magnitude': float, 'creationdate': 'datetime'}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"

        elif pardic['id'] not in [ox[0] for ox in
                                  self.execute_sql('SELECT id FROM object;')]:
            return -1, "ERROR: no object with the id!"
        keys.remove('id')
        # remove any keys that are invalid or not allowed to be updated
        for key in reversed(keys):
            if key not in param_types.keys():
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_update_sql(pardic, keys, 'object')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_program sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_program sql command failed with " \
                       "a ProgrammingError!"
        return pardic['id'], "Program updated, columns " + str(keys)[1:-1]

    def get_from_program(self, values, where_dict=None, compare_dict=None):
        """
        select values from `program`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options: ['id', 'designator', 'name', 'group_id', 'PI',
                                  'time_allocated', 'priority', 'inidate',
                                  'enddate', 'color']

        Returns:
            list of tuples containing the values for each program
            matching the criteria

            empty list if no programs match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'name': str, 'designator': str,
                          'group_id': int, 'PI': str, 'color': str,
                          'time_allocated': 'timedelta', 'priority': float,
                          'inidate': 'datetime', 'enddate': 'datetime'}

        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'program')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: get_from_program sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: get_from_program sql command failed with " \
                       "a ProgrammingError!"
        return results

    def add_allocation(self, pardic):
        """
        creates a new allocation entry

        Args:
            pardic (dict):
                required:
                    'program_id' (int) (id of associated program)
                    'designator' (str) (e.g. 2016A-I0001)
                optional:
                    'time_allocated' (datetime.timedelta object or
                                      float/int seconds)
                    'time_spent' (datetime.timedelta object or
                                  float/int seconds)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)

        Returns:
            (-1, "ERROR...") if it failed to add

            (id (long), "Allocation added") if the program is added successfully
        """
        param_types = {'id': int, 'program_id': int, 'designator': str,
                       'time_allocated': 'timedelta', 'time_spent': 'timedelta',
                       'inidate': 'datetime', 'enddate': 'datetime',
                       'color': str, 'active': bool}
        aid = _id_from_time()
        pardic['id'] = aid
        keys = list(pardic.keys())

        if 'program_id' in keys:
            if pardic['program_id'] not in [obj[0] for obj in
                                            self.execute_sql(
                                                'SELECT id FROM program;')]:
                return -1, "ERROR: no program with that id exists!"

        for key in ['program_id', 'designator']:
            if key not in keys:  # check for required key
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(keys):  # remove any extraneous keys
            if key not in ['id', 'program_id', 'designator', 'time_allocated',
                           'time_spent', 'inidate', 'enddate',
                           'color', 'active']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        program_sql = _generate_insert_sql(pardic, keys, 'allocation')
        try:
            self.execute_sql(program_sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_program sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_program sql command failed with " \
                       "a ProgrammingError!"
        return aid, "Allocation added"

    def update_allocation(self, pardic):
        """
        updates an allocation entry

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'time_allocated' (datetime.timedelta object or
                                      float/int seconds)
                    'time_spent' (datetime.timedelta object or
                                  float/int seconds)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)
                    'active' (boolean)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Allocation updated, columns 'column_names'")
             if the allocation is updated successfully
        """
        param_types = {'id': int, 'time_allocated': 'timedelta',
                       'time_spent': 'timedelta', 'inidate': 'datetime',
                       'enddate': 'datetime', 'color': str, 'active': bool}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"

        elif pardic['id'] not in [ax[0] for ax in
                                  self.execute_sql(
                                      'SELECT id FROM allocation;')]:
            return -1, "ERROR: no allocation with the id!"
        keys.remove('id')
        # remove any keys that are invalid or not allowed to be updated
        for key in reversed(keys):
            if key not in ['time_allocated', 'time_spent', 'inidate',
                           'enddate', 'color', 'active']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        print(pardic, keys)

        sql = _generate_update_sql(pardic, keys, 'allocation')

        print(sql)

        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_allocation sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_allocation sql command failed with " \
                       "a ProgrammingError!"
        return pardic['id'], "Allocation updated, columns " + str(keys)[1:-1]

    def get_from_allocation(self, values, where_dict=None, compare_dict=None):
        """
        select values from `allocation`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options: ['id', 'designator', 'program_id, 'time_spent',
                                  'time_allocated', 'inidate', 'enddate',
                                  'color', 'active']

        Returns:
            list of tuples containing the values for each program
             matching the criteria

            empty list if no programs match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'program_id': int, 'designator': str,
                          'time_spent': 'timedelta', 'color': str,
                          'time_allocated': 'timedelta', 'inidate': 'datetime',
                          'enddate': 'datetime', 'active': bool}

        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'allocation')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: get_from_allocation sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: get_from_allocation sql command failed with " \
                       "a ProgrammingError!"
        return results

    def add_object(self, pardic):
        """
        Creates a new object

        Args:
            pardic (dict):
                required:
                    'name' (str),
                    'typedesig' (str),
                required for a fixed object:
                    'ra' (float) ra in degrees,
                    'dec' (float) dec in degrees
                optional:
                    'magnitude' (float) (preferably 'r' filter),
                    'epoch' (str),
                    'iauname' (str),
                    'marshal_id' (int/long)

                'typedesig' should be one of:
                    'f' (fixed),
                    'v' (periodic fixed),
                    'P' (built-in planet or satellite name),
                    'e' (heliocentric elliptical),
                    'h' (heliocentric hyperbolic),
                    'p' (heliocentric parabolic),
                    'E' (geocentric elliptical)

        Returns:
            (-1, "ERROR...") if it failed to add

            (id (long), "Object added") if the object is added successfully
        """
        param_types = {'id': int, 'name': str, 'typedesig': str, 'ra': float,
                       'dec': float, 'epoch': str, 'iauname': str,
                       'marshal_id': int, 'magnitude': float}
        oid = _id_from_time()
        pardic['id'] = oid
        obj_keys = list(pardic.keys())
        if 'marshal_id' in obj_keys:
            if pardic['marshal_id'] in [obj[0] for obj in
                                        self.execute_sql(
                                            'SELECT marshal_id FROM object')]:
                return -1, "ERROR: object exists!"
        if 'epoch' not in obj_keys:
            obj_keys.append('epoch')
            pardic['epoch'] = 'J2000'

        # check if 'name' and 'typedesig' are provided
        for key in ['name', 'typedesig']:
            if key not in obj_keys:
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(obj_keys):  # remove any extraneous keys
            if key not in ['id', 'name', 'typedesig', 'ra', 'dec', 'epoch',
                           'marshal_id', 'iauname', 'magnitude']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(obj_keys, pardic, param_types)
        if type_check:
            return -1, type_check

        if pardic['typedesig'] == 'f':
            for key in ['ra', 'dec']:
                if key not in obj_keys:
                    return -1, "ERROR: %s not provided!" % (key,)
            dup = self.execute_sql("SELECT id, name FROM object "
                                   "WHERE q3c_radial_query(ra, dec,"
                                   " '%s', '%s', .000278)"
                                   % (pardic['ra'], pardic['dec']))
            if dup:  # if there is already an object within an arcsecond
                # check for same coords, same name
                if pardic['name'] in [dx[1] for dx in dup]:
                    lis = [dx[1] for dx in dup]
                    idx = lis.index(pardic['name'])
                    return (-1, "ERROR: The object '%s' is already "
                                "in the database with id %s"
                            % (pardic['name'], dup[idx][0]))
                print("there is already an object within 1 arcsec "
                      "of given coordinates with id: %s, name: %s"
                      % (dup[0][0], dup[0][1]))

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_object sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_object sql command failed with " \
                           "a ProgrammingError!"
            return oid, "Fixed object added"
        elif pardic['typedesig'] == 'v':
            for key in ['ra', 'dec']:
                if key not in obj_keys:
                    return -1, "ERROR: %s not provided!" % (key,)
            dup = self.execute_sql("SELECT id, name FROM object "
                                   "WHERE q3c_radial_query(ra, dec, "
                                   "'%s', '%s', .000278)"
                                   % (pardic['ra'], pardic['dec']))
            if dup:  # if there is already an object within an arcsecond
                # check for same coords, same name
                if pardic['name'] in [dx[1] for dx in dup]:
                    lis = [dx[1] for dx in dup]
                    idx = lis.index(pardic['name'])
                    return (-1, "ERROR: The object '%s' is already in the "
                                "database with id %s"
                            % (pardic['name'], dup[idx][0]))
                print("there is already an object within 1 arcsec "
                      "of given coordinates with id: %s, name: %s"
                      % (dup[0][0], dup[0][1]))

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_object sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_object sql command failed with " \
                           "a ProgrammingError!"
            return oid, "Fixed object added, input period data into table " \
                        "`periodic`"

        elif pardic['typedesig'] in ['h', 'E', 'e', 'p']:
            function_dict = {'e': 'add_elliptical_heliocentric',
                             'h': 'add_hyperbolic_heliocentric',
                             'p': 'add_parabolic_heliocentric',
                             'E': 'add_earth_satellite'}

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_object sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_object sql command failed with " \
                           "a ProgrammingError!"
            return oid, "Non-fixed object added, orbit parameters " \
                        "can be added with `%s`" % function_dict[pardic['typedesig']]
        elif pardic['typedesig'] == 'P':
            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_object sql command failed with " \
                           "an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_object sql command failed with " \
                           "a ProgrammingError!"
            return oid, "Non-fixed object added"
        else:
            return -1, "ERROR: typedesig provided was invalid, " \
                       "it must be one of ['f', 'h', 'E', 'e', 'p', 'P']!"

    # TODO: write update_object

    def get_from_object(self, values, where_dict=None, compare_dict=None):
        """
        select values from `objects`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options: ['id', 'marshal_id', 'name', 'iauname', 'ra',
                                  'dec', 'typedesig', 'epoch']

        Returns:
            list of tuples containing the values
            for each user matching the criteria

            empty list if no objects match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'marshal_id': int, 'name': str,
                          'iauname': str, 'ra': float, 'dec': float,
                          'typedesig': str, 'epoch': float, 'magnitude': float}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'object')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def get_objects_near(self, ra, dec, radius, values=None):
        """
        get object entries with coordinates within a radius

        Args:
            ra (float or int): ra in degrees of the object
            dec (float or int): dec in degrees of the object
            radius (float or int): radius in arcseconds
            values (list): which parameters to return
                    options: 'id', 'marshal_id', 'name', 'iauname', 'ra', 'dec',
                     'typedesig', 'epoch'
                    defaults to ['id', 'name', 'epoch']

        Returns:
            list of tuples [(id, name, epoch)] for each object in the radius

            [] if there are no objects found

            (-1, "ERROR...") if there is an issue
        """
        if values is None:
            values = ['id', 'name', 'epoch']
        for value in values[::-1]:
            if value not in ['id', 'marshal_id', 'name', 'iauname', 'ra', 'dec',
                             'typedesig', 'epoch']:
                return -1, "ERROR: %s is an invalid column name!" % (value,)
        if not (isinstance(ra, float) or isinstance(ra, int)):
            return -1, "ERROR: parameter ra must be of type 'float' " \
                       "or type 'int'!"
        if not (isinstance(dec, float) or isinstance(dec, int)):
            return -1, "ERROR: parameter dec must be of type 'float' " \
                       "or type 'int'!"
        if not (isinstance(radius, float) or isinstance(radius, int)):
            return -1, "ERROR: parameter radius must be of type 'float' " \
                       "or type 'int'!"

        objects = self.execute_sql("SELECT id, name, epoch FROM object WHERE "
                                   "q3c_radial_query(ra, dec, '%s', '%s', '%s')"
                                   % (ra, dec, .000278 * radius))
        return objects

    def get_object_id_from_name(self, object_name):
        """
        finds the id of an object given its name or part of its name

        Args:
            object_name (str): part of the name of an object

        Returns:
            list of tuples [(id, full name)] for each object with
            a matching name

            [] if no matching object name is found
        """
        object_name = object_name.lower()
        sql = "SELECT id, name FROM object WHERE LOWER(name) LIKE '%s%s%s'" \
              % ('%', object_name, '%')
        obj = self.execute_sql(sql)
        return obj

    def add_elliptical_heliocentric(self, orbit_params):
        """
        Adds the orbit parameters for an elliptical heliocentric orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id' (int/long),
                    'inclination' (float),
                    'longascnode_O' (float) (lon. of ascending node),
                    'perihelion_o' (float) (arg. of perihelion),
                    'a' (float) (mean distance AU),
                    'n' (float) (mean daily motion deg/day),
                    'e' (float) (eccentricity),
                    'M' (float) (mean anomaly),
                    'mjdepoch' (float) (epoch, time of 'M'),
                    'D': (float) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components
                                  of magnitude model)
                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'inclination': float,
                       'longascnode_O': float, 'perihelion_o': float,
                       'a': float, 'n': float, 'e': float, 'M': float,
                       'mjdepoch': float, 'D': float, 'M1': float,
                       'M2': float, 's': float}
        eid = _id_from_time()
        orbit_params['id'] = eid
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['inclination', 'longascnode_O', 'perihelion_o', 'a', 'n',
                    'e', 'M', 'mjdepoch', 'D', 'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(orb_keys):
            if key not in ['id', 'inclination', 'longascnode_O', 'perihelion_o',
                           'a', 'n', 'e', 'M', 'mjdepoch', 'D', 'M1', 'M2', 's',
                           'object_id']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return -1, type_check

        orb_sql = _generate_insert_sql(orbit_params, orb_keys,
                                       'elliptical_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_elliptical_orbit sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_elliptical_orbit sql command failed with " \
                       "a ProgrammingError!"
        return eid, "Elliptical heliocentric orbit added"

    def get_from_elliptical_heliocentric(self, values, where_dict=None,
                                         compare_dict=None):
        """
        select values from `elliptical_heliocentric`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'object_id' (int/long),
                'inclination' (float),
                'longascnode_O' (float) (lon. of ascending node),
                'perihelion_o' (float) (arg. of perihelion),
                'a' (float) (mean distance AU),
                'n' (float) (mean daily motion deg/day),
                'e' (float) (eccentricity),
                'M' (float) (mean anomaly),
                'mjdepoch' (float) (epoch, time of 'M'),
                'D': (float) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching
            the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'object_id': int, 'inclination': float,
                          'longascnode_O': float, 'perihelion_o': float,
                          'a': float, 'n': float, 'e': float, 'M': float,
                          'mjdepoch': float, 'D': float, 'M1': float,
                          'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'elliptical_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_hyperbolic_heliocentric(self, orbit_params):
        """
        Adds the orbit parameters for an hyperbolic heliocentric orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id' (int/long),
                    'T' ('year-month-day' or `~astropy.time.Time` object),
                    'inclination' (float),
                    'longascnode_O' (float) (lon. of ascending node),
                    'perihelion_o' (float) (arg. of perihelion),
                    'e' (float) (eccentricity),
                    'q' (float) (perihelion distance AU),
                    'D': (float) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components
                                  of magnitude model)

                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'T': 'date', 'e': float,
                       'inclination': float, 'longascnode_O': float,
                       'perihelion_o': float, 'q': float, 'D': float,
                       'M1': float, 'M2': float, 's': float}
        hid = _id_from_time()
        orbit_params['id'] = hid
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'e',
                    'q', 'D', 'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(orb_keys):
            if key not in ['id', 'T', 'inclination', 'longascnode_O',
                           'perihelion_o', 'e', 'q', 'D', 'M1', 'M2', 's',
                           'object_id']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return -1, type_check

        orb_sql = _generate_insert_sql(orbit_params, orb_keys,
                                       'hyperbolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_hyperbolic_orbit sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_hyperbolic_orbit sql command failed with " \
                       "a ProgrammingError!"
        return hid, "Hyperbolic heliocentric orbit added"

    def get_from_hyperbolic_heliocentric(self, values, where_dict=None,
                                         compare_dict=None):
        """
        select values from `hyperbolic_heliocentric`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'object_id' (int/long),
                'T' ('year-month-day' or `~astropy.time.Time` object),
                'inclination' (float),
                'longascnode_O' (float) (lon. of ascending node),
                'perihelion_o' (float) (arg. of perihelion),
                'e' (float) (eccentricity),
                'q' (float) (perihelion distance AU),
                'D': (float) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit
            matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'object_id': int, 'T': 'date', 'e': float,
                          'inclination': float, 'longascnode_O': float,
                          'perihelion_o': float, 'q': float, 'D': float,
                          'M1': float, 'M2': float, 's': float,
                          'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'hyperbolic_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_parabolic_heliocentric(self, orbit_params):
        """
        Adds the orbit parameters for an parabolic heliocentric orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id'(int),
                    'T' (date str),
                    'inclination' (float),
                    'perihelion_o' (float) (arg. of perihelion),
                    'q' (float) (perihelion distance),
                    'longascnode_O' (float) (lon. of ascending node),
                    'D': (float) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components
                                  of magnitude model)

                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'T': 'date',
                       'inclination': float, 'longascnode_O': float,
                       'perihelion_o': float, 'q': float, 'D': float,
                       'M1': float, 'M2': float, 's': float}
        pid = _id_from_time()
        orbit_params['id'] = pid
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'q',
                    'D', 'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(orb_keys):
            if key not in ['id', 'T', 'inclination', 'longascnode_O',
                           'perihelion_o', 'q', 'D', 'M1', 'M2', 's',
                           'object_id']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return -1, type_check

        orb_sql = _generate_insert_sql(orbit_params, orb_keys,
                                       'parabolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_parabolic_orbit sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_parabolic_orbit sql command failed with " \
                       "a ProgrammingError!"
        return pid, "Parabolic heliocentric orbit added"

    def get_from_parabolic_heliocentric(self, values, where_dict=None,
                                        compare_dict=None):
        """
        select values from `parabolic_heliocentric`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'object_id'(int),
                'T' (date str),
                'inclination' (float),
                'perihelion_o' (float) (arg. of perihelion),
                'q' (float) (perihelion distance),
                'longascnode_O' (float) (lon. of ascending node),
                'D': (float) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit
            matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'object_id': int, 'T': 'date', 'inclination': float,
                          'longascnode_O': float, 'perihelion_o': float,
                          'q': float, 'D': float, 'M1': float, 'M2': float,
                          's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'parabolic_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_earth_satellite(self, orbit_params):
        """
        Adds the orbit parameters for an Earth satellite orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id' (int/long),
                    'T' ('year-month-day') (epoch of other fields),
                    'inclination' (float),
                    'ra' (float) (ra of ascending node),
                    'e' (float) (eccentricity),
                    'pedigree' (float) (arg. of pedigree),
                    'M' (float) (mean anomaly),
                    'n' (float) (mean motion, revs/day),
                    'decay' (float) (orbit decay rate, rev/day^2),
                    'reforbit' (int) (integral reference orbit number at epoch),

                optional:
                    'drag' (float) (drag coefficient, 1/(Earth radii))

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'T': 'date', 'e': float,
                       'inclination': float, 'ra': float, 'pedigree': float,
                       'M': float, 'n': float, 'decay': float, 'reforbit': int,
                       'drag': float}
        eid = _id_from_time()
        orbit_params['id'] = eid
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                    'decay', 'reforbit', 'object_id']:
            if key not in orb_keys:
                return -1, "ERROR: %s not provided!" % (key,)
        for key in reversed(orb_keys):
            if key not in ['id', 'T', 'inclination', 'ra', 'e', 'pedigree',
                           'M', 'n', 'decay', 'reforbit', 'drag', 'object_id']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return -1, type_check

        orb_sql = _generate_insert_sql(orbit_params, orb_keys,
                                       'earth_satellite')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_earth_satellite_orbit sql command failed " \
                       "with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_earth_satellite_orbit sql command failed " \
                       "with a ProgrammingError!"
        return eid, "Earth satellite orbit added"

    def get_from_earth_satellite(self, values, where_dict=None,
                                 compare_dict=None):
        """
        select values from `earth_satellite`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'object_id' (int/long),
                'T' ('year-month-day') (epoch of other fields),
                'inclination' (float),
                'ra' (float) (ra of ascending node),
                'e' (float) (eccentricity),
                'pedigree' (float) (arg. of pedigree),
                'M' (float) (mean anomaly),
                'n' (float) (mean motion, revs/day),
                'decay' (float) (orbit decay rate, rev/day^2),
                'reforbit' (int) (integral reference orbit number at epoch),
                'drag' (float) (drag coefficient, 1/(Earth radii))

        Returns:
            list of tuples containing the values for each orbit
            matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'object_id': int, 'T': 'date', 'e': float,
                          'inclination': float, 'ra': float, 'pedigree': float,
                          'M': float, 'n': float, 'decay': float,
                          'reforbit': int, 'drag': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'earth_satellite')
        print(sql)
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_periodic(self, pardic):
        """
        Creates a new entry in the periodic table with the parameters
        specified in the dictionary.

        Args:
            pardic:
                required:
                    'object_id' (int/long)
                    'phasedays' (float) (period in days)
                    one of:
                        'mjd0' (float) (mjd of period 0)
                        'phi' (float) (phase at jd=2,400,000 as in Sesar+2017)

        Returns:
            (-1: "ERROR...") if there is an issue

            (id (long), "Periodic information added")
             if it completes successfully
        """
        param_types = {'id': int, 'object_id': int, 'mjd0': float,
                       'phi': float, 'phasedays': float}
        pid = _id_from_time()
        pardic['id'] = pid
        # TODO: which parameters are required? test
        keys = list(pardic.keys())

        # check we have all STRICTLY required keys
        for key in ['object_id', 'phasedays']:
            if key not in keys:
                return -1, "ERROR: %s not provided!" % (key,)

        # check there's at least one of ('mjd0', 'phi')
        if 'mjd0' in keys and 'phi' in keys:
            # if there are both, check they don't conflict

            # phi at jd=2,400,000 (ie mjd=-0.5) calculated using mjd0
            phi_mjd0 = ((-0.5 - pardic['mjd0']) / pardic['phasedays']) % 1
            # the 0.5 parts are to account for eg phi1=0.999, phi2=0.0001
            phase_diff = (phi_mjd0 - pardic['phi'] + 0.5) % 1 - 0.5

            if abs(phase_diff) > 0.03:
                return -1, "ERROR: mjd0 {mjd0} and phi {phi} are " \
                           "overdetermined and inconsistent".format(**pardic)
            else:
                print("mjd0 and phi differ by just %f periods" % phase_diff)
                del pardic['phi']  # we really shouldn't overdetermine the db

        elif 'mjd0' not in keys and 'phi' not in keys:
            # we can't be missing both
            return -1, "ERROR: mjd0 and phi both not provided!"
            
        # check we have no extraneous keys
        for key in reversed(keys):  # why are we reversing this?
            if key not in ['id', 'object_id', 'mjd0', 'phi', 'phasedays']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'periodic')
        print(sql, type_check)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_periodic sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_periodic sql command failed with " \
                       "a ProgrammingError!"
        return pid, "Periodic infomration added"

    def get_from_periodic(self, values, where_dict=None, compare_dict=None):
        """
        select values from `periodic`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'object_id' (int/long)
                'mjd0' (float)
                'phasedays' (float)

        Returns:
            list of tuples containing the values for each orbit
            matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'object_id': int, 'mjd0': float, 'phasedays': float,
                          'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'periodic')
        print(sql)
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def _add_planet_satellite_orbit(self, orbit_params):
        raise NotImplementedError
        # TODO: actually implement? maybe let higher level do this
        # TODO: query ??? for if it already exists

    def add_request(self, pardic):
        """
        Add a request

        Args:
            pardic (dict):
                required:
                    'object_id' (int/long),
                    'user_id' (int/long),
                    'allocation_id' (int/long),
                    'exptime' (str '{spec_duration, phot_duration}'),
                    'priority' (float),
                    'inidate' ('year-month-day') (start of observing window),
                    'enddate' ('year-month-day') (end of observing window),
                    'nexposures' or 'obs_seq' (below),

                optional:
                    'seq_repeats' (int) (defaults to 1),
                    'marshal_id' (int/long),
                    'maxairmass' (float) (max allowable airmass
                                          for observation, default 2.5),
                    'cadence' (float) (time between periods),
                    'phasesamples' (float) (how many samples in a period),
                    'sampletolerance' (float) (how much tolerance in when
                                               the samples can be taken),
                    'nexposures' (str '{# of ifu, # of u, # of g, # of r,
                                        # of i}'),
                    'obs_seq' (str e.g. '{3g, 3r, 1i, 1ifu, 2i}' for 3 of g,
                               then 3 of r, then 1 i, 1 ifu, 1 i),
                    'max_fwhm' (float),
                    'min_moon_dist' (float) (default 45 deg?),
                    'max_moon_illum' (float) (fraction),
                    'max_cloud_cover' (float) (fraction)

                Note:
                    the numbers in 'obs_seq' must be single digit,
                    spec/phot_duration should be duration per exp

        Returns:
            (-1, "ERROR...") if there is an issue with the input

            (id (long), "Request added") if there are no errors
        """
        # TODO: get a better description of cadence/phasesamples/sampletolerance
        param_types = {'id': int, 'object_id': int, 'user_id': int,
                       'allocation_id': int, 'exptime': str,
                       'priority': float, 'external_id': int, 'shareid': int,
                       'inidate': 'date', 'enddate': 'date', 'marshal_id': int,
                       'maxairmass': float, 'cadence': float,
                       'phasesamples': float, 'sampletolerance': float,
                       'nexposures': str, 'obs_seq': str, 'max_fwhm': float,
                       'min_moon_dist': float, 'max_moon_illum': float,
                       'max_cloud_cover': float, 'seq_repeats': int,
                       'seq_completed': int, 'status': str}
        rid = _id_from_time()
        pardic['id'] = rid
        requests = self.execute_sql(
            "SELECT object_id, allocation_id FROM request "
            "WHERE status != 'EXPIRED';")
        # check allocation_id, issue warning if it is a repeat, but allow
        if (pardic['object_id'], pardic['allocation_id']) in requests:
            print("program %s has already requested object: %s"
                  % (pardic['allocation_id'], pardic['object_id']))
            # TODO: add to output
        if pardic['object_id'] not in [obj[0] for obj in
                                       self.execute_sql(
                                           'SELECT id FROM object;')]:
            return -1, "ERROR: object does not exist!"
        if pardic['user_id'] not in [user[0] for user in
                                     self.execute_sql('SELECT id FROM users;')]:
            return -1, "ERROR: user does not exist!"
        if 'seq_repeats' not in pardic.keys():
            pardic['seq_repeats'] = 1

        if 'obs_seq' in pardic.keys():
            nexpo = pardic['obs_seq'][1:-1].split(',')
            nexposure = [0, 0, 0, 0, 0]
            for entry in nexpo:
                if entry[1:] == 'ifu':
                    nexposure[0] += int(entry[0])
                elif entry[1:] == 'u':
                    nexposure[1] += int(entry[0])
                elif entry[1:] == 'g':
                    nexposure[2] += int(entry[0])
                elif entry[1:] == 'r':
                    nexposure[3] += int(entry[0])
                elif entry[1:] == 'i':
                    nexposure[4] += int(entry[0])
            # nexposure = [n * pardic['seq_repeats'] for n in nexposure]
            # make sure that 'nexposures' and 'obs_seq' are consistent
            # if 'nexposures' in pardic.keys():
            #    if not '{
            #   %s, %s, %s, %s, %s}' % tuple(nexposure) == pardic['nexposures']:
            #        return (-1,
            #                "ERROR: nexposures and obs_seq are inconsistent!")
            # else:
            #    pardic['nexposures'] = '{
            #           %s, %s, %s, %s, %s}' % tuple(nexposure)

        elif not ('nexposures' in pardic.keys() or 'obs_seq' in pardic.keys()):
            return -1, "ERROR: nexposures or obs_seq is required!"

        keys = list(pardic.keys())
        default_params = ['object_id', 'user_id', 'allocation_id', 'exptime',
                          'priority', 'inidate', 'enddate']
        # check that all required values are provided
        for param in default_params:
            if param not in keys:
                return -1, "ERROR: %s not in dictionary!" % (param,)
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['id', 'object_id', 'user_id', 'allocation_id',
                           'exptime', 'priority', 'status', 'inidate',
                           'enddate', 'marshal_id', 'maxairmass', 'cadence',
                           'seq_completed', 'phasesamples', 'sampletolerance',
                           'nexposures', 'obs_seq', 'seq_repeats', 'max_fwhm',
                           'min_moon_dist', 'max_moon_illum', 'max_cloud_cover',
                           'external_id', 'shareid']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'request')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_request sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_request sql command failed with " \
                       "a ProgrammingError!"
        return rid, "Request added"

    def update_request(self, pardic):
        """
        Updates the request table with the parameters from the dictionary.

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'status' (str),
                    'maxairmass' (float),
                    'max_fwhm' (float),
                    'min_moon_dist' (float) (degrees),
                    'max_moon_illum' (float) (fraction),
                    'max_cloud_cover' (float) (fraction),
                    'priority' (float),
                    'seq_completed' (int),
                    'inidate' ('year-month-day'),
                    'enddate' ('year-month-day'),
                    'last_obs_jd' (float)
                Note: 'status' can be 'PENDING', 'ACTIVE', 'COMPLETED',
                                        'REDUCED', 'CANCELED', or 'EXPIRED'

        Returns:
            (-1, "ERROR...") if there was an issue with the updating

            (id, "Requests updated, columns 'column_names'")
                if the update was successful
        """
        param_types = {'id': int, 'object_id': int, 'user_id': int,
                       'allocation_id': int, 'exptime': str, 'priority': float,
                       'inidate': 'date', 'enddate': 'date', 'marshal_id': int,
                       'maxairmass': float, 'cadence': float,
                       'phasesamples': float, 'sampletolerance': float,
                       'nexposures': str, 'obs_seq': str, 'max_fwhm': float,
                       'min_moon_dist': float, 'max_moon_illum': float,
                       'max_cloud_cover': float, 'seq_repeats': int,
                       'seq_completed': int, 'status': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: no id provided!"
        elif not (isinstance(pardic['id'], int) or isinstance(pardic['id'],
                                                              long)):
            return -1, "ERROR: parameter id must be of type 'int'!"
        if pardic['id'] not in [rx[0] for rx in
                                self.execute_sql('SELECT id FROM request;')]:
            return -1, "ERROR: request does not exist!"
        keys.remove('id')
        if 'status' in keys:
            if pardic['status'] not in ['PENDING', 'ACTIVE', 'COMPLETED',
                                        'CANCELED', 'OBSERVED', 'EXPIRED']:
                return -1, "ERROR: %s is an invalid status value!" \
                       % (pardic['status'],)
        # remove any keys that are invalid or not allowed to be updated
        for key in reversed(keys):
            if key not in ['id', 'object_id', 'user_id', 'allocation_id',
                           'exptime', 'priority', 'status', 'inidate',
                           'enddate', 'marshal_id', 'maxairmass', 'cadence',
                           'seq_completed', 'phasesamples', 'sampletolerance',
                           'nexposures', 'obs_seq', 'seq_repeats', 'max_fwhm',
                           'min_moon_dist', 'max_moon_illum',
                           'max_cloud_cover']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_update_sql(pardic, keys, 'request', True)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_request sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_request sql command failed with " \
                       "a ProgrammingError!"

        return pardic['id'], "Requests updated, columns " + str(keys)[1:-1]

    def get_from_request(self, values, where_dict=None, compare_dict=None):
        """
        select values from `request`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'object_id' (int/long),
                'user_id' (int/long),
                'allocation_id' (int/long),
                'exptime' (str),
                'priority' (float),
                'inidate' ('year-month-day'),
                'enddate' ('year-month-day'),
                'marshal_id' (int/long),
                'maxairmass' (float),
                'cadence' (float),
                'phasesamples' (float),
                'sampletolerance' (float),
                'filters' (str),
                'nexposures' (str),
                'obs_seq' (str),
                'status' (str),
                'external_id' (int),
                'shareid' (int),
                'creationdate' ('year-month-day'),
                'lastmodified' ('year-month-day'),
                'max_fwhm' (float),
                'min_moon_dist' (float) (default 45 deg?),
                'max_moon_illum' (float) (fraction),
                'max_cloud_cover' (float) (fraction),
                'seq_repeats' (int),
                'seq_completed' (int),
                'last_obs_jd' (float)

        Returns:
            list of tuples containing the values for each request
            matching the criteria,

            empty list if no requests match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'object_id': int, 'user_id': int,
                          'allocation_id': int, 'exptime': str, 'status': str,
                          'external_id': int, 'shareid': int, 'priority': float,
                          'inidate': 'date', 'enddate': 'date',
                          'marshal_id': int, 'maxairmass': float,
                          'cadence': float, 'phasesamples': float,
                          'sampletolerance': float, 'filters': str,
                          'nexposures': str, 'obs_seq': str,
                          'creationdate': 'date', 'lastmodified': 'date',
                          'seq_repeats': int, 'seq_completed': int,
                          'last_obs_jd': float, 'max_fwhm': float,
                          'min_moon_dist': float, 'max_moon_illum': float,
                          'max_cloud_cover': float}
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'request')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def expire_requests(self, update_growth=True, update_fritz=True,
                        send_alerts=False, testing=False, interactive=False):
        """
        Updates the request table. For all the active requests that were
            not completed, and had an expiry date before than NOW(),
            are marked as "EXPIRED".

        Returns:
            (N, "Requests expired")
        """
        # tests written
        # 1. Get a list of request that are set to be expired.
        ret = self.execute_sql("SELECT u.id, o.name, u.email, r.marshal_id, "
                               "r.external_id, r.id "
                               "FROM request r "
                               "INNER JOIN object o ON (object_id = o.id) "
                               "INNER JOIN users u ON (user_id = u.id) "
                               "WHERE r.enddate < NOW() "
                               "AND (r.status = 'PENDING' "
                               "OR r.status = 'ACTIVE');")

        n_expired = len(ret)
        print("Found %d expiring requests" % n_expired)

        if testing:
            return n_expired, "Testing"

        for items in ret:
            # Reset DB status in request table
            sql = "UPDATE request SET status='EXPIRED', lastmodified=NOW() " \
                  "WHERE id=%d;" % items[5]
            self.execute_sql(sql)
            print("Request %d status set to EXPIRED in database" % items[5])

            if items[3] < 0:
                print("Calib request %s expired." % items[1])

            # Which marshal are we from?
            elif items[4] == 2:    # Fritz request
                print("Expiring Fritz request %d for target %s" % (items[3],
                                                                   items[1]))
                if send_alerts:
                    print("Sending email to %s" % items[2])
                    self.send_email_by_request(
                        to=items[2],
                        subject='Fritz request for target %s '
                                'has expired' % items[1],
                        template='expired_request',
                        template_dict={
                            'object_name': items[1],
                            'marshal_url': fritz_view_source_url}
                    )
                if update_fritz:
                    print("Updating Fritz marshal")
                    res = update_status_request("Expired", items[3], 'fritz')
                    print(res)

            else:               # Growth request
                print("Expiring Growth request %d for target %s" % (items[3],
                                                                    items[1]))
                if send_alerts:
                    print("Sending email to %s" % items[2])
                    self.send_email_by_request(
                        to=items[2],
                        subject='Growth request for target %s '
                                'has expired' % items[1],
                        template='expired_request',
                        template_dict={
                            'object_name': items[1],
                            'marshal_url': growth_view_source_url}
                    )
                if update_growth:
                    from growth import growth
                    # if the entry is greater than 1000
                    # then it should have come from the GROWTH MARSHAL

                    if items[3] and items[3] > 1000:
                        print("Updating Growth marshal")
                        res = growth.update_request(request_id=items[3],
                                                    output_dir='/scr/rsw/',
                                                    status='EXPIRED')
                        print(res)
            if interactive:
                q = input("next: ")
                if "Q" in q.upper():
                    break

        # sql = "UPDATE request SET status='EXPIRED', lastmodified=NOW() " \
        #      "WHERE enddate < NOW() AND status ='PENDING';"
        # self.execute_sql(sql)

        return n_expired, "Requests expired"

    def add_observation(self, header_dict):
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
        header_types = {'id': int, 'object_id': int, 'request_id': int,
                        'mjd': float, 'airmass': float, 'airmass_end': float,
                        'parang': float, 'parang_end': float, 'exptime': float,
                        'fitsfile': str, 'lst': str, 'ra': float, 'dec': float,
                        'tel_az': float, 'tel_el': float, 'tel_pa': float,
                        'ra_off': float, 'dec_off': float, 'imtype': str,
                        'camera': str, 'filter': str,
                        'time_elapsed': float}

        required_keys = ['object_id', 'request_id', 'mjd', 'airmass', 'exptime',
                         'fitsfile', 'lst', 'ra', 'dec', 'tel_az', 'tel_el',
                         'tel_pa', 'ra_off', 'dec_off']

        new_observation_id = _id_from_time()
        header_dict['id'] = new_observation_id

        header_keys = list(header_dict.keys())

        # Test for required keys
        for key in required_keys:
            if key not in header_keys:
                return -1, "ERROR: %s not provided!" % (key,)
        # Test for valid keys
        for key in reversed(header_keys):
            if key not in header_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(header_keys, header_dict, header_types)
        if type_check:
            return -1, type_check
        # Check if observation already exists in DB
        observation_id = self.get_from_observation(['id'],
                                                   {'fitsfile':
                                                    header_dict['fitsfile']},
                                                   {'fitsfile': '~'})
        if observation_id:
            return -1, "ERROR: %s already in database with id %d!" % \
                   (header_dict['fitsfile'], int(observation_id[0][0]))
        else:
            sql = _generate_insert_sql(header_dict, header_keys, 'observation')
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return -1, "ERROR(exc): adding observation sql command " \
                           "failed with an IntegrityError!"
            except psycopg2.IntegrityError:
                return -1, "ERROR(psycopg2): adding observation sql command " \
                           "failed with an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: adding observation sql command failed " \
                           "with a ProgrammingError!"
            except psycopg2.errors.NumericValueOutOfRange:
                return -1, "ERROR: adding observation sql command failed " \
                           "with a NumericValueOutOfRange error!"
            return new_observation_id, "Observation added"

    def update_observation(self, pardic):
        """

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
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

            (id, "Observation updated, columns 'column_names'")
                    if it completed successfully
        """
        param_types = {'id': int, 'object_id': int, 'request_id': int,
                       'mjd': float, 'airmass': float, 'airmass_end': float,
                       'parang': float, 'parang_end': float, 'exptime': float,
                       'fitsfile': str, 'lst': str, 'ra': float, 'dec': float,
                       'tel_az': float, 'tel_el': float, 'tel_pa': float,
                       'ra_off': float, 'dec_off': float, 'imtype': str,
                       'camera': str, 'filter': str,
                       'time_elapsed': float}

        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"
        elif pardic['id'] not in [ox[0] for ox in self.execute_sql(
                'SELECT id FROM observation;')]:
            return -1, "ERROR: observation does not exist!"
        keys.remove('id')

        for key in reversed(keys):
            if key not in param_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_update_sql(pardic, keys, 'observation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR(exc): adding observation sql command " \
                       "failed with an IntegrityError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): adding observation sql command " \
                       "failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_observation sql command failed with a " \
                       "ProgrammingError!"
        return pardic['id'], "Observation updated, columns " + str(keys)[1:-1]

    def get_from_observation(self, values, where_dict=None, compare_dict=None):
        """
        select values from `observation`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
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
                'airmass_end' (float),
                'parang' (float),
                'parang_end' (float),
                'ra_off' (float),
                'dec_off' (float),
                'imtype' (str),
                'time_elapsed' (datetime.timedelta object or
                               float/int seconds),
                'filter' (str),
                    options - 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b', 'NA'
                'camera' (str)

        Returns:
            list of tuples containing the values for each observation
            matching the criteria

            empty list if no observations match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'object_id': int, 'request_id': int,
                          'mjd': float, 'airmass': float, 'airmass_end': float,
                          'parang': float, 'parang_end': float,
                          'exptime': float, 'fitsfile': str, 'lst': str,
                          'ra': float, 'dec': float, 'tel_az': float,
                          'tel_el': float, 'tel_pa': float, 'ra_off': float,
                          'dec_off': float, 'imtype': str, 'camera': str,
                          'filter': str, 'time_elapsed': float}
        # checks type
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'observation')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_telescope_stats(self, tel_stats):
        """
        Adds telescope stats associated with an observation

        Args:
            tel_stats:
                required:
                    'observation_id' (int/long),
                    'date' ('year-month-day'),
                    'dome_status' (str),
                    'in_temp' (float),
                    'in_humidity' (float),
                    'in_dew' (float),
                    'out_temp' (float),
                    'out_humidity' (float),
                    'out_dew' (float),
                    'wind_dir' (float),
                    'wsp_cur' (float),
                    'wsp_avg' (float),
                    'mir_temp' (float),
                    'top_air' (float),
                    'pri_temp' (float),
                    'sec_temp' (float),
                    'flo_temp' (float),
                    'bot_temp' (float),
                    'mid_temp' (float),
                    'top_temp' (float)

        Returns:
            (-1, "ERROR...") if an issue occurs

            (id (long), "Telescope stats added") if successful
        """
        telstat_types = {'id': int, 'date': 'date', 'dome_status': str,
                         'in_temp': float, 'in_humidity': float,
                         'in_dew': float, 'out_temp': float,
                         'out_humidity': float, 'out_dew': float,
                         'wind_dir': float, 'wsp_cur': float, 'wsp_avg': float,
                         'mir_temp': float, 'top_air': float, 'pri_temp': float,
                         'sec_temp': float, 'flo_temp': float,
                         'bot_temp': float, 'mid_temp': float,
                         'top_temp': float, 'observation_id': int}
        stat_keys = list(tel_stats.keys())
        for key in reversed(stat_keys):
            if key not in ['id', 'date', 'dome_status', 'in_temp',
                           'in_humidity', 'in_dew', 'out_temp', 'out_humidity',
                           'out_dew', 'wind_dir', 'wsp_cur', 'wsp_avg',
                           'mir_temp', 'top_air', 'pri_temp', 'sec_temp',
                           'flo_temp', 'bot_temp', 'mid_temp', 'top_temp',
                           'observation_id']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(stat_keys, tel_stats, telstat_types)
        if type_check:
            return -1, type_check

        stat_sql = _generate_insert_sql(tel_stats, stat_keys, 'telescope_stats')
        try:
            self.execute_sql(stat_sql)
        except exc.IntegrityError:
            return -1, "ERROR: adding tel_stats sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: adding tel_stats sql command failed with " \
                       "a ProgrammingError!"

        return id, "Telescope stats added"

    # TODO: write update_observation() and update_telescope_stats()

    def get_from_telescope_stats(self, values, where_dict=None,
                                 compare_dict=None):
        """
        select values from `telescope_stats`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'observation_id' (int/long),
                'date' ('year-month-day'),
                'dome_status' (str),
                'in_temp' (float),
                'in_humidity' (float),
                'in_dew' (float),
                'out_temp' (float),
                'out_humidity' (float),
                'out_dew' (float),
                'wind_dir' (float),
                'wsp_cur' (float),
                'wsp_avg' (float),
                'mir_temp' (float),
                'top_air' (float),
                'pri_temp' (float),
                'sec_temp' (float),
                'flo_temp' (float),
                'bot_temp' (float),
                'mid_temp' (float),
                'top_temp' (float)

        Returns:
            list of tuples containing the values for telescope stats
            matching the criteria

            empty list if no telescope_stats entries match
            the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'date': 'date', 'dome_status': str, 'in_temp': float,
                          'in_humidity': float, 'in_dew': float,
                          'out_temp': float, 'out_humidity': float,
                          'out_dew': float, 'wind_dir': float,
                          'wsp_cur': float, 'wsp_avg': float, 'mir_temp': float,
                          'top_air': float, 'pri_temp': float,
                          'sec_temp': float, 'flo_temp': float,
                          'bot_temp': float, 'mid_temp': float,
                          'top_temp': float, 'observation_id': int, 'id': int}
        # checks type and if the sql generation returned an error
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'telescope_stats')
        if sql[0] == 'E':
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    # TODO: separate add_phot/spec and update_phot/spec
    def add_phot(self, pardic):
        """
        Adds reduced photometry or, if the observation already has photometry,
        updates it

        Args:
            pardic (dict):
                required:
                    'phot_calib_id' (int/long),
                    'observation_id' (int/long),
                    'astrometry' ('true' or 'false'),
                    'filter' (str),
                    'reducedfile' (abspath str),
                    'sexfile' (abspath str),
                    'maskfile' (abspath str),
                    'pipeline' (str),
                    'marshal_phot_id' (int/long)

        Returns:
            (-1, "ERROR...") if there was an issue

            (id (long), "Photometry added") if the photometry was
                    added successfully

            (id (long), "Photometry updated for observation_id ...")
                if the photometry existed and was updated
        """
        param_types = {'id': int, 'observation_id': int, 'astrometry': 'bool',
                       'filter': str, 'reducedfile': str, 'sexfile': str,
                       'maskfile': str, 'pipeline': str, 'marshal_phot_id': int,
                       'phot_calib_id': int}
        pid = _id_from_time()
        pardic['id'] = pid
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return -1, "ERROR: observation_id not provided!"
        phot_id = self.get_from_phot(
            ['id'], {'observation_id': pardic['observation_id']})
        # if there is already an entry for that observation, update instead
        if phot_id:
            if phot_id[0] == -1:
                return phot_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['astrometry', 'filter', 'reducedfile', 'sexfile',
                               'maskfile', 'pipeline', 'marshal_phot_id',
                               'phot_calib_id']:
                    return -1, "ERROR: %s is an invalid key!" % (key,)
            pardic['id'] = phot_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return -1, type_check

            update_sql = _generate_insert_sql(pardic, keys, 'spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_reduced_photometry update sql command " \
                           "failed with an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_reduced_photometry update sql command " \
                           "failed with a ProgrammingError!"
            return (phot_id[0][0], "Photometry updated for observation_id "
                                   "%s, columns " % (pardic['observation_id'],)
                    + str(keys)[1:-1])

        obs = self.get_from_observation(['fitsfile'],
                                        {'id': pardic['observation_id']})
        if not obs:
            return -1, "ERROR: no observation with the observation_id"
        elif obs[0] == -1:
            return obs
        else:
            pass
            # TODO: generate the filter here?

        phot_calib = self.get_from_phot_calib(['id'],
                                              {'id': pardic['phot_calib_id']})
        if not phot_calib:
            return -1, "ERROR: no phot_calib with the phot_calib_id"
        elif phot_calib[0] == -1:
            return phot_calib

        for key in ['observation_id', 'astrometry', 'filter', 'reducedfile',
                    'sexfile', 'maskfile', 'pipeline',
                    'phot_calib_id']:  # include 'marshal_phot_id'?
            if key not in keys:
                return -1, "ERROR: %s not provided!" % (key,)

        for key in reversed(keys):  # remove any invalid keys
            if key not in ['id', 'observation_id', 'astrometry', 'filter',
                           'reducedfile', 'sexfile', 'maskfile', 'pipeline',
                           'marshal_phot_id', 'phot_calib_id']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_reduced_photometry sql command failed " \
                       "with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_reduced_photometry sql command failed " \
                       "with a ProgrammingError!"

    def get_from_phot(self, values, where_dict=None, compare_dict=None):
        """
        select values from `phot`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'phot_calib_id' (int/long),
                'observation_id' (int/long),
                'astrometry' ('true' or 'false'),
                'filter' (str),
                'reducedfile' (abspath str),
                'sexfile' (abspath str),
                'maskfile' (abspath str),
                'pipeline' (str),
                'marshal_phot_id' (int/long)

        Returns:
            list of tuples containing the values for phot entries
            matching the criteria

            empty list if no phot entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'observation_id': int, 'astrometry': 'bool',
                          'filter': str, 'reducedfile': str, 'sexfile': str,
                          'maskfile': str, 'pipeline': str,
                          'marshal_phot_id': int, 'phot_calib_id': int,
                          'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'phot')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_spec(self, pardic, update=False):
        """
        Adds the reduced spectrum or, if the observation already has a spectrum,
            updates it

        Args:
            pardic (dict):
                required:
                    'observation_id' (int/long),
                    'spec_calib_id' (int/long),
                    'asciifile' (abspath str),
                    'quality' (int),
                optional:
                    'fitsfile' (abspath str),
                    'npyfile' (abspath str),
                    'imgset' (str),
                    'cubefile' (abspath str),
                    'standardfile' (abspath str),
                    'marshal_id' (int/long),
                    'skysub' ('true' or 'false'),
                    'extract_x' (float),
                    'extract_y' (float),
                    'extract_pa' (float),
                    'extract_a' (float),
                    'extract_b' (float),
                    'ad_red' (float),
                    'ad_blue' (float),
                    'prlltc' (float),
                    'flexure_x_corr_nm' (float),
                    'flexure_y_corr_pix' (float),
                    'reducer' (str),
                    'fwhm' (float),
                    'background' (float),
                    'line_fwhm' (float),
                    'pos_ok' (bool),
                    'pos_x_spax' (float),
                    'pos_y_spax' (float),
                    'psf_model' (str),
                    'psf_fwhm' (float),
                    'psf_ell' (float),
                    'psf_adr_pa' (float),
                    'psf_adr_z' (float),
                    'psf_adr_c2' (float),
                    'fluxcal' (bool),
                    'fluxcalfile', (abspath str)
                    'extr_type' (str)
            update (bool): set to true to update an existing entry

        Returns:
            (-1, "ERROR...") if there was an issue

            (id (long), "Spectrum added")  if the spectrum
                                            was added successfully

            (id (long), "Spectrum updated for observation_id ...")
                        if the spectrum existed and was updated
        """
        # TODO: add the following columns once created in db table 'spec'
        # 'cube_id': int
        param_types = {'id': int, 'spec_calib_id': int, 'observation_id': int,
                       'asciifile': str, 'npyfile': str, 'fitsfile': str,
                       'imgset': str, 'quality': int, 'cubefile': str,
                       'standardfile': str, 'marshal_spec_id': int,
                       'skysub': bool, 'extract_x': float, 'extract_y': float,
                       'extract_pa': float, 'extract_a': float,
                       'extract_b': float, 'ad_red': float, 'ad_blue': float,
                       'prlltc': float, 'reducer': str, 'airmass': float,
                       'atmcorr': float,
                       'fwhm': float, 'background': float, 'line_fwhm': float,
                       'pos_ok': bool, 'srcpos': str, 'pos_x_spax': float,
                       'pos_y_spax': float, 'psf_model': str, 'psf_fwhm': float,
                       'psf_ell': float, 'psf_adr_pa': float,
                       'psf_adr_z': float, 'psf_adr_c2': float, 'fluxcal': bool,
                       'fluxcalfile': str, 'extr_type': str, 'cube_id': int}

        required_keys = ('spec_calib_id', 'observation_id', 'fitsfile',
                         'quality')

        # Are we updating?
        if update:
            # Must have ID
            if 'id' in pardic:
                spec_id = pardic['id']
            else:
                return -1, "ERROR: must include id if updating record"
            # Verify spec_id
            if spec_id is None or spec_id <= 0:
                return spec_id, "ERROR: something went wrong getting spec id!"
            # Validate input keys
            keys = list(pardic.keys())
            for key in reversed(keys):  # TODO: test the updating
                if key not in param_types:
                    return -1, "ERROR: %s is an invalid key!" % (key,)
            # Update existing record with input values
            update_sql = _generate_update_sql(pardic, keys, 'spec')
            # print(update_sql)
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_spec sql command failed with an " \
                           "IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_spec sql command failed with a " \
                           "ProgrammingError!"
            return (spec_id,
                    "Spectrum updated for spec_id %s, columns " %
                    (pardic['id'],)
                    + str(keys)[1:-1])
        # New entry in spec table
        else:
            new_spec_id = _id_from_time()
            if not new_spec_id:
                new_spec_id = _id_from_time()
            if not new_spec_id:
                return -1, "ERROR: bad id from time"
            pardic['id'] = new_spec_id

            # get pardic keys
            keys = list(pardic.keys())

            # Test if observation id is valid
            obs_id = self.get_from_observation(['id'],
                                               {'id': pardic['observation_id']})
            if not obs_id:
                return -1, "ERROR: no observation exists with the given id!"
            elif obs_id[0] == -1:
                return obs_id, "ERROR: something went wrong " \
                               "getting observation id!"
            # Test for required keys
            for key in required_keys:
                if key not in keys:
                    return -1, "ERROR: %s not provided!" % (key,)
            # Test for valid keys
            for key in reversed(keys):
                if key not in param_types:
                    return -1, "ERROR: %s is an invalid key!" % (key,)
            # Check input parameter types
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return -1, type_check
            # Generate SQL command
            sql = _generate_insert_sql(pardic, keys, 'spec')
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_spec sql command failed with an " \
                           "IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_spec sql command failed with a " \
                           "ProgrammingError!"
            return new_spec_id, "Spectrum added with columns " + str(keys)[1:-1]

    def get_from_spec(self, values, where_dict=None, compare_dict=None):
        """
        select values from `spec`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'observation_id' (int/long),
                'spec_calib_id' (int/long),
                'asciifile' (abspath str),
                'quality' (int),
                'fitsfile' (abspath str),
                'npyfile' (abspath str),
                'imgset' (str),
                'cubefile' (abspath str),
                'standardfile' (abspath str),
                'marshal_id' (int/long),
                'skysub' ('true' or 'false'),
                'extract_x' (float),
                'extract_y' (float),
                'extract_pa' (float),
                'extract_a' (float),
                'extract_b' (float),
                'ad_red' (float),
                'ad_blue' (float),
                'prlltc' (float),
                'flexure_x_corr_nm' (float),
                'flexure_y_corr_pix' (float),
                'reducer' (str),
                'fwhm' (float),
                'background' (float),
                'line_fwhm' (float),
                'pos_ok' (bool),
                'pos_x_spax' (float),
                'pos_y_spax' (float),
                'psf_model' (str),
                'psf_fwhm' (float),
                'psf_ell' (float),
                'psf_adr_pa' (float),
                'psf_adr_z' (float),
                'psf_adr_c2' (float),
                'fluxcal' (bool),
                'fluxcalfile', (abspath str)
                'extr_type' (str)

        Returns:
            list of tuples containing the values for spectra matching
                the criteria

            empty list if no spec entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'spec_calib_id': int,
                          'observation_id': int, 'asciifile': str,
                          'npyfile': str, 'fitsfile': str, 'imgset': str,
                          'quality': int, 'cubefile': str, 'standardfile': str,
                          'marshal_spec_id': int, 'skysub': bool,
                          'extract_x': float, 'extract_y': float,
                          'extract_pa': float, 'extract_a': float,
                          'extract_b': float, 'ad_red': float, 'ad_blue': float,
                          'prlltc': float, 'flexure_x_corr_nm': float,
                          'flexure_y_corr_pix': float, 'reducer': str,
                          'fwhm': float, 'background': float,
                          'line_fwhm': float, 'pos_ok': bool, 'srcpos': str,
                          'pos_x_spax': float, 'pos_y_spax': float,
                          'psf_model': str, 'psf_fwhm': float, 'psf_ell': float,
                          'psf_adr_pa': float, 'psf_adr_z': float,
                          'psf_adr_c2': float, 'fluxcal': bool,
                          'fluxcalfile': str, 'extr_type': str}

        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'spec')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_metrics_phot(self, pardic):
        """
        Creates an entry in metrics_phot or updates an existing metrics entry

        Args:
            pardic (dict):
                required:
                    'phot_id' (int/long),
                    'fwhm' (float),
                    'background' (float),
                    'zp' (float),
                    'zperr' (float),
                    'ellipticity' (float),
                    'nsources' (int)

        Returns:
            (-1, "ERROR...") if there was an issue

            (id (long), "Photometry metrics updated for phot_id ...")
                            if it updated existing metrics

            (id (long), "Photometry metrics added") if the metrics
                            were added successfully
        """
        param_types = {'id': int, 'phot_id': int, 'fwhm': float,
                       'background': float, 'zp': float, 'zperr': float,
                       'ellipticity': float, 'nsources': int}
        mid = _id_from_time()
        pardic['id'] = mid
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'phot_id' not in keys:
            return -1, "ERROR: phot_id not provided!"
        metric_id = self.get_from_metrics_phot(['id'],
                                               {'phot_id': pardic['phot_id']})
        # if there is already an entry for that observation, update instead
        if metric_id:
            if metric_id[0] == -1:
                return metric_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fwhm', 'background', 'zp', 'zperr',
                               'ellipticity', 'nsources']:
                    return -1, "ERROR: %s is an invalid key!" % (key,)
            pardic['id'] = metric_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return -1, type_check

            update_sql = _generate_insert_sql(pardic, keys, 'metrics_phot')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_metrics_phot update sql command " \
                           "failed with an IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_metrics_phot update sql command " \
                           "failed with a ProgrammingError!"
            return (metric_id[0][0], "Photometry metrics updated for phot_"
                                     "id %s, columns " % (pardic['phot_id'],)
                    + str(keys)[1:-1])
        ph_id = self.get_from_phot(['id'], {'id': pardic['phot_id']})
        if not ph_id:
            return -1, "ERROR: no photometry exists with the given id!"
        elif ph_id[0] == -1:
            return ph_id

        for key in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity',
                    'nsources']:  # phot_id already tested
            if key not in keys:
                return -1, "ERROR: %s not provided!" % (key,)

        for key in reversed(keys):
            if key not in ['id', 'phot_id', 'fwhm', 'background', 'zp', 'zperr',
                           'ellipticity', 'nsources']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'metrics_phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_metrics_phot sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_metrics_phot sql command failed with " \
                       "a ProgrammingError!"
        return mid, "Photometry metrics added"

    def get_from_metrics_phot(self, values, where_dict=None, compare_dict=None):
        """
        select values from `metrics_phot`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'phot_id' (int/long),
                'fwhm' (float),
                'background' (float),
                'zp' (float),
                'zperr' (float),
                'ellipticity' (float),
                'nsources' (int)

        Returns:
            list of tuples containing the values for metrics
            matching the criteria

            empty list if no metrics_phot entries match
            the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'phot_id': int, 'fwhm': float, 'background': float,
                          'zp': float, 'zperr': float, 'ellipticity': float,
                          'nsources': int, 'id': int}

        # checks type and if the sql generation returned an error
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'metrics_phot')
        if sql[0] == 'E':
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_phot_calib(self, pardic):
        """
        Creates a new object in the phot calib table with the parameters
        specified in the dictionary.

        Args:
            pardic:
                required:
                    'bias' (abspath str)
                    'flat' (abspath str)

        Returns:
            (-1: "ERROR...") if there is an issue

            (id (long), "Photometry calibration added")
                                if it completes successfully
        """
        param_types = {'id': int, 'bias': str, 'flat': str}
        cid = _id_from_time()
        pardic['id'] = cid
        keys = list(pardic.keys())

        for key in ['bias', 'flat']:  # phot_id already tested
            if key not in keys:
                return -1, "ERROR: %s not provided!" % (key,)

        for key in reversed(keys):
            if key not in ['id', 'bias', 'flat']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'phot_calib')
        print(sql, type_check)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_phot_calib sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_phot_calib sql command failed with " \
                       "a ProgrammingError!"
        return cid, "Photometry calibration added"

    def update_phot_calib(self, pardic):
        """
        updates a phot_calib entry

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'bias' (abspath str)
                    'flat' (abspath str)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Phot_calib updated, columns 'column_names'")
                            if the entry is updated successfully
        """
        param_types = {'id': int, 'bias': str, 'flat': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"

        elif pardic['id'] not in [cx[0] for cx in
                                  self.execute_sql(
                                      'SELECT id FROM phot_calib;')]:
            return -1, "ERROR: no phot_calib entry with the id!"
        keys.remove('id')
        # remove any keys that are invalid or not allowed to be updated
        for key in reversed(keys):
            if key not in ['bias', 'flat']:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_update_sql(pardic, keys, 'phot_calib')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_phot_calib sql command failed with " \
                       "an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_phot_calib sql command failed with " \
                       "a ProgrammingError!"
        return pardic['id'], "Phot_calib updated, columns " + str(keys)[1:-1]

    def get_from_phot_calib(self, values, where_dict=None, compare_dict=None):
        """
        select values from `phot_calib`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long)
                'bias' (abspath str)
                'flat' (abspath str)
        Returns:
            list of tuples containing the values for calibration
            matching the criteria

            empty list if no phot_calib entries
            match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'bias': str, 'flat': str}

        # checks type and if the sql generation returned an error
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'phot_calib')
        if sql[0] == 'E':
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_spec_calib(self, pardic):
        """
        Creates a new object in the spec calib table with the parameters
        specified in the dictionary.

        Args:
            pardic:
                required:
                    'flat' (abspath str)
                    'hexagrid' (abspath str)
                    'tracematch' (abspath str)
                    'tracematch_withmasks' (abspath str)
                    'wavesolution' (abspath str)
                    'utdate' (YYYYMMDD date str)

                optional:
                    'dome_master' (abspath str)
                    'bias_slow_master' (abspath str)
                    'bias fast_master' (abspath str)
                    'cosmic_filter' (bool)
                    'drpver' (float)
                    'hg_master' (abspath str)
                    'cd_master' (abspath str)
                    'xe_master' (abspath str)
                    'avg_rms' (abspath str)
                    'min_rms' (abspath str)
                    'max_rms' (abspath str)
                    'dispersionmap' (abspath str)
                    'flatmap' (abspath str)
                    'nspaxels' (int)
                    'wave_rms_avg' (float)
                    'wave_rms_min' (float)
                    'wave_rms_max' (float)
                    'width_rms_avg' (float)
                    'width_rms_min' (float)
                    'width_rms_max' (float)

        Returns:
            (-1: "ERROR...") if there is an issue

            (id (long), "Spectrum calibration added") if successfully added
        """
        param_types = {'id': int, 'dome_master': str, 'bias_slow_master': str,
                       'bias_fast_master': str, 'flat': str,
                       'cosmic_filter': bool, 'drpver': str, 'hg_master': str,
                       'cd_master': str, 'xe_master': str, 'avg_rms': str,
                       'min_rms': str, 'max_rms': str, 'hexagrid': str,
                       'tracematch': str, 'tracematch_withmasks': str,
                       'wavesolution': str, 'dispersionmap': str,
                       'flatmap': str, 'nspaxels': int, 'wave_rms_avg': float,
                       'wave_rms_min': float, 'wave_rms_max': float,
                       'width_rms_avg': float, 'width_rms_min': float,
                       'width_rms_max': float, 'utdate': str}

        required_keys = ['flat', 'hexagrid', 'tracematch',
                         'tracematch_withmasks', 'wavesolution', 'utdate']

        new_spec_calib_id = _id_from_time()
        pardic['id'] = new_spec_calib_id

        keys = list(pardic.keys())

        # Test for required keys
        for rkey in required_keys:
            if rkey not in keys:
                return -1, "ERROR: %s not provided!" % (rkey,)

        # Test if spec_calib already has an entry on this utdate
        spec_calib_id = self.get_from_spec_calib(['id'],
                                                 {'utdate': pardic['utdate']})
        # Already exists
        if spec_calib_id:
            # Verify spec_calib_id
            if spec_calib_id[0] == -1:
                return spec_calib_id, "ERROR: something went wrong " \
                                      "getting spec calib id!"
            # Validate input keys
            for key in reversed(keys):
                if key not in param_types:
                    return -1, "ERROR: %s is an invalid key!" % (key,)
            # Check input parameter types
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return -1, type_check
            # Update existing record with input values
            pardic['id'] = spec_calib_id[0][0]
            update_sql = _generate_update_sql(pardic, keys, 'spec_calib')
            print(update_sql)
            # Attempt to update record
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_spec sql command failed with an " \
                           "IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_spec sql command failed with a " \
                           "ProgrammingError!"
            except psycopg2.IntegrityError:
                return -1, "ERROR(psycopg2): updating spec_calib sql command " \
                           "failed with an IntegrityError!"
            return (spec_calib_id[0][0],
                    "Spec calib updated for spec_calib_id %s, columns " %
                    (pardic['id'],) + str(keys)[1:-1])
        # New entry
        else:
            # Validate input keys
            for key in reversed(keys):
                if key not in param_types:
                    return -1, "ERROR: %s is an invalid key!" % (key,)
            # Check input parameter types
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return -1, type_check
            # Generate SQL command
            sql = _generate_insert_sql(pardic, keys, 'spec_calib')
            print(sql)
            # Attempt to insert record
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return -1, "ERROR: add_spec_calib sql command failed with an " \
                           "IntegrityError!"
            except exc.ProgrammingError:
                return -1, "ERROR: add_spec_calib sql command failed with a " \
                           "ProgrammingError!"
            except psycopg2.IntegrityError:
                return -1, "ERROR(psycopg2): inserting spec_calib sql " \
                           "command failed with an IntegrityError!"
            return new_spec_calib_id, "Spectrum calibration added"

    def update_spec_calib(self, pardic):
        """
        updates a spec_calib entry

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'flat' (abspath str)
                    'hexagrid' (abspath str)
                    'tracematch' (abspath str)
                    'tracematch_withmasks' (abspath str)
                    'wavesolution' (abspath str)
                    'utdate' (YYYYMMDD date str)
                    'dome_master' (abspath str)
                    'bias_slow_master' (abspath str)
                    'bias fast_master' (abspath str)
                    'cosmic_filter' (bool)
                    'drpver' (float)
                    'hg_master' (abspath str)
                    'cd_master' (abspath str)
                    'xe_master' (abspath str)
                    'avg_rms' (abspath str)
                    'min_rms' (abspath str)
                    'max_rms' (abspath str)
                    'dispersionmap' (abspath str)
                    'flatmap' (abspath str)
                    'nspaxels' (int)
                    'wave_rms_avg' (float)
                    'wave_rms_min' (float)
                    'wave_rms_max' (float)
                    'width_rms_avg' (float)
                    'width_rms_min' (float)
                    'width_rms_max' (float)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Spec_calib updated, columns 'column_names'")
                        if the entry is updated successfully
        """
        param_types = {'id': int, 'dome_master': str, 'bias_slow_master': str,
                       'bias_fast_master': str, 'flat': str,
                       'cosmic_filter': bool, 'drpver': str, 'hg_master': str,
                       'cd_master': str, 'xe_master': str, 'avg_rms': str,
                       'min_rms': str, 'max_rms': str, 'hexagrid': str,
                       'tracematch': str, 'tracematch_withmasks': str,
                       'wavesolution': str, 'dispersionmap': str,
                       'flatmap': str, 'nspaxels': int, 'wave_rms_avg': float,
                       'wave_rms_min': float, 'wave_rms_max': float,
                       'width_rms_avg': float, 'width_rms_min': float,
                       'width_rms_max': float, 'utdate': str}

        keys = list(pardic.keys())
        # Check if id included
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"
        # Test if id exists
        elif pardic['id'] not in [sx[0] for sx in self.execute_sql(
                'SELECT id FROM spec_calib;')]:
            return -1, "ERROR: no spec_calib entry with the id!"
        keys.remove('id')
        # Check if input keys are valid
        for key in reversed(keys):
            if key not in param_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        # Do we have any keys left?
        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        # Verify types
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check
        # Make update sql command
        sql = _generate_update_sql(pardic, keys, 'spec_calib')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_spec_calib sql command failed with an " \
                       "IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_spec_calib sql command failed with a " \
                       "ProgrammingError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): updating spec_calib sql command " \
                       "failed with an IntegrityError!"
        return pardic['id'], "Spec_calib updated, columns " + str(keys)[1:-1]

    def get_from_spec_calib(self, values, where_dict=None, compare_dict=None):
        """
        select values from `spec_calib`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long)
                'flat' (abspath str)
                'hexagrid' (abspath str)
                'tracematch' (abspath str)
                'tracematch_withmasks' (abspath str)
                'wavesolution' (abspath str)
                'utdate' (YYYYMMDD date str)
                'dome_master' (abspath str)
                'bias_slow_master' (abspath str)
                'bias fast_master' (abspath str)
                'cosmic_filter' (bool)
                'drpver' (float)
                'hg_master' (abspath str)
                'cd_master' (abspath str)
                'xe_master' (abspath str)
                'avg_rms' (abspath str)
                'min_rms' (abspath str)
                'max_rms' (abspath str)
                'dispersionmap' (abspath str)
                'flatmap' (abspath str)
                'nspaxels' (int)
                'wave_rms_avg' (float)
                'wave_rms_min' (float)
                'wave_rms_max' (float)
                'width_rms_avg' (float)
                'width_rms_min' (float)
                'width_rms_max' (float)
        Returns:
            list of tuples containing the values for calibration matching
                the criteria

            empty list if no spec_calib entries match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {'id': int, 'dome_master': str,
                          'bias_slow_master': str, 'bias_fast_master': str,
                          'flat': str, 'cosmic_filter': bool, 'drpver': str,
                          'hg_master': str, 'cd_master': str, 'xe_master': str,
                          'avg_rms': str, 'min_rms': str, 'max_rms': str,
                          'hexagrid': str, 'tracematch': str,
                          'tracematch_withmasks': str, 'wavesolution': str,
                          'dispersionmap': str, 'flatmap': str, 'nspaxels': int,
                          'wave_rms_avg': float, 'wave_rms_min': float,
                          'wave_rms_max': float, 'width_rms_avg': float,
                          'width_rms_min': float, 'width_rms_max': float,
                          'utdate': str}
        # checks types
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'spec_calib')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_cube(self, pardic):
        """
        Creates a new object in the cube table.

        Args:
            pardic (dict):
                required:
                    'observation_id' (int)
                    'ccd_x_flex_corr' (bool)
                    'ccd_x_flex_px' (float)
                    'ccd_y_flex_corr' (bool)
                    'ccd_y_flex_px' (float)
                    'atm_corr' (bool)
                    'atm_source' (str)
                    'atm_mean_corr' (float)
                    'spec_calib_id' (int)

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Cube added")
        """
        param_types = {
            'id': int, 'observation_id': int,
            'ccd_x_flex_corr': bool, 'ccd_x_flex_px': float,
            'ccd_y_flex_corr': bool, 'ccd_y_flex_px': float,
            'atm_corr': bool, 'atm_source': str, 'atm_mean_corr': float,
            'spec_calib_id': int
        }
        required_keys = ['observation_id', 'spec_calib_id',
                         'ccd_x_flex_corr', 'ccd_x_flex_px',
                         'ccd_y_flex_corr', 'ccd_y_flex_px',
                         'atm_corr', 'atm_source', 'atm_mean_corr',
                         'spec_calib_id']

        new_cube_id = _id_from_time()
        pardic['id'] = new_cube_id

        keys = list(pardic.keys())

        # Test for required keys
        for rkey in required_keys:
            if rkey not in keys:
                return -1, "ERROR: %s not provided!" % (rkey,)
        # Validate input keys
        for key in reversed(keys):
            if key not in param_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        # Check input parameter types
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'cube')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_cube sql command failed with an " \
                       "IntegrityError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): add_cube sql command failed with an " \
                       "IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_cube sql command failed with a " \
                       "ProgrammingError!"
        return new_cube_id, "Cube added"

    def update_cube(self, pardic):
        """
        Updates a cube entry.

        Args:
            pardic (dict):
                required:
                    'id' (int)
                optional:
                    'observation_id' (int)
                    'ccd_x_flex_corr' (bool)
                    'ccd_x_flex_px' (float)
                    'ccd_y_flex_corr' (bool)
                    'ccd_y_flex_px' (float)
                    'atm_corr' (bool)
                    'atm_source' (str)
                    'atm_mean_corr' (float)
                    'spec_calib_id' (int)

        Returns:
            (-1, "ERROR...") it failed to update

            (id (long), "Cube updated, columns 'column_names'")
                        if the entry is updated successfully
        """
        param_types = {
            'id': int, 'observation_id': int,
            'ccd_x_flex_corr': bool, 'ccd_x_flex_px': float,
            'ccd_y_flex_corr': bool, 'ccd_y_flex_px': float,
            'atm_corr': bool, 'atm_source': str, 'atm_mean_corr': float,
            'spec_calib_id': int
        }
        keys = list(pardic.keys())
        # Check if id included
        if 'id' not in keys:
            return -1, "ERROR: id not provided!"
        # Test if id exists
        elif pardic['id'] not in [cx[0] for cx in self.execute_sql(
                'SELECT id FROM cube;')]:
            return -1, "ERROR: no cube entry with the id!"
        keys.remove('id')
        # Check if input keys are valid
        for key in reversed(keys):
            if key not in param_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)
        # Check input parameter types
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check
        # Generate update sql command
        sql = _generate_update_sql(pardic, keys, 'cube')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_cube sql command failed with an " \
                       "IntegrityError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): add_cube sql command failed with an " \
                       "IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_cube sql command failed with a " \
                       "ProgrammingError!"
        return pardic['id'], "Cube updated, columns " + str(keys)[1:-1]

    def get_from_cube(self, values, where_dict=None, compare_dict=None):
        """
        select values from `flexure`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long)
                'observation_id' (int)
                'ccd_x_flex_corr' (bool)
                'ccd_x_flex_px' (float)
                'ccd_y_flex_corr' (bool)
                'ccd_y_flex_px' (float)
                'atm_corr' (bool)
                'atm_source' (str)
                'atm_mean_corr' (float)
                'spec_calib_id' (int)

        Returns:
            list of tuples containing the values for cube matching the criteria

            empty list if no cube entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {
            'id': int, 'observation_id': int,
            'ccd_x_flex_corr': bool, 'ccd_x_flex_px': float,
            'ccd_y_flex_corr': bool, 'ccd_y_flex_px': float,
            'atm_corr': bool, 'atm_source': str, 'atm_mean_corr': float,
            'spec_calib_id': int
        }
        # checks types
        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'cube')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def add_classification(self, pardic):
        """
        Creates a classification object attached to the reduced spectrum.

        Args:
            pardic (dict):
                required:
                    'spec_id' (int/long),
                    'object_id' (int/long),
                    'classification' (str),
                    'redshift' (float),
                    'redshift_err' (float),
                    'classifier' (str),
                    'score' (float)
                optional:
                    'auto' (bool)
                    'phase' (float),
                    'phase_err' (float),
                    'score_err' (float),
                    'score_type' (str),
                    'class_source' (str),
                    'class_template' (str)
        Returns:
            (-1, "ERROR...") if there is an issue

            (id (long), 'Classification added") if it was successful
        """
        param_types = {'id': int, 'spec_id': int, 'object_id': int,
                       'classification': str, 'redshift': float, 'auto': bool,
                       'redshift_err': float, 'classifier': str,
                       'score': float, 'score_err': float,
                       'phase': float, 'phase_err': float, 'score_type': str,
                       'class_source': str, 'class_template': str}
        required_keys = ['spec_id', 'object_id', 'classification', 'redshift',
                         'redshift_err', 'classifier', 'score']
        # Get id number
        class_id = _id_from_time()
        pardic['id'] = class_id

        keys = list(pardic.keys())

        # Do we have all required keys?
        for key in required_keys:
            if key not in keys:
                return -1, "ERROR: %s not provided!" % (key,)

        # Has this classification already been entered?
        classified = self.get_from_classification(
            ['classification', 'redshift', 'redshift_err', 'spec_id'],
            {'spec_id': pardic['spec_id'], 'classifier': pardic['classifier']})
        if classified:
            if classified[0] == -1:
                return classified
            return (-1, "ERROR: entry exists for that spectrum and classifier"
                        " with classification %s, redshift %s, "
                        "redshift_err %s. Use `update_classification` if "
                        "necessary."
                    % (classified[0][0], classified[0][1], classified[0][2]))

        # Make sure all keys are valid
        for key in reversed(keys):
            if key not in param_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)

        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sql = _generate_insert_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: add_classification sql command failed with an " \
                       "IntegrityError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): add_classification sql command " \
                       "failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: add_classification sql command failed with a " \
                       "ProgrammingError!"
        return class_id, "Classification added"

    def update_classification(self, pardic):
        """
        Update a classification

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                     OR
                     'spec_id' (int/long),
                     'classifier' (str)
                optional:
                    'object_it' (int/long),
                    'classification' (str),
                    'redshift' (float),
                    'redshift_err' (float),
                    'auto' (bool),
                    'phase' (float),
                    'phase_err' (float),
                    'score' (float),
                    'score_err' (float),
                    'redshift' (float),
                    'redshift_err' (float),
                    'classifier' (str),
                    'score' (float),
                    'score_type' (str),
                    'class_source' (str),
                    'class_template' (str),
                Note: this function will not modify id, spec_id,
                classifier or object_id

        Returns:
            (-1, "ERROR...") if there is an issue

            (id, "Classification updated, columns 'column_names'") if it was
            successful
        """
        param_types = {'id': int, 'spec_id': int, 'object_id': int,
                       'classification': str, 'redshift': float, 'auto': bool,
                       'redshift_err': float, 'classifier': str,
                       'score': float, 'score_err': float,
                       'phase': float, 'phase_err': float, 'score_type': str,
                       'class_source': str, 'class_template': str}

        keys = list(pardic.keys())

        # Check id
        if 'id' in keys:
            class_id = pardic['id']
            id_classifier = self.get_from_classification(
                ['spec_id', 'classifier'], {'id': pardic['id']})
            if not id_classifier:
                return -1, "ERROR: no classification entry with the given id!"
            elif id_classifier[0] == -1:
                return id_classifier
            # Must match classifier and spec_id
            if 'classifier' in keys:
                if not pardic['classifier'] == id_classifier[0][1]:
                    return -1, "ERROR: classifier provided does not match " \
                               "classification id!"
            if 'spec_id' in keys:
                if not pardic['spec_id'] == id_classifier[0][0]:
                    return -1, "ERROR: spec_id provided does not match " \
                               "classification id!"
            keys.remove('id')
        elif 'spec_id' in keys and 'classifier' in keys:
            class_id = self.get_from_classification(
                ['id'], {'spec_id': pardic['spec_id'],
                         'classifier': pardic['classifier']})
            if not class_id:
                return -1, "ERROR: no classification entry with the given " \
                           "spec_id and classifier!"
            elif class_id[0] == -1:
                return class_id
            pardic['id'] = class_id[0][0]
        else:
            return -1, "ERROR: needs id or both spec_id and classifier"

        # Make sure all keys are valid
        for key in reversed(keys):
            if key not in param_types:
                return -1, "ERROR: %s is an invalid key!" % (key,)

        if len(keys) == 0:
            return -1, "ERROR: no parameters given to update!"
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return -1, type_check

        sp_id = self.get_from_spec(['id'], {'id': pardic['spec_id']})
        if not sp_id:
            return -1, "ERROR: no spectrum exists with the given spec_id!"
        elif sp_id[0] == -1:
            return sp_id
        obj_id = self.get_from_object(['id'], {'id': pardic['object_id']})
        if not obj_id:
            return -1, "ERROR: no object exists with the given object_id!"
        elif obj_id[0] == -1:
            return obj_id

        sql = _generate_update_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: update_classification sql command failed with " \
                       "an IntegrityError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): add_classification sql command " \
                       "failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: update_classification sql command failed with " \
                       "a ProgrammingError!"
        return class_id, "Classification updated, columns " + str(keys)[1:-1]

    def get_from_classification(self, values, where_dict=None,
                                compare_dict=None):
        """
        select values from `classification`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            compare_dict (dict): default is {}
                'param': 'inequality' (i.e. '>', '<', '>=', '<=', '<>', '!='))
                if no inequality is provided, '=' is assumed
            values/keys options:
                'id' (int/long),
                'spec_id' (int/long),
                'classifier' (str),
                'classification' (str),
                'redshift' (float),
                'redshift_err' (float),
                'auto' (bool),
                'phase' (float),
                'phase_err' (float),
                'score' (float),
                'score_err' (float),
                'score_type' (str),
                'class_source' (str),
                'class_template' (str)

        Returns:
            list of tuples containing the values for classifications
            matching the criteria

            empty list if no classification entries match the
            ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        if where_dict is None:
            where_dict = {}
        if compare_dict is None:
            compare_dict = {}
        allowed_params = {
                'id': int, 'spec_id': int, 'object_id': int,
                'classification': str, 'redshift': float, 'auto': bool,
                'redshift_err': float, 'classifier': str,
                'score': float, 'score_err': float,
                'phase': float, 'phase_err': float,
                'score_type': str, 'class_source': str,
                'class_template': str
        }

        sql = _generate_select_sql(values, where_dict, allowed_params,
                                   compare_dict, 'classification')
        if sql[0] == 'E':  # if the sql generation returned an error
            return -1, sql

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return -1, "ERROR: sql command failed with an IntegrityError!"
        except psycopg2.IntegrityError:
            return -1, "ERROR(psycopg2): add_classification sql command " \
                       "failed with an IntegrityError!"
        except exc.ProgrammingError:
            return -1, "ERROR: sql command failed with a ProgrammingError!"
        return results

    def update_allocation_time(self, allocation_id):
        """
        Searches for all the completed requests under that allocation and
        updates the time take to execute it.
        Updates the field" time_spent

        Args:
            allocation_id (int): int
                Id of the allocation to update.
        """

        exptime = self.execute_sql("""SELECT exptime 
                                    FROM request 
                                    WHERE status ='COMPLETED' 
                                    AND allocation_id = %d;""" % allocation_id)

        if len(exptime) == 0:
            ret, msg = self.update_allocation(
                {"id": allocation_id, "time_spent": timedelta(days=0.0)})
        else:
            reqsum = [np.sum(e) for e in exptime]
            time_spent = float(np.sum(reqsum))

            ret, msg = self.update_allocation(
                {"id": allocation_id, "time_spent": time_spent})
        print(ret, msg)

    def get_allocation_spent_time(self, allocation_id, inidate, enddate):
        """
        Searches for all the completed requests under that allocation and
        updates the time take to execute it.
        Updates the field" time_spent

        Args:
            allocation_id (int): int
                Id of the allocation to update.

            inidate: datetime
                Starting time to compute the spent time for the allocation.

            enddate: datetime
                End of the time to compute the spent time for the allocation.
        """

        exptime = self.execute_sql("""SELECT exptime 
                                    FROM request 
                                    WHERE status ='COMPLETED' 
                                    AND allocation_id = %d
                                    AND lastmodified>DATE('%s') 
                                    AND lastmodified<DATE('%s');"""
                                   % (allocation_id, inidate, enddate))

        if len(exptime) == 0:
            time_spent = 0
        else:
            reqsum = [np.sum(e) for e in exptime]
            time_spent = np.sum(reqsum)

        return time_spent

    def update_all_allocations(self):
        """
        Updates all the active allocation's time.
        """
        alloc = self.get_from_allocation(["id"], {"active": True})
        al = list(set(alloc))
        for a in al:
            self.update_allocation_time(int(a[0]))

    def send_email_by_request(self, requestid=None, to=None, subject=None,
                              template='', template_dict=None):
        """
        Allow 
        :param requestid:
        :param to:
        :param subject:
        :param template:
        :param template_dict:
        :return: 
        """
        msg = EmailMessage()
        # 1. First start by making sure the request id is valid and that we can
        #    extract a user email from the request.  Ignore the email address
        #    in the request in the case where there is a "to" address given
        if to:
            msg['To'] = to
        else:
            try:
                userid = self.get_from_request(
                    values=['user_id'], where_dict={'id': requestid})[0][0]
            except IndexError:
                return "Unable to retrieve userid from database"

            try:
                user_email = self.get_from_users(values=['email'],
                                                 where_dict={'id': userid})
            except IndexError:
                return "Unable to retrieve user_email from database"

            msg['To'] = user_email[0]

        msg['Subject'] = subject

        # ToDo: Setup auto email user
        msg['From'] = 'No_reply_sedm_robot@astro.caltech.edu'

        # 2. Choose the email template to use
        template_file = '%s/email_templates/%s.txt' % (SITE_ROOT,
                                                       template.lower())
        print(template_file)
        if os.path.exists(template_file):
            with open(template_file) as f:
                email_template = f.read()

            email_content = email_template.format(**template_dict)

            msg.set_content(email_content)

            print(msg)

        else:
            return "Template does not exist"

        # 3. Send the message
        # Send the message via local SMTP server.
        with smtplib.SMTP('smtp-server.astro.caltech.edu') as s:
            s.send_message(msg)


def _data_type_check(keys, pardic, value_types):
    """
    make sure the values given match the data types required in the database

    Args:
        keys (list): keys to be tested
        pardic (dict): keys and keys' values
        value_types (dict): keys and types the values should be
                            (e.g. {'ra': float})

        types are: str, float, int, 'date', 'datetime', 'bool'

    Returns:
        "ERROR..." if a value was of the wrong type

        None if all matched
    """
    for key in keys:
        if value_types[key] == str:
            if pardic[key] is not None:
                pass  # all values are added to the queries as strings anyway
            else:
                return "ERROR: %s must not be None!" % (key,)
        elif value_types[key] == float:
            try:
                pardic[key] = float(pardic[key])
            except (ValueError, TypeError):
                return "ERROR: %s must be of %s!" % (key, str(float)[1:-1])
        elif value_types[key] == int:
            try:
                pardic[key] = long(pardic[key])
            except (ValueError, TypeError):
                return "ERROR: %s must be of %s or %s!" % (key, str(int)[1:-1],
                                                           str(long)[1:-1])
        elif value_types[key] == 'date':
            try:
                pardic[key] = str(Time(pardic[key])).split(' ')[0]
            except (ValueError, TypeError):
                return "ERROR: %s must be of the format 'year-month-day'!" \
                       % (key,)
        elif value_types[key] == 'datetime':
            try:
                pardic[key] = str(Time(pardic[key]))
            except (ValueError, TypeError):
                return "ERROR: %s must be of " \
                       "the format 'year-month-day hour:minute:second'!" \
                       % (key,)
        elif value_types[key] == 'timedelta':
            if not isinstance(pardic[key], timedelta):
                try:
                    pardic[key] = timedelta(0, pardic[key])
                except (TypeError, ValueError):
                    return "ERROR: %s must be a datetime.timedelta object " \
                           "or a float(seconds)!" % (key,)
        elif value_types[key] == 'bool':
            if not (pardic[key] == 'false' or pardic[key] == 'true'):
                return "ERROR: %s must be either 'true' or 'false'" % (key,)
        elif not isinstance(pardic[key], value_types[key]):
            return "ERROR: parameter %s must be of %s!" % (
                key, str(value_types[key])[1:-1])
    return None


def _generate_select_sql(values, where_dict, allowed_params, compare_dict,
                         table):
    """
    generate the sql for a select query

    Args:
        values (list): list of names of values
            list of values to return
        where_dict (dict):
            {'param':'value',...} adds WHERE param='value'...
        allowed_params (dict):
            the parameters of the table to be queried and their type
             {'param':type,...}
        compare_dict (dict): default is {}
            {'param': 'inequality',...} (i.e. '>', '<', '>=', '<=', '<>', '!='))
            if no inequality is provided, '=' is assumed
        table (str):
            name of the table

    Returns:
        str (sql select query)

        "ERROR..." if they type_check fails
    """
    for value in reversed(values):
        if value not in allowed_params.keys():
            return "ERROR: %s is an invalid column name!" % (value,)
    where_keys = list(where_dict.keys())
    for param in reversed(where_keys):
        if param not in allowed_params:
            return "ERROR: requested condition on nonexistent column '%s'!" % (
                param,)
    type_check = _data_type_check(where_keys, where_dict, allowed_params)
    if type_check:
        return type_check
    if not values:
        return "ERROR: no valid values requested!"

    values_str = ''
    for value in values:
        values_str += value + ', '
    values_str = values_str[:-2]

    sql = "SELECT %s FROM %s" % (values_str, table)
    if where_keys:
        sql += " WHERE"
        for key in where_keys:
            sql += " %s %s '%s' AND" % (key, compare_dict.get(key, '='),
                                        where_dict[key])
        sql = sql[:-4] + ";"
    return sql


def _generate_insert_sql(pardic, param_list, table):
    """
    generate the sql for an insert command

    Args:
        pardic (dict): (same as given to calling function)
        param_list (list): list of names of parameters to insert
        table (str): name of table

    Returns:
        sql string
    """
    columns = "("
    values = "("
    for param in param_list:
        if not pardic[param] is None:  # make sure that there is a value
            columns += (param + ", ")
            values += "'%s', " % (pardic[param],)
    columns = columns[:-2] + ')'
    values = values[:-2] + ')'
    sql = "INSERT INTO %s %s VALUES %s;" % (table, columns, values)
    return sql


def _generate_update_sql(pardic, param_list, table, lastmodified=False):
    """
    generate the sql for an update command

    Args:
        pardic (dict): (same as given to upper function) (must contain 'id')
        param_list (list): list of names of parameters to update
        table (str): name of table
        lastmodified (bool): if the table has a lastmodified column

    Returns:
        sql string
    """
    sql = "UPDATE %s SET " % (table,)
    for param in param_list:
        if not pardic[param] is None:  # it may be a key with nothing in it
            sql += " %s = '%s'," % (param, pardic[param])
    if lastmodified:
        sql += " lastmodified = 'NOW()' "
    else:
        sql = sql[:-1]
    sql += " WHERE id = %s;" % (pardic['id'],)
    return sql


def _id_from_time():
    """Generate an id from the current time of format YYYYMMDDHHMMSSsss"""
    time = Time.now()
    tid = time.iso
    tid = tid.replace(
        '-', '').replace(' ', '').replace(':', '').replace('.', '')
    return long(tid)


if __name__ == "__main__":

    sedmdb = SedmDB()
    print(sedmdb.expire_requests(send_alerts=True))
    print(sedmdb.update_all_allocations())

    # How to add a request by hand (better to use web form)
    #
    # obsdict = {
    #    'object_id': 20190811030012746,
    #    'user_id': 189,
    #    'allocation_id': 20180131224646741,
    #    'exptime': '{2400}',
    #    'priority': 3.05,
    #    'inidate': '2019-08-09',
    #    'enddate': '2019-08-12',
    #    'obs_seq': '{1ifu}'
    # }
    # print(x.add_request(obsdict))
