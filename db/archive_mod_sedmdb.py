import sqlalchemy.pool as pool
from sqlalchemy import exc, create_engine
import psycopg2
import numpy as np
import subprocess
import warnings
from astropy.time import Time
from datetime import timedelta
from werkzeug.security import generate_password_hash


# Singleton/SingletonPattern.py

class SedmDB:
    class __SedmDB:
        def __init__(self, dbname, host):
            """
            Creates the instance of db connections.
            Needs the username as a parameter.
            The password for SedmDB must be stored in ~/.pgpass
            """
            cmd = "cat ~/.pgpass | grep sedmdbtest | grep -v '#' | awk -F ':' '{print $4}' | head -n 1"

            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            self.user_sedmdb = p.stdout.read().replace('\n', '')
            self.dbname = dbname
            self.host = host

            self.pool_sedmdb = pool.QueuePool(self.__getSedmDBConn__, max_overflow=10, pool_size=2, recycle=True)

        def __str__(self):
            return repr(self)

        def __getSedmDBConn__(self):
            """
            Creates the connection to SedmDB.
            """
            sedmdbcon = psycopg2.connect(host=self.host, port="5432", dbname=self.dbname,
                                         user=self.user_sedmdb)
            return sedmdbcon

    instance = None

    def __init__(self, dbname='sedmdbtest', host='localhost'):
        """
        Makes sure only one instance is created.
        """
        if not SedmDB.instance:
            SedmDB.instance = SedmDB.__SedmDB(dbname, host)

        self.sso_objects = None

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def execute_sql(self, sql):
        """
        Runs the SedmDB sql query in a safe way through the DBManager.

        Returns the object with the results.
        """
        conn = self.pool_sedmdb.connect()
        cursor = conn.cursor()
        try:
            cursor.execute(sql)
        except exc.DBAPIError, e:
            # an exception is raised, Connection is invalidated.
            if e.connection_invalidated:
                print "Connection was invalidated!"
        if 'SELECT' in sql[:8]:
            print cursor
            obj = cursor.fetchall()
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
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())
        if 'username' not in keys:
            return (-1, "ERROR: no username provided!")
        # check for duplicate username
        usernames = [user[0] for user in self.execute_sql('SELECT username FROM users')]
        if pardic['username'] in usernames:
            return (-1, "ERROR: user with that username exists!")
        if 'password' in keys:
            pardic['password'] = generate_password_hash(pardic['password'])
        for key in reversed(keys):  # remove group keys and any other bad keys
            if key not in ['id', 'username', 'name', 'email', 'password']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        sql = _generate_insert_sql(pardic, keys, 'users')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_user sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_user sql command failed with a ProgrammingError!")
        return (id, "User added")

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

            (id (long), "User updated, columns 'column_names'") if the user is updated successfully
        """
        param_types = {'id': int, 'name': str, 'email': str, 'password': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: id not provided!")

        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM users;')]:
            return (-1, "ERROR: no user with the id!")
        if 'password' in keys:
            pardic['password'] = generate_password_hash(pardic['password'])

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['name', 'email', 'password']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'users')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_user sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_user sql command failed with a ProgrammingError!")
        return (pardic['id'], "User updated, columns " + str(keys)[1:-1])

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
            (-1, "ERROR...") if there was an issue (user doesn't exist, not enough information in pardic)

            (0, "User removed") if the removal was successful
        """
        if 'username' in pardic.keys():
            user_id = self.get_from_users(['id'], {'username': pardic['username']})
            if user_id:
                if user_id[0] == -1:  # if get_from_users failed
                    return user_id
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'" % (user_id[0][0],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (user_id[0][0],))
                return (0, "User removed")
            else:
                return (-1, "ERROR: no user with that username!")
        elif 'id' in pardic.keys():
            if not (isinstance(pardic['id'], int) or isinstance(pardic['id'], long)):
                return (-1, "ERROR: id must be of type 'int'")
            if pardic['id'] in [x[0] for x in self.execute_sql('SELECT id FROM users;')]:
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'" % (pardic['id'],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (pardic['id'],))
                return (0, "User removed")
            else:
                return (-1, "ERROR: no user with that id!")
        else:
            return (-1, "ERROR: username or id required!")

    def get_from_users(self, values, where_dict={}, compare_dict={}):
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
            list of tuples containing the values for each user matching the criteria
            
            empty list if no users match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'username': str, 'name': str, 'email': str, 'password': str}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'users')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)
        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in designator

        Args:
            pardic (dict):
                required:
                    'designator' (str)

        Returns:
            (-1, "ERROR...") if no designator was provided or there is already a group with it

            (id (long), "Group added") if the adding was successful
        """
        id = _id_from_time()
        pardic['id'] = id
        if 'designator' not in pardic.keys():
            return (-1, 'ERROR: no group designator provided!')
        groups = [des[0] for des in self.execute_sql('SELECT designator FROM groups;')]
        if pardic['designator'] not in groups:
            sql = ("INSERT INTO groups (id, designator) VALUES ('%s', '%s')" % (pardic['id'], pardic['designator']))
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_group sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_group sql command failed with a ProgrammingError!")
            return (id, "Group added")
        else:
            return (-1, "ERROR: group exists!")

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
        if user not in [user_id[0] for user_id in self.execute_sql('SELECT id FROM users')]:
            return (-1, "ERROR: user does not exist!")
        if group not in [group_id[0] for group_id in self.execute_sql('SELECT id FROM groups')]:
            return (-1, "ERROR: group does not exist!")
        usergroups = self.execute_sql('SELECT user_id, group_id FROM usergroups')
        if (user, group) in usergroups:
            return (-1, "ERROR: user already in group!")
        else:
            sql = "INSERT INTO usergroups (user_id, group_id) VALUES ('%s','%s')" % (user, group)
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_to_group sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_to_group sql command failed with a ProgrammingError!")
            return (user, "User added to group")

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
        if user not in [user_id[0] for user_id in self.execute_sql('SELECT id FROM users')]:
            return (-1, "ERROR: user does not exist!")
        if group not in [group_id[0] for group_id in self.execute_sql('SELECT id FROM groups')]:
            return (-1, "ERROR: group does not exist!")
        usergroups = self.execute_sql('SELECT user_id, group_id FROM usergroups')
        if (user, group) not in usergroups:
            return (-1, "ERROR: user not in group!")
        else:
            sql = "DELETE FROM usergroups WHERE user_id='%s' AND group_id='%s'" % (user, group)
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return (-1, "ERROR: remove_from_group sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: remove_from_group sql command failed with a ProgrammingError!")
            return (0, "User removed from group")

    def get_from_usergroups(self, values, where_dict={}, compare_dict={}):
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
            list of tuples containing the values for each usergroup matching the criteria

            empty list if no usergroups match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'user_id': int, 'group_id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'usergroups')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: get_from_usergroups sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: get_from_usergroups sql command failed with a ProgrammingError!")
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
                    'PI' (str)
                optional:
                    'time_allocated' (datetime.timedelta object or float/int seconds)
                    'priority' (float)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)

        Returns:
            (-1, "ERROR...") if it failed to add

            (id (long), "Program added") if the program is added successfully
        """
        param_types = {'id': int, 'name': str, 'designator': str, 'group_id': int, 'PI': str, 'color': str,
                       'time_allocated': 'timedelta', 'priority': float, 'inidate': 'datetime', 'enddate': 'datetime'}
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())

        if 'designator' in keys:
            if pardic['designator'] in [obj[0] for obj in self.execute_sql('SELECT designator FROM program')]:
                return (-1, "ERROR: a program with that designator already exists!")

        for key in ['designator', 'name', 'group_id', 'PI']:  # check for required keys
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(keys):  # remove any extraneous keys
            if key not in ['id', 'name', 'designator', 'group_id', 'PI', 'color'
                           'time_allocated', 'priority', 'inidate', 'enddate']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        program_sql = _generate_insert_sql(pardic, keys, 'program')
        try:
            self.execute_sql(program_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_program sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_program sql command failed with a ProgrammingError!")
        return (id, "Program added")

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
                    'time_allocated' (datetime.timedelta object or float/int seconds)
                    'priority' (float)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Program updated, columns 'column_names'") if the program is updated successfully
        """
        param_types = {'id': int, 'time_allocated': 'timedelta', 'name': str, 'PI': str, 'priority': float,
                       'inidate': 'datetime', 'enddate': 'datetime', 'color': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: id not provided!")

        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM program;')]:
            return (-1, "ERROR: no program with the id!")

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['time_allocated', 'name', 'PI', 'priority', 'inidate', 'enddate', 'color']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'program')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_program sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_program sql command failed with a ProgrammingError!")
        return (pardic['id'], "Program updated, columns " + str(keys)[1:-1])

    def get_from_program(self, values, where_dict={}, compare_dict={}):
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
            values/keys options: ['id', 'designator', 'name', 'group_id', 'PI', 'time_allocated',
                                  'priority', 'inidate', 'enddate', 'color']

        Returns:
            list of tuples containing the values for each program matching the criteria

            empty list if no programs match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'name': str, 'designator': str, 'group_id': int, 'PI': str, 'color': str,
                          'time_allocated': 'timedelta', 'priority': float,
                          'inidate': 'datetime', 'enddate': 'datetime'}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'program')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: get_from_program sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: get_from_program sql command failed with a ProgrammingError!")
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
                    'time_allocated' (datetime.timedelta object or float/int seconds)
                    'time_spent' (datetime.timedelta object or float/int seconds)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)

        Returns:
            (-1, "ERROR...") if it failed to add

            (id (long), "Allocation added") if the program is added successfully
        """
        param_types = {'id': int, 'program_id': int, 'designator': str, 'time_allocated': 'timedelta', 'time_spent': 'timedelta',
                       'inidate': 'datetime', 'enddate': 'datetime', 'color': str, 'active': bool}
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())

        if 'program_id' in keys:
            if pardic['program_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM program;')]:
                return (-1, "ERROR: no program with that id exists!")

        for key in ['program_id', 'designator']:
            if key not in keys:  # check for required key
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(keys):  # remove any extraneous keys
            if key not in ['id', 'program_id', 'designator', 'time_allocated', 'time_spent', 'inidate', 'enddate', 'color', 'active']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        program_sql = _generate_insert_sql(pardic, keys, 'allocation')
        try:
            self.execute_sql(program_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_program sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_program sql command failed with a ProgrammingError!")
        return (id, "Allocation added")

    def update_allocation(self, pardic):
        """
        updates an allocation entry

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'time_allocated' (datetime.timedelta object or float/int seconds)
                    'time_spent' (datetime.timedelta object or float/int seconds)
                    'inidate' ('year-month-day hour:minute:second')
                    'enddate' ('year-month-day hour:minute:second')
                    'color' (hex color code)
                    'active' (boolean)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Allocation updated, columns 'column_names'") if the allocation is updated successfully
        """
        param_types = {'id': int, 'time_allocated': 'timedelta', 'time_spent': 'timedelta',
                       'inidate': 'datetime', 'enddate': 'datetime', 'color': str, 'active': bool}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: id not provided!")

        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM allocation;')]:
            return (-1, "ERROR: no allocation with the id!")

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['time_allocated', 'time_spent', 'inidate', 'enddate', 'color', 'active']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'allocation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_allocation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_allocation sql command failed with a ProgrammingError!")
        return (pardic['id'], "Allocation updated, columns " + str(keys)[1:-1])

    def get_from_allocation(self, values, where_dict={}, compare_dict={}):
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
            values/keys options: ['id', 'designator', 'program_id, 'time_spent', 'time_allocated', 'inidate', 'enddate', 'color', 'active']

        Returns:
            list of tuples containing the values for each program matching the criteria

            empty list if no programs match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'program_id': int, 'designator': str, 'time_spent': 'timedelta', 'color': str,
                          'time_allocated': 'timedelta', 'inidate': 'datetime', 'enddate': 'datetime', 'active': bool}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'allocation')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: get_from_allocation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: get_from_allocation sql command failed with a ProgrammingError!")
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
                    'f' (fixed), 'v' (periodic fixed), 'P' (built-in planet or satellite name), 'e' (heliocentric elliptical),
                    'h' (heliocentric hyperbolic), 'p' (heliocentric parabolic), 'E' (geocentric elliptical)

        Returns:
            (-1, "ERROR...") if it failed to add

            (id (long), "Object added") if the object is added successfully
        """
        param_types = {'id': int, 'name': str, 'typedesig': str, 'ra': float, 'dec': float, 'epoch': str,
                       'iauname': str, 'marshal_id': int, 'magnitude': float}
        id = _id_from_time()
        pardic['id'] = id
        obj_keys = list(pardic.keys())
        if 'marshal_id' in obj_keys:
            if pardic['marshal_id'] in [obj[0] for obj in self.execute_sql('SELECT marshal_id FROM object')]:
                return (-1, "ERROR: object exists!")

        for key in ['name', 'typedesig']:  # check if 'name' and 'typedesig' are provided
            if key not in obj_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(obj_keys):  # remove any extraneous keys
            if key not in ['id', 'name', 'typedesig', 'ra', 'dec', 'epoch', 'marshal_id', 'iauname', 'magnitude']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(obj_keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        pardic['name'] = pardic['name'].lower()  # make all of the names the same format for consistant searching
        if pardic['typedesig'] == 'f':
            for key in ['ra', 'dec']:
                if key not in obj_keys:
                    return (-1, "ERROR: %s not provided!" % (key,))
            dup = self.execute_sql("SELECT id, name FROM object WHERE q3c_radial_query(ra, dec, '%s', '%s', .000278)"
                                   % (pardic['ra'], pardic['dec']))
            if dup:  # if there is already an object within an arcsecond
                if pardic['name'] in [x[1] for x in dup]:  # check for same coords, same name
                    lis = [x[1] for x in dup]
                    idx = lis.index(pardic['name'])
                    return(-1, "ERROR: The object '%s' is already in the database with id %s"
                               % (pardic['name'], dup[idx][0]))
                print("there is already an object within 1 arcsec of given coordinates with "
                      "id: %s, name: %s" % (dup[0][0], dup[0][1]))

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")
            return (id, "Fixed object added")
        elif pardic['typedesig'] == 'v':
            for key in ['ra', 'dec']:
                if key not in obj_keys:
                    return (-1, "ERROR: %s not provided!" % (key,))
            dup = self.execute_sql("SELECT id, name FROM object WHERE q3c_radial_query(ra, dec, '%s', '%s', .000278)"
                                   % (pardic['ra'], pardic['dec']))
            if dup:  # if there is already an object within an arcsecond
                if pardic['name'] in [x[1] for x in dup]:  # check for same coords, same name
                    lis = [x[1] for x in dup]
                    idx = lis.index(pardic['name'])
                    return(-1, "ERROR: The object '%s' is already in the database with id %s"
                               % (pardic['name'], dup[idx][0])) 
                print("there is already an object within 1 arcsec of given coordinates with "
                      "id: %s, name: %s" % (dup[0][0], dup[0][1]))

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")
            return (id, "Fixed object added, input period data into table `periodic`")

        elif pardic['typedesig'] in ['h', 'E', 'e', 'p']:
            function_dict = {'e': 'add_elliptical_heliocentric', 'h': 'add_hyperbolic_heliocentric',
                             'p': 'add_parabolic_heliocentric', 'E': 'add_earth_satellite'}

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")
            return (id, "Non-fixed object added, orbit parameters can be added with `%s`"
                    % (function_dict[pardic['typedesig']],))
        elif pardic['typedesig'] == 'P':
            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")
            return (id, "Non-fixed object added")
        else:
            return (-1, "ERROR: typedesig provided was invalid, it must be one of ['f', 'h', 'E', 'e', 'p', 'P']!")

    # TODO: write update_object

    def get_from_object(self, values, where_dict={}, compare_dict={}):
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
            values/keys options: ['id', 'marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch']

        Returns:
            list of tuples containing the values for each user matching the criteria

            empty list if no objects match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'marshal_id': int, 'name': str, 'iauname': str, 'ra': float, 'dec': float,
                          'typedesig': str, 'epoch': float}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'object')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def get_objects_near(self, ra, dec, radius, values=['id', 'name', 'epoch']):
        """
        get object entries with coordinates within a radius

        Args:
            ra (float or int): ra in degrees of the object
            dec (float or int): dec in degrees of the object
            radius (float or int): radius in arcseconds
            values (list): which parameters to return
                    options: 'id', 'marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch'
                    defaults to ['id', 'name', 'epoch']

        Returns:
            list of tuples [(id, name, epoch)] for each object in the radius

            [] if there are no objects found

            (-1, "ERROR...") if there is an issue
        """
        for value in values[::-1]:
            if value not in ['id', 'marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch']:
                return (-1, "ERROR: %s is an invalid column name!" % (value,))
        if not (isinstance(ra, float) or isinstance(ra, int)):
            return (-1, "ERROR: parameter ra must be of type 'float' or type 'int'!")
        if not (isinstance(dec, float) or isinstance(dec, int)):
            return (-1, "ERROR: parameter dec must be of type 'float' or type 'int'!")
        if not (isinstance(radius, float) or isinstance(radius, int)):
            return (-1, "ERROR: parameter radius must be of type 'float' or type 'int'!")

        objects = self.execute_sql("SELECT id, name, epoch FROM object WHERE "
                                   "q3c_radial_query(ra, dec, '%s', '%s', '%s')" % (ra, dec, .000278*radius))
        return objects

    def get_object_id_from_name(self, object_name):
        """
        finds the id of an object given its name or part of its name

        Args:
            object_name (str): part of the name of an object

        Returns:
            list of tuples [(id, full name)] for each object with a matching name

            [] if no matching object name is found
        """
        object_name = object_name.lower()
        sql = "SELECT id, name FROM object WHERE name LIKE '%s%s%s'" % ('%', object_name, '%')
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
                    'mjdepoch' (int/long) (epoch, time of 'M'),
                    'D' (int/long) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components of magnitude model)
                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'inclination': float, 'longascnode_O': float, 'perihelion_o': float,
                       'a': float, 'n': float, 'e': float, 'M': float, 'mjdepoch': int, 'D': int, 'M1': float,
                       'M2': float, 's': float}
        id = _id_from_time()
        orbit_params['id'] = id
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['inclination', 'longascnode_O', 'perihelion_o', 'a', 'n', 'e',
                    'M', 'mjdepoch', 'D', 'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['id', 'inclination', 'longascnode_O', 'perihelion_o', 'a', 'n', 'e',
                           'M', 'mjdepoch', 'D', 'M1', 'M2', 's', 'object_id']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'elliptical_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_elliptical_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_elliptical_orbit sql command failed with a ProgrammingError!")
        return (id, "Elliptical heliocentric orbit added")

    def get_from_elliptical_heliocentric(self, values, where_dict={}, compare_dict={}):
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
                'mjdepoch' (int) (epoch, time of 'M'),
                'D' (int) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'object_id': int, 'inclination': float, 'longascnode_O': float, 'perihelion_o': float,
                          'a': float, 'n': float, 'e': float, 'M': float, 'mjdepoch': int, 'D': int, 'M1': float,
                          'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'elliptical_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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
                    'D' (int) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components of magnitude model)

                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'longascnode_O': float,
                       'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float}
        id = _id_from_time()
        orbit_params['id'] = id
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'e', 'q', 'D',
                    'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['id', 'T', 'inclination', 'longascnode_O', 'perihelion_o', 'e', 'q', 'D',
                           'M1', 'M2', 's', 'object_id']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'hyperbolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_hyperbolic_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_hyperbolic_orbit sql command failed with a ProgrammingError!")
        return (id, "Hyperbolic heliocentric orbit added")

    def get_from_hyperbolic_heliocentric(self, values, where_dict={}, compare_dict={}):
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
                'D' (int) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'longascnode_O': float,
                          'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'hyperbolic_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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
                    'D' (int) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components of magnitude model)

                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'id': int, 'object_id': int, 'T': 'date', 'inclination': float, 'longascnode_O': float,
                       'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float}
        id = _id_from_time()
        orbit_params['id'] = id
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'q', 'D',
                    'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['id', 'T', 'inclination', 'longascnode_O', 'perihelion_o', 'q', 'D',
                           'M1', 'M2', 's', 'object_id']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'parabolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_parabolic_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_parabolic_orbit sql command failed with a ProgrammingError!")
        return (id, "Parabolic heliocentric orbit added")

    def get_from_parabolic_heliocentric(self, values, where_dict={}, compare_dict={}):
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
                'D' (int) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'object_id': int, 'T': 'date', 'inclination': float, 'longascnode_O': float,
                          'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'parabolic_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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
        param_types = {'id': int, 'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'ra': float,
                       'pedigree': float, 'M': float, 'n': float, 'decay': float, 'reforbit': int, 'drag': float}
        id = _id_from_time()
        orbit_params['id'] = id
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                    'decay', 'reforbit', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['id', 'T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                           'decay', 'reforbit', 'drag', 'object_id']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'earth_satellite')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_earth_satellite_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_earth_satellite_orbit sql command failed with a ProgrammingError!")
        return (id, "Earth satellite orbit added")

    def get_from_earth_satellite(self, values, where_dict={}, compare_dict={}):
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
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'ra': float,
                          'pedigree': float, 'M': float, 'n': float, 'decay': float, 'reforbit': int,
                          'drag': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'earth_satellite')
        print sql
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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
                    'maxairmass' (float) (max allowable airmass for observation, default 2.5),
                    'cadence' (float) (time between periods),
                    'phasesamples' (float) (how many samples in a period),
                    'sampletolerance' (float) (how much tolerance in when the samples can be taken),
                    'nexposures' (str '{# of ifu, # of u, # of g, # of r, # of i}'),
                    'obs_seq' (str e.g. '{3g, 3r, 1i, 1ifu, 2i}' for 3 of g, then 3 of r, then 1 i, 1 ifu, 1 i),
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
        param_types = {'id': int, 'object_id': int, 'user_id': int, 'allocation_id': int, 'exptime': str, 'priority': float,
                       'inidate': 'date', 'enddate': 'date', 'marshal_id': int, 'maxairmass': float, 'cadence': float,
                       'phasesamples': float, 'sampletolerance': float, 'nexposures': str, 'obs_seq': str,
                       'max_fwhm': float, 'min_moon_dist': float, 'max_moon_illum': float, 'max_cloud_cover': float,
                       'seq_repeats': int}
        id = _id_from_time()
        pardic['id'] = id
        requests = self.execute_sql("SELECT object_id, allocation_id FROM request WHERE status != 'EXPIRED';")
        # check allocation_id, issue warning if it is a repeat, but allow
        if (pardic['object_id'], pardic['allocation_id']) in requests:
            print "program %s has already requested object: %s" % (pardic['allocation_id'], pardic['object_id'])
            # TODO: add to output
        if pardic['object_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM object;')]:
            return (-1, "ERROR: object does not exist!")
        if pardic['user_id'] not in [user[0] for user in self.execute_sql('SELECT id FROM users;')]:
            return (-1, "ERROR: user does not exist!")
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
            nexposure = [n*pardic['seq_repeats'] for n in nexposure]
            # make sure that 'nexposures' and 'obs_seq' are consistent
            if 'nexposures' in pardic.keys():
                if not '{%s, %s, %s, %s, %s}' % tuple(nexposure) == pardic['nexposures']:
                    return (-1, "ERROR: nexposures and obs_seq are inconsistent!")
            else:
                pardic['nexposures'] = '{%s, %s, %s, %s, %s}' % tuple(nexposure)

        elif not ('nexposures' in pardic.keys() or 'obs_seq' in pardic.keys()):
            return (-1, "ERROR: nexposures or obs_seq is required!")

        keys = list(pardic.keys())
        default_params = ['object_id', 'user_id', 'allocation_id', 'exptime', 'priority',
                          'inidate', 'enddate']
        for param in default_params:  # check that all required values are provided
            if param not in keys:
                return (-1, "ERROR: %s not in dictionary!" % (param,))
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['id', 'object_id', 'user_id', 'allocation_id', 'exptime', 'priority',
                           'inidate', 'enddate', 'marshal_id', 'maxairmass', 'cadence',
                           'phasesamples', 'sampletolerance', 'nexposures', 'obs_seq', 'seq_repeats',
                           'max_fwhm', 'min_moon_dist', 'max_moon_illum', 'max_cloud_cover']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'request')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_request sql command failed with a ProgrammingError!")
        return (id, "Request added")

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
                Note: 'status' can be 'PENDING', 'ACTIVE', 'COMPLETED', 'REDUCED', 'CANCELED', or 'EXPIRED'

        Returns:
            (-1, "ERROR...") if there was an issue with the updating

            (id, "Requests updated, columns 'column_names'") if the update was successful
        """
        param_types = {'id': int, 'status': str, 'maxairmass': float, 'priority': float,
                       'inidate': 'date', 'enddate': 'date', 'seq_completed': int, 'last_obs_jd': float,
                       'max_fwhm': float, 'min_moon_dist': float, 'max_moon_illum': float, 'max_cloud_cover': float}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: no id provided!")
        elif not (isinstance(pardic['id'], int) or isinstance(pardic['id'], long)):
            return (-1, "ERROR: parameter id must be of type 'int'!")
        if pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM request;')]:
            return (-1, "ERROR: request does not exist!")
        if 'status' in keys:
            if pardic['status'] not in ['PENDING', 'ACTIVE', 'COMPLETED', 'CANCELED', 'EXPIRED']:
                return (-1, "ERROR: %s is an invalid status value!" % (pardic['status'],))
        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['status', 'maxairmass', 'priority', 'inidate', 'enddate', 'seq_completed',
                           'max_fwhm', 'min_moon_dist', 'max_moon_illum', 'max_cloud_cover', 'last_obs_jd']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'request', True)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_request sql command failed with a ProgrammingError!")

        return (pardic['id'], "Requests updated, columns " + str(keys)[1:-1])

    def get_from_request(self, values, where_dict={}, compare_dict={}):
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
                'obs_seq' (str)
                'status' (str),
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
            list of tuples containing the values for each request matching the criteria,

            empty list if no requests match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'object_id': int, 'user_id': int, 'allocation_id': int, 'exptime': str, 'status': str,
                          'priority': float, 'inidate': 'date', 'enddate': 'date', 'marshal_id': int,
                          'maxairmass': float, 'cadence': float, 'phasesamples': float, 'sampletolerance': float,
                          'filters': str, 'nexposures': str, 'obs_seq': str, 'creationdate': 'date',
                          'lastmodified': 'date', 'seq_repeats': int, 'seq_completed': int, 'last_obs_jd': float,
                          'max_fwhm': float, 'min_moon_dist': float, 'max_moon_illum': float, 'max_cloud_cover': float}
        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'request')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def expire_requests(self):
        """
        Updates the request table. For all the active requests that were not completed,
            and had an expiry date before than NOW(), are marked as "EXPIRED".

        Returns:
            (0, "Requests expired")
        """
        # TODO: move to logic layer? (requires sql "knowledge")
        # TODO: make it more discerning of other statuses
        # tests written
        sql = "UPDATE request SET status='EXPIRED' WHERE enddate < 'NOW()' AND status != 'COMPLETED';"
        self.execute_sql(sql)
        return (0, "Requests expired")

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
                    'tel_ra' (str),
                    'tel_dec' (str),
                    'tel_az' (float),
                    'tel_el' (float),
                    'tel_pa' (float),
                    'ra_off' (float),
                    'dec_off' (float),
                optional:
                    'imtype' (str),
                    'time_elapsed' (datetime.timedelta object or float/int seconds),
                    'filter' (str),
                        options: 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b'
                    'camera' (str)

        Returns:
            (-1, "ERROR...") if there was an issue

            (id (long), "Observation added") if it completed successfully
        """
        header_types = {'id': int, 'object_id': int, 'request_id': int, 'mjd': float, 'airmass': float,
                        'exptime': float, 'fitsfile': str, 'lst': str, 'ra': float, 'dec': float, 'tel_ra': str,
                        'tel_dec': str, 'tel_az': float, 'tel_el': float, 'tel_pa': float, 'ra_off': float,
                        'dec_off': float, 'imtype': str, 'camera': str, 'filter': str, 'time_elapsed': 'timedelta'}
        id = _id_from_time()
        if 'id' not in header_dict.keys():        
            header_dict['id'] = id

        header_keys = list(header_dict.keys())
        if 'filter' in header_keys:
            if header_dict['filter'] not in ['u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b']:  # check the filter is valid
                return (-1, "ERROR: invalid filter given!")

        for key in ['object_id', 'request_id', 'mjd', 'airmass', 'exptime', 'fitsfile', 'lst',
                    'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az', 'tel_el', 'tel_pa', 'ra_off', 'dec_off']:
            if key not in header_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(header_keys):
            if key not in ['id', 'object_id', 'request_id', 'mjd', 'airmass', 'exptime',
                           'fitsfile', 'imtype', 'lst', 'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az',
                           'tel_el', 'tel_pa', 'ra_off', 'dec_off', 'camera', 'filter', 'time_elapsed']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(header_keys, header_dict, header_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(header_dict, header_keys, 'observation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: adding observation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: adding observation sql command failed with a ProgrammingError!")
        return (id, "Observation added")

    def update_observation(self, pardic):
        """

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'mjd' (float),
                    'airmass' (float),
                    'exptime' (float),
                    'fitsfile' (abspath str),
                    'lst' (str),
                    'ra' (float),
                    'dec' (float),
                    'tel_ra' (str),
                    'tel_dec' (str),
                    'tel_az' (float),
                    'tel_el' (float),
                    'tel_pa' (float),
                    'ra_off' (float),
                    'dec_off' (float),
                    'imtype' (str),
                    'camera' (str),
                    'time_elapsed (datetime.timedelta object or float/int seconds),
                    'filter' (str)
                        options: 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b'

        Returns:
            (-1, "ERROR...") if there was an issue

            (id, "Observation updated, columns 'column_names'") if it completed successfully
        """
        param_types = {'id': int, 'mjd': float, 'airmass': float, 'filter': str,
                       'exptime': float, 'fitsfile': str, 'lst': str, 'ra': float, 'dec': float, 'tel_ra': str,
                       'tel_dec': str, 'tel_az': float, 'tel_el': float, 'tel_pa': float, 'ra_off': float,
                       'dec_off': float, 'imtype': str, 'camera': str, 'time_elapsed': 'timedelta'}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: id not provided!")
        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM observation;')]:
            return (-1, "ERROR: observation does not exist!")

        if 'filter' in keys:
            if pardic['filter'] not in ['u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b']:  # check the filter is valid
                return (-1, "ERROR: invalid filter given!")

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['mjd', 'airmass', 'exptime', 'fitsfile', 'lst', 'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az',
                            'tel_el', 'tel_pa', 'ra_off', 'dec_off', 'imtype', 'camera', 'filter', 'time_elapsed']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'observation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_observation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_observation sql command failed with a ProgrammingError!")
        return (pardic['id'], "Observation updated, columns " + str(keys)[1:-1])

    def get_from_observation(self, values, where_dict={}, compare_dict={}):
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
                'tel_ra' (str),
                'tel_dec' (str),
                'tel_az' (float),
                'tel_el' (float),
                'tel_pa' (float),
                'ra_off' (float),
                'dec_off' (float),
                'imtype' (str),
                'camera' (str),
                'time_elapsed (datetime.timedelta object or float/int seconds),
                'filter' (str)
                    options: 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b'

        Returns:
            list of tuples containing the values for each observation matching the criteria

            empty list if no observations match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'object_id': int, 'request_id': int, 'mjd': float, 'airmass': float,
                          'exptime': float, 'fitsfile': str, 'lst': str, 'ra': float, 'dec': float, 'tel_ra': str,
                          'tel_dec': str, 'tel_az': float, 'tel_el': float, 'tel_pa': float, 'ra_off': float,
                          'dec_off': float, 'imtype': str, 'camera': str, 'id': int, 'filter': str,
                          'time_elapsed': 'timedelta'}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'observation')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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
        telstat_types = {'id': int, 'date': 'date', 'dome_status': str, 'in_temp': float, 'in_humidity': float, 'in_dew': float,
                         'out_temp': float, 'out_humidity': float, 'out_dew': float, 'wind_dir': float,
                         'wsp_cur': float, 'wsp_avg': float, 'mir_temp': float, 'top_air': float, 'pri_temp': float,
                         'sec_temp': float, 'flo_temp': float, 'bot_temp': float, 'mid_temp': float, 'top_temp': float,
                         'observation_id': int}
        stat_keys = list(tel_stats.keys())
        for key in reversed(stat_keys):
            if key not in ['id', 'date', 'dome_status', 'in_temp', 'in_humidity', 'in_dew', 'out_temp', 'out_humidity',
                           'out_dew', 'wind_dir', 'wsp_cur', 'wsp_avg', 'mir_temp', 'top_air', 'pri_temp', 'sec_temp',
                           'flo_temp', 'bot_temp', 'mid_temp', 'top_temp', 'observation_id']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(stat_keys, tel_stats, telstat_types)
        if type_check:
            return (-1, type_check)

        stat_sql = _generate_insert_sql(tel_stats, stat_keys, 'telescope_stats')
        try:
            self.execute_sql(stat_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: adding tel_stats sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: adding tel_stats sql command failed with a ProgrammingError!")

        return (id, "Telescope stats added")

    # TODO: write update_observation() and update_telescope_stats()

    def get_from_telescope_stats(self, values, where_dict={}, compare_dict={}):
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
            list of tuples containing the values for telescope stats matching the criteria

            empty list if no telescope_stats entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'date': 'date', 'dome_status': str, 'in_temp': float, 'in_humidity': float, 'in_dew': float,
                          'out_temp': float, 'out_humidity': float, 'out_dew': float, 'wind_dir': float,
                          'wsp_cur': float, 'wsp_avg': float, 'mir_temp': float, 'top_air': float, 'pri_temp': float,
                          'sec_temp': float, 'flo_temp': float, 'bot_temp': float, 'mid_temp': float, 'top_temp': float,
                          'observation_id': int, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'telescope_stats')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    # TODO: separate add_phot/spec and update_phot/spec
    def add_phot(self, pardic):
        """
        Adds reduced photometry or, if the observation already has photometry, updates it

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

            (id (long), "Photometry added") if the photometry was added successfully

            (id (long), "Photometry updated for observation_id ...") if the photometry existed and was updated
        """
        param_types = {'id': int, 'observation_id': int, 'astrometry': 'bool', 'filter': str, 'reducedfile': str, 'sexfile': str,
                       'maskfile': str, 'pipeline': str, 'marshal_phot_id': int, 'phot_calib_id': int}
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return (-1, "ERROR: observation_id not provided!")
        phot_id = self.get_from_phot(['id'], {'observation_id': pardic['observation_id']})
        if phot_id:  # if there is already an entry for that observation, update instead
            if phot_id[0] == -1:
                return phot_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['astrometry', 'filter', 'reducedfile', 'sexfile',
                               'maskfile', 'pipeline', 'marshal_phot_id', 'phot_calib_id']:
                    return (-1, "ERROR: %s is an invalid key!" % (key,))
            pardic['id'] = phot_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return (-1, type_check)

            update_sql = _generate_insert_sql(pardic, keys, 'spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_reduced_photometry update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_reduced_photometry update sql command failed with a ProgrammingError!")
            return (phot_id[0][0], "Photometry updated for observation_id %s, columns " % (pardic['observation_id'],)
                                                                                + str(keys)[1:-1])

        obs = self.get_from_observation(['fitsfile'], {'id': pardic['observation_id']})
        if not obs:
            return (-1, "ERROR: no observation with the observation_id")
        elif obs[0] == -1:
            return obs
        else:
            pass
            # TODO: generate the filter here?

        phot_calib = self.get_from_phot_calib(['id'], {'id': pardic['phot_calib_id']})
        if not phot_calib:
            return (-1, "ERROR: no phot_calib with the phot_calib_id")
        elif phot_calib[0] == -1:
            return phot_calib

        for key in ['observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile',
                    'maskfile', 'pipeline', 'phot_calib_id']:  # include 'marshal_phot_id'?
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):  # remove any invalid keys
            if key not in ['id', 'observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile',
                           'maskfile', 'pipeline', 'marshal_phot_id', 'phot_calib_id']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with a ProgrammingError!")

    def get_from_phot(self, values, where_dict={}, compare_dict={}):
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
            list of tuples containing the values for phot entries matching the criteria

            empty list if no phot entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'observation_id': int, 'astrometry': 'bool', 'filter': str, 'reducedfile': str, 'sexfile': str,
                          'maskfile': str, 'pipeline': str, 'marshal_phot_id': int, 'phot_calib_id': int, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'phot')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_spec(self, pardic):
        """
        Adds the reduced spectrum or, if the observation already has a spectrum, updates it

        Args:
            pardic (dict):
                required:
                    'observation_id' (int/long),
                    'fitsfile' (abspath str),
                    'npyfile' (abspath str),
                    'asciifile' (abspath str),
                    'imgset' (str),
                    'quality' (int),
                    'cubefile' (abspath str),
                    'standardfile' (abspath str),
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
                    'reducer' (float),
                    'fwhm' (float),
                    'background' (float),
                    'line_fwhm' (float)

        Returns:
            (-1, "ERROR...") if there was an issue

            (id (long), "Spectrum added")  if the spectrum was added successfully

            (id (long), "Spectrum updated for observation_id ...") if the spectrum existed and was updated
        """
        param_types = {'id': int, 'observation_id': int,
                       'imgset': str, 'quality': int, 'cubefile': str, 'standardfile': str, 'skysub': 'bool',
                       'extract_x': float, 'extract_y': float, 'extract_pa': float, 'extract_a': float,
                       'extract_b': float, 'ad_red': float, 'ad_blue': float, 'prlltc': float,
                       'flexure_x_corr_nm': float, 'flexure_y_corr_pix': float, 'reducer': float, 'fwhm': float,
                       'background': float, 'line_fwhm': float, 'fitsfile': str, 'npyfile': str, 'asciifile': str}
        id = _id_from_time()
        pardic['id'] = id
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return (-1, "ERROR: observation_id not provided!")
        spec_id = self.get_from_spec(['id'], {'observation_id': pardic['observation_id']})
        if spec_id:  # if there is already an entry for that observation, update instead
            if spec_id[0] == -1:
                return spec_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fitsfile', 'npyfile', 'asciifile', 'imgset', 'quality',
                               'cubefile', 'standardfile', 'marshal_spec_id', 'skysub', 'extract_x', 'extract_y',
                               'extract_pa', 'extract_a', 'extract_b', 'ad_red', 'ad_blue', 'prlltc',
                               'flexure_x_corr_nm', 'flexure_y_corr_pix', 'reducer', 'fwhm', 'background', 'line_fwhm']:
                    return (-1, "ERROR: %s is an invalid key!" % (key,))
            pardic['id'] = spec_id[0][0]
            update_sql = _generate_update_sql(pardic, keys, 'spec')
            print update_sql
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_reduced_spectrum update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_reduced_spectrum update sql command failed with a ProgrammingError!")
            return (spec_id[0][0], "Spectrum updated for observation_id %s, columns " % (pardic['observation_id'],)
                                                                                + str(keys)[1:-1])
        obs_id = self.get_from_observation(['id'], {'id': pardic['observation_id']})
        if not obs_id:
            return (-1, "ERROR: no observation exists with the given id!")
        elif obs_id[0] == -1:
            return obs_id

        for key in ['observation_id', 'fitsfile', 'imgset', 'quality', 'cubefile',
                    'standardfile', 'skysub']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['id', 'observation_id', 'fitsfile', 'npyfile', 'asciifile', 'imgset', 'quality',
                           'cubefile', 'standardfile', 'marshal_spec_id', 'skysub', 'extract_x', 'extract_y',
                           'extract_pa', 'extract_a', 'extract_b', 'ad_red', 'ad_blue', 'prlltc',
                           'flexure_x_corr_nm', 'flexure_y_corr_pix', 'reducer', 'fwhm', 'background', 'line_fwhm']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'spec')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with a ProgrammingError!")

    def get_from_spec(self, values, where_dict={}, compare_dict={}):
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
                'fitsfile' (abspath str),
                'npyfile' (abspath str),
                'asciifile' (abspath str),
                'imgset' (str),
                'quality' (int),
                'cubefile' (abspath str),
                'standardfile' (abspath str),
                'skysub' ('true' or 'false')
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
                'reducer' (float),
                'fwhm' (float),
                'background' (float),
                'line_fwhm' (float)

        Returns:
            list of tuples containing the values for spectra matching the criteria

            empty list if no spec entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'observation_id': int,  'fitsfile': str, 'npyfile': str, 'asciifile': str,
                          'imgset': str, 'quality': int, 'cubefile': str, 'standardfile': str,
                          'skysub': 'bool', 'id': int,
                          'extract_x': float, 'extract_y': float, 'extract_pa': float, 'extract_a': float,
                          'extract_b': float, 'ad_red': float, 'ad_blue': float, 'prlltc': float,
                          'flexure_x_corr_nm': float, 'flexure_y_corr_pix': float, 'reducer': float, 'fwhm': float,
                          'background': float, 'line_fwhm': float}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'spec')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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

            (id (long), "Photometry metrics updated for phot_id ...") if it updated existing metrics

            (id (long), "Photometry metrics added") if the metrics were added successfully
        """
        param_types = {'id': int, 'phot_id': int, 'fwhm': float, 'background': float, 'zp': float,
                       'zperr': float, 'ellipticity': float, 'nsources': int}
        id = _id_from_time()
        pardic['id'] = id
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'phot_id' not in keys:
            return (-1, "ERROR: phot_id not provided!")
        metric_id = self.get_from_metrics_phot(['id'], {'phot_id': pardic['phot_id']})
        if metric_id:  # if there is already an entry for that observation, update instead
            if metric_id[0] == -1:
                return metric_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:
                    return (-1, "ERROR: %s is an invalid key!" % (key,))
            pardic['id'] = metric_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return (-1, type_check)

            update_sql = _generate_insert_sql(pardic, keys, 'metrics_phot')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_metrics_phot update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_metrics_phot update sql command failed with a ProgrammingError!")
            return (metric_id[0][0], "Photometry metrics updated for phot_id %s, columns " % (pardic['phot_id'],)
                                                                                    + str(keys)[1:-1])
        ph_id = self.get_from_phot(['id'], {'id': pardic['phot_id']})
        if not ph_id:
            return (-1, "ERROR: no photometry exists with the given id!")
        elif ph_id[0] == -1:
            return ph_id

        for key in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:  # phot_id already tested
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['id', 'phot_id', 'fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'metrics_phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_metrics_phot sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_metrics_phot sql command failed with a ProgrammingError!")
        return (id, "Photometry metrics added")

    def get_from_metrics_phot(self, values, where_dict):
        """
        select values from `metrics_phot`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
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
            list of tuples containing the values for metrics matching the criteria

            empty list if no metrics_phot entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'phot_id': int, 'fwhm': float, 'background': float, 'zp': float,
                          'zperr': float, 'ellipticity': float, 'nsources': int, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'metrics_phot')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_phot_calib(self, pardic):
        """
        Creates a new object in the phot calib table with the parameters specified in the dictionary.

        Args:
            pardic:
                required:
                    'bias' (abspath str)
                    'flat' (abspath str)

        Returns:
            (-1: "ERROR...") if there is an issue

            (id (long), "Photometry calibration added") if it completes successfully
        """
        param_types = {'id': int, 'bias': str, 'flat': str}
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())

        for key in ['bias', 'flat']:  # phot_id already tested
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['id', 'bias', 'flat']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'phot_calib')
        print sql, type_check
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_phot_calib sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_phot_calib sql command failed with a ProgrammingError!")
        return (id, "Photometry calibration added")

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

            (id (long), "Phot_calib updated, columns 'column_names'") if the entry is updated successfully
        """
        param_types = {'id': int, 'bias': str, 'flat': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: id not provided!")

        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM phot_calib;')]:
            return (-1, "ERROR: no phot_calib entry with the id!")

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['bias', 'flat']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'phot_calib')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_phot_calib sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_phot_calib sql command failed with a ProgrammingError!")
        return (pardic['id'], "Phot_calib updated, columns " + str(keys)[1:-1])

    def get_from_phot_calib(self, values, where_dict={}, compare_dict={}):
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
            list of tuples containing the values for calibration matching the criteria

            empty list if no phot_calib entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'bias': str, 'flat': str}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'phot_calib')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_spec_calib(self, pardic):
        """
        Creates a new object in the spec calib table with the parameters specified in the dictionary.

        Args:
            pardic:
                required:
                    'bias' (abspath str)
                    'flat' (abspath str)
                optional:
                    'dome' (abspath str)
                    'cosmic_filter' (bool)
                    'drpver' (float)
                    'Hg_master' (abspath str)
                    'Cd_master' (abspath str)
                    'Xe_master' (abspath str)
                    'avg_rms' (abspath str)
                    'min_rms' (abspath str)
                    'max_rms' (abspath str)

        Returns:
            (-1: "ERROR...") if there is an issue

            (id (long), "Spectrum calibration added") if it completes successfully
        """
        param_types = {'id': int, 'dome': str, 'bias': str, 'flat': str, 'cosmic_filter': bool, 'drpver': float,
                       'Hg_master': str, 'Cd_master': str, 'Xe_master': str, 'avg_rms': str, 'min_rms': str,
                       'max_rms': str}
        id = _id_from_time()
        pardic['id'] = id
        # TODO: which parameters are required? test
        keys = list(pardic.keys())

        for key in ['bias', 'flat']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['id', 'dome', 'bias', 'flat', 'cosmic_filter', 'drpver', 'Hg_master', 'Cd_master',
                           'Xe_master', 'avg_rms', 'min_rms', 'max_rms']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'spec_calib')
        print sql, type_check
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_spec_calib sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_spec_calib sql command failed with a ProgrammingError!")
        return (id, "Spectrum calibration added")

    def update_spec_calib(self, pardic):
        """
        updates a spec_calib entry

        Args:
            pardic (dict):
                required:
                    'id' (int/long)
                optional:
                    'bias' (abspath str)
                    'flat' (abspath str)
                    'dome' (abspath str)
                    'cosmic_filter' (bool)
                    'drpver' (float)
                    'Hg_master' (abspath str)
                    'Cd_master' (abspath str)
                    'Xe_master' (abspath str)
                    'avg_rms' (abspath str)
                    'min_rms' (abspath str)
                    'max_rms' (abspath str)

        Returns:
            (-1, "ERROR...") if it failed to update

            (id (long), "Spec_calib updated, columns 'column_names'") if the entry is updated successfully
        """
        param_types = {'id': int, 'dome': str, 'bias': str, 'flat': str, 'cosmic_filter': bool, 'drpver': float,
                       'Hg_master': str, 'Cd_master': str, 'Xe_master': str, 'avg_rms': str, 'min_rms': str,
                       'max_rms': str}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: id not provided!")

        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM spec_calib;')]:
            return (-1, "ERROR: no spec_calib entry with the id!")

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['dome', 'bias', 'flat', 'cosmic_filter', 'drpver', 'Hg_master', 'Cd_master',
                           'Xe_master', 'avg_rms', 'min_rms', 'max_rms']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'spec_calib')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_spec_calib sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_spec_calib sql command failed with a ProgrammingError!")
        return (pardic['id'], "Spec_calib updated, columns " + str(keys)[1:-1])

    def get_from_spec_calib(self, values, where_dict={}, compare_dict={}):
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
                'dome' (abspath str)
                'bias' (abspath str)
                'flat' (abspath str)
                'cosmic_filter' (bool)
                'drpver' (float)
                'Hg_master' (abspath str)
                'Cd_master' (abspath str)
                'Xe_master' (abspath str)
                'avg_rms' (abspath str)
                'min_rms' (abspath str)
                'max_rms' (abspath str)
        Returns:
            list of tuples containing the values for calibration matching the criteria

            empty list if no spec_calib entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'dome': str, 'bias': str, 'flat': str, 'cosmic_filter': bool, 'drpver': float,
                          'Hg_master': str, 'Cd_master': str, 'Xe_master': str, 'avg_rms': str, 'min_rms': str,
                          'max_rms': str}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'spec_calib')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_flexure(self, pardic):
        """
        Creates a new object in the flexure table.

        Args:
            pardic (dict):
                required:
                    'rms' (float),
                    'spec_id_1' (int),
                    'spec_id_2' (int),
                    'timestamp1' ('year-month-day hour:minute:second'),
                    'timestamp2' ('year-month-day hour:minute:second')

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Flexure added")
        """
        param_types = {'id': int, 'rms': float, 'spec_id_1': int, 'spec_id_2': int,
                       'timestamp1': 'datetime', 'timestamp2': 'datetime'}
        # TODO: find out what the 'rms' is here
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())
        for key in ['rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(keys):
            if key not in ['id', 'rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))
        sp1_id = self.get_from_spec(['id'], {'id': pardic['spec_id_1']})
        if not sp1_id:
            return (-1, "ERROR: no spectrum exists with the given spec_id_1!")
        elif sp1_id[0] == -1:
            return sp1_id
        sp2_id = self.get_from_spec(['id'], {'id': pardic['spec_id_2']})
        if not sp2_id:
            return (-1, "ERROR: no spectrum exists with the given spec_id_2!")
        elif sp2_id[0] == -1:
            return sp2_id

        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'flexure')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_flexure sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_flexure sql command failed with a ProgrammingError!")
        return (id, "Flexure added")

    def get_from_flexure(self, values, where_dict={}, compare_dict={}):
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
                'id' (int/long),
                'rms' (float),
                'spec_id_1' (int),
                'spec_id_2' (int),
                'timestamp1' ('year-month-day hour:minute:second'),
                'timestamp2' ('year-month-day hour:minute:second')

        Returns:
            list of tuples containing the values for flexure matching the criteria

            empty list if no flexure entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'rms': float, 'spec_id_1': int, 'spec_id_2': int,
                          'timestamp1': 'datetime', 'timestamp2': 'datetime'}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'flexure')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
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
                    'phase' (float),
                    'phase_err' (float)
        Returns:
            (-1, "ERROR...") if there is an issue

            (id (long), 'Classification added") if it was successful
        """
        param_types = {'id': int, 'spec_id': int, 'object_id': int, 'classification': str, 'redshift': float,
                       'redshift_err': float, 'classifier': str, 'score': float, 'phase': float, 'phase_err': float}
        # TODO: clean up the required parameters, test
        id = _id_from_time()
        pardic['id'] = id
        keys = list(pardic.keys())
        for key in ['spec_id', 'object_id', 'classification', 'redshift', 'redshift_err', 'classifier', 'score']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        classified = self.get_from_classification(['classification', 'redshift', 'redshift_err'],
                                     {'spec_id': pardic['spec_id'], 'classifier': pardic['classifier']})
        if classified:
            if classified[0] == -1:
                return classified
            return (-1, "ERROR: entry exists for that spectrum and classifier with classification %s, redshift %s, "
                        "redshift_err %s. Use `update_classification` if necessary."
                        % (classified[0][0], classified[0][1], classified[0][2]))
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['id', 'spec_id', 'object_id', 'classification', 'redshift', 'redshift_err', 'classifier', 'score',
                           'phase', 'phase_err']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))

        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_classification sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_classification sql command failed with a ProgrammingError!")
        return (id, "Classification added")

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
                    'classification' (str),
                    'redshift' (float),
                    'redshift_err' (float),
                    'phase' (float),
                    'phase_err' (float),
                    'score' (float)
                Note: this function will not modify id, spec_id, classifier or object_id

        Returns:
            (-1, "ERROR...") if there is an issue

            (id, "Classification updated, columns 'column_names'") if it was successful
        """
        param_types = {'id': int, 'spec_id': int, 'classifier': str, 'classification': str, 'redshift': float,
                       'redshift_err': float, 'phase': float, 'phase_err': float, 'score': float}
        keys = list(pardic.keys())
        if 'id' in keys:
            id_classifier = self.get_from_classification(['spec_id', 'classifier'], {'id': pardic['id']})
            if not id_classifier:
                return (-1, "ERROR: no classification entry with the given id!")
            elif id_classifier[0] == -1:
                return id_classifier

            if 'classifier' in keys:
                if not pardic['classifier'] == id_classifier[0][1]:
                    return (-1, "ERROR: classifier provided does not match classification id!")
            if 'spec_id' in keys:
                if not pardic['spec_id'] == id_classifier[0][0]:
                    return (-1, "ERROR: spec_id provided does not match classification id!")
        elif 'spec_id' in keys and 'classifier' in keys:
            id = self.get_from_classification(['id'], {'spec_id': pardic['spec_id'],
                                                       'classifier': pardic['classifier']})
            if not id:
                return (-1, "ERROR: no classification entry with the given spec_id and classifier!")
            elif id[0] == -1:
                return id
            pardic['id'] = id[0][0]
        else:
            return (-1, "ERROR: needs id or both spec_id and classifier")

        for key in reversed(keys):  # remove 'object_id', 'classifier' and any invalid keys
            if key not in ['id', 'classification', 'redshift', 'redshift_err', 'phase', 'phase_err', 'score']:
                return (-1, "ERROR: %s is an invalid key!" % (key,))

        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sp_id = self.get_from_spec(['id'], {'id': pardic['spec_id']})
        if not sp_id:
            return (-1, "ERROR: no spectrum exists with the given spec_id!")
        elif sp_id[0] == -1:
            return sp_id
        obj_id = self.get_from_object(['id'], {'id': pardic['object_id']})
        if not obj_id:
            return (-1, "ERROR: no object exists with the given object_id!")
        elif obj_id[0] == -1:
            return obj_id

        sql = _generate_update_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_classification sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_classification sql command failed with a ProgrammingError!")
        return (pardic['id'], "Classification updated, columns " + str(keys)[1:-1])

    def get_from_classification(self, values, where_dict={}, compare_dict={}):
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
                'phase' (float),
                'phase_err' (float),
                'score' (float)

        Returns:
            list of tuples containing the values for classifications matching the criteria

            empty list if no classification entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'spec_id': int, 'object_id': int, 'classification': str, 'redshift': float, 'id': int,
                          'redshift_err': float, 'classifier': str, 'score': float, 'phase': float, 'phase_err': float}

        sql = _generate_select_sql(values, where_dict, allowed_params, compare_dict, 'classification')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results


def _data_type_check(keys, pardic, value_types):
    """
    make sure the values given match the data types required in the database

    Args:
        keys (list): keys to be tested
        pardic (dict): keys and keys' values
        value_types (dict): keys and types the values should be (e.g. {'ra': float})

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
                return "ERROR: %s must be of %s or %s!" % (key, str(int)[1:-1], str(long)[1:-1])
        elif value_types[key] == 'date':
            try:
                pardic[key] = str(Time(pardic[key])).split(' ')[0]
            except (ValueError, TypeError):
                return "ERROR: %s must be of the format 'year-month-day'!" % (key,)
        elif value_types[key] == 'datetime':
            try:
                pardic[key] = str(Time(pardic[key]))
            except (ValueError, TypeError):
                return "ERROR: %s must be of the format 'year-month-day hour:minute:second'!" % (key,)
        elif value_types[key] == 'timedelta':
            if not isinstance(pardic[key], timedelta):
                try:
                    pardic[key] = timedelta(0, pardic[key])
                except (TypeError, ValueError):
                    return "ERROR: %s must be a datetime.timedelta object or a float(seconds)!" % (key,)
        elif value_types[key] == 'bool':
            if not (pardic[key] == 'false' or pardic[key] == 'true'):
                return "ERROR: %s must be either 'true' or 'false'" % (key,)
        elif not isinstance(pardic[key], value_types[key]):
            return "ERROR: parameter %s must be of %s!" % (key, str(value_types[key])[1:-1])
    return None


def _generate_select_sql(values, where_dict, allowed_params, compare_dict, table):
    """
    generate the sql for a select query

    Args:
        values (list): list of names of values
            list of values to return
        where_dict (dict):
            {'param':'value',...} adds WHERE param='value'...
        allowed_params (dict):
            the parameters of the table to be queried and their type {'param':type,...}
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
            return (-1, "ERROR: %s is an invalid column name!" % (value,))
    where_keys = list(where_dict.keys())
    for param in reversed(where_keys):
        if param not in allowed_params:
            return "ERROR: requested condition on nonexistent column '%s'!" % (param,)
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
            sql += " %s %s '%s' AND" % (key, compare_dict.get(key, '='), where_dict[key])
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
        if pardic[param]:  # make sure that there is a value
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
    sql = "UPDATE %s SET" % (table,)
    for param in param_list:
        if pardic[param]:  # it may be a key with nothing in it
            sql += " %s = '%s'," % (param, pardic[param])
    if lastmodified:
        sql += " lastmodified = 'NOW()'"
    else:
        sql = sql[:-1]
    sql += " WHERE id = %s;" % (pardic['id'],)
    return sql


def _id_from_time():
    """Generate an id from the current time of format YYYYMMDDHHMMSSsss"""
    time = Time.now()
    id = time.iso.translate(None, '- :.')
    return long(id)
