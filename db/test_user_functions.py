import SedmDb

# default entries (shouldn't be changed by functions)
# table: users; 'id': 1, 'username': 'default_user', 'name': 'apo', 'email': 'kiwi'
# table: gropus; 'id': 1, 'designator': 'default group'
# table: groups; 'id': 2, 'designator': 'second default group'
db = SedmDb.SedmDB()


# test add_group
def test_add_group():  # add_to_group is tested below
    # test unsuccessful add_group
    assert db.add_group({'designator': 'default group'}) == (-1, "ERROR: group exists!")
    assert db.add_group({'id': 500}) == (-1, "ERROR: no group designator provided!")
    # TODO: complete add_group functions?


# testing user/usergroup manipulation
def test_user_manipulation():
    user_dict = {'username': 'test_user', 'name': 'nmo', 'email': 'dodo'}
    # test a successful add_user
    added = db.add_user(user_dict)
        # check that the function returned a positive and that it was successful
    assert added[0] == 0
    assert [user[0] for user in db.execute_sql("SELECT username FROM users WHERE username='test_user';")]
    # test get_from_users (no WHERE, WHERE, WHERE with no results)
    assert (db.execute_sql('SELECT username, id, name FROM users') ==
            db.get_from_users(['username', 'id', 'name'], {}))
    assert (db.execute_sql("SELECT username, id, name FROM users WHERE name='nmo'") ==
            db.get_from_users(['username', 'id', 'name'], {'name': 'nmo'}))
    assert (db.execute_sql("SELECT username, id, name FROM users WHERE username='aaa'") ==
            db.get_from_users(['username', 'id', 'name'], {'username': 'aaa'}))
    # test a successful remove_user
    removed = db.remove_user({'username': 'test_user'})
    assert removed[0] == 0
    assert not db.get_from_users(['username'], {'username': 'test_user'})
    # test unsuccessful add_user
    add_again = db.add_user({'username': 'default_user', 'name': 'apo', 'email': 'kiwi'})
    assert add_again[0] == -1
    no_username = db.add_user({'name': 'apo'})
    assert no_username[0] == -1
    # test unsuccessful remove_user
    assert db.remove_user({'username': 'unused_name'}) == (-1, "ERROR: no user with that username!")
    assert db.remove_user({'id': 0}) == (-1, "ERROR: no user with that id!")
    assert db.remove_user({'name': 'neo', 'email': '3po'}) == (-1, "ERROR: username or id required!")
    assert db.remove_user({}) == (-1, "ERROR: username or id required!")
    # test successful add_to_group
    assert db.add_usergroup(1, 2)[0] == 0
    assert db.add_usergroup(1, 1)[0] == 0  # add to a second group
    usergroups = db.execute_sql("SELECT user_id, group_id FROM usergroups")
    assert (1, 2) in usergroups
    assert (1, 1) in usergroups
    # test get_from_usergroups
    assert usergroups == db.get_from_usergroups(['user_id', 'group_id'], {})
    assert not db.get_from_usergroups(['user_id'], {'user_id': 0})
    assert db.get_from_usergroups(['id'], {}) == (-1, "ERROR: no valid values requested!")
    # test unsuccessful add_to_group
    assert db.add_usergroup(1, 2) == (-1, "ERROR: user already in group!")
    assert db.add_usergroup(1, 0) == (-1, "ERROR: group does not exist!")
    assert db.add_usergroup(0, 1) == (-1, "ERROR: user does not exist!")
    # test successful remove_from_group
    assert db.remove_from_group(1, 2)[0] == 0
    assert db.remove_from_group(1, 1)[0] == 0
    assert not db.execute_sql("SELECT user_id, group_id FROM usergroups WHERE user_id=1")
    # test unsuccessful remove_from_group
    assert db.remove_from_group(0, 1) == (-1, "ERROR: user does not exist!")
    assert db.remove_from_group(1, 0) == (-1, "ERROR: group does not exist!")
    assert db.remove_from_group(1, 1) == (-1, "ERROR: user not in group!")

if __name__ == '__main__':
    test_add_group()
    test_user_manipulation()

