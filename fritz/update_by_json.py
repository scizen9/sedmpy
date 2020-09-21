import requests
import json


def get_marshal_id(request_id):
    """
    The request_id can be found in the header of SEDm fits files as 'req_id'
    :param request_id: 
    :return: 
    """
    url = 'http://127.0.0.1:5000/get_marshal_id'
    payload = {'request_id': request_id}

    headers = {'content-type': 'application/json'}
    json_data = json.dumps(payload)
    response = requests.post(url, data=json_data, headers=headers)

    ret = json.loads(response.text)
    print(ret)
    if 'error' in ret:
        return False
    else:
        return ret['marshal_id']


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


def update_request(status, marshal_id, instrument_id=65,
                   user='', pwd='', save_as=''):
    """
    Function to update the status of any request as long as it has
    not been deleted. The new status will show up on the status section
    of the request on the fritz marshal.

    :param filename:
    :param status:
    :param instrument_id:
    :param request_id:
    :param output_dir:
    :param save:
    :return:
    """
    fritz_base_url = 'http://skipper.caltech.edu:8080/cgi-bin/fritz/'
    fritz_stat_url = fritz_base_url + 'update_followup_status.cgi'

    # 1. Make sure that we have the required fields
    if not marshal_id:
        print("Can't update a request without the request id")
        return False

    if not save_as:
        save_as = "%s.txt" % (str(marshal_id))

    output_file = save_as

    # 2. Create the new status dictionary
    status_config = {'instrument_id': instrument_id,
                     'request_id': marshal_id,
                     'new_status': status}

    # 3. Write and read in the json file to memory
    json_file = open(write_json_file(status_config, output_file), 'r')

    # 4. Send the request, close the file, and save if needed
    ret = requests.post(fritz_stat_url, auth=(user, pwd),
                        files={'jsonfile': json_file})
    json_file.close()

    # 5. Print the request response for the user
    print(ret)

    return True


if __name__ == "__main__":
    marshal_id = get_marshal_id(request_id=20181022160114460)
    update_request('Phot Completed', marshal_id, instrument_id=65, user='',
                   pwd='', save_as='')
