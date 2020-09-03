import os
import json
import glob

SITE_ROOT = os.path.abspath(os.path.dirname(__file__)+'/../..')

with open(os.path.join(SITE_ROOT, 'config', 'marshals.json')) as data_file:
    params = json.load(data_file)

print(params)


def update_status_request(status, request_id, marshal_name,
                          update_type='jsonfile', save=True,
                          output_file='', testing=False):
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
    :param testing:
    :return:
    """
    instrument_id = ""

    # 1. Get the instrument id
    if marshal_name.lower() not in params['marshals']:
        return {"iserror": True, "msg": "Marshal: %s not found in config file"}

    if save:
        if not output_file:
            output_file = os.path.join(marshal_name.lower(), "_", request_id,
                                       ".json")

            if os.path.exists(output_file):
                files = sorted(glob.glob("*_%s_*"))
                if not files:
                    output_file = os.path.join(marshal_name.lower(), "_",
                                               request_id, "_", "1", ".json")
                else:
                    last_file_count = files[-1].split('_')[-1].replace('.json', '')
                    last_file_count = int(last_file_count) + 1
                    output_file = os.path.join(marshal_name.lower(), "_",
                                               request_id, "_",
                                               str(last_file_count),
                                               ".json")

    # 2. Create the new status dictionary
    status_config = {'instrument_id': instrument_id,
                     'request_id': request_id,
                     'new_status': status}


