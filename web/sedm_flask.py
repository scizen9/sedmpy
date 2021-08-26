from flask import Flask, render_template, request, jsonify, redirect, url_for, flash, send_from_directory
import sys

sys.path.append('/scr2/sedm/sedmpy/')
import flask_login
from web.forms import *
import json
import os
import web.model as model
# from werkzeug.datastructures import ImmutableMultiDict
from bokeh.resources import INLINE
from marshals import watcher
resources = INLINE
js_resources = resources.render_js()
css_resources = resources.render_css()

# config
SECRET_KEY = 'secret'
USERNAME = 'admin'
PASSWORD = 'default'

app = Flask(__name__)
app.config.from_object(__name__)
login_manager = flask_login.LoginManager()
login_manager.init_app(app)

config = model.get_config_paths()


class User(flask_login.UserMixin):
    pass


@login_manager.user_loader
def load_user(user_id):
    users = model.get_from_users(user_id)
    if not users:
        return None
    user = User()
    user.id = users[0][0]
    user.name = users[0][1]
    flask_login.login_user(user)  # , remember=True)
    return user


@app.route('/')
def home():
    if flask_login.current_user.is_authenticated:
        sedm_dict = model.get_homepage(flask_login.current_user.id,
                                       flask_login.current_user.name)

        return render_template('sedm.html', sedm_dict=sedm_dict)
    else:
        return redirect('login')


@app.route('/request', methods=['GET', 'POST'])
def requests():
    form = AddFixedRequest()

    # 1. If the request method is of type post then we expect this to be a
    #    submission.
    if request.method == 'POST':
        content = request.form
        # print(flask_login.current_user)
        if not flask_login.current_user.is_authenticated:
            user_id = content['user_id']
        else:
            user_id = flask_login.current_user.id
        req_dict, form = model.process_request_form(content, form, user_id)
        req_dict['message'] = req_dict['message'].replace("--", "<br>")

        return render_template('request.html', req_dict=req_dict, form=form)

    # 2. We should also be able to handle request that are of type json for
    #    automated ingestion.
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=False)

    req_dict, form = model.get_request_page(flask_login.current_user.id,
                                            form,
                                            content=content)

    return render_template('request.html', req_dict=req_dict, form=form)


@app.route('/add_csv', methods=['GET', 'POST'])
def add_csv():
    form = AddCSVRequest()

    # 1. If the request method is of type post then we expect this to be a
    #    submission.
    if request.method == 'POST':
        content = request.form
        req_dict, form = model.process_add_csv(content, form,
                                               flask_login.current_user.id)

        req_dict['message'] = req_dict['message'].replace("--", "<br>")

        return render_template('add_csv.html', req_dict=req_dict, form=form)

    # 2. We should also be able to handle request that are of type json for
    #    automated ingestion.
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=False)

    # 3. Generate the webpage dictionary if no POST or JSON request
    req_dict, form = model.get_add_csv(flask_login.current_user.id, form,
                                       content=content)

    return render_template('add_csv.html', req_dict=req_dict, form=form)


@app.route('/get_rc_redux_product', methods=['GET', 'POST'])
# @flask_login.login_required
def get_rc_redux_product():
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        obsdate = datetime.datetime.utcnow().strftime("%Y-%m-%d")
        return jsonify(model.get_rc_redux_products(obsdate))
    return jsonify(model.get_rc_redux_products(**content))


@app.route('/data_access/<path:instrument>', methods=['GET'])
# @flask_login.login_required
def data_access(instrument):
    if not flask_login.current_user.is_authenticated:
        return redirect('login')

    # 1. If the request method is of type post then we expect this to be a
    #    submission.
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=False)

    content['user_id'] = flask_login.current_user.id
    content['camera_type'] = instrument.lower()
    if instrument.lower() == 'ifu':
        out = model.get_science_products(**content)
        return render_template('view_data.html', sedm_dict=out)
    else:
        out = model.get_rc_redux_products(**content)
        # print(out)
        return render_template('view_data_redux.html', sedm_dict=out)


@app.route('/data_r/<path:filename>')
# @flask_login.login_required
def data_static_r(filename):
    # print(filename, 'this is the filename in data')
    base_path = model.base_dir
    # print(os.path.join(base_path, filename), "this is the full path")
    _p, _f = os.path.split(os.path.join(base_path, filename))
    return send_from_directory(_p, _f)


@app.route('/data/<path:filename>')
# @flask_login.login_required
def data_static(filename):
    """
     Get files from the archive
    :param filename:
    :return:
    """
    _p, _f = os.path.split(filename)
    if _f.startswith('finder') and 'ACQ' in _f:
        if os.path.exists(os.path.join(config['path']['path_archive'], _p,
                                       'finders', _f)):
            return send_from_directory(
                os.path.join(config['path']['path_archive'], _p, 'finders'), _f)
        else:
            return send_from_directory(os.path.join(config['path']['path_phot'],
                                                    _p, 'finders'), _f)
    elif _f.startswith('rc') or _f.startswith('finder') or 'ACQ' in _f:
        # print("In rc path")
        # print(_p, _f)
        if _f.startswith('rc'):
            if 'redux' in _p:
                return send_from_directory(
                    os.path.join(config['path']['path_phot'], _p), _f)
            else:
                # print("USING THE HERE")
                base_obspath = os.path.join(config['path']['path_redux_phot'],
                                            _p, 'pngraw')
                pathlist = ['acquisition', 'bias', 'dome', 'focus',
                            'guider', 'science', 'twilight']
                for i in pathlist:
                    test_path = os.path.join(base_obspath, i)
                    if os.path.exists(os.path.join(test_path, _f)):
                        return send_from_directory(test_path, _f)
        else:
            return send_from_directory(os.path.join(config['path']['path_phot'],
                                                    _p), _f)
    else:
        return send_from_directory(os.path.join(config['path']['path_archive'],
                                                _p), _f)


@app.route('/visibility')
def active_visibility():
    sedm_dict = model.get_active_visibility(flask_login.current_user.id)
    sedm_dict['js_resources'] = INLINE.render_js()  # TODO
    sedm_dict['css_resources'] = INLINE.render_css()

    return render_template('visibility.html', sedm_dict=sedm_dict)


@app.route('/weather_stats', methods=['GET', 'POST'])
def weather_stats():
    # 1. If the request method is of type post then we expect this to be a
    #    submission.
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=True)

    out = model.get_weather_stats(**content)
    out['js_resources'] = INLINE.render_js()
    out['css_resources'] = INLINE.render_css()
    return render_template('weather_stats.html', sedm_dict=out)


@app.route('/add_growth', methods=['GET', 'POST'])
def add_growth():
    x = request.files['jsonfile'].read()
    if x:
        content = json.loads(x)
        output = open(
            '/scr2/sedm/sedmpy/web/static/request_%s.txt' %
            datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S.%f"), 'w')
        data = json.dumps(content)
        output.write(data)
        output.close()
        return 'Content-type: text/html\n <title>Test CGI</title>'
    else:
        print('Not a json file')
        return 'ERROR'


@app.route('/add_request', methods=['GET', 'POST'])
def add_request():
    origin_url = request.environ['REMOTE_ADDR']

    # First check for the type of request
    if request.data:
        # Assume data string in json format
        content = json.loads(request.data)
    elif request.is_json:
        content = json.loads(request.get_json())
    elif 'jsonfile' in request.files:
        data = request.files['jsonfile'].read()
        content = json.loads(data)
    else:
        return 'Content-type: text/html\n ' \
               '<title>Invalid Formatting of Request</title>'

    # Next add the origins url to determine where the request came from
    content['origins_url'] = origin_url
    print(content)
    try:
        watcher.process_new_request(content, isfile=False)
    except Exception as e:
        print(str(e))

    # Now write the request to file
    output = open(
        'new_request_%s.txt'
        % datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S.%f"), 'w')

    data = json.dumps(content)
    output.write(data)
    output.close()

    return 'Content-type: text/html\n <title>Accepted Request</title>'


@app.route('/add_fritz', methods=['GET', 'POST'])
def add_fritz():
    print(request.data)
    print(request.form)
    origin_url = request.url
    if request.data:
        content = json.loads(request.data)
        content['origins_url'] = origin_url
        # output = open('fritz_request_%s.txt' %
        # datetime.datetime.utcnow().strftime(
        #    "%Y%m%d_%H_%M_%S.%f"), 'w')
        output = open(
            '/scr2/sedm/sedmpy/web/static/fritz_request_%s.txt'
            % datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S.%f"), 'w')

        data = json.dumps(content)
        output.write(data)
        output.close()
        return 'Content-type: text/html\n <title>Accepted Fritz CGI</title>'
    if request.is_json:
        content = json.loads(request.get_json())
        content['origins_url'] = origin_url
        output = open(
            '/scr2/sedm/sedmpy/web/static/fritz_request_%s.txt'
            % datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S.%f"), 'w')
        data = json.dumps(content)
        output.write(data)
        output.close()
        return 'Content-type: text/html\n <title>Accepted Fritz CGI</title>'
    if request.form:
        content = request.form.to_dict(flat=True)
        content['origins_url'] = origin_url
        output = open(
            '/scr2/sedm/sedmpy/web/static/fritz_request_%s.txt'
            % datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S.%f"), 'w')
        data = json.dumps(content)
        output.write(data)
        output.close()
        return 'Content-type: text/html\n <title>Accepted Fritz CGI</title>'

    if 'jsonfile' in request.files:
        x = request.files['jsonfile'].read()
    else:
        return "No json file"
    if x:
        content = json.loads(x)
        content['origins_url'] = origin_url
        output = open(
            '/scr2/sedm/sedmpy/web/static/fritz_request_%s.txt'
            % datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S.%f"), 'w')
        data = json.dumps(content)
        output.write(data)
        output.close()
        return 'Content-type: text/html\n <title>Accepted Fritz CGI</title>'
    else:
        print('Not a json file')
        return 'ERROR'


@app.route('/get_marshal_id', methods=['GET', 'POST'])
def get_marhsal_id():
    # 1. If the request method is of type post then we expect this to be a
    #    submission.

    if request.is_json:
        if isinstance(request.get_json(), dict):
            content = request.get_json()
        else:
            content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=True)

    out = model.get_marshal_id(**content)
    return jsonify(out)


@app.route('/get_user_observations', methods=['GET', 'POST'])
def get_user_observations():
    # 1. If the request method is of type post then we expect this to be a
    #    submission.

    if request.is_json:
        if isinstance(request.get_json(), dict):
            content = request.get_json()
        else:
            content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=True)

    out = model.get_user_observations(**content)
    return jsonify(out)


@app.route('/objects', methods=['GET', 'POST'])
def objects():
    form = FindObject()

    if request.method == 'POST':
        object_name = form.object_name.data
        user_id = flask_login.current_user.id
        out = model.get_object(object_name, user_id)
        return render_template('objects2.html', sedm_dict=out, form=form)

    return render_template('objects2.html', sedm_dict={}, form=form)


@app.route('/project_stats', methods=['GET', 'POST'])
def project_stats():
    # form = AddFixedRequest()
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=True)

    out = model.get_project_stats(content, user_id=flask_login.current_user.id)
    out['js_resources'] = INLINE.render_js()
    out['css_resources'] = INLINE.render_css()
    return render_template('project_stats.html', sedm_dict=out)


@app.route('/scheduler', methods=['GET', 'POST'])
def scheduler():
    out = model.get_schedule()
    return render_template('scheduler.html', sedm_dict=out)


@app.route('/search/get_objects', methods=['GET', 'POST'])
def get_object():
    if request.is_json:
        content = json.loads(request.get_json())
    elif request.method == 'POST':
        content = request.form.to_dict(flat=False)
    else:
        content = request.args.to_dict(flat=False)

    out = model.get_object_info(**content)
    return jsonify(out)


@app.route('/search/get_object_values', methods=['GET', 'POST'])
def get_object_values():
    if request.is_json:
        content = json.loads(request.get_json())
    elif request.method == 'POST':
        content = request.form.to_dict(flat=False)
    else:
        content = request.args.to_dict(flat=False)

    out = model.get_object_values(content['id'])

    return jsonify(out)


@app.route('/status/get_update', methods=['GET', 'POST'])
def update():
    out = model.get_status()

    return jsonify(out)


@app.route('/monitor', methods=['GET', 'POST'])
def monitor():
    out = model.get_obstimes()
    return render_template('monitor.html', sedm_dict=out)


@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        # Login and validate the user.index
        # user should be an instance of your `User` class
        username = form.username.data
        password = form.password.data

        ret = model.check_login(username, password)
        if ret[0]:
            user = User()
            user.id, username = [ret[1], username]
            flask_login.login_user(user)
            flash("Logged in as %s" % username)
            return redirect(url_for('home'))
        else:
            return render_template('login.html', message=ret[1], form=form)

    return render_template('login.html', form=form)


@app.route('/passchange', methods=['GET', 'POST'])
def login_change():
    form = PassChangeForm()

    # 1. If the request method is of type post then we expect this to be a
    #    submission.
    if request.method == 'POST':
        out = model.password_change(form, flask_login.current_user.id)
        if out['message'] == 'Password Changed!':
            return redirect(url_for('login'))
        else:
            return render_template('change_pass.html', sedm_dict=out, form=form)

    return render_template('change_pass.html', sedm_dict={'message': ''},
                           form=form)


@app.route("/logout")
@flask_login.login_required
def logout():
    flask_login.logout_user()
    return redirect(url_for('login'))


if __name__ == '__main__':
    app.run(debug=True)
