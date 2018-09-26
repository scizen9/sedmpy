from flask import Flask, render_template, request, jsonify, redirect, url_for, flash, send_from_directory
import sys
sys.path.append('/scr2/sedm/sedmpy/')
import flask_login
from web.forms import *
import json
import os
import web.model as model

from bokeh.resources import INLINE

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
    flask_login.login_user(user) # , remember=True)
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
        req_dict, form = model.process_request_form(content, form,
                                                    flask_login.current_user.id)
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


@app.route('/data_access/<path:instrument>', methods=['GET'])
@flask_login.login_required
def data_access(instrument):
    # 1. If the request method is of type post then we expect this to be a
    #    submission.
    if request.is_json:
        content = json.loads(request.get_json())
    else:
        content = request.args.to_dict(flat=False)

    content['user_id'] = flask_login.current_user.id
    content['camera_type'] = instrument.lower()
    out = model.get_science_products(**content)
    print(out)
    return render_template('view_data.html', sedm_dict=out)


@app.route('/data/<path:filename>')
@flask_login.login_required
def data_static(filename):
    '''
     Get files from the archive
    :param filename:
    :return:
    '''
    _p, _f = os.path.split(filename)

    if _f.startswith('finder') and 'ACQ' in _f:
        return send_from_directory(os.path.join(config['path']['path_phot'], _p, 'finders'), _f)
    elif _f.startswith('rc') or _f.startswith('finder') or 'ACQ' in _f:
        return send_from_directory(os.path.join(config['path']['path_phot'], _p), _f)
    else:
        return send_from_directory(os.path.join(config['path']['path_archive'], _p), _f)


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


@app.route('/objects', methods=['GET', 'POST'])
def objects():
    #form = AddFixedRequest()
    return render_template('sedm_base.html')

@app.route('/project_stats', methods=['GET', 'POST'])
def project_stats():
    #form = AddFixedRequest()
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
    table = model.get_schedule()
    return render_template('scheduler.html', sedm_dict={'schedulerTable': table})


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

    return render_template('login.html',  form=form)


@app.route("/logout")
@flask_login.login_required
def logout():
    flask_login.logout_user()
    return redirect(url_for('login'))


if __name__ == '__main__':
    app.run(debug=True)
