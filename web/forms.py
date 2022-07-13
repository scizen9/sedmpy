from flask_wtf import FlaskForm
from wtforms import fields, validators
import datetime

SECRET_KEY = 'secret'


class LoginForm(FlaskForm):
    username = fields.StringField('username')
    password = fields.PasswordField('password')


class PassChangeForm(FlaskForm):
    password = fields.PasswordField(
        'Old Password', validators=[validators.input_required()])
    pass_new = fields.PasswordField(
        'New Password', validators=[
            validators.input_required(),
            validators.EqualTo('pass_conf', message='Passwords must match')])
    pass_conf = fields.PasswordField('Confirm New Password',
                                     validators=[validators.input_required()])


class AddCSVRequest(FlaskForm):
    """
    Class to handle csv request
    """

    submit = fields.SubmitField("Submit Request")
    textbox = fields.TextAreaField("Add CSV Request")


class FindObject(FlaskForm):
    object_name = fields.StringField("Enter object name")


class AddFixedRequest(FlaskForm):
    # object data
    request_id = fields.HiddenField('Hidden')
    object_id = fields.HiddenField('Hidden')
    marshal_id = fields.HiddenField('Hidden', default=-1)
    user_id = fields.HiddenField('Hidden')
    allocation_id = fields.HiddenField('Hidden')
    obj_name = fields.StringField('Object name (Names should be limited to 26 '
                                  'Characters and have no blank spaces or '
                                  'special characters such as ["*_.], '
                                  'dashes are acceptable)',
                                  [validators.input_required()])
    status = fields.SelectField('Status',
                                coerce=str,
                                choices=[('DEFAULT', "--------"),
                                         ('EXPIRED', 'EXPIRED'),
                                         ('COMPLETED', 'COMPLETED'),
                                         ('ACTIVE', 'ACTIVE'),
                                         ('CANCELED', 'CANCELED'),
                                         ('PENDING', 'PENDING')])
    obj_ra = fields.StringField('RA (HH:MM:SS or Degrees)')
    obj_dec = fields.StringField('DEC (DD:MM:SS or Degrees)')
    obj_epoch = fields.StringField('EPOCH', default='2000')
    obj_mag = fields.FloatField('Magnitude')
    allocation = fields.SelectField('Allocation')

    priority = fields.FloatField('Priority', default=.99)
    filters_op = fields.SelectField('Filters',
                                    coerce=str,
                                    choices=[(' 1, 1, 1, 1}', 'u-g-r-i'),
                                             (' 0, 1, 1, 1}', 'g-r-i'),
                                             (' 0, 1, 0, 0}', 'g'),
                                             (' 0, 0, 1, 0}', 'r'),
                                             ('0, 0, 0, 1}', 'i'),
                                             ('0, 0, 0, 0}', '-')])
    do_r = fields.BooleanField('r')
    do_g = fields.BooleanField('g')
    do_i = fields.BooleanField('i')
    do_u = fields.BooleanField('u')
    r_exptime = fields.IntegerField('Exptime (seconds)', default=0)
    g_exptime = fields.IntegerField('Exptime (seconds)', default=0)
    i_exptime = fields.IntegerField('Exptime (seconds)', default=0)
    u_exptime = fields.IntegerField('Exptime (seconds)', default=0)
    r_repeats = fields.IntegerField('# Repeats', default=1)
    g_repeats = fields.IntegerField('# Repeats', default=1)
    i_repeats = fields.IntegerField('# Repeats', default=1)
    u_repeats = fields.IntegerField('# Repeats', default=1)

    rc_use_mag = fields.BooleanField(
        'Use object magnitude to set exposure time')
    seq_repeats = fields.IntegerField('# of repeated sequences', default=1)
    seq_completed = fields.IntegerField('# of sequences already completed',
                                        default=0)
    ifu = fields.BooleanField('Take IFU exposure')
    rc = fields.BooleanField('Take RC exposure')
    ifu_use_mag = fields.BooleanField(
        'Use object magnitude to set exposure time')
    ifu_exptime = fields.IntegerField(
        'Enter total IFU exposure time in seconds')
    #ab = fields.BooleanField('Select For AB pair')

    cadence = fields.FloatField('Cadence', default=None)

    min_moon_dist = fields.FloatField(
        ' Min moon dist. (deg)',
        validators=(validators.Optional(),
                    validators.number_range(0., 180.)), default=30)

    max_moon_illum = fields.FloatField(
        ' Max moon phase (0. to 1.)',
        validators=(validators.Optional(),
                    validators.number_range(0., 1.)), default=1)
    maxairmass = fields.FloatField(
        ' Max airmass',
        validators=(validators.Optional(),
                    validators.number_range(1, 5)), default=2.5)
    max_cloud_cover = fields.FloatField(
        ' Max cloud cover (0. to 1.)',
        validators=(validators.Optional(),
                    validators.number_range(0., 1.)), default=1)
    max_fwhm = fields.FloatField(
        ' Max FWHM (asec)',
        validators=(validators.Optional(),
                    validators.number_range(0., 10.)), default=10)
    phasesamples = fields.FloatField(' Samples per period', default=None)
    sampletolerance = fields.FloatField(' Samples tolerance', default=None)
    inidate = fields.DateField(' Start date (Y-m-d)',
                               validators=[validators.input_required()],
                               format='%Y-%m-%d',
                               default=datetime.datetime.utcnow())

    enddate = fields.DateField(' End date (Y-m-d)',
                               validators=[validators.input_required()],
                               format='%Y-%m-%d',
                               default=(datetime.datetime.utcnow() +
                                        datetime.timedelta(days=3)))
    creationdate = fields.Label('creationdate', 'Creation date')
    lastmodified = fields.Label('lastmodified', 'Last modified')
    submit_req = fields.SubmitField('Submit request')
