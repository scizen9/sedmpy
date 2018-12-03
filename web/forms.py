from flask_wtf import FlaskForm
from wtforms import fields, validators
import datetime

SECRET_KEY = 'secret'


class LoginForm(FlaskForm):
    username = fields.StringField('username')
    password = fields.PasswordField('password')


class PassChangeForm(FlaskForm):
    password = fields.PasswordField('Old Password',
                                    validators=[validators.input_required()])
    pass_new = fields.PasswordField('New Password',
                                    validators=[validators.input_required(),
                                                validators.EqualTo('pass_conf',
                                                                   message='Passwords must match')] )
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
    obj_name = fields.StringField('Object Name (Names should be limited to 26 '
                                  'Characters and have no blank spaces are '
                                  'special characters such as ["*_.] '
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
    allocation = fields.SelectField('allocation')

    priority = fields.FloatField('priority', default=.99)
    filters_op = fields.SelectField('filters',
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

    rc_use_mag = fields.BooleanField('Use object magnitude to set exposure time')
    seq_repeats = fields.IntegerField('# of Repeated Sequences', default=1)
    seq_completed = fields.IntegerField('# of Sequences Already Completed', default=0)
    ifu = fields.BooleanField('Take IFU Exposure')
    rc = fields.BooleanField('Take RC Exposure')
    ifu_use_mag = fields.BooleanField('Use object magnitude to set exposure time')
    ifu_exptime = fields.IntegerField('Enter Total IFU Exposure Time in seconds')
    ab = fields.BooleanField('Select For AB pair')

    cadence = fields.FloatField('cadence', default=None)

    min_moon_dist = fields.FloatField('minimum moon distance (degrees)',
                                      validators=(validators.Optional(),
                                                  validators.number_range(0., 180.)),
                                      default=30)

    max_moon_illum = fields.FloatField('Maximum Moon Illumination (fractional 0 to 1)',
                                       validators=(validators.Optional(),
                                                   validators.number_range(0., 1.)),
                                       default=1)
    maxairmass = fields.FloatField('Maximum Airmass',
                                   validators=(validators.Optional(),
                                               validators.number_range(1, 5)),
                                   default=2.5)
    max_cloud_cover = fields.FloatField('Maximum Cloud Cover (fractional)',
                                        validators=(validators.Optional(),
                                                    validators.number_range(0., 1.)),
                                        default=1)
    max_fwhm = fields.FloatField('Maximum FWHM',
                                        validators=(validators.Optional(),
                                                    validators.number_range(0., 10.)),
                                        default=10)
    phasesamples = fields.FloatField('samples per period', default=None)
    sampletolerance = fields.FloatField('samples tolerance', default=None)
    inidate = fields.DateField('start date (Y-m-d)',
                               validators=[validators.input_required()],
                               format='%Y-%m-%d',
                               default=datetime.datetime.utcnow())

    enddate = fields.DateField('end date (Y-m-d)',
                               validators=[validators.input_required()],
                               format='%Y-%m-%d',
                               default=(datetime.datetime.utcnow() +
                                        datetime.timedelta(days=3)))
    creationdate = fields.Label('creationdate', 'Creation Date')
    lastmodified = fields.Label('lastmodified', 'Last Modified')
    last_obs_jd = fields.Label('last_obs_jd', 'Last observation')
    submit_req = fields.SubmitField('Submit Request')
