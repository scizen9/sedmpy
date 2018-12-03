import os
import subprocess
import shutil
from astropy.io.ascii import SExtractor

SITE_ROOT = os.path.realpath(os.path.dirname(__file__))
config_dir = os.path.join(SITE_ROOT, 'sexConfig')


def run_sex(image, outdir="", default_to_image_dir=True, overwrite=True):
    """

    :param image:
    :param mask:
    :param cosmics:
    :param overwrite:
    :return:
    """

    if not os.path.exists(image):
        return False

    if not outdir and default_to_image_dir:
        outdir = os.path.dirname(image)
        outdir = os.path.join(outdir, 'sextractor')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    sexcat = os.path.join(outdir,
                          os.path.basename(image).replace(".fits", ".sex"))

    if not overwrite and os.path.exists(sexcat):
        print('Sextractor already run on %s' % image)
        return sexcat

    cmd = "sex -c %s/daofind.sex %s" % (config_dir, image)
    subprocess.call(cmd, stdout=subprocess.DEVNULL, shell=True)

    shutil.move("image.sex", sexcat)

    return sexcat


def analyze_sex_by_pandas(catalog, filter_data=True, return_type='dict'):
    """
    
    :param catalog: 
    :param return_type: 
    :return: 
    """

    sex = SExtractor()
    data = sex.read(catalog)
    h = data.to_pandas()
    if filter_data:
        df = h.loc[h['FLAGS'] == 0]
        df = df[(df['X_IMAGE'] < 885) | (df['X_IMAGE'] > 1540)]
        df = df[(df['Y_IMAGE'] < 850) | (df['Y_IMAGE'] > 1125)]
        filterd_df = df[(df['ELLIPTICITY'] < .3)]

        return filterd_df