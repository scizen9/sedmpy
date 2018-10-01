import os
import time
import subprocess
import shutil


def run_sex(image, outdir="", overwrite=False):
    """

    :param image:
    :param mask:
    :param cosmics:
    :param overwrite:
    :return:
    """

    start = time.time()
    sexcat = os.path.join(outdir,
                            os.path.basename(image).replace(".fits", ".sex"))

    if not overwrite and os.path.exists(sexcat):
        print('Sextractor already run on %s' % image)
        return sexcat


    cmd = "sex -c daofind.sex %s" % image
    ret = subprocess.call(cmd, shell=True)
    print(ret)

    shutil.move("image.sex", sexcat)

    return sexcat