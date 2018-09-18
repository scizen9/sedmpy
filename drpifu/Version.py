import subprocess
import os
import datetime


def ifu_drp_version():
    try:
        cwd = os.getcwd()
        home = os.environ["HOME"]
        os.chdir(home+"/sedmpy/drpifu")
        ver = subprocess.check_output(["git", "describe", "--always"],
                                      universal_newlines=True).strip()
        os.chdir(cwd)
    except:
        ver = "%s" % datetime.datetime.now()

    return str(ver)
