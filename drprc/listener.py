# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 22:22:15 2016

@author: nadiablago
"""
import os
import time
try:
    import recenter_ifu
except ImportError:
    import drprc.recenter_ifu as recenter_ifu
import socket
try:
    import sextractor
except ImportError:
    import drprc.sextractor as sextractor
import sys
try:
    import sao
except ImportError:
    import drprc.sao as sao
import subprocess
import datetime
import logging
import numpy as np

from configparser import ConfigParser
import codecs

import matplotlib
matplotlib.use("Agg", warn=False)

parser = ConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.read_file(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')

_host = parser.get('listener', 'host')
_port = parser.getint('listener', 'port')
_alivefile = parser.get('listener', 'alivefile')
_rsa = parser.get('listener', 'rsa')


def start_listening_loop():
    """
    Start accepting connections from pylos.
    """

    # Set address
    ip = _host
    port = _port
    
    # Log into a file
    log_format = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
    root_dir = _logpath

    now = datetime.datetime.utcnow()
    timestamp = datetime.datetime.isoformat(now)
    timestamp = timestamp.split("T")[0]
    logging.basicConfig(
        format=log_format,
        filename=os.path.join(root_dir, "listener_{0}.log".format(timestamp)),
        level=logging.INFO)
    logger = logging.getLogger('listener')
            
    # bind socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.bind((ip, port))
    except socket.error:
        # If there is already a listener there, it will generate an error.
        # In this case we just return.
        return
    s.listen(10)
    
    logger.info("Starting a new listener. Hello World!")
    
    # create continous while loop to listen for request
    # Exit the loop at 11:00AM, as the new day will start.
    while True:
        connection, caddress = s.accept()
        time.sleep(1)
        cmd = "touch %s" % _alivefile
        subprocess.call(cmd, shell=True)

        while True:

            subprocess.call(cmd, shell=True)
            
            data = connection.recv(2048)
            logger.info("Incoming command: %s " % data)
            
            now = datetime.datetime.utcnow()
            timestamp = datetime.datetime.isoformat(now)
            timestamp = timestamp.split("T")[0]
            logging.basicConfig(
                format=log_format,
                filename=os.path.join(root_dir,
                                      "listener_{0}.log".format(timestamp)),
                level=logging.INFO)
            logger = logging.getLogger('listener')
    
            command = data.split(b',')[0]
            
            if b"FOCUSIFU" in command:
                logger.info("Finding the best focus for the IFU.")
                lfiles = data.split(b",")[1:]
                for i, lf in enumerate(lfiles):
                    lf = lf.replace(b"raw", b"phot")
                    lfiles[i] = lf
                # Wait until the image is available.
                while not os.path.isfile(lfiles[-1]):
                    time.sleep(0.5)
                focus, sigma = sextractor.get_focus_ifu(lfiles, plot=True)
                connection.sendall(b"%.2f,%.2f\n" % (focus, sigma))
                logger.info("Selected focus ifu: %.2f,%.2f\n" % (focus, sigma))
            elif b"FOCUS" in command:
                logger.info("Finding the best focus for RC.")
                lfiles = data.split(b",")[1:]
                for i, lf in enumerate(lfiles):
                    lf = lf.replace(b"raw", b"phot")
                    lfiles[i] = lf
                # Wait until the image is available.
                while not os.path.isfile(lfiles[-1]):
                    time.sleep(0.5)
                focus, sigma = sextractor.get_focus(lfiles, plot=True)
                connection.sendall(b"%.2f,%.2f\n" % (focus, sigma))
                logger.info("Selected focus: %.2f,%.2f\n" % (focus, sigma))
            elif b"SAO" in command:
                try:
                    logger.info("Looking for a nice SAO star.")
                    name, ra, dec = sao.get_sao()
                    connection.sendall(b"%d,%s,%s,%s\n" % (0, name, ra, dec))
                    logger.info("Found star. Returning: %d,%s,%s,%s\n" %
                                (0, name, ra, dec))
                except Exception as e:
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
                    connection.sendall(
                        b"%d,%s,%s,%s\n" %
                        (-1, "null", (datetime.datetime.utcnow()[4]+3)*15, 40))

            elif b"STATS" in command:
                logger.info("Finding the statistics for the image.")
                lfile = data.split(b",")[1]
                lfile = lfile.strip()
                lfile = lfile.replace(b"raw", b"phot")
                while not os.path.isfile(lfile):
                    time.sleep(0.5)
                try:
                    nsources, fwhm, ellipticity, bkg = \
                        sextractor.get_image_pars(lfile)
                    connection.sendall(
                        b"%d,%d,%d,%.3f,%.3f\n" %
                        (0, nsources, fwhm, ellipticity, bkg))
                except Exception as e:
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
                    connection.sendall(b"%d,%d,%d,%d,%.3f\n" % (-1, 0, 0, 0, 0))
                    
            elif b"FWHM" in command:
                logger.info("Finding the FWHM for the last 3 images.")
                mydir = timestamp.replace("-", "")
                statsfile = os.path.join("%s/%s/stats" % (_photpath, mydir),
                                         "stats_%s.log" % mydir)
                if os.path.isfile(statsfile):
                    stats = np.genfromtxt(statsfile, dtype=None, delimiter=",")
                    if len(stats) > 3:
                        connection.sendall(b"%d,%.2f\n" %
                                           (0, np.median(stats["f4"][-4:])))
                        logger.info("Sent (%d,%.2f)" %
                                    (0, np.median(stats["f4"][-4:])))
                    else:
                        connection.sendall(b"%d,%.2f\n" %
                                           (0, np.median(stats["f4"])))
                        logger.info("Sent (%d,%.2f)" %
                                    (0, np.median(stats["f4"])))

                else:
                    logger.error("Stats file %s does not exists!" % statsfile)
                    connection.sendall(b"%d,%.2f\n" % (-1, 0))

            elif b"OFFSET" in command:
                try:
                    is_ab = False
                    astrometry = 0
                    image = data.split(b",")[2].rstrip()
                    logger.info(
                        "Get Offsets AB=%s for image %s. Astrometry active=%s" %
                        (is_ab, image, astrometry))
                    # if(astrometry==0):
                    astrofile = os.path.basename(image)
                    date = astrofile.split("_")[0].replace(b"rc", b"")
                    astrofile = astrofile.replace("rc",
                                                  "a_rc").replace(".new",
                                                                  ".fits")
                    endpath = "%s/%s/%s" % (_photpath, date, astrofile)
                    if not os.path.isdir(os.path.dirname(endpath)):
                        os.makedirs(os.path.dirname(endpath))
                    os.system(
                        'scp %s developer@p200-guider.palomar.caltech.edu:%s %s'
                        % (_rsa, image, endpath))
                    astro = False
                    # Only run astrometry if image is unavailable.
                    if not os.path.isfile(endpath):
                        astro = True
                        logger.error(
                            "Astrometry resolved image %s could not be copied "
                            "into %s. Setting astrometry to True." %
                            (image, endpath))
                        endpath = os.path.join(os.path.dirname(endpath),
                                               os.path.basename(image))
                    res = recenter_ifu.main(endpath, astro=astro, doplot=True)
                    retcode = res[0]
                    offsets = res
                    logger.info("OFFSETS %s  Return code %d" %
                                (offsets, retcode))
                    connection.sendall(b"%d,%s,%s\n" %
                                       (retcode, offsets[1], offsets[2]))
                except Exception as e:
                    logger.error("Error occurred when processing command  " +
                                 str(data))
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
                    connection.sendall(b"%d,%s,%s\n" % (-1, 0, 0))
            elif b"PING" in command:
                try:
                    logger.info("Received a PING command. Responding...")
                    connection.sendall(b"PONG\n")
                    logger.info("PONG sent.\n")
                except Exception as e:
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
            else:
                logger.error("Unknown command: %s" % command)
                break
        # Exit and restart the listener if the program is running for more
        # than 12 hours and it is later than 10AM, so we don't disrupt the
        # night scheduler.
        '''if (datetime.datetime.utcnow()-now).total_seconds() > 12*3600. and 
        (datetime.datetime.utcnow()).hour>17:
            logger.info( "Programe is running for more than 12h and now is 
            later than 10AM. Restarting")
            sys.exit(0)'''


if __name__ == '__main__':
    """
    If the program has not been running since 1 minute,
    then the main funciton relaunches the main loop.
    Otherwise it does nothing.
    """
    # If the file was modified last time more than 60s ago,
    # relaunch the listener.
    modified = datetime.datetime.strptime(
        time.ctime(os.path.getmtime(_alivefile)), "%a %b  %d %H:%M:%S %Y")
    t_now = datetime.datetime.now()
    
    if (t_now - modified).seconds > 10:
        start_listening_loop()
