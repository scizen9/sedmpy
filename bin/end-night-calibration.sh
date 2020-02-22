#!/bin/bash

source ~/.bashrc

echo "Starting the end of night calibration"

python /scr2/sedmdrp/sedmpy/drprc/rcred.py 1>/tmp/sedmdrp_rcred_info 2> /tmp/sedmdrp_rcred_error
