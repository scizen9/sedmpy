#!/bin/bash

source ~/.bashrc

echo "Starting the end of night calibration"

python /scr2/sedmdrp/sedmpy/drprc/rcred.py 1>/tmp/sedmdrp_rcred_info 2> /tmp/sedmdrp_rcred_error
#python /scr2/nblago/kpy/SEDMrph/rcred.py 1>/tmp/rcred_info 2>/tmp/rcred_error
#python /scr2/nblago/kpy/SEDMrph/zeropoint.py 1>/tmp/zeropoint_info 2>/tmp/zeropoint_error  
#python /scr2/nblago/kpy/SEDMrph/app_phot.py 1>/tmp/apphot_info 2>/tmp/apphot_error

