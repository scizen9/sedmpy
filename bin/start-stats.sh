#!/bin/bash

source ~/.bashrc

echo "Starting the stats compilation for the weather"
python /scr2/sedmdrp/sedmpy/drprc/stats.py 1>/tmp/sedmdrp_stats_error 2>>/tmp/sedmdrp_stats_error
#python /scr2/nblago/kpy/SEDMrph/stats.py  1>/tmp/stats_error 2>>/tmp/stats_error

