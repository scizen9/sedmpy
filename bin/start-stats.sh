#!/bin/bash

source ~/.bashrc

echo "Starting the stats compilation for the weather"
python /data/sedmdrp/sedmpy/drprc/stats.py 1>/tmp/sedmdrp_stats_error 2>>/tmp/sedmdrp_stats_error
