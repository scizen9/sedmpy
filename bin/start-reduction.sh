#!/bin/bash

source ~/.bashrc

echo "Starting the daily reduction"
python /data/sedmdrp/sedmpy/drprc/realtimered.py -n -l 1>/tmp/sedmdrp_reduction_info 2>/tmp/sedmdrp_reduction_error
