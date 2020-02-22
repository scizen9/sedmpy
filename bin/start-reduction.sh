#!/bin/bash

source ~/.bashrc

echo "Starting the daily reduction"
python /scr2/sedmdrp/sedmpy/drprc/realtimered.py 1>/tmp/sedmdrp_reduction_info 2>/tmp/sedmdrp_reduction_error
