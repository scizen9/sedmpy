#!/bin/bash

source ~/.bashrc

echo "Starting the daily reduction"
python /scr2/sedmdrp/sedmpy/drprc/realtimered.py -f -o 1>/tmp/reduction_info 2>/tmp/reduction_error
