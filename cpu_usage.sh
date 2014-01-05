#!/bin/bash 

## If we have given a process name to monitor

while [ 1 ]; do
    top -bd .10 -n 1 | head -15 >> $HOME/.cpu_usage
    sleep 60;
done
