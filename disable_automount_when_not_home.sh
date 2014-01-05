#!/bin/bash

ping=`ping -c 1 badabing.local 2>/dev/null | wc | cut -d ' ' -f 7`;
running=`ps -ef | grep automount | grep sbin | wc -l`

## If badabing is not accessible
if [ $ping == 0 ]; then
    if [ $running == 1 ] ; then
    	echo "Stopping autofs...";
	service autofs stop;
    else
	echo "Autofs is not running. Exiting";
    fi
## If it is
else
    if [ $running == 0 ]; then
	echo "Starting autofs";
	service autofs start;
    else
	echo "Autofs is already running. Exiting";
	exit;
    fi
fi
