#!/usr/bin/env bash

##############################################
# Hack to fix laptop screen going blank.     #
# This will switch it off and it should	     #
# be reinitialized when the script continues #
##############################################
xrandr --output DP-3 --off 



if ( xrandr | grep VGA | grep -w connected >/dev/null )
then
    notify-send "Extending desktop to VGA screen"
    xrandr --output DP-3 --auto --output VGA-0 --auto --right-of DP-3 --primary
else
    if(xrandr | grep DP-2 | grep connected >/dev/null )
    then
	notify-send "Extending desktop to DisplayPort screen"
	xrandr --output DP-3 --auto --output DP-2 --auto --right-of DP-3 --primary
    else
	notify-send "No known screens found"
    fi
    
fi
