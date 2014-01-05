#!/bin/bash


layout=`xset -q | grep "LED mask:" | awk '{ print $10 }'`


if (( $layout == 00000000 )); then
    layout=us
fi

if (( $layout == 00001000 )); then
    layout=es
fi

if (( $layout == 00000000 )); then
    layout=us
fi

