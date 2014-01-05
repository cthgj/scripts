#!/bin/bash

status=`xinput  list-props "AlpsPS/2 ALPS DualPoint TouchPad" | grep "Device Enabled" | gawk '{print $NF}'`;

if (( $status==1 )); then
    `xinput -set-int-prop "AlpsPS/2 ALPS DualPoint TouchPad" "Device Enabled" 8 0`
else
    `xinput -set-int-prop "AlpsPS/2 ALPS DualPoint TouchPad" "Device Enabled" 8 1`
fi
