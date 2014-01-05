#!/bin/sh
while true; do
    watch -g dpkg --get-selections
    echo CHANGED
    sleep 1
done
