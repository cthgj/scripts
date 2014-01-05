#!/usr/bin/env bash

while true; do
#  aticonfig --odgt --adapter="$gpu" | 
 cat ~/foo/a |  awk -v t="$1" '(/Sensor:/ && $(NF-1) < t ){exit(1)}' || echo reboot 
   sleep 1
done

