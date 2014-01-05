#!/bin/bash
modprobe i8k

min_fanspeed=3

while true; do
   fanspeeds=(`i8kctl fan`)
    if [[ ${fanspeeds[1]} -lt 3 ]]; then
        i8kctl fan - $min_fanspeed >& /dev/null
    fi
done
