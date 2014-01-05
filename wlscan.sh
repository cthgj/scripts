#!/bin/bash

#for n in [0..$1];
let n=10
if(($1))
then
    let n=$1;
fi

while (($n >0))
do
    /sbin/iwlist eth1 scan | color -c red  -l freebosc
    let n=n-1
    sleep 5;
done
