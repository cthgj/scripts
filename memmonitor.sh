#!/bin/bash 

mem=0;
## If we have given a process name to monitor
while [ 1 ]; do
    pid=`pgrep $1`
    size=$(ps  -o size  -p "$pid" --no-headers); 
    bytes=$((1024*1024*1024))
    gigs=$(($size/$bytes))
    echo "$size $bytes $gigs"
    exit
    if [[ $aa > 0 ]] 
    then
    size=`top -cbd .10 -n 1 | grep -w $1 | grep -v grep | gawk '(\$NF=="$1"){print \$7}' | perl -ne 's/(\d+)\w+/\$1/; print;'`
	echo "YES $size top -cbd .10 -n 1 | grep -w $1 | grep -v grep | gawk '(\$NF=="$1"){print \$7}' | perl -ne 's/(\d+)\w+/\$1/; print;'"
    if [[ $size > $mem ]]
       then 
	mem=$size;
	echo "Max so far : $mem";
    fi
    sleep 5 ;
    else
	echo "Greatest memory use for $f was : $mem";
	exit
    fi
done
