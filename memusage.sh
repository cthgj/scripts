#!/bin/bash 

## If we have given a process name to monitor
if [[ $1 ]];     then

    while [ 1 ]; do
#	pidno=`ps xo pid,comm  | grep -v "ps -ao" | grep -v emacs | grep -v grep | grep -v $0 | grep $1 | gawk '{print \$1}' | head -1`
	pidno=`pgrep $1`
	mem=`top -bd .10 -p $pidno -n 1  | grep $pidno | gawk '{print \$10}' | sed 's/\..*//'`
	echo $mem;
# if ! [[ "$yournumber" =~ ^[0-9]+$ ]] ; then
#    exec >&2; echo "error: Not a number"; exit 1
# fi
	## If the second parameter passed is a number
	## take it as the mem limit
	if [[ "$2" =~ ^[0-9]+$ ]]; then
	    limit=$2;
	else
	    limit=70;
	fi

	if [ $mem -gt $limit ]
	then
	    kill $pidno;
	    exit
	fi
	sleep 5;
    done

## else just monitor/kill the one that uses the most memory
## when total memory is almost full
else
    while [ 1 ]; do
	pidno=`top -bd .10 -n 1 | grep "^[ 0-9]" | sort -nk 10 | tail -1 | gawk '{print $1}'`;
	name=`top -bd .10 -n 1 | grep "^[0-9]" | sort -nk 10 | tail -1 | gawk '{print $NF}'`;
	memfree=`free | grep "\-/+" |gawk '{print $4}'`; #top -bd .10 -n 1 | rgrep "^[0-9]" | sort -nk 10 | tail -1 | gawk '{print $10}'`
	if [ $memfree -lt  50000 ]
	then              
	    echo "Killing $pidno"
	    kill $pidno && echo -n "Killed $name ($pidno): " >> $HOME/.memusage.er; date >> $HOME/.memusage.er; 
	fi
    done
fi                   
