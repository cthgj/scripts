#!/bin/sh 


while [ 1 ]; do
    date=`date +%D" "%H:%M:%S`
    t1=`inxi -F | grep Sensors | cut -f 7 -d " " | perl -ne '/(\d+)/; print "\$1\n"'`
    t2=`inxi -F | grep Sensors | cut -f 11 -d " " | perl -ne '/(\d+)C/; print "\$1\n"'`
    speed=`inxi -F | grep Speeds: | cut -f 15,16 -d " "`
    temp=$(($t1+$t2))
    if [ "$temp" -ge 120 ];
    then
        echo "====== $date, CPU ($speed): $t1 GPU : $t2 ========" >> ~/.temperatures
        top -cbd .10 -n 1 |  gawk 'NR>6' | head >> ~/.temperatures
        echo "\n===============================================================\n">> ~/.temperatures
        sleep 5 ;

    fi
done
