#!/bin/bash
#ip=`wget -qO - http://cfaj.freeshell.org/ipaddr.cgi`
echo "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbB"
screens=`disper -l | grep -c display`
card=`/sbin/ifconfig | grep -B 1 Bcast | head -1 | gawk '{print $1}'`




if [ $screens == '1' ]
then
    conkyrc=~/.conkyrc_one_screen
else
    conkyrc=~/.conkyrc_two_screens
#
    conkyrc=~/.conkyrc_left_screen
fi

if [ $1 ]
then
    conkyrc=$1;
fi
echo "cc $conkyrc";
echo "sed 's/eth./$card/g' $conkyrc > ~/silly_old_tempfile"
sed "s/eth./$card/g" $conkyrc  > ~/silly_old_tempfile
#cp ~/silly_old_tempfile ~/.conkyrc;

killall -9 conky
#sleep 20
echo "conky -c ~/silly_old_tempfile "
perl ~/scripts/weather.pl mar c l > .conky_location
conky -c ~/silly_old_tempfile &

