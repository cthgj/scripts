#!/bin/bash


ip=`wget -qO - http://cfaj.freeshell.org/ipaddr.cgi`
if [ $ip == '78.227.40.153' ]
then
    wakeonlan 00:1f:c6:cb:2f:20
else
    wakeonlan -i 78.227.40.153 00:1f:c6:cb:2f:20
fi
