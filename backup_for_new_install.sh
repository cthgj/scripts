#!/bin/bash

date=`date +%F`
mkdir $date;
sudo cp -rv /var/spool/ /var/www /etc /boot /VirtualBox $date;
ls / > $date/list_root
