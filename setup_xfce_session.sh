#!/bin/bash
cd /home/cchapple
setxkbmap -option grp:switch,grp:alt_shift_toggle us,es,gr
gnome-do &
start_conky.sh ~/.conkyrc&
compiz --replace&
gtk-window-decorator --replace&
killall tint2; 
tint2&
