#!/bin/bash

if [ ! "$(pidof conky)" ]
then
   #echo "the variable X is not the empty string"
   echo "Conky is starting..."
   feh --bg-scale `grep 'wallpaper=' ~/.kde/share/config/plasma-desktop-appletsrc | tail --lines=1| sed 's/wallpaper=//'`


   /usr/bin/conky -c $1 &
else
    echo "Killing conky..."
    killall conky  
    sleep 3
    echo "Conky is starting..."
    feh --bg-scale `grep 'wallpaper=' ~/.kde/share/config/plasma-desktop-appletsrc | tail --lines=1| sed 's/wallpaper=//'`

/usr/bin/conky -c $1 &    

  
fi
