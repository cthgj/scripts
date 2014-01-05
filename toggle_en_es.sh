#!/bin/bash 

current=`setxkbmap -query | grep layout | fold -s2 | tail -n 1`
if [ "$current" == 'us' ]
then
    setxkbmap -layout es
else
    setxkbmap -layout us
fi
