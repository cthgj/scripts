#!/bin/bash 

current=`setxkbmap -query | grep layout | fold -s2 | tail -n 1`
if [ "$current" == 'ua' ]
then
    setxkbmap -layout gr
else
    setxkbmap -layout ua
fi
