#!/bin/bash

ma="m";
ka="k";
name="cchapple";
if [ $1 = "$ma" ]; then
    name="mmariotti";
fi
echo $name;

ssh $name@apolo.crg.es