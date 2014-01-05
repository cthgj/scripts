#!/bin/bash
function usage () {
    me=`basename $0`;
    echo -e "USAGE:\n\t$me [q] script.cyto \n\n\tThe -q flag will cause cytoscape to export the network view as a png in the current folder and exit.\n"
    exit;
}
if [[ $1 =~ "-h" ]]
then
    usage
fi
if (( $# < 1 ))
then
    usage
fi
quit=0;
if [[ $1 =~ "-" ]]
then
    if (( $# < 2 ))
    then
	usage
    fi
    if [[ $1 -ne "-q" ]]
    then
	usage
    fi
    quit=1
    script=$2;
else
    script=$1;
fi

tmp_script=/tmp/`echo $RANDOM`;
cat $script > $tmp_script

if [[ $quit == 1 ]]
then
    png=`echo $script | sed 's/.cyto//'`.png
    echo "network view export file=\"./$png\" type=png zoom=5" >> $tmp_script
    echo "quit" >> $tmp_script
fi

$HOME/Cytoscape_v2.8.3/cytoscape.sh -S $tmp_script &

