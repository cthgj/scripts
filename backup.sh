#!/bin/bash

# This script will backup the desired dir or file to the speciofied location 
# (local or remote). Default is to backup ~/research either to badabing or 
# petitbonum, depending on ip. 

# USAGE: backup.sh user@dest.com:/path/to/backup

EXCLUDES="$HOME/.rsync_excludes"


echo "STARTING : $(date)" 1>&2

function usage () {
    me=`basename $0`;
    echo -e "USAGE:\n\t$me destination source \n\n\tSource can be any directory or file.\n\tDestination can either be a local directory, or a remote location in the following format:\n\n\t\tusername@server:/destination/\ne.g.\n\t$me cchapple@10.1.1.54:/ptitbbackup/cchapple/backup/ /home/terdon/research/ \n\n\tIf source is ommitted, ~/research is assumed, if destination is omitted then either\n\tpetitbonum or badabing is chosen depending on my IP address\n"

echo "EXAMPLES:"
echo -e "backup.sh                 : will backup $HOME/research to badabing or petitbonum depending on ip"
echo -e "backup.sh badabing        : will backup $HOME/research to badabing:/mnt/disk2/backup/research/marseille"
echo -e "backup.sh badabing ./kaka : will backup the file/dir ./kaka to badabing:/mnt/disk2/backup/"
echo "backup.sh cchapple@apolo.crg.es:/users/rg/cchapple/backup ./kaka"


echo "Runs (eg) : rsync -h --progress --stats -r -tgo -p -l -D --update -e ssh $HOME/research cchapple@10.1.1.54:/cshared/cchapple/backup/"
    exit;
}
if [[ $1 =~ "-h" ]]
then
    usage
fi

AGENT="ssh-agent -s"
if [ ! -d $HOME/.ssh/agent ]; then
        mkdir -p $HOME/.ssh/agent
fi
#
# Start an agent if there isn't one running already.
#
pid=`ps -u$LOGNAME | grep ssh-age | awk '{print $1}'`
if [ -z "$pid" ]; then
        $AGENT | grep -v echo > $HOME/.ssh/agent/$HOST & pid=$!
        sleep 1 # Let it fork and stuff
fi

dest='petitbonum';

echo "=================================" >&2
date >&2
echo "=================================" >&2

## If we have at least 1 argument
if ((  $# >= 1 ))
then
    ## If we have exactly 2 args, source is the 2nd
    if (( $# == 2 ))
	then
	dest=$1;
	sources=(${sources[@]} $2);
    ## If we have >2 args, there is something wrong
    elif (( $# > 2 ))
    then
	usage;
    ## If we have one, that is the source
    else
	sources=(${sources[@]} $1);
	echo "aaaa ${sources[@]} :: $dest";
    fi
    
 
## If we have no args, default behaviour
elif (( $# == 0 ))
then
    sources[0]="$HOME/scripts";
    sources[1]="$HOME/research";
    sources[2]="$HOME/doc";
    sources[3]="$HOME/.[a-zA-z]*";
    sources[4]="/var/www";
    sources[5]="/etc";
    sources[6]="/boot";
    echo "Sources are: ${sources[@]}"
fi

ip=`/sbin/ifconfig | grep  'inet addr:' | grep -v 127.0.0.1 | cut -d: -f2| awk '{ print $1}'`;
if [[ $ip =~  '10.1.1' ]]
then
    echo -e "\n\t Backing up to jolitorax \n";
    dest="cchapple@10.1.1.26:/cshared/cchapple/backup";
elif [[ $ip =~ '192.168' ]]
then
    echo -e "\n\tBacking up to badabing\n";
    dest="lacoloc@badabing:/home/lacoloc/backup/daily"
    ssh lacoloc@badabing rm /home/lacoloc/backup/daily/backup_finished.txt
    
else
    echo "unknown network";
    exit;
fi


## Backup!
for i in "${sources[@]}"
do
    if [[ "$dest" =~ badabing ]]
    then
	cmd="rsync -h --progress --stats -zr -tgo -p -l -D  --update -e ssh --exclude-from="$EXCLUDES" $i $dest"
	echo "$cmd"
	eval $cmd
    else
	cmd="rsync -h --progress --stats -r -tgo -p -l -D --update -e ssh $i $dest;"
	echo "$cmd"
	eval $cmd
    fi
done

if [[ "$dest" =~ badabing ]]; then
    ssh lacoloc@badabing touch /home/lacoloc/backup/daily/backup_finished.txt
fi

echo "" >&2;
exit;
