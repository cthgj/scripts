#!/bin/bash 

if [[ $1 =~ "-h"  ||  $1 =~ "help" ]];
then
    cat <<EOF
PASOS: primero copias TODOS los datos a la carpeta hoy. Depspues, en el terminal:
       cd hoy
       juli_raton_patontos.sh
EOF
 exit;
fi

basedir=`pwd`;
for n in $(lowest_dir.pl $basedir);
do 
    cd $n;
    if [ -e todos.txt ]
    then 
	rm todos.txt
    fi
    if [ -e sumas.txt ]
    then 
	rm sumas.txt
    fi
    echo -n "Doing $n...";
    juli_raton.pl *;
    cd $basedir;
    echo "Done";
done
