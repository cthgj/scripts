#!/bin/bash 


PREFIX="";
[ $1 ] && PREFIX=$1;

## Create the PREFIX dir if it does not exist
mkdir -p $PREFIX

species=("Human" "Saccharomyces cerevisiae" "Mus musculus" "Drosophila melanogaster" "Caenorhabditis elegans");
files=(human fungi rodents invertebrates invertebrates);
names=(human yeast mouse fly worm);

date=`date`;
echo "======$date=====" >>  $PREFIX/log;
ip=`/sbin/ifconfig | grep  'inet addr:' | grep -v 127.0.0.1 | cut -d: -f2| awk '{ print $1}'`;

if [[ $1 =~ 'test' ]]
then
for (( i=0;i<5;i++)); 
do
    echo "zcat $PREFIX/uniprot_sprot_${files[$i]}.dat.gz | extract_species_from_flat.pl -s \"${species[$i]}\" > $PREFIX/${names[$i]}.flat"
echo "zcat $PREFIX/uniprot_trembl_${files[$i]}.dat.gz | extract_species_from_flat.pl -s \"${species[$i]}\" >> $PREFIX/${names[$i]}.flat"
done
exit;
fi


## Only run this if we are at the TAGC
if [[ $ip =~  '10.1.1' ]]
then
    cd $PREFIX;
    for (( i=0;i<5;i++)); 
    do
	time0=`stat -c %Y uniprot_trembl_${files[$i]}.dat.gz 2>/dev/null `;
	time1=`stat -c %Y uniprot_sprot_${files[$i]}.dat.gz 2>/dev/null`;
	wget -nv -N ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_${files[$i]}.dat.gz ;
	wget -nv -N ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_${files[$i]}.dat.gz ;
	time2=`stat -c %Y uniprot_trembl_${files[$i]}.dat.gz `;
	time3=`stat -c %Y uniprot_sprot_${files[$i]}.dat.gz `;
	update=0;

	if [ "$time0" != "$time2" ]
	    then
	    update=1;
	fi
	if [ "$time1" != "$time3" ]
	    then
	    update=1;
	fi

	if [ $update == 1 ]
	then
	   echo  "${files[$i]} has changed, updating ${names[$i]}.flat...";
	   zcat uniprot_sprot_${files[$i]}.dat.gz | extract_species_from_flat.pl -s "${species[$i]}" > ${names[$i]}.flat
	   zcat uniprot_trembl_${files[$i]}.dat.gz | extract_species_from_flat.pl -s "${species[$i]}" >> ${names[$i]}.flat;
	   echo "Done";
	   echo "${species[$i]} WAS updated" >>  log
	   
	else
	    echo "${species[$i]} NOT updated" >> log
	fi
    done
    
else
    echo "Sorry, this will only run from the TAGC. This is what it would have done though:";
    for (( i=0;i<5;i++)); 
    do
	echo "wget -nv -N ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_${files[$i]}.dat.gz"
	echo "wget -nv -N ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_${files[$i]}.dat.gz"
	echo "zcat $PREFIX/uniprot_sprot_${files[$i]}.dat.gz | extract_species_from_flat.pl -s \"${species[$i]}\" > ${names[$i]}.flat"
	echo "zcat $PREFIX/uniprot_trembl_${files[$i]}.dat.gz | extract_species_from_flat.pl -s \"${species[$i]}\" >> ${names[$i]}.flat"

    done
fi


