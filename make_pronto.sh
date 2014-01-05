#!/usr/bin/env bash

###########################################################
# IMPORTANT: This script assumes passwordless ssh access  #
# to the database server and passwordless access (through #
# a ~/.my.cnf file) for the $dbuser.			  #
###########################################################


########################
# Initialize variables #
########################
## database name
database="pronto_"$(date +%b%d%Y)

## database host
host="10.1.3.30"
## database user
dbuser=root

## date string
date=$(date +%b%d-%Y)

## local data directory
local="local"
## data directory, default is ./data, good idea to make that a 
## link to $local/data
data=$(pwd)/data
## uniprot data directory
uniprot="$data/uniprot"
## Administrator's email
admin=cchapple@tagc.univ-mrs.fr

## logfile
LOG=$(mktemp)
echo $LOG >&2

## List of species names
species=(human fly mouse yeast worm)

#########################################################
# Function to facilitate printing of commands being run #
#########################################################
runthis(){
    ok=0;
    echo "$@" >> $LOG
    eval "$@" 2>> $LOG || ok=1;
    if [ "$ok" = 1 ]; then
    	echo "FAILED: $@" >> $LOG;
    	exit;
    fi
}


####################################################
# Create the local directory if it does not exist  #
# It is recommended that this directory be created #
# manually and made to point to a local drive to   #
# avoid extended I/O operations over the network   #
####################################################
mkdir -p $local

#####################################################
# Create the uniprot directory if it does not exist #
#####################################################
mkdir -p $uniprot

#####################################################
# Create the $data directory if it does not exist   #
#####################################################
mkdir -p $data


###################
# Get annotations #
###################
echo -e "\n\n ANNOTATIONS \n\n"  >> $LOG;
DownloadDate=$(date +'%b %d %Y');
for s in ${species[@]}; do
    dir=$(echo $s | tr '[:lower:]' '[:upper:]')
    rm $data/gene_association.$s $data/gene_association.goa_$s.gz $data/gene_association.goa_$s 2>/dev/null
    runthis "wget -nv ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/$dir/gene_association.goa_$s.gz -O $data/gene_association.goa_$s.gz >/dev/null && gunzip $data/gene_association.goa_$s.gz" 
    ###############################################
    # Build the uniAC<=>uniID map for the species #
    ###############################################
    runthis "gawk -F'\t' 'BEGIN{OFS=\"\t\"}!/!/{print \$2,\$11}' $data/gene_association.goa_$s | sed 's/|.*//' > $data/$s.ac2id"
    runthis "rm $data/$s.map; ln -s  $s.ac2id $data/$s.map"
    runthis "rm $data/gene_association.$s; ln -s gene_association.goa_$s $data/gene_association.$s"
done

echo -e "\n\n ONTOLOGIES \n\n"  >> $LOG;

#################################
# Get the GO.terms_alt_ids file #
#################################
runthis "wget -nv http://www.geneontology.org/doc/GO.terms_alt_ids -O $data/GO.terms_alt_ids"

###########################
# Get the ontology itself #
###########################
runthis "wget -nv http://geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo -O $data/gene_ontology.$(date +%b%d-%Y).obo"
runthis "ln -s gene_ontology.$(date +%b%d-%Y).obo $data/gene_ontology.obo"

######################################################################
# For some reason, the ontology can contain terms that are absent in #
# GO.terms_alt_ids, and that causes the calculate_prob scripts	     #
# problems. Find them and add them to GO.terms_alt_ids.		     #
######################################################################
perl -e 'my %k; open(A,"$ARGV[0]"); 
         while(<A>){
          chomp; @a=split(/\t/); $k{$a[0]}=$a[$#a]; 
          my @b=split(/\s+/,$a[1]); map{next unless /^GO:\d+$/; $k{$_}=$a[$#a];}@b;
         } close(A); 
         my ($go, $o); my %oo=("biological_process"=>"P", "cellular_component"=>"C", "molecular_function"=>"F"); 
         open(B,"$ARGV[1]");
         while(<>){
          chomp; 
          if(/^id:\s*(GO:\d+)/){$go=$1;} 
          elsif(/^name:\s*(.+?)\s*$/){$desc=$1}
          elsif(/^namespace:\s*(.+?)\s*$/){
           $o=$oo{$1}; 
           unless(defined($k{$go})){
            print "$go\t\t$desc\t$o\t\n"
           }
          }
         }' $data/GO.terms_alt_ids $data/gene_ontology.obo >> $data/GO.terms_alt_ids


#######################
# Build the genealogy #
#######################
echo -e "\n\n GENEALOGY \n\n"  >> $LOG;
runthis "OntoGenealogy.pl $data/gene_ontology.$(date +%b%d-%Y).obo"
rm "$data/all_three.genealogy";
runthis "cat $data/molecular_function.$(date +%b%d-%Y).genealogy $data/biological_process.$(date +%b%d-%Y).genealogy $data/cellular_component.$(date +%b%d-%Y).genealogy > $data/all_three_$(date +%b%d-%Y).genealogy && ln -s  all_three_$(date +%b%d-%Y).genealogy  $data/all_three.genealogy"

###################################################
# Get updated uniprot flat files for each species #
###################################################
echo -e "\n\n UPDATE UNIPROT \n\n"  >> $LOG;
runthis "update_uniprot.sh $uniprot"

################################################################
# Build the networks. The script below will query the various  #
# interactome databases, download their interactions and build #
# a new network for each species. This can take a while.       #
################################################################
echo -e "\n\n PSIQCUICK \n\n"  >> $LOG;

for s in ${species[@]}; do
    runthis "nice psicquic_wrapper.pl -FBvl --binary --net-file $s.$date.gr  --psi2id $s.$date.psi2id --flat-file $uniprot/$s.flat --missing-psi $s.$date.missing_psi -p $s.$date.psi -s $s 2>$s.$date.er"
    ####################################
    # Make a link to each new file    #
    ####################################
    rm $s.gr 2>/dev/null; ln -s $s.$date.nr.gr $s.gr
    rm $s.missing_psi;    ln -s $s.$date.missing_psi $s.missing_psi 
    rm $s.psi2id;         ln -s $s.$date.psi2id $s.psi2id;  
    rm $s.psi;            ln -s $s.$date.psi $s.psi;  
done



Create and run the scripts that will be run to generate the probs #

echo -e "\n\n PROBABILITIES \n\n"  >> $LOG;
for s in ${species[@]}; do 
    echo launching $s >> $LOG;
    runthis "make_species_prob_script.sh $s $data > $s.sh; chmod 755 $s.sh"
    runthis "make_species_prob_script.sh $s inter $data > $s.inter.sh; chmod 755 $s.inter.sh" 
done
for s in ${species[@]}; do 
    runthis "./$s.sh && ./$s.inter.sh &"
    error_files=(${error_files[@]} $s.error $s.inter.error);
done

################################################################
# The scripts ($s.sh) created above will print "All done" into #
# $s.error and $s.inter.error. Look for that to check whether  #
# they have finished.					       #
################################################################
while [[ $(cat "${error_files[@]}" | grep -c "All done") -lt "${#error_files[@]}" ]]; do 
    sleep 5; 
done;

###########################################
# Calculate the number of prots annotated #
# to each term.				  #
###########################################
for s in ${species[@]}; do 
    runthis "count_annotations_c1.pl -va $data/gene_association.goa_$s -g $data/all_three.genealogy > $data/$s.counts"
    runthis "count_interaction_gos.pl -va $data/gene_association.goa_\"$s\" -g $data/all_three.genealogy -n $s.gr -m $data/$s.ac2id  > $data/inter_\"$s\".counts"

done

###################################
# Calculate term precision values #
###################################
tt=$(mktemp)
runthis "awk '/^GO/{print \$1}' $data/GO.terms_alt_ids > $tt"
runthis "term_precision.pl -g $data/all_three.genealogy -a $data/GO.terms_alt_ids $tt | sort | uniq > $data/precision.all"




###################################################################################
######################## Setup the PrOnto database ################################
###################################################################################

echo -e "\n\n DATABASE \n\n"  >> $LOG;

###################################################
# This Perl script creates a new, empty database  #
# called $database on server $host		  #
###################################################
runthis "make_prontoDB.pl $database $host"

#################################################
# Disable database indexing to speed up loading #
#################################################
for s in ${species[@]}; do
    runthis "ssh $dbuser@$host myisamchk --silent --keys-used=0  -rq /data/database/$database/$s" 
    runthis "ssh $dbuser@$host myisamchk --silent --keys-used=0  -rq /data/database/$database/inter_$s "
done

##################################################
# Copy the probability flat files to the DB host #
##################################################
files=""
for s in ${species[@]}; do files="$files $data/$s.prob $data/$s.inter.prob"; done
runthis "rsync -zhv $files $dbuser@$host:/data/database/$database/ >> $LOG"

######################################
# Copy the network files to the host #
######################################
files=""
for s in ${species[@]}; do files="$files $s.gr"; done
runthis "rsync -Lzhv $files $dbuser@$host:/var/www_8080/PrOnto/data/ >> $LOG"

############################################
# Load the probabilities into the database #
############################################
for s in ${species[@]}; do 
    runthis "mysql -h 10.1.3.30 -u $dbuser -D $database -e \"LOAD DATA INFILE '$s.prob' INTO TABLE $s IGNORE 1 LINES\"" 
    runthis "mysql -h 10.1.3.30 -u $dbuser -D $database -e \"LOAD DATA INFILE '$s.inter.prob' INTO TABLE inter_$s IGNORE 1 LINES;\""
done

################################
# Rebuild the database indeces #
################################
for s in ${species[@]}; do
    runthis "ssh $dbuser@$host myisamchk --silent -rq /data/database/$database/$s" 
    runthis "ssh $dbuser@$host myisamchk --silent -rq /data/database/$database/inter_$s"
done


############################################
# Populate the precision and counts tables #
############################################
files=""
for s in ${species[@]}; do files="$files $data/$s.counts $data/inter_\"$s\".counts"; done
runthis "rsync -zhv $files $dbuser@$host:/data/database/$database/"
runthis "rsync -zhv $data/precision.all $dbuser@$host:/data/database/$database/"

for s in ${species[@]}; do 
    runthis "mysql -h 10.1.3.30 -u $dbuser -D $database -e \"LOAD DATA INFILE '$s.counts' INTO TABLE counts_$s\"" 
    runthis "mysql -h 10.1.3.30 -u $dbuser -D $database -e \"LOAD DATA INFILE 'inter_$s.counts' INTO TABLE inter_counts_$s;\""
done

    runthis "mysql -h 10.1.3.30 -u $dbuser -D $database -e \"LOAD DATA INFILE 'precision.all' INTO TABLE precision;\""


###############################################
# Make sure everything was inserted correctly #
###############################################
for s in ${species[@]}; do
    lines=$(wc -l "$data/$s.prob" | cut -d ' ' -f 1)
    dblines=$(mysqlshow -h $host -u $dbuser --status pronto_Nov082013 $s | grep Fix | gawk -F'[ |]*' '{print $6}')
    ilines=$(wc -l "$data/$s.inter.prob" | cut -d ' ' -f 1)
    idblines=$(mysqlshow -h $host -u $dbuser --status pronto_Nov082013 inter_$s | grep Fix | gawk -F'[ |]*' '{print $6}')
    ###############################################################
    # The first lines of the input file is he header which is	  #
    # ignored. Therefore, the number of lines in the database's   #
    # table should be one less than the number of lines in the 	  #
    # prob file. If this is not the case, somethinf went wrong.   #
    ###############################################################
    if ! [[ $((lines - dblines)) = 1 &&  $((ilines - idblines)) = 1 ]]; then
	echo "ERROR while updating PrOnto for $s on $date" | sendmail $admin
	echo -e "ERROR for $s:\n Prob file has $lines lines, DB has $dblines\nInterP file has $ilines lines, DB has $idblines" >> $LOG
	exit;
    fi
done;


###################################################
# If all tables were OK in the previous step,	  #
# remove the old database and link to the new one #
###################################################

########################################################
# Link the database "pronto" to the new version	       #
########################################################
runthis "ssh $dbuser@$host rm /data/database/pronto"
runthis "ssh $dbuser@$host ln -s /data/database/$database /data/database/pronto" 

#################################
# Delete the .prob files        #
#################################
runthis "ssh $dbuser@$host rm /data/database/$database/*.prob" 

##################################
# Delete the old database	 #
##################################
runthis "mysql -h 10.1.3.30 -u $dbuser -e 'drop schema $database'"


###################################
# Collect the stats for each onto #
###################################
files="";
ifiles="";
for s in ${species[@]}; do
    # runthis "pronto_stats.pl $data/GO.terms_alt_ids $data/$s.prob > $data/$s.stats"
    # runthis "pronto_stats.pl $data/GO.terms_alt_ids $data/$s.inter.prob > $data/$s.inter.stats"
    files="$files $data/$s.stats";
    ifiles="$ifiles $data/$s.inter.stats";
done
AnnotStats=$(runthis "make_stats_table.pl 'Annotation Probabilities' $files" );
InterStats=$(runthis "make_stats_table.pl 'Interaction Probabilities' $ifiles" );

################################################
# Update the about page with the download date #
################################################
runthis "ssh root@10.1.3.30 \"perl -i -pe 's/<span id=\\\"updated\\\">.+?<\/span>/<span id=\\\"updated\\\"> $DownloadDate<\/span>/' /var/www_8080/PoGo/cgi-bin/pogo.php\""

########################################
# Update the about page with the stats #
########################################
runthis "ssh root@10.1.3.30 \"update_pronto_about.pl '/var/www_8080/PoGo/cgi-bin/pogo.php' '$AnnotStats' '$InterStats' > /tmp/a && mv /tmp/a /var/www_8080/PoGo/cgi-bin/pogo.php\"" 

#############################
# Update the downloads page #
#############################
runthis "ssh root@10.1.3.30 \"update_pronto_netstats.pl '/var/www_8080/PoGo/cgi-bin/pogo.php' > /tmp/a && mv /tmp/a /var/www_8080/PoGo/cgi-bin/pogo.php\"" ;



##################################
# Email the administrator	 #
##################################
echo "PrOnto updated succesfully on $date" | sendmail $admin



#################################
# Restart the MySQL service     #
#################################
ssh $dbuser@$host service mysql restart >>$LOG 2>&1 


