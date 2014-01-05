#parse_egosecis.sh

# This script parses SECISearch output of EGO files to return the names of the clusters with >1 and > 2 secis. 

filename=$1

name=${filename%.secis}


gawk -F"-" '/>/{print $1}'< $filename | sed 's/>//' | uniq -cd | sort -nr > clusters_found
perl -ne '/\s(\d{1,2})\s/; if($1>2){print;}' clusters_found > cluster_morethan2

