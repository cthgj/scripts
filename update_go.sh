#!/bin/bash 

date=`date`;
echo "======$date=====" >>  /home/terdon/research/GO/latest/log;

cd /home/terdon/research/GO/latest/;
wget http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz



