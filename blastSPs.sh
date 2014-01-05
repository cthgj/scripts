#!/bin/bash
## Run with m to do masking

if [ $1 = "m" ]
    then
    mask=1;
    SEQ=$2;
    seqfile=$SEQ.masked
else 
    SEQ=$1
    mask=0
    seqfile=$SEQ
fi
PARAM="/home/ug/cchapple/tetraodon.param.3.No_SECIS.cDNA_full_length.param";
echo $mask
echo $SEQ
DB="/disk2/scratch/cchapple/DBs/Human_26_SPs_2_Cys_1tetra.fa";
outfile="$SEQ"_SPs.out;

if [ $mask = 1 ]
    then
    pb blastall -F F  -e 0.1 -p blastx -i $SEQ -d $DB -o $outfile 
    alignthingie.pl -Svg qa $outfile > "$SEQ"_SPs.aln
    cat $SEQ | gffmask.pl $outfile.gff > $seqfile
fi

test -d gff ||  mkdir gff
test -d fasta || mkdir fasta
test -r fasta/idlist && rm fasta/idlist
test -r $SEQ.geneid && rm $SEQ.geneid
test -r $SEQ.geneid.error && rm $SEQ.geneid.error


SECISearch.pl -vIscGp n -o gff $seqfile > $seqfile.gff
c=1;
for name in $(gawk '{print $1}' $seqfile.gff); 
  do 
  gawk NR==$c $seqfile.gff > gff/$name.gff
  echo $name >> fasta/idlist
  let c=c+1
  
done


cd fasta
retrieveseqs.pl -is ../$seqfile idlist 
cd ../

echo $PWD

for name in $(/bin/ls -1 fasta/*fa | sed 's/.fa//' | sed 's/fasta\///');
  do
echo "running geneid...parameter : $PARAM"

~scaste/geneid/bin/geneid -vP $PARAM -R gff/$name.gff fasta/$name.fa >> $SEQ.geneid 2>> $SEQ.geneid.error;

done
