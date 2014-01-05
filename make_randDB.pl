#!/usr/bin/perl -w
use DBI();
use List::Util 'shuffle';
use Getopt::Std;
my ($dbh,$dsn);
my ($icands,$Bcands)=(0,0);
my %opts;
getopts('p:i:r:R:s:n:g:',\%opts);
my $prob_file=$opts{p}||die("Need a prob file, -p\n");
my $inter_prob_file=$opts{i}||die("Need an inter prob file, -1\n");
my $rand_file=$opts{r}||die("Need a rand file name, -r\n");#$ARGV[2] . ".randfile";
my $inter_rand_file=$opts{R}||die("Need an inter rand file name, -R\n");#$ARGV[3] . ".randfile";
my $rand_file_name=`basename $rand_file`;
chomp($rand_file_name);
my $inter_rand_file_name=`basename $inter_rand_file`;
chomp($inter_rand_file_name);
my $database="rand";
my $table=$opts{s}||"human";
my $go_pairs = $opts{g}||die("Need a list of GO pairs of interest, -g\n");
my $inter_table="inter_" . $table;
my $user="root";
print STDERR "Counting rows ($prob_file)...";
my $ROWS=`wc -l $prob_file | gawk '{print \$1}'`;
chomp($ROWS);
print STDERR "$ROWS\nCounting rows ($inter_prob_file)...";
my $inter_ROWS=`wc -l $inter_prob_file | gawk '{print \$1}'`;
my $number_of_pairs=$opts{n}||die("Need a number of pairs to select\n");
chomp($inter_ROWS);
print STDERR "$inter_ROWS\n";

###############
# SUBROUTINES #
###############
sub retry{
    my $c=shift;
    print STDERR "$c\n";
    my $ok=1;
    my $n=0;
    while ($ok!=0) {
    	 $ok=system("$c 1>&2 2>/dev/null");
    	 $n++;
	 $ok!=0 && do {print STDERR "Retrying $c ($n)\n"};
	 sleep 1;
	 if ($ok > 10){
	     system("echo 'Something is stuck, looped $c 10 times' | sendmail 'karolos.chapple\@gmail.com'");
	 }

    }
}

###################################################################
# Make the various prob files that we will be needing so the long #
# R step only needs to be run once				  #
###################################################################


open(R,">./rscript.R")||die "Could not open rscript.R for writing: $!\n";

## See http://stackoverflow.com/q/8273313/1081936 for the sample function
## missing dep on jolitorax, had to compile a newer R version locally
    print R<<EndOf;
df<-read.delim("$prob_file", nrows=$ROWS,header=F, skip=1,colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
pp<-read.table("pairs")
idf<-read.delim("$inter_prob_file", nrows=$inter_ROWS,header=F, skip=1,colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
library(kimisc)
EndOf

for (my $i=1; $i<=100; $i++) {
    my $out_file=$rand_file . $i;
    my $inter_out_file=$inter_rand_file . $i;

print R<<EndOf;
df1<- transform(sample(df,size=nrow(pp)),V1=pp\$V1)
write.table(df1,"$out_file", quote=F, sep='\\t', col=FALSE, row=FALSE)
df2<- transform(sample(idf,size=nrow(pp)),V1=pp\$V1)
write.table(df1,"$inter_out_file", quote=F, sep='\\t', col=FALSE, row=FALSE)
EndOf
    
}
close(R);
#########
# Run R #
#########
my $r=`cat ./rscript.R`;
print STDERR "$r\n------------\n";
print STDERR "~/R-2.15.1/bin/Rscript ./rscript.R:\n";
`nice ~/R-2.15.1/bin/Rscript ./rscript.R 2>rscript.er`;


####################################################
# Now create the database 100 times and run moonGO #
####################################################

for (my $i=1; $i<=100; $i++) {
    print STDERR "---------------- $i ---------------------\n";

    #######################
    # Repopulate database #
    #######################
    my $cmd="ssh $user\@10.1.3.30 \"service mysql stop\"";
    
    #######################################################
    # Sometimes it times out, try again until we connect. #
    #######################################################
    retry($cmd);

    $cmd="ssh $user\@10.1.3.30 \"myisamchk --keys-used=0 -rq /data/database/$database/$table\"";
    retry($cmd);
    $cmd="ssh $user\@10.1.3.30 \"myisamchk --keys-used=0 -rq /data/database/$database/$inter_table\"";
    retry($cmd);
    
    $cmd="ssh $user\@10.1.3.30 \"service mysql start\"";
    retry($cmd);
    
    ## Copy rand file to server
    $cmd="scp $rand_file$i $user\@10.1.3.30:/data/database/$database/$rand_file_name";
    retry($cmd);
    $cmd="scp $inter_rand_file$i $user\@10.1.3.30:/data/database/$database/$rand_file_name";
    retry($cmd);

    $cmd="ssh $user\@10.1.3.30 \"mysqlimport --local --delete -pyenapas $database /data/database/$database/$rand_file_name\"";
    retry($cmd);
    $cmd="ssh $user\@10.1.3.30 \"mysqlimport --local --delete -pyenapas $database /data/database/$database/$inter_rand_file_name\"";
    retry($cmd);
    
    $cmd="ssh $user\@10.1.3.30 \"service mysql stop\"";
    retry($cmd);

    $cmd="ssh $user\@10.1.3.30 \"myisamchk -rq /data/database/$database/$table\"";
    retry($cmd);
    $cmd="ssh $user\@10.1.3.30 \"myisamchk -rq /data/database/$database/$inter_table\"";
    retry($cmd);

    $cmd="ssh $user\@10.1.3.30 service mysql start";
    retry($cmd);
   
    # print STDERR "ssh $user\@10.1.3.30 \"myisamchk --keys-used=0 -rq /data/database/$database/$table\"\n";
    # `ssh $user\@10.1.3.30 "myisamchk --keys-used=0 -rq /data/database/$database/$table"`;
    # print STDERR "ssh $user\@10.1.3.30 \"mysqlimport --local --delete -pyenapas rand /data/database/rand/cchapple/$rand_file\"\n";
    # `ssh $user\@10.1.3.30 "mysqlimport --local --delete -pyenapas rand /data/database/rand/$rand_file"`;
    # print STDERR "ssh $user\@10.1.3.30 \"myisamchk -rq /data/database/$database/$table\"\n";
    # `ssh $user\@10.1.3.30 "myisamchk -rq /data/database/$database/$table"`;
    # print STDERR "ssh $user\@10.1.3.30 service mysql restart\n";
    # `ssh $user\@10.1.3.30 service mysql restart`;
    
    my $rows=`ssh $user\@10.1.3.30 "mysqlshow --status $database -p'yenapas'" | grep -w $table | gawk '{print \$10}'`; 
    chomp($rows);
    my $inter_rows=`ssh $user\@10.1.3.30 "mysqlshow --status $database -p'yenapas'" | grep -w $inter_table | gawk '{print \$10}'`; 
    chomp($inter_rows);
    print STDERR "$table has $rows rows\n$inter_table has $inter_rows\n";
    

    ##############
    # Run MoonGO #
    ##############
     # print STDERR "nice moonGO.pl -A $database  -vm B $table.gr 2>&1 | grep \"candidates\" | tail -n 1 |  gawk '{print \$1}'\n";
     # my $Br=`nice moonGO.pl -A $database -vm B $table.gr 2>&1 | grep "candidates" | tail -n 1 |  gawk '{print \$1}'`;
    my $com="nice moonGO.pl -A $database  -m i $table.gr > human.rand_$i.i.moon && grep \"candidates\"  human.rand_$i.i.moon | tail -n 1 |  gawk '{print \$1}'";
    print STDERR "$com\n";
    `$com`;
    $com="nice moonGO.pl -A $database  -im i -L human.rand_$i.i.moon $table.gr > human.rand_$i.comb.i.moon && gawk '\$1!~/#/{print \$1}' human.rand_$i.comb.i.moon | sort | uniq | wc -l";
    print STDERR "$com\n";
    my $in=`$com`;
    chomp($in);
    $icands+=$in;
     # chomp($Br);
     # $Bcands+=$Br;
#     print STDOUT "$i : Inter.: $in\tBridge: $Br\n";
     print STDOUT "$i : $in\n";
    
}
print "---------\n";
print "Intersection\t" . $icands/100 . "\n";
#print "Bridge\t" . $Bcands/100 . "\n";
exit(0);

    ###################
    # Randomize probs #
    ###################
    # my (%map,%k);
    # open(A, "$prob_file")||die("Could not open $prob_file for reading: $!\n");
    # while (<A>) { ## read 
    # 	chomp;
    # 	if ($.==1 && /^#/) {
    # 	    next;
    # 	}
    # 	print STDERR "*" if $. % 10000000==0;
    # 	/^(.+?)\t(.+)/; 
    # 	$k{$1}=$2;
    # }
    # close(A); 
    # print STDERR "\nShuffling progress: ";
    # my @rn=shuffle(keys(%k)); 
    # my $k=0;
    # foreach (keys(%k)) {
    # 	print STDERR "*" if $k % 10000000==0;
    # 	$map{$_}=$k{$rn[$k]} ;
    # 	$k++;
	
    # }
#    for (my $k=0; $k<=$#nn; $k++) {
    # 	print STDERR "*" if $k % 10000000==0;
    # 	$map{$nn[$k]}=$k{$rn[$k]} 
    # }
    # print STDERR "\nWriting...\n";
    # open(B,">$rand_file")||die("Could not open $rand_file for writing: $!\n");
    # map{print B "$_\t$map{$_}\n"}keys(%map);
    # close(B);
