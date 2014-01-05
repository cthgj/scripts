#!/usr/bin/perl -w
use DBI();
use List::Util 'shuffle';
my ($dbh,$dsn);
my ($icands,$Bcands);
my $prob_file=$ARGV[0];
my $rand_file="/root/data/database/rand/randfile";
for (my $i=1; $i<=100; $i++) {

    ###################
    # Randomize probs #
    ###################
    my %k; 
    my %map;
    open(A, "$prob_file");
    while (<A>) {
	chomp;
	/^(.+?)\t(.+)/; 
	$k{$1}=$2;
    }
    close(A); 
    my @rn=shuffle(keys(%k)); 
    my @nn=keys(%k); 
    for (my $infile=0; $i<=$#nn; $i++) {
	$map{$nn[$i]}=$k{$rn[$i]} 
    }
    open(B,">$rand_file");
    map{print B "$_\t$map{$_}\n"}keys(%map);
    close(B);

    ############################
    # Declare database details #
    ############################
    #my $host="10.1.3.30";
    my $host="127.0.0.1";
    my $port=3306;
    my $database="rand";
    my $user="root";
    my $pw="yenapas";
    my $table="$ARGV[2]";    

    $dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
    $dbh=DBI->connect($dsn,$user,$pw);

    ## Truncate table
    my $query="TRUNCATE TABLE $table";
    my $result=$dbh->prepare("$query");       
    $result->execute;      

    ## Repopulate
    `myisamchk --keys-used=0 -rq /data/database/$database/$table`;
    `mysql -u root -p'yenapas' -e "LOAD DATA INFILE '$rand_file' INTO TABLE $table IGNORE 1 LINES; " rand`;
    `myisamchk -rq /data/database/$database/$table`;

    ##############
    # Run MoonGO #
    ##############
    my $Br=`moonGO.rand.pl -vm B human.gr 2>&1 | grep "candidates" | tail -n 1 |  gawk '{print \$1}'`;
    my $in=`moonGO.rand.pl -vm i human.gr 2>&1 | grep "candidates" | tail -n 1 |  gawk '{print \$1}'`;
    $icands+=$in;
    $Bcands+=$Br
    print "$i\t$in\t$Br\n";
}
print "---------\n";
print "Intersection\t" . $icands/100 . "\n";
print "Bridge\t" . $Bcands/100 . "\n";
