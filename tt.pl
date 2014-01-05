#!/usr/bin/env perl 
use DBI();
my $db_host="127.0.0.1";
my $db_port=3307;

my $host=$db_host;
my $TESTS=93257;
my $user="root";
my $pw="yenapas";
my $table="human";    
my $database="pogo";
my $port=$db_port;

if ($ARGV[0]) {
    $port=3306;
    $host="10.1.3.30";
}

$dsn = "DBI:mysql:database=$database;host=$host;port=$port;"; 
##Select database 
$dbh=DBI->connect($dsn,$user,$pw);

my $gopair="GO:0016070_GO:0050896";
my $query="SELECT P_low from $table where gopair='$gopair'";  
	my $result=$dbh->prepare("$query");                        
	$result->execute;                                          
	my $ref = $result->fetchrow_hashref();
print "Unco: ",$ref->{'P_low'},"\nCorr: ";
print $ref->{'P_low'} * $TESTS,"\n";
