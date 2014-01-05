#!/usr/bin/env perl 


use DBI();
use 5.10.0;
my $port=3306;
my $host=$ARGV[1];
my $user="root";
my $pw="yenapas";

usage() unless $ARGV[1];

my @date=localtime();
my $day=sprintf("%02d", $date[3]);
my $month=sprintf("%02d", $date[4]+1);
my $year = sprintf("%02d", $date[5] % 100);
my $database=$ARGV[0] || die ("Need a database name as an argument\n");


my $dsn = "DBI:mysql:host=$host;port=$port;"; 
my $dbh=DBI->connect($dsn,$user,$pw);

###########################
# Create the new database #
###########################
$dbh->do("CREATE DATABASE IF NOT EXISTS $database");

#####################
# Create the tables #
#####################
$dbh->do("CREATE TABLE IF NOT EXISTS $database.human (
  gopair CHAR(25) NOT NULL,
  go1 CHAR(12) NOT NULL,
  go2 CHAR(12) NOT NULL,
  tot_onto INTEGER NOT NULL,
  go1_num INTEGER NOT NULL,
  go2_num INTEGER NOT NULL,
  overlap INTEGER NOT NULL,
  P_low DOUBLE NOT NULL,
  P_high DOUBLE NOT NULL,
  low_holm DOUBLE NOT NULL,
  low_hoch DOUBLE NOT NULL,
  low_BH DOUBLE NOT NULL,
  low_BY DOUBLE NOT NULL,
  high_holm DOUBLE NOT NULL,
  high_hoch DOUBLE NOT NULL,
  high_BH DOUBLE NOT NULL,
  high_BY DOUBLE NOT NULL,
  KEY (go1),
  KEY (go2),
  PRIMARY KEY  (gopair)
);");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.mouse like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.fly like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.worm like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.yeast like $database.human;");

#############################
# Create interactome tables #
#############################
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_human like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_mouse like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_fly like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_worm like $database.human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_yeast like $database.human;");

############################################
# Create the count tables. These hold the  #
# number of prots annotated per GO.	   #
############################################
$dbh->do("CREATE TABLE IF NOT EXISTS $database.counts_human (go CHAR(12) NOT NULL, onto CHAR(2) NOT NULL, count INTEGER NOT NULL, PRIMARY KEY  (go, onto));");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.counts_mouse like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.counts_fly like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.counts_worm like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.counts_yeast like $database.counts_human;");

$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_counts_human like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_counts_mouse like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_counts_fly like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_counts_worm like $database.counts_human;");
$dbh->do("CREATE TABLE IF NOT EXISTS $database.inter_counts_yeast like $database.counts_human;");


##################################################
# Create the precision table. It will hold the   #
# precision values of each GO term.		 #
##################################################
$dbh->do("CREATE TABLE IF NOT EXISTS $database.precision (
  go CHAR(12) NOT NULL,
  prec INTEGER NOT NULL,
  PRIMARY KEY  (go))
");



##########################
# Disconnect from the DB #
##########################
$dbh->disconnect();


################################
# Print a simple usage message #
################################
sub usage{
    $0=~/.+\/(.+)/;
    my $name = $1;
    print STDERR <<EndOfHelp;

USAGE:  
    $name [database name] [database server]

	$name will create a new, empty PrOnto database on $host. 
        The script needs a single argument which will be the database's name.
EXAMPLE:

    $name testDB 192.167.12.1

EndOfHelp
exit;
}
