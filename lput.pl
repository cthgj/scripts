#!/usr/bin/perl 
if(!$ARGV[0] || $ARGV[0] eq "help" || $ARGV[0] eq "-h" || $ARGV[0] eq "--help"){
    print "Usage: lput <HOST> <SOURCE FILE> <DESTINATION>\n";
    exit;
}
my ($scpopts,$source,$dest,$host);
my $user="cchapple";

$scpopts="-r";
$host=$ARGV[0];
$source = $ARGV[1];

if($host eq "cumin")      {$host="10.1.1.105"}
elsif($host eq "biproc")  {$scpopts.=" -P 24222"; $host="10.1.1.53"}
elsif($host eq "badabing"){$user="juka";        $host="badabing.dyndns.org"}
elsif($host eq "docpad")  {$user="terdon";      $host="docpad.local"}
elsif($host eq "badabing.local"){$user="juka";        $host="badabing.local"}
elsif($host eq "petit")  {$host="10.1.1.54"}
elsif($host eq "petitbonum")  {$host="10.1.1.54"}
elsif($host eq "goudurix")  {$host="10.1.1.48"}
elsif($host eq "db")  {$host="10.1.3.30"; $user="root"}
elsif($host eq "cleopatre")  {$scpopts.=" -P 24222"; $host="139.124.66.43"; $dest="/cobelix/chapple/";}
#my $source = $ARGV[0] || &usage("USAGE: dget SOURCE DEST\nNeed a source file");
#my $dest = $ARGV[1] || &usage("USAGE: dget SOURCE DEST\nNeed a destination");
#my $host="10.1.1.105";

$dest .= $ARGV[2];

print STDERR "scp $scpopts $source $user\@$host:$dest\n";

system("scp $scpopts $source $user\@$host:$dest");


sub usage{
print "$_[0]\n";
exit(0);
}
