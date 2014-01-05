#!/usr/bin/perl 
if($ARGV[0] eq "help" || $ARGV[0] eq "-h" || $ARGV[0] eq "--help"){
    print "Usage: lget [scp options ]<HOST> <SOURCE FILE> <DESTINATION>\n";
    exit;
}
my ($scpopts,$source,$dest,$host);
my $user="cchapple";
if($ARGV[0]=~/^-/){
    $scpopts=$ARGV[0];
    $host=$ARGV[1];
    $source = $ARGV[2];
    $dest = $ARGV[3];
}
else{
    $scpopts=undef;
    $host=$ARGV[0];
    $source = $ARGV[1];
    $dest = $ARGV[2];
}


if($host eq "cumin")      {$host="10.1.1.105"}
elsif($host eq "biproc")  {$scpopts="-P 24222" unless defined $scpopts ; $host="10.1.1.53"}
elsif($host eq "badabing"){$user="juka";        $host="badabing.dyndns.org"}
elsif($host eq "badabing.local"){$user="juka";        $host="badabing.local"}
elsif($host eq "docpad")  {$user="terdon";      $host="docpad.local"}
elsif($host eq "petit")  {$host="10.1.1.54"}
elsif($host eq "petitbonum")  {$host="10.1.1.54"}
elsif($host eq "cleopatre")  {$scpopts.=" -P 24222"; $host="139.124.66.43"}
#my $source = $ARGV[0] || &usage("USAGE: dget SOURCE DEST\nNeed a source file");
#my $dest = $ARGV[1] || &usage("USAGE: dget SOURCE DEST\nNeed a destination");
#my $host="10.1.1.105";

print STDERR "scp $scpopts $user\@$host:$source $dest\n";

system("scp $scpopts $user\@$host:$source $dest");

