#!/usr/bin/perl -w
use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Std;
my %opts;
getopts('hv', \%opts);
my $verbose=$opts{v}||undef;
my $list = $ARGV[0]; # File containg list of UniProt identifiers.

my $base = 'http://www.uniprot.org';
my $tool = 'batch';
my $agent = LWP::UserAgent->new;
push @{$agent->requests_redirectable}, 'POST';
my @names;
open(F,"$list")||die("cannot open $list: $!\n");
while(<F>){
    chomp;
    push @names,$_;
}
close(F);
my $query="http://pir15.georgetown.edu/cgi-bin/idmapping_http_client?from=ACC&to=ID&ids=";
my $n=0;
my $k=scalar(@names);
my $c=1;
my $response;
foreach my $name (@names){
    # open(N,">$name")||die();
    # print N "$name\n";
    # close(N);
    
#print "pir15.georgetown.edu/cgi-bin/idmapping_http_client?from=ACC&to=ID&ids="; my $c=0; while(<>){chomp; print "$_+" if $c<500; if ($c==500){print "$_\npir15.georgetown.edu/cgi-bin/idmapping_http_client?from=ACC&to=ID&ids=" ; $c=0;}; $c++} print "\n"'
    $c=0 if $c%500==0;
    if($c==0) {
	my $try_again=1;
	$query=~s/\+$//;
	while($try_again==1){
	    $response=&get_response($query);
	    if ($response->is_success){$try_again=0}
	    else {print STDERR "retrying...\n";}
	}
	if($response->is_success){
	    $n++;
	    printf(STDERR "$n of $k\r") if $verbose;
	    my $out=$response->content;
	    $out=~s/^\s*$//gm;
	    $out=~s/^([^\s]+)\r\n/$1\tdeleted\r\n/gm;
	    print STDOUT "$out\n"; 
	}
	else{
	    die 'Failed, got ' . $response->status_line .
		' for ' . $response->request->uri . "\n";
	}
        $query="http://pir15.georgetown.edu/cgi-bin/idmapping_http_client?from=ACC&to=ID&ids=";
    }
    else{
	$query.="$name+";
    }
    $c++;
}
my $try_again=1;
while($try_again==1){
    $response=&get_response($query);
    if ($response->is_success){$try_again=0}
    else {print STDERR "retrying...\n";}
}
if($response->is_success){
    $n++;
    printf(STDERR "$n of $k\r") if $verbose;
    my $out=$response->content;
    $out=~s/^\s*$//gm;
    $out=~s/^([^\s]+)\r\n/$1\tdeleted\r\n/gm;
    print STDOUT "$out\n"; 
}
else{
    die 'Failed, got ' . $response->status_line .
	' for ' . $response->request->uri . "\n";
}

sub get_response{
    my $url=shift;
    my $response = $agent->get($url);

    return($response);

}
