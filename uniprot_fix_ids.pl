#!/usr/bin/perl -w
use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Std;
my %opts;
getopts('hv', \%opts);
my $verbose=$opts{v}||undef;

my $base = 'http://www.uniprot.org';
my $tool = 'batch';
my $agent = LWP::UserAgent->new;
push @{$agent->requests_redirectable}, 'POST';

####################################
# Collect the IDs we want to match #
####################################
my @names;
foreach my $ar(@ARGV){
    if (-e $ar){
	open(A, $ar);
	while(<A>){
	    chomp;
	    push @names,$_;
	}
    }
    else{push @names,$ar};
}

my $n=0;
my $k=scalar(@names);
foreach my $name (@names){
    # open(N,">$name")||die();
    # print N "$name\n";
    # close(N);
    my $response;
    my $try_again=1;
    while($try_again==1){
	
	$response=&get_response($name);
	if ($response->is_success){$try_again=0}
	else {print STDERR "retrying...\n" if $verbose;}
    }
    if($response->is_success){
	$n++;
	printf(STDERR "$n of $k\r") if $verbose;
	if(length($response->content)==0){
	    print  "$name\t$name\tdeleted\n"; 
	}
	else{
	    $response->content=~/(>.+?)\s/;
	    my $a=$1;
	    $a=~s/>..\|//;
	    my @ids=split(/\|/,$a);
	    my $ko=shift(@ids);
	    print STDOUT "$name\t$ko\t@ids\n"; 
	}
    }
    else{
	
	    die 'Failed, got ' . $response->status_line .
	    ' for ' . $response->request->uri . "\n";
    }
    unlink($name);
}
sub get_response{
    my $response = $agent->post("$base/$tool/",
				[ 'id' => $_[0],
				  'format' => 'fasta',
				],
				'Content_Type' => 'form-data');
    while (my $wait = $response->header('Retry-After')) {
	print STDERR "Waiting ($wait)...\n";
	sleep $wait;
	$response = $agent->get($response->base);
    }
    return($response);

}
