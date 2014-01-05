#!/usr/bin/perl -w
use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Std;
my %opts;
getopts('vhmf:F:t:', \%opts);
my $format=$opts{f}||'fasta';
my $find_missing=$opts{m}||undef;
my $list = $ARGV[0]; # File containg list of UniProt identifiers.
our $verbose=$opts{v}||undef;
my $base = 'http://www.uniprot.org';
my $tool = 'batch';
my $agent = LWP::UserAgent->new;
push @{$agent->requests_redirectable}, 'POST';

my $response;
my $try_again=1;
while($try_again==1){
    $response=&get_response;
    if(defined($response->is_success)){$try_again=0}
    else{print STDERR "Failed, retrying...\n"; print "-";}
}

if ($find_missing){
    if($response->is_success){
       if(length($response->content)==0){
	   print  "$list\tdeleted\n"; 
       }
    else{
	my @ids=($response->content =~ />.*?\|(.+)\|/g);
	print  "$list\t@ids\n"; 
	#print $response->content;
    }
   }

   else{
	die 'Failed, got ' . $response->status_line .
	' for ' . $response->request->uri . "\n";
   }
}
else{
    $response->is_success ?
	print $response->content :
	die 'Failed, got ' . $response->status_line .
	' for ' . $response->request->uri . "\n";
}


sub get_response{
    my $response = $agent->post("$base/$tool/",
				[ 'file' => [$list],
				  'format' => $format,
				],
				'Content_Type' => 'form-data');
    while (my $wait = $response->header('Retry-After')) {
	print STDERR "Waiting ($wait)...\n";
	sleep $wait;
	$response = $agent->get($response->base);
    }
    return($response);

}
