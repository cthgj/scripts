#!/usr/bin/perl
# Description: This script runs the NetGlycate 1.0.ws0 Web Service. It reads a FASTA file from STDIN and produces predictions in a simple table.
# Author: Edita Bartaseviciute
# Email: edita@cbs.dtu.dk
# Version: 1.0 ws0
# Date: 2009-07-10
# usage: perl netglycate.pl [-g] < example.fsa

use strict;
# include standard XML::Compile helper functions (used to initiate WSDL proxys)
require "xml-compile.pl"; # downloadable from the same site as this script

#taking option from a command line
my $graph = $ARGV[0];

# create proxy to NetGlycate Web Service
my $netglycate = WSDL2proxy ( 'http://www.cbs.dtu.dk/ws/NetGlycate/NetGlycate_1_0_ws0.wsdl' );

# append schema definitions
$netglycate = appendSchemas ( $netglycate , 
	"http://www.cbs.dtu.dk/ws/common/ws_common_1_0b.xsd" ,
	"http://www.cbs.dtu.dk/ws/NetGlycate/ws_netglycate_1_0_ws0.xsd"
);
# create hash of operations from proxy
my %ops = addOperations ( $netglycate ) ;

# Get sequence in fasta format from STDIN
my @fasta;
my $entry = -1;

while (<STDIN>) {
	if (/^>(.*)/) {
		my ($id , $comment) = split (" ",$1);
		$entry++;
		$fasta[$entry]->{id} = $id;
		$fasta[$entry]->{comment} = $comment if defined $comment;
	} elsif (/^([A-Za-z]+)/) {
		$fasta[$entry]->{seq} .= $1;
	}
}

# Create sequence for request
my @sequence;
for ( my $i = 0 ; $i < scalar ( @fasta ) ; $i ++ ) {
	push @sequence , { id => $fasta[$i]->{id} , comment => $fasta[$i]->{comment} , seq => $fasta[$i]->{seq} };
}

# Do the request
my $response;
if ($graph) {	
	$response = $ops{runService}->(
	 parameters => {
	  parameters => {
	   graph => 'required',
	   sequencedata => {sequence => [@sequence]} } });
}
else {
	$response = $ops{runService}->(
	 parameters => {
	  parameters => {
	   sequencedata => {sequence => [@sequence]} } });
}

# uncomment the two following lines to inspect the structure of $response
#use Data::Dumper;
#print Dumper($response);   

#get job id which can be used to get the results later
my $jobid;

if ( ! defined ( $response->{parameters}->{queueentry}) ) {
	die "error obtaining jobid\n";
} else {
	$jobid = $response->{parameters}->{queueentry}->{jobid};
	print STDERR "# waiting for job $jobid";
	my $status = "UNKNOWN";;
	# poll the queue
	while ( $status =~ /ACTIVE|RUNNING|QUEUED|WAITING|PENDING|UNKNOWN/ ) {
		my $response = $ops{pollQueue}->( job => { job  => { jobid => $jobid }   }) ;
		$status = $response->{queueentry}->{queueentry}->{status};
		print STDERR  ". [$status]\r";
	}
	die "\nunexpected job status '$status'\n" unless $status eq "FINISHED";
	print STDERR "\n# job has finished\n";
}
# when the job is done, fetch the result
$response = $ops{fetchResult}->(job => { jobid => $jobid });

# uncomment the two following lines to inspect the structure of $response
#use Data::Dumper;
#print Dumper($response); 

#printing the results (suitable for one sequence)
foreach my $ann (@{$response->{parameters}->{anndata}->{ann}}) {
	my $sequence = $sequence[0]->{seq};
	my $length = length($sequence);
	print ">$ann->{sequence}->{id}\t$length amino acids\n#\n",
	      "# netglycate-1.0 prediction results\n#\n",
	      "# Sequence\t\t\t#    Score\t\tAnswer\n",
   	      "# -----------------------------------------------------------\n";

	foreach my $annrecord (@{$ann->{annrecords}->{annrecord}}) {
		my $pos = sprintf ("%4s", $annrecord->{pos});
		my $score = sprintf ("%6.3f", $annrecord->{score}[0]->{value});
		my $comment;
		if (defined $annrecord->{comment}) {
			$comment = "YES";
		}
		else {
			$comment = ".";
		}
		
		print "# $ann->{sequence}->{id}\t\t     $pos   $score   $annrecord->{feature}     $comment\n";
					
	}	
	print "#\n";
}
#returning the files with graphs if option -g was used
if ($graph){
	foreach my $image (@{$response->{parameters}->{image}}) {
		open ("OUT", ">", $image->{comment});
		my $decoded_image = decode_base64($image->{content});
		print OUT $decoded_image;
		close OUT;
	}
}
