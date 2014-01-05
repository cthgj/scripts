#!/usr/bin/perl
$|=1;
my $VERSION = "1.02"; my $AUTHOR = "Prilusky"; my $YEAR = "2003";
use LWP::Simple;

# simple example to retrieve only the foldindex value

my $c=0;
open(A, "$ARGV[0]");
my $id;
my $seq;

while(<A>){
    chomp;
    next if /^\s*$/;
    if(/>(.+?)\s/){
	$c>0 && do {
	    print STDERR "$c\r";
	    my $index=&getFoldIndex($seq);
	    print "$id\t$index\n";	 
	};
	$id=$1;
	$seq="";

	$c++;
    }
    else {
	$seq.=$_;
    }    
}
print STDERR "$c \r";
my $index=&getFoldIndex($seq);
print "$id\t$index\n";	 

sub getFoldIndex {
  my ($aa) = @_;
  $aa =~ s/\W//g;
  my $content = get("http://bioportal.weizmann.ac.il/fldbin/findex?m=xml&sq=$aa");
  my ($findex) = $content =~ /<findex>([\-\.\d]+)<\/findex>/;
  return $findex;
}
