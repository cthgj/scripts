#!/usr/bin/perl -w

### This script will take a file and add to each ID the start nt
### and strand of that ID's element. Scirpt must be modified depending on input file

die("Need a file as input\n") unless $ARGV[0];



&secis if ($ARGV[0] =~ /\.secis/);
&gff   if ($ARGV[0] =~ /\.gff/);


sub secis
{
    open(FILE,"$ARGV[0]");
    while(<FILE>)
    {
	if (/^>/)
	{
	    /^>(.*?)\s\[(\d+)\s-(.*?)$/;
	    my $id = $1 ._ . $2;
	    my $start = $2;
	    chomp $start;
	    my $rest = $3;
	    if (/strand/)
	    {
		$strand = '-'
		}
	    else
	    {
		$strand = '+';
	    }		
	    print STDOUT ">$id" ."_$strand [$start -$rest\n";
	    
	}
	else
	{
	    print STDOUT "$_";
	}
	
    }
}

sub gff
{
    open(FILE,"$ARGV[0]");
    while(<FILE>)
    {
	/^.*?\t(.*?)$/;
	my $rest = $1;	    
	/^(.*?)\t.*?\t.*?\t(\d*?)\t(\d*?)\t.*?\t(.)/;
	my $id = $1;
	my $start = $2;
	my $end = $3;
	my $strand = $4;
	print STDOUT "$id" ."_$start" ."_$strand\t$rest\n" if $start < $end;
	print STDOUT "$id" ."_$end" ."_$strand\t$rest\n" if $start > $end;
    }
}

	
