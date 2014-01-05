#!/usr/bin/perl -w


### This script should be run on output off Add_Cys_pos2queryID.awk or similar scripts. 
### It will return a list of peptides supported by >5 ests and those ests.
### Needs the following lines in the input file
###
### Query= <query_name>_<[*/C]position>
### Sbjct= <subject_name> <[*/C]position> [+/-]


my $query;
my @sbjct;
my %pairs;
my $prevquery = 'hello';
my $subject;
my $line;

my $q;
my $s; 
my $que = 0; ### counter to define first query


my $n = 0;

while(<>)
{

    next if /^\s*$/;   ### skip blank lines
    
    next unless (/Query=/ || /Sbjct=/);  ### keep only subject and query id lines

     &get_line();  ### returns array with "@lines[n] = query_name subject_name"
}

my $a;

my $sca = scalar(@lines);



foreach my $line (@lines)  ### foreach query-subject pair..
{
   
    $line =~ /^(.*?)\s(.*?)$/;
    $query = $1;
    $subject = $2;

    if ($query eq $prevquery)   ### if same as last query
    {
	push @sbjct, $subject;
	next;
    }
    else
    {
	
	$pairs{$prevquery} = [@sbjct] if scalar(@sbjct) >= 5; ### if query supported by >4 subjects asdd to hash
	$prevquery = $query;
	@sbjct = ();   ### set @sbjct = 0 for next set of query-sbjct pairs
	push @sbjct, $subject; ### start over...
	

    }
    next;
}  

my @keys = keys(%pairs);

map {print STDOUT "$_ @{$pairs{$_}}\n";} @keys;

sub get_line
{
    

    if ((/Query=/) && ($que == 0))  ### if this is the first query
    { 
	/=\s(.*?)\s/; ### get query name
	$q=$1;
	$prevquery = $q; 
	$que = 1;
	
    } 
    elsif((/Query=/) && ($que == 1)) ### if this is not the first query
    {
	/=\s(.*?)\s/; ### get query name
	$q=$1;
	
	
    }
    else
    {
	
	/=\s(.*?)\s(.*?)\s(.)/; 
	$s= $1 . '_' . $2 . '_'. $3;  ### concatenate subject name_featureposition_strand

	$line = $q .' ' .$s;  ### concatenate query and sbjct lines into one
	push @lines, $line;
	    
    }
    
    return @lines;
    
}



