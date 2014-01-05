#!/usr/local/bin/perl -w

#  Basic structure : %query{query_name} = %hash{subject_name} = @list of lists of hsps
#
#
use strict;

my $SCOREMIN=0;
my $SCOREMAX=10000;
my $EVALUE=1;
my $IDENTITY=1;
my $PAT="*"; # aa in the query seq.
my $bit=0;
my $MATCH="*"; 
my $c = 0;
my $counter = 0;
my $score = undef;
my($align,$query_ID,$identities,$letters,$subject_ID,@subjects,$id,$bits,$frame,%HSP,%subj,%hits,%scores,@scores);
my $j=0;

while(<>)
{ 
    
    next if /^\s+$/;
    
### Get query name and make hash{q_name}= hash{subjects}
### 
    if (/^Query=/)
    {
	
	if($query_ID)  ## If this is not the first query
	{
	    $subjects{$subject_ID} = [@sbjcts];
	    $query{$query_ID} = [%subjects]; 
	    @subjects = ();
	}
	/Query=\s(.*)/;
	$query_ID = $1;
	
	next;
    }

### Get query length
    if (/letters\)/)
    {
	/(\d+)/;
	$letters = $1;
	next;
    }

### Get Identities line
    if (/Identities/)
    {
	$identities = $_;
	chomp $identities;
	/^.*?\((\d{1,2}%)\)/;
	$id = $1;  ## Actual identities %age
	next;
    }

### Get subject ID
    if (/^>/)
    {
	$counter=0;
	$subject_ID = $_;
	chomp $subject_ID;
	
	next;	    
    }	
     
### get score line
    if(/Score\s=/)
    {
	
	




	if ($score)
	{
	   
	    #$hits{$query_ID} =  @{ $HSP{$subject_ID} };
	    my $k = $subject_ID . $counter;
	    @{ $HSP{$k} } =  @{ $HSP{$subject_ID} };
	    $hits{$query_ID} =  @{ $HSP{$k} } ;   ### hits{query_id} contains the HSP of current query/subject pair
	  
	    #print STDOUT "a : @{$HSP{$subject_ID}}";  
	}
	$counter++; print STDOUT "counter : $counter\n";
	$score = $_;
	chomp $score;
	/^.*?=\s(.*?)\sbits/;
	$bits = $1;  ### get bit score
	push @scores, $score;
	


	push @{ $HSP{$subject_ID} } , $score;
	push @{ $HSP{$subject_ID} } ,"\n";
	print STDOUT "c : $c\n";



	unless ($c == 0)
	{
	 #   print STDOUT "align : $align\n"; 
	    my @a = split(/\d+/,$align);
	    map {print STDOUT "$_";} @a;
	}	
#	$scores{$query_ID}{$subject_ID} = $score;

	next;
    }
   
### Get frame 

    if (/Frame/)
    {
	/Frame\s=\s(..)/;
	$frame = $1; 
	next;
    }

### get alignment

    if((/^Query:/ || /^Sbjct:/ || (/^\s+/ && $_ !~ /[az\d]/)) && ($_ !~ /Score/))
    { 
	push @sbs, $_;
	
	

#	push @{ $HSP{$subject_ID} } , $_ ;
#	push @{ $HSP{$subject_ID} } ,"\n" if /^Sbjct:/;
#	my $i = $align;
#	$align = &get_lines(' ');	  
#	$align = $align . $i;
#	$c++; 
	 
    }
   
}

my @queries = keys(%hits);   ### @queries : list of query_IDs
my @subs = keys(%HSP);   ### @subjects : list of subject_IDs



sub output
{

    foreach my $q (@queries)
    {
	
	foreach my $s (@subs)
	{ 
	    
	    print STDOUT "#################################################\n\n";
	    print STDOUT "Query= $q\n\n\n";
	    print STDOUT "$s\n\n\n";
#	$subj{$query_ID} = [@subjects];
	    
#	print STDOUT "aaaaa : @{$subj{$query_ID}}\n"; die();
	    foreach my $su (@{$subj{$q}})
	    {
#	    print STDOUT "su : $su\n";
	    }
	    
	    foreach my $joe (@{$HSP{$s}})
	    {
		print STDOUT "$joe";
	    }
	}
    }
    
}

sub get_lines {
	my $spc = $_[0];  ### spc : spacer
	my $tmp;
	# local($tmp);
	while (<>) {
	    
		last if /^\s*$/;
	#	print STDOUT "$_";
	#	chop;
		s/^\s*/$spc/;
		$tmp .= $_;
	    };
	$tmp;
}



#######################################################  NOTES  ########################################################  
# Try getting all HSP lines into one variable and split to print                                                       
#                  
#                  
#                  
#                  
#######################################################  NOTES  ########################################################  
