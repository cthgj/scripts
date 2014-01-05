#!/usr/bin/perl


# use Switch;
use feature "switch";

#------------------------------------------------------------------------
#
# This script computes the genealogy of all terms of an ontology
# You can also compute the genealogy by using a different root
# e.g. use the root "development" rather than "biological process"

#-------------------------------------------------------------------------

if ($#ARGV == -1)
{
    ERROR();
}
else
{
    $ontofile=$ARGV[0];
}
if ($#ARGV == 1) {$newroot = $ARGV[1];}	    


open(ONTOFILE,$ontofile) || die "cannot open $ontofile\n";
while ($line=<ONTOFILE>)
{

    chomp $line;
    
    given ($line)
    {
	
	when ( /\[Term\]/ )             
	{
	    $is_a_term=1;
	}
	when ( /\[Typedef\]/)           
	{
	    $is_a_term=0;
	}  
	when ( /default-namespace: / )
	{
	    $namespace=$line;
	    $namespace =~ s/default-namespace: //;
	    @NS=($namespace);
	}
	when ( /^id: /)
	{
	    $id=$line;       
	    $id =~ s/^id: //;
	    $namespace{$id}=$namespace;
	    @{$ALT_ID{$id}}=($id);
	}
	when ( /^alt_id: /)
	{
	    $alt_id=$line;       
	    $alt_id =~ s/^alt_id: //;
	    push(@{$ALT_ID{$id}},$alt_id);
	}
	when ( /^name:/)
	{
	    $name=$line;
	    $name =~ s/^name: //;
	    $term{$id} = $name;
	}
	when ( /^namespace:/)
	{
	    $namespace = $line;
	    $namespace =~ s/^namespace: //;
	    $namespace{$id}=$namespace;
	    if (not grep(/$namespace/,@NS))
	    {
		push(@NS,$namespace);
	    }
	}
	when ( /^is_a: /)
	{
	    ($parent_id) = $line =~ /([a-zA-Z_]+:[0-9]+)/;
	    push(@{$children{$parent_id}},$id);
	    push(@{$parents{$id}},$parent_id);
	}
	when ( /^relationship:\spart_of/)
	{
#	    $line =~ s/^relationship:\spart_of //;
	    ($parent_id) = $line =~ /([a-zA-Z_]+:[0-9]+)/;
	    ## Some part_of relationships involve different ontologies
	    ## e.g. GO:0042910(MF) part_of GO:0042908(BP). Do not count
	    ## these as parents.
#	    if($namespace{$parent_id} eq $namespace{$id}){
		push(@{$children{$parent_id}},$id);
                push(@{$parents{$id}},$parent_id);
#	    }
	    
	}
	when ( /is_obsolete:\strue/)
	{
	    $obsolete{$id}=1;
	    $n{$namespace}--;
	}
	when ("")
	{
	    if ((not $already{$id}) and ($is_a_term) and (not ($obsolete{$id}))  )
	    {
		$idisterm{$id}=1;
		push(@allids,$id);
		push(@{$terms{$namespace}},$id);
		$n{$namespace}++;
		$already{$id}=1;
	    }
	}
	
    }
}

foreach $ns (@NS)
{
    
    $ontofile=~/.+\.(.+?).obo/;
    my $outfile="$ns.$1.genealogy";

    ## for psi-mi ontology, the default-namespace field is annoying
    ## this fixes it
    if($outfile=~/\//){
	next;
    }
    my $outdir=`dirname $ontofile`;
    chomp($outdir);
    my $outfile="$outdir/$outfile";
    chomp($outfile);
    open($ns."_FH","> $outfile") || die "Cannot open file $outfile :: $ontofile";
}

foreach $id (@allids)
{

#    print STDERR "id=$id\n";
#    print "$id\t$goterm{$id}\n";
#    print "children : ",@{$children{$id}},"\n";

    $name{$id} =~ s/\s/_/g;
    $ns=$namespace{$id};
	
#	    print "id=$id ns=$ns\n";
# 	%leaf_already_visited=();
# 	$n_leaves=NumberLeaves($id);
# 	$prec_leaves = log($n_leaves/$nleaves{$myontology})/log($minleaves);
	
    %parent_already_visited=();
    
    
    @ALLPARENTS=();
    @genealogy = FETCHPARENTS($id);
    foreach $i (@{$ALT_ID{$id}})
    {
	
	PrintOut($ns,($i,@genealogy[1..$#genealogy]));
    }
}


foreach $ns (@NS)
{
    close($ns."_FH");
}


#--------------------------------------------------------------------------------
#    SUBROUTINES
#--------------------------------------------------------------------------------



sub ERROR
{
    die "Usage: ./OntoGenealogy <ontology file (obo format)>\n";
}


sub FETCHPARENTS
{
    my $id = shift @_;
    my @parents = @{$parents{$id}};

    if (not $parent_already_visited{$id})
    {
	push(@ALLPARENTS,$id);
	$parent_already_visited{$id}=1;
    }
    if ($#parents == -1)
    {
	return @ALLPARENTS;
    }
    else
    {
	foreach $parent (@parents)
	{
	    FETCHPARENTS($parent);
	}
	return @ALLPARENTS;
    }
}
	
sub PrintOut
{
    my $ns = shift @_;
    my $FH=$ns."_FH";
    my $el = shift @_;
    print $FH "$el";
#    if($el eq "GO:0008559"){print "b\n";print "aaa $_ " if /GO:0042908/}@genealogy;

    foreach my $el (@_)
    {
	print $FH "\t$el";
    }
    print $FH "\n";
}
