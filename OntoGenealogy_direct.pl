#!/usr/bin/perl


use Switch;

#######################################################################
# This script will return the first two levels of ancestors for each  #
# GO term. That is the direct ancestors and their direct ancestors    #
#######################################################################

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
    switch ($line)
    {
	case /\[Term\]/              
	{
	    $is_a_term=1;
	}
	case /\[Typedef\]/           
	{
	    $is_a_term=0;
	}  
	case /default-namespace/
	{
	    $namespace=$line;
	    $namespace =~ s/default-namespace: //;
	    @NS=($namespace);
	}
	case /^id/
	{
	    $id=$line;       
	    $id =~ s/^id: //;
	    $namespace{$id}=$namespace;
	    @{$ALT_ID{$id}}=($id);
	}
	case /^alt_id/
	{
	    $alt_id=$line;       
	    $alt_id =~ s/^alt_id: //;
	    push(@{$ALT_ID{$id}},$alt_id);
	}
	case /^name. /
	{
	    $name=$line;
	    $name =~ s/^name: //;
	    $term{$id} = $name;
	}
	case /^namespace/
	{
	    $namespace = $line;
	    $namespace =~ s/^namespace: //;
	    $namespace{$id}=$namespace;
	    if (not grep(/$namespace/,@NS))
	    {
		push(@NS,$namespace);
	    }
	}
	case /^is_a/
	{
#	    $line =~ s/^is_a: //;
	    ($parent_id) = $line =~ /([a-zA-Z_]+:[0-9]+)/;
	    push(@{$children{$parent_id}},$id);
	    push(@{$parents{$id}},$parent_id);
	}
	case /^relationship: part_of/
	{
#	    $line =~ s/^relationship: part_of //;
	    ($parent_id) = $line =~ /([a-zA-Z_]+:[0-9]+)/;
	    push(@{$children{$parent_id}},$id);
	    push(@{$parents{$id}},$parent_id);
	    
	}
	case /is_obsolete: true/
	{
	    $obsolete{$id}=1;
	    $n{$namespace}--;
	}
	case ""
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
    $outfile=$ns.".direct.geneology";
    open($ns."_FH","> $outfile") || die "Cannot open file $outfile";
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
    @geneology = @{$parents{$id}};

    map{push @geneology, @{$parents{$_}} }@{$parents{$id}};
    foreach $i (@{$ALT_ID{$id}})
    {
	PrintOut($ns,($i,@geneology));
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
    die "Usage: ./OntoGeneology <ontology file (obo format)>\n";
}


sub FETCHPARENTS
{
    my $id = shift @_;
    my @Parents = @{$parents{$id}};
    print "$id :: @parents\n";
    if (not $parent_already_visited{$id})
    {
	push(@Parents,$id);
	$parent_already_visited{$id}=1;
    }
    if ($#Parents == -1)
    {
	return @Parents;
    }
    else
    {
	foreach $parent (@Parents)
	{
	    push @Parents, @{$parents{$parent}};
	}
	 return @Parents;
    }
}
	
sub PrintOut
{
    my $ns = shift @_;
    my $FH=$ns."_FH";
    my $el = shift @_;

    print $FH "$el";
    foreach my $el (@_)
    {
	print $FH "\t$el";
    }
    print $FH "\n";
}
