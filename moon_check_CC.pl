#!/usr/bin/perl -w
use strict;
use Getopt::Std;


my (%inter_prob,%proteins,%ancestors,%opts,%terms_to_GOs,%offspring,%synonyms,%net_names,%candidates);
getopts('vChn:g:o:a:s:',\%opts);
my $gaf_annotations_file=$opts{a};
my $network_file=$opts{n}||die("Need a network file\n");
my $subonto=$opts{o}||"C";
my $interactome_probs=$opts{i}||'/home/terdon/research/testing/new/data/human.CC.inter.prob';
#my $no_anc=$opts{C}||1;
my $verbose=$opts{v}||undef;
my $synfile=$opts{s}||undef;
my $geneology_file=$opts{g};
my $go_terms_file='/home/terdon/research/testing/new/data/GO.terms_alt_ids';
my $have_already_read_terms_file=0;
my $prob_file_read=0;
$synonyms{LOADED}=0;
while(<>){
    chomp;
    # /^(.+?)\s.+?:\s+(.+?)\s/;
    # my $bait=&get_name($1);
    # my $target=&get_name($2);
    # $net_names{$bait}=$1;
    # $net_names{$target}=$2;
    # $candidates{$bait}{INTERACTORS}{$target}++;
    my $name=&get_name($_);
     $candidates{$name}++;
}
&load_network;
&load_ancestors();
&parse_gaf_file();


bait:foreach my $bait (keys(%candidates)){
    #print "b : $bait\n";
    target:foreach my $target (keys(%{$proteins{$bait}{INTERACTORS}})){
	my $cand=0;
	#print "t : $target\n";
	bgo:foreach my $b_go (keys(%{$proteins{$bait}{GOs}})){
	    	#print "bg : $b_go\n";
	    ## next target if it shares a GO with bait
	    if (defined($proteins{$target}{GOs}{$b_go})){
		$cand=0;
		next target 
	    }
	    tgo:foreach my $t_go (keys(%{$proteins{$target}{GOs}})){
		#print "tg : $t_go\n";
		## If b_go and t_go are related, skip
		if((defined($offspring{$b_go}{$t_go}))||
		   (defined($offspring{$t_go}{$b_go}))){
		    $cand=0;
		    next target;
		}
		my $P=&interactome_probability($b_go,$t_go);
		print "p:$b_go,$t_go $P \n" unless $P==1;
		my @b = sort {$b lt $a} ($b_go,$t_go);
		my $gopair=join("_",@b);
		print "$gopair\n";
		if($P<=0.9999943){
		    $cand=1;
		}
		else{
		    $cand=0;
		    next target;
		}
	    }
	}
	my (@a,@b);
	map{push @a, &terms_to_GOs($_,1)}keys(%{$proteins{$bait}{GOs}});
	map{push @b, &terms_to_GOs($_,1)}keys(%{$proteins{$target}{GOs}});
	if ($cand==1){	
	    print "$net_names{$bait} ";
	    map{print "$_,"}@a;
	    print "\t$net_names{$target} ";
	    map{print ", $_"}@b;
	    print "\n";
	}
    }

}


###############################################################
sub load_network {
    print STDERR "Loading Network...\n" if $verbose;
    open(NET,"$network_file")|| die("Cannot open $network_file : $!\n");
    while(<NET>){
	next if /^\d+$/;
	next if /^!/;
	next if /^\s*$/;
	my ($bait1,$target1)=split(/\s+/,$_);
	my $bait=&get_name($bait1,"load_network");
	my $target=&get_name($target1,"load_network");
	$net_names{$bait}=$bait1;
	$net_names{$target}=$target1;
 	## %proteins will hold all the info on all the 
	## network's proteins
	$proteins{$bait}{INTERACTORS}{$target}++;
	$proteins{$target}{INTERACTORS}{$bait}++;
    }
    close(NET);
}

############################################################
sub parse_gaf_file{    
    print STDERR "Parsing GAF file...\n" if $verbose;
    open(A,"$gaf_annotations_file")|| die("Cannot open $gaf_annotations_file:$!\n");
    while(<A>){
	next if /^!/;
	chomp;
	my @tmparray=split(/\t/,$_);
	next unless $tmparray[8] eq $subonto;
	next if $tmparray[3] !~/^$/;
	my $name=&get_name($tmparray[1],"parse_gaf_file");
	## If this protein is a candidate
	if(defined($proteins{$name})){
	    $proteins{$name}{GOs}{$tmparray[4]}="DIR";
	    ## Inherit ancestor annotations
#	    map{$proteins{$name}{GOs}{$_}="ANC"}@{$ancestors{$tmparray[4]}} unless $no_anc;
	}
    }
    close(A);
}

############################################################
sub load_ancestors{
    open (ANC,"$geneology_file")|| die("cannot open $geneology_file:$!\n");
    while(<ANC>){
	chomp;
	my @terms=split(/\t/);
	my %seen;
	@terms = grep { ! $seen{$_} ++ } @terms;
	my $child=shift(@terms);
	$ancestors{$child}=[@terms];
	map{$offspring{$_}{$child}++}@terms;
    }
    close(ANC);
}


#################################################
sub get_name{
    my $name=$_[0];
    my $old_name=$name;
    my $called_by=$_[1]||"null";
    my $reverse=$_[2]||undef;
    if ($synonyms{LOADED}==0){
	if(-e $synfile){
	    open(S,$synfile);
	    while(<S>){
		if(/^(.+)\t(.+)/){
		    $synonyms{NAME}{$1}=$2;
		}	
	    }
	    close(S);
	}
	else{die("Could not open $synfile\n" );
	 }	
 	$synonyms{LOADED}=1;
    }
    if(defined($synonyms{NAME}{$name})){
	$name=$synonyms{NAME}{$name} ;
    }
    return $name;
    
}
############################################################
sub terms_to_GOs{
    my $term=shift;
    my $mode=shift; ## 0 will return GO:xxx, 1 will return term name
#    unless($term eq 'biological_process'){die("crapiola : -$term-\n");}
    $term=~s/_/ /g;
    if($have_already_read_terms_file==0){
	open(T,"$go_terms_file")|| die("Cannot open terms file : $!\n");
	while(my $line=<T>){
	    next if $line=~/^\!/; 
	    chomp($line);
            ## For some reason, latest GO terms file had
	    ## biological_process instead of biological process.
	    ## So, deal with that:
	    $line=~s/_/ /g; 

	    my @aa=split(/\t/, $line);
	    my @a=@aa;
	    my @terms=($aa[0]);
	    if($line=~/obs$/){pop(@a);}
	    if($aa[1] =~/GO:/){
		push @terms,split(/\s+/,$aa[1]);
	    }
	    if($a[$#a] ne $subonto){
		$terms_to_GOs{TERMS}{$a[$#a-1]}="BAD";
		map{$terms_to_GOs{GOs}{$_}="BAD";}@terms;
	    }
	    else{
		$terms_to_GOs{TERMS}{$a[$#a-1]}=$a[0];
		map{$terms_to_GOs{GOs}{$_}=$a[$#a-1];}@terms;
	    }
	}
	close(T);
	$have_already_read_terms_file=1;
    }

#    &debug("term : $term, id:$terms_to_GOs{TERMS}{$term}, id:$terms_to_GOs{TERMS}{$term} " );
#    print STDERR "term : $term, id:$terms_to_GOs{TERMS}{$term} \n";
    $mode==0 ? 
	return($terms_to_GOs{TERMS}{$term}) :
	return($terms_to_GOs{GOs}{$term}) ;
    
}


############################################################
sub interactome_probability{
    my $go1=shift;
    my $go2=shift;
    my @b = sort {$b lt $a} ($go1,$go2);
    my $gopair=join("_",@b);
    return($inter_prob{$gopair})  if defined($inter_prob{$gopair});
    if (($prob_file_read==1) && (not defined($inter_prob{$gopair}))){
	$inter_prob{$gopair}=1;
	return($inter_prob{$gopair});
    }
    if ($go1 eq $go2){  
	$inter_prob{$gopair}=1;
	return($inter_prob{$gopair});
    }
    ## If one term is an ancestor of the other prob=1;
    if((defined($offspring{$go1}{$go2}))||
       (defined($offspring{$go2}{$go1}))){
	$inter_prob{$gopair}=1; 
	return($inter_prob{$gopair});	
    }
    if ((-e "$interactome_probs") && ($prob_file_read==0)){
	open(A,"$interactome_probs")||die("Cannot open $interactome_probs : $!\n");
	print STDERR "\nReading $interactome_probs...";
	while(<A>){
	    chomp;
	    print STDERR "." if $. % 100000 ==0;
	    my @tt=split(/\t/);
	    my $tempgo=$tt[0];
	    $inter_prob{$tempgo}=$tt[6];
	}
	close(A);
	$prob_file_read=1;
	print STDERR "Done\n";
    }
   
    defined($inter_prob{$gopair}) ? 
	return($inter_prob{$gopair}):
	return(1);
}
