#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my (%all_terms,%ancestors,%stats,%gos,%go_count,%opts,%counted);
getopts('pv',\%opts);
my $verbose=$opts{v}||undef;
my $print_to_file=$opts{p}||undef; ## do not run rscript, just print it
my $subonto="P";
my $gaf_file=$ARGV[0]||"gene_association.goa_human";
my $geneology_file=$ARGV[1]||"biological_process.genealogy";
my $gopair_file=$ARGV[2]||"gopairs";
#my $stats_file=$opts{s}||undef;
my $prot_name="Bob";



open(A,"$geneology_file")||die("cannot open $geneology_file: $!\n");
while(<A>){
    chomp;
    my @terms=split(/\t/);
#    map {$all_terms{$_}=1}@terms;
    my $child=shift(@terms);
    @{$ancestors{$child}}=@terms;
}
close(A);
my $silivar="!";
open(A,"$gaf_file")||die("cannot open $gaf_file: $!\n");
while(<A>){
    next if /^$silivar/; ## silly thing cause emacs screws up the syntax highlihting when /!/;
    chomp;
    my @tmparray=split(/\t/);
    next unless $tmparray[8] eq $subonto;
    next unless $tmparray[3]=~/^$/;  ## skip the NOT annotations
    my $nn=$tmparray[2]; ## $nn == protein name
    $gos{$nn}{$tmparray[4]}++; ## $tmparray[4] == GO
}
close(A);
my @prots=keys(%gos);
my $total=scalar(@prots); ## total is the number of annotated proteins
my $c=0;

## Inherit parent terms
foreach my $prot (@prots) {
    $c++;
   # printf(STDERR "[$c, of $total]\r");
    my %protsgos;
    foreach my $goterm (keys(%{$gos{$prot}})){
	$go_count{$goterm}++;
	$protsgos{$goterm}++;
	map{
	    $gos{$prot}{$_}++; 
	    $go_count{$_}++ unless defined($protsgos{$_});
	}@{$ancestors{$goterm}};
    }
## @GOs is now all the terms of $prot
    my @GOs=keys(%{$gos{$prot}});
    #print "gg : @GOs\n";
    for (my $n=0; $n<=$#GOs;$n++){
	my $jojo=$n+1;
	for (my $k=$n+1; $k<scalar(@GOs);$k++){
	    #print "aa : $n:$k " . scalar(@GOs) ." dd\n";
	    my $go1=$GOs[$n];
	    my $go2=$GOs[$k];
	    $go1=~/GO:(\d+)/;
	    my $l=$1;
	    $go2=~/GO:(\d+)/;
	    my $k=$1;
	    my $GOpair;
	    $k>$l ? ($GOpair=$go1 . "xxx" . $go2) : ($GOpair=$go2 . "xxx" . $go1);
	    $stats{$GOpair}++;
	}##end for my $k
    }##end for my $n
}##end foreach prot
print STDERR "\n-------------------------------------------------------------------------\n";

## OK, now deal with all go pairs
open(A,"$gopair_file")||die("cannot open $gopair_file:$!\n");
while(<A>){
    printf(STDERR "$.\r") if $. % 100000==0;
    chomp;
    my $GOpair=$_;
    my ($go1,$go2)=split(/xxx/);
    if(defined($stats{$go1 . "xxx" . $go2})){
	$GOpair=$go1 . "xxx" . $go2;
    }
    elsif(defined($stats{$go2 . "xxx" . $go1})){
	$GOpair=$go2 . "xxx" . $go1;
    }
    else{
	next;	
    }
    $go_count{$go1}=0 unless defined($go_count{$go1});
    $go_count{$go2}=0 unless defined($go_count{$go2});

    if(defined($stats{$GOpair})){
## R stuff
## $stats{$GOpair} : number of times these 2 GOs coincide
## $go_count{$go1} : number of prots with go1 
## $total-$go_count{$go1} : 
## $go_count{$go2} : number of prots with go2




	print STDOUT "hhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = FALSE)\n";
	print STDOUT "lhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = TRUE)\n";
	print STDOUT "hypers=c(\"$go1\",\"$go2\",hhyper,lhyper)\n";
	print STDOUT "write(hypers,file=\"/ptitbbackup/cchapple_temp/hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
    } 
    # if ($c % 100000==0){
    # 	$k++;
    # 	print STDERR "\n Sending R\n";
    # 	system("R CMD BATCH  /ptitbbackup/cchapple_temp/rscript");
    # 	print STDERR "\n$k"."st R sent";
    # 	open(L,">>/ptitbbackup/cchapple_temp/prog_log");
    # 	print L "\n$k"."st R sent\n";
    # }

}
close(A);



# my @GOs=keys(%all_terms);
# $total=scalar(@GOs);
# for (my $n=0; $n<=$#GOs;$n++){
#     my $jojo=$n+1;
#     for (my $k=$n; $k<$#GOs;$k++){
# 	printf(STDERR "[$jojo/$k, of $total]\r");
	
# 	next if $GOs[$n] eq $GOs[$k];
# 	my $go1=$GOs[$n];
# 	my $go2=$GOs[$k]; 
# 	my $GOpair=$go1 . "xxx" . $go2;
# 	if(defined($stats{$go1 . "xxx" . $go2})){
# 	    $GOpair=$go1 . "xxx" . $go2;
# 	}
# 	elsif(defined($stats{$go2 . "xxx" . $go1})){
# 	    $GOpair=$go2 . "xxx" . $go1;
# 	}
# 	else{
# 	    $stats{$go1 . "xxx" . $go2}=$stats{$go2 . "xxx" . $go1}=0;
# 	    $GOpair=$go1 . "xxx" . $go2;	
# 	}
# 	## Skip those terms that are not present in the annotations
# 	unless(defined($go_count{$go1}) and defined( $go_count{$go2})){
# 	    $stats{$GOpair}=undef;
# 	    next;
# 	}
# 	print STDOUT "hhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = FALSE)\n";
# 	print STDOUT "lhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = TRUE)\n";
# 	print STDOUT "hypers=c(\"$go1\",\"$go2\",hhyper,lhyper)\n";
# 	print STDOUT "write(hypers,file=\".hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
# 	$stats{$GOpair}=undef;
	
#     }
#     $go_count{$GOs[$n]}=undef;
#     shift(@GOs);
    
# }

###########################################################3


#  while(<A>){
# 	next if $.<3;
# 	next unless /\|/;
# 	print STDERR "." if $. %10000 == 0 && $verbose;
# 	print STDERR "[$.]\n"  if $. %1000000 == 0 && $verbose;
# 	s/^\s*//;
# 	chomp;
# 	my @tt=split(/\s+\|\s*/);
# 	next if $tt[1] eq $tt[2];
# 	$go_count{$tt[1]}++ unless defined($counted{$prot_name}{$tt[1]});
# 	$go_count{$tt[2]}++ unless defined($counted{$prot_name}{$tt[2]});
# 	$counted{$prot_name}{$tt[1]}=$counted{$prot_name}{$tt[2]}=1;
# 	$prot_name=$tt[0];
# 	my ($go1,$go2)=($tt[1],$tt[2]);
# 	my $GOpair=$tt[1] . "xxx" . $tt[2];
# 	$stats{$GOpair}++;

#     }
#     open(S, ">stats");
#     map{print S "$_ : $stats{$_}\n"}keys(%stats);
#     close(S);

#     open(G,">gocount");
#     map{print G "$_ : $go_count{$_}\n"}keys(%go_count);
#     close(G);
#     close(A);
#     print STDERR "\n" if $verbose;
#     %stats=();
#     %go_count=();
#     exit(0);



# my $c=0;
#     unlink(".rscript") if (-e ".rscript");
#     unlink(".hyper") if (-e ".hyper");

# my @GOs=keys(%go_count);
# my $total=scalar(@GOs);

# if($print_to_file){
#     for (my $n=0; $n<$total;$n++){
# 	my $jojo=$n+1;
# 	for (my $k=$n; $k<$total;$k++){
# 	    printf(STDERR "[$jojo/$k, of $total]\r");
# 	    next if $GOs[$n] eq $GOs[$k];
# 	    my $go1=$GOs[$n];
# 	    my $go2=$GOs[$k];
# 	    my $GOpair=$go1 . "xxx" . $go2;
# 	    if(defined($stats{$go1 . "xxx" . $go2})){
# 		$GOpair=$go1 . "xxx" . $go2;
# 	    }
# 	    elsif(defined($stats{$go2 . "xxx" . $go1})){
# 		$GOpair=$go2 . "xxx" . $go1;
# 	    }
# 	    else{
# 		$stats{$go1 . "xxx" . $go2}=$stats{$go2 . "xxx" . $go1}=0;
# 		$GOpair=$go1 . "xxx" . $go2;	
# 	    }
# 	    print STDOUT "hhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = FALSE)\n";
# 	    print STDOUT "lhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = TRUE)\n";
# 	    print STDOUT "hypers=c(\"$go1\",\"$go2\",hhyper,lhyper)\n";
# 	    print STDOUT "write(hypers,file=\".hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
# 	    $stats{$GOpair}=undef;
# 	}
#     }
# }
# else{
#     for (my $n=0; $n<$total;$n++){
# 	my $jojo=$n+1;
# 	for (my $k=$n; $k<$total;$k++){
# 	    printf(STDERR "[$jojo/$k, of $total]\r");
# 	    next if $GOs[$n] eq $GOs[$k];
# 	    my $go1=$GOs[$n];
# 	    my $go2=$GOs[$k];
# 	    my $GOpair=$go1 . "xxx" . $go2;
# 	    if(defined($stats{$go1 . "xxx" . $go2})){
# 		$GOpair=$go1 . "xxx" . $go2;
# 	    }
# 	    elsif(defined($stats{$go2 . "xxx" . $go1})){
# 		$GOpair=$go2 . "xxx" . $go1;
# 	    }
# 	    else{
# 		$stats{$go1 . "xxx" . $go2}=$stats{$go2 . "xxx" . $go1}=0;
# 		$GOpair=$go1 . "xxx" . $go2;	
# 	    }

# 	    open(R,">>.rscript");
# 	    print R  "hhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = FALSE)\n";
# 	    print R "lhyper=phyper($stats{$GOpair},$go_count{$go1}," . ($total-$go_count{$go1}) .",$go_count{$go2},lower.tail = TRUE)\n";
# 	    print R "hypers=c(\"$go1\",\"$go2\",hhyper,lhyper)\n";
# 	    print R "write(hypers,file=\".hyper\",ncolumns=4,append=TRUE,sep=\"\t\")\n";
# 	    close(R);

# 	    if($n % 1000 == 0){
# 		system("R CMD BATCH  .rscript");
# 		unlink(".rscript");
# 	    }
# 	}## end for my $k;
#     }
#     %stats=();

#     open(A,".hyper")||die("hell\n");
#     while(<A>){
# 	next if /^>/;
# 	print;
#     }
# }
