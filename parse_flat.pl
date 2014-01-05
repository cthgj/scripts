#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Term::ANSIColor; 

my (%names,%opts,%acc,%ids);
getopts('hvmco:ia', \%opts);
my $id="";
my $map=$opts{m}||undef;
my $color=$opts{c}||undef;
my $ids=$opts{i}||undef;
my $accessions=$opts{a}||undef;
my $org=$opts{o}||undef;

print STDERR "CAREFULL unfinished script...\n";



### Extract desired organisms from large flat file
if($org){
    my $want=0;
    if($org eq "mus"){$org=10090}
    elsif($org eq "hs"){$org=9606}
    elsif($org eq "dro"){$org=7227}
    elsif($org eq "ele"){$org=6239}
    elsif($org eq "scc"){$org=4932}
    else{die("-o must be a ncbi tax id or one of mus,hs,dro,scc or ele\n") unless $org=~/^\d+$/};
    if($ARGV[0] =~ /\.gz$/) 
    {
	open(F,"zcat $ARGV[0] |")||die("Need a flat file : $!\n");
    } 
    else{
	open(F,"$ARGV[0]")||die("Need a flat file :$!\n");
    }
    my @lines;
    my $o;
    while(<F>){
	my @acc_list=();
	$.==1 && do{
	    die("This is not a flat file\n") unless /^ID/;
	};
	chomp;
	push @lines,$_;
	if(/^ID\s+(.+?)\s/){
	    $id=$1;
	}
	elsif (/^OX.+=(\d+)/){
	    $o=$1;
	    $o == $org ? ($want=1) : ($want=0);
#	    print "aaaaaa $org : $1  :$want\n";
	}
	elsif(/^\/\//){
	    if ($want==1){
		map{print "$_\n"}@lines;
	    }
		@lines=();
		$want=0;
	}
	
	
    }

}













if($map){
    open(A,"$ARGV[0]")||die("Need a list of names\n");
    while(<A>){
	chomp;
	unless (defined($accessions) || defined($ids)){
	    /_/ ? ($ids=1) : ($accessions=1) ;
	}
	#if($accessions){die("File must contain accessions NOT ids. For ids use -i\n") unless /\d/;}
	$names{$_}++;
    }
    close(A);
    open(F,"$ARGV[1]")||die("Need a flat file\n");
    while(<F>){
	my @acc_list=();
	$.==1 && do{
	    die("This is not a flat file\n") unless /^ID/;
	};
	if(/^ID\s+(.+?)\s/){
	    $id=$1;
	}
	elsif (/^AC\s+(.+)/){
	    @acc_list=split(/;\s*/,$1);
#	    print "aa : @acc_list\n";
	#    die("Damn uniprot: $acc_list[0]\n") if defined($acc{$acc_list[0]});
	    map{ push @{$acc{$_}},$acc_list[0] }@acc_list;
	    $acc{$id}=[@acc_list];
	    $ids{$acc_list[0]}=$id;
	    
	}
    }
## If the file of names contains accessions
    if($accessions){
	foreach my $name (keys(%names)){
	    if(defined($acc{$name})){
		if($color){
		    if($name ne ${$acc{$name}}[0] ){
			my $a="$name\t$ids{$acc{$name}[0]}\t@{$acc{$name}}\n";
			$a=~s/($name)/color("bold blue").$name.color("reset")/ige;
			$a=~s/(${$acc{$name}}[0])/color("bold red").${$acc{$name}}[0].color("reset")/ige;
			print "$a";
			
		    }
#die("$name : $ids{$name}:@{$acc{$name}}\n") unless defined($ids{$name});
		    else{print "$name\t$ids{$acc{$name}[0]}\t@{$acc{$name}}\n";}
		}
		else{print "$name\t$ids{$acc{$name}[0]}\t@{$acc{$name}}\n";}
	    }
	    else{print "$name\tdeleted   \t$name\n";}
	}
    }
## If the file of names contains ids. This does exactly the sam as the one before at the moment.
    else{
	foreach my $name (keys(%names)){
	    if(defined($acc{$name})){
		if($color){
		    if($name ne ${$acc{$name}}[0] ){
			my $a="$name\t$ids{$acc{$name}[0]}\t@{$acc{$name}}\n";
			$a=~s/($name)/color("bold blue").$name.color("reset")/ige;
			$a=~s/(${$acc{$name}}[0])/color("bold red").${$acc{$name}}[0].color("reset")/ige;
			print "$a";
			
		    }
#die("$name : $ids{$name}:@{$acc{$name}}\n") unless defined($ids{$name});
		    else{print "$name\t$ids{$acc{$name}[0]}\t@{$acc{$name}}\n";}
		}
		else{print "$name\t$ids{$acc{$name}[0]}\t@{$acc{$name}}\n";}
	    }
	    else{print "$name\tdeleted   \t$name\n";}
	}
    }
}


