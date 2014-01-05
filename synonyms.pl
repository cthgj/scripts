#!/usr/bin/perl
#
use Getopt::Std;
my %opts;
getopts('dihns:',\%opts) || do { print "Invalid option, try 'synonyms.pl -h' for more information\n"; exit(1); };
&usage() && exit(1) if $opts{h};

my $print_dbsym=$opts{d}||undef;
my $print_dbID=$opts{i}||undef;
my $species=$opts{s}||undef;
my $print_species_name=$opts{n}||undef;

print STDERR "No species name specified (-s), using HUMAN\n" if $print_species_name;

my $nl="\n";
if ($^O =~ /x$/i){$nl="\n"}
else{$nl="\r\n"}

my %synonyms;
$synonyms{LOADED}=0;
if(-e  $ARGV[1]){
    open(A,"$ARGV[1]");
    while(<A>){
	next if /^\s+$/;
	next if /^\d+\s*$/;
	chomp;
	s/\r//g;

	if($print_species_name){
	    print "$_\t";
	    my $aa=get_name($_);
	    print "${$aa}[0]\n";  
	}
	else{
	    my $aa=get_name($_);
	    print "${$aa}[0]\t";
	    shift(@{$aa});
	    print "@{$aa}$nl";
	}
    }
}
else{
    my @names=split(/,/,$ARGV[1]);
    map{
	my $aa=get_name($_);
	print "${$aa}[0]\t";
	shift(@{$aa});
	print "@{$aa}$nl";
    }@names;
}

sub get_name{
    my $name=shift;
    my @names=(); 
    if ($synonyms{LOADED}==0){
	open(B,"$ARGV[0]")|| die("Cannot open $ARGV[0]:$!\n");
	while(<B>){
	    my $silivar="!";
	    next if /^$silivar/; ## silly thing cause emacs scres up the syntax highlihting when /!/;
	    my @tmparray=split(/\t/);
	    my @a=split(/\|/,$tmparray[10]);
	    push @a,$tmparray[1], $tmparray[2] ;
	    if($print_dbsym){
		map{
		    $synonyms{$_}[0]=$tmparray[2]; 
		    $synonyms{$tmparray[2]}[0]=$tmparray[2];
		}@a;
	    }
	    elsif($print_dbID){
		map{
		    $synonyms{$_}[0]=$tmparray[1]; 
		    $synonyms{$tmparray[1]}[0]=$tmparray[1];
		}@a;
	    }
	    else{  
		@{$synonyms{$tmparray[2]}}=@a;
		map{
		    $synonyms{$_}=\@a;
		}@a;
	    }
	    $synonyms{LOADED}=1;
	}	
	
    }
    close(B);
    if(defined($synonyms{$name})){
	@names=@{$synonyms{$name}};
	my @aa;
	if($species){
	    map{
		if (/_$species/i){ 
		    push @aa,$_ ;
		}
	    }@{$synonyms{$name}};
	}
	scalar(@aa)>0 ? (unshift(@names,@aa)) : (unshift(@names,$name));
    }
#    else{$names[0]=$name; }
    else{unshift(@names,$name,$name);}
    return \@names;
}

sub usage{

print <<EndOfHelp

USAGE:  
	synonyms.pl [options] <GAF FILE> <INPUT NAMES>

synonyms.pl will take a GAF format file  (e.g. gene_association.goa_human)
and a list of names (one per line) as input and will return the synonyms
found for each name given in the GAF file. Unless you run it with the 
"-d" option, it will return ALL synonyms found in the GAF file.

OPTIONS:

    -h : Print this help
    -d : Only print the "DB Object Symbol", eg ARRB1_HUMAN => ARRB1
    -i : Only print the "DB Object ID", eg ARRB1_HUMAN => P49407
    -n : Print the name passed and the _SPECIES format name it was
         changed to ONLY.
    -s : Give a species name to make the primary protein name be of the 
         format NAME_SPECIES, 
EndOfHelp


}



