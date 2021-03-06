#!/usr/bin/perl
#
# 
# $Id: anysearch.pl,v 1.2 2003/08/27 13:31:41 cchapple Exp cchapple $
#

use strict;
use Getopt::Std;




#################################################################
##           Define locations of external programs             ##
#################################################################
my $scan1GB = "/usr/local/molbio/bin/scan_for_matches";  # "/home/ug/cchapple/bin/scan_1GB";
my $RNAfold = "/usr/local/molbio/bin/RNAfold"; #/home/ug/cchapple/bin/RNAfold13"; 
my $ghostvw = "/usr/bin/gs";
my $bdir    = "/home/ug/cchapple/share/SECISearch";


#################################################################
##     Get cmd-line options and declare some variables         ##
#################################################################

my %opts;
getopts('r:p:f:e:E:o:GRHITPYOBShvtlgsndci',\%opts) || do { &usage; exit(1); };

### Input Options

my $file    = $opts{f};         # sequence file
my $pat     = $opts{p};         # -p specifies the pattern for scan_for_matches
my %e1      = (pat_standard_I_II => -10, pat_non_standard_I_II => -10, pat_twilight_I_II => -10 ,pat_Sep20 => -3.7,  pat_dm2 => -3.7, pat_dm => -3.7, pat_s => -8.5, pat_c => -3.7);
my %e2      = (pat_standard_I_II => -20, pat_non_standard_I_II => -20, pat_twilight_I_II => -20, pat_Sep20 => -12.6, pat_dm2 => -12.6, pat_dm => -12.6, pat_s => -15, pat_c => -12.6);




### Output Options

my $I       = $opts{I} || 0;    # -I : bypass the imager2 function
my $P       = $opts{P} || 0;    # -A : bypass scan_for_matches 
my $H       = $opts{H} || 0;    # -H : produce HTML output
my $v       = $opts{v} || 0;    # Verbose
my $t       = $opts{t} || 0;    # Return <FILENAME>.hits
my $en      = $opts{n} || 0;    # Don't unlink energy param file
my $T       = $opts{T} || 0;    # Print no of SECIS found
my $gen     = $opts{G} || 0;    # Program runs in general mode, not SECISearch.
my $g       = $opts{g} || 0;    # Produce gff output
my $s       = $opts{s} || 0;    # Do not return <FILENAME>.secis
my $R       = $opts{R} || 0;    # Don't run RNAfold 
my $rules   = $opts{r} || 0;    # Get rules to generate constraints
my $output  = $opts{o} || 'fs'; # what is printed at the STDOUT

my $h       = $opts{h} || 0;    # cry for help
my $debug   = $opts{d} || 0;    #run with debugging info

### Scanning Options

my $c;              # -c scan_for_matches will scan complementary strand
my $Y_flag  = 1;    # S,Y,O,B are filters for specific SECIS features
my $O_flag  = 1;    #
my $B_flag  = 1;    #
my $S_flag  = 1;    # 

# NOTE: See sub get_params for the rest




#################################################################
##               Define some global variables                  ##
#################################################################

my $dir = $ENV{PWD} || './';   # get name of current directory
my $hits;
my $filename;
my $patname;
my $pngdir;
my $e1;
my $e2;
my $totalhits = 0;
my $secis;
my $mode;

my $units;
my $unit8;
my $lastunit;
my $unit7;

### anysearch variables

my %hash;
my @rule;
my @num;

####################################################################
##                         Define Variables                       ##
####################################################################

my ($alpha, $beta, $ci, $cmd, $curpair, $currentseq, $currentstruct, $energy, $fontsize,
    $fonttype, $fullstructen, $gam, $hn, $index, $j, $k, $n, $name, $offx, $offy, 
    $outfilename, $pi, $posbase1, $posbase2, $r, $rotangle, $scale, $start, $stop, 
    $structmod, $substruct, $upstemen, $x1, $x2, $xcenter, $xf, $xstart, $xstop, 
    $y1, $y2, $ycenter, $yf, $ystart, $ystop, $i);

my (@_a, @bold, @circle, @comment, @energy, @fullstructen, @left, @name, @nonbold, @pairof,
    @parts, @png, @seq, @start, @stop, @texta, @upstemen, @x, @xf, @xlocal, @xnew,
    @y, @yf, @ylocal, @ynew);

 
my @a = ();
my @b = ();

##############################    

    
 



###########################################################################################
###############                   Core program START                      #################
###########################################################################################

&help();  
                      &debug("help finished\n");
&get_params();
                      &debug("get_param finished\n");
&get_rules();
                      &debug("get_rules finished\n");
&prepare_ground();
                      &debug("prepare_ground finished\n");
&scan_for_matches();
                      &debug("scan_for_matches finished\n");
&skip_rnafold();
                      &debug("skip_rnafold finished\n");

if ($mode eq 'secis') 
{
    &thermodynamic_eval();
                      &debug("thermodynamic_eval finished\n");
}
elsif ($mode eq 'gen')
{
    &thermodynamic_eval_gen();
                      &debug("thermodynamic_eval_gen finished\n");
}
###########################################################################################
###############                      Core program END                     #################
###########################################################################################

                       &debug("main loop finished");
&secis_output();
                       &debug("secis_output finished\n");
&patscan_hits();
                       &debug("patscan_hits finished\n");
&standard_out();
                       &debug("standard_out finished\n");
&HTML_output();
                       &debug("HTML_output finished\n");
&clean_up();    
                       &debug("clean_up finished\n");
exit(0);
    



###########################################################################################
###############                         Program END                       #################
###########################################################################################






############################################################################################
################                        FUNCTIONS                         ##################
############################################################################################



    
   
####################################################################
##        Collect necessary variables and parameters              ##
####################################################################

sub get_params
{

### Determine which mode the program is running in
   
    
   
    if ($rules || $gen) {$gen = 1; $secis = 0}
    else {$secis = 1; $gen = 0;}
    
### Set mode-dependant scanning options

    if ($opts{Y} || $gen){$Y_flag = 0}
    if ($opts{O} || $gen){$O_flag = 0}
    if ($opts{B} || $gen){$B_flag = 0}
    if ($opts{S} || $gen){$S_flag = 0}


### Allow input seq from STDIN
    
    if (!$opts{f})
    {
	my @inputseq = <STDIN>;
	open (INFILE,">$dir/inseq.fa");
	print INFILE "@inputseq";
	close(INFILE);
	
	$file = "inseq.fa";
    }
    
### Get user defined pattern or use default one
    
    if ($opts{p})
    {
	
	if ($opts{p} eq 'n')               ### non-standard
	{
	    $pat = "$bdir/pat_non_standard_I_II"; 
	}
	
	elsif ($opts{p} eq 's')            ### standard
	{
	    $pat = "$bdir/pat_standard_I_II"; 
	}
	elsif ($opts{p} eq 't')           ### The twilight zone
	{
	    $pat = "$bdir/pat_twilight_I_II"; 
	}
	elsif (-r "$pat")
	{
	    $pat = $pat;
	    $gen = 1;
	}
	elsif ($opts{p} !~ /\//)
	{
	    $pat = "$bdir/$opts{p}";
	}
	elsif ($opts{p} =~ /\//)
	{
	    $pat = "$opts{p}";
	    $gen = 1;
	}
    }
    else 
    {
	print STDERR "*** no pattern specified, or specified pattern not found, using default pattern\n";
	$pat = "$bdir/pat_standard_I_II";
    }
    
    
###  Juggle variable names around to accomodate directory structure   
    
    
    if ($file =~ m/\//)          #m/\/((\w*?)|(\w*?\..*?))$/)
    {
	$file =~ m/\/([^\/]*)$/g;   #m/(\w*?\.*\w*?\.*)$/g;
	$filename = $1;
    }
    else 
    {
	$filename = $file;
    }
   
    if ($pat =~ m/\/((\w*)|(\w*\.\w*))$/g)
    {
	$patname = $1;
    }
    else 
    {
	$patname = $pat;
    }
    
### Allow for user-defined energy cutoff values     
    
    
    
    if ($opts{e}){$e1 = $opts{e};} else {$e1 = $e1{$patname} || -10}   
    if ($opts{E}){$e2 = $opts{E};} else {$e2 = $e2{$patname} || -20} 
    
    
###  Have scan_for_matches search complementary strand      

    if ($opts{c}){$c = '-c';}    

### Create png directory
    if ($I)
    {
	mkdir("$dir/$filename\.$patname\_png", 0777) unless (-d "$dir/$filename\.$patname\_png");
	$pngdir = "$dir/$filename\.$patname\_png";
    }
###
    


    &debug("##################");
    if ($gen){$mode = 'general'}
    elsif($secis){$mode = 'secis'}
    &debug("Mode is $mode");
    &debug("\$dir is $dir");
    &debug( "\$filename is $filename");
    &debug( "\$file is $file");
    &debug( "\$e1 is $e1");
    &debug( "\$e2 is $e2");
    &debug( "\$pat is $pat");
    &debug( "\$patname is $patname");
    &debug( "\$rules is $rules");
    &debug( "\$H is $H");
#    &debug( "\$ is $");
#    &debug( "\$ is $");
    &debug("##################\n");
    
}

##################################################################
##      Get pattern rules to generate RNAfold constraints       ##
##################################################################

sub get_rules
{
    if ($rules)
    {
	open (RULES, "$rules") || die "cannot open $rules: $!";
	while (<RULES>)
	{
	    
	    m/(\d{1,2})\s*(\(|\)|\<|\>|x|\.)$/og;
	    my $num  = $1;
	    my $rule = $2;
	    $hash{$num} = $rule;
	    push @num, $num;
	    push @rule, $rule;
	    $rule=0;
	    $num=0;
	    
	}

	close(RULES);
    }
}


##################################################################
##        Delete files from a possible previous instance         ##
##################################################################


sub prepare_ground
{
    if (($g) || ($output eq 'gff')) {unlink("$dir/$filename\.$patname\_gff")}
    unlink("$dir/$filename\.$patname\_en");
}


##################################################################
##    Running scan_for_matches with choosen pattern             ##
##################################################################


sub scan_for_matches
{
   if ($P) 
    {
	if ($v || $debug)
	{
	    print STDERR "*** ignoring scan_for_matches\n";
	}
	$hits = "$file";
    }
    else
    {      

	if ($v || $debug)
	{
	    print STDERR "*** running scan_for_matches, please wait...\n";
	}
	&debug("### PatScan Command: $scan1GB $c $pat < $file > $dir/$filename\.$patname\_hits ###\n");
	system("$scan1GB $c $pat < $file > $dir/$filename\.$patname\_hits") || die "cannot run scan_for_matches: $!, perhaps the input file is empty...";
	
	$hits = "$dir/$filename\.$patname\_hits";
    }
}    

####################################################################
##                        Ignore RNAfold                          ##
####################################################################

sub skip_rnafold
{
    if ($R)
    {
	open (PWP, "$dir/$file\.$patname\_hits");
	while (<PWP>)
	{
	    if (/^>/)
	    {
		print; 
	    }
	    else
	    {
		s/(\w)\s(\w)/$1$2/g;
		print;
	    }
	}
	
	if ($v)
	{
	    print STDERR "\n*** WARNING: These sequences have NOT been evaluated thermodynamically\n\n";
	}
	exit();
    }
}


######################################################################################################
##                           Evaluate hits thermodynamically (SECIS)                                ##
######################################################################################################

sub thermodynamic_eval
{
    if ($v || $debug)
    {
	print STDERR "*** scan_for_matches finished, evaluating hits thermodynamically...\n";
    }
    
    unless ($I)  ### Make png dir
    {
	mkdir("$dir/$filename\.$patname\_png", 0777) unless (-d "$dir/$filename\.$patname\_png");
	$pngdir = "$dir/$filename\.$patname\_png";
    }
    
    my $seqname;
    open (HITS, "$hits") || die("cannot open $hits: $!");

    my @hitcount = <HITS>;
    close(HITS);
    my $hitcount = @hitcount;
    if ($hitcount == 0)
    {
	die("\n*** Input file is empty!!!\n") if $P;
	die("\n*** No matching sequences found.\n") if !$P;
    }
    open (HITS, "$hits") || die("cannot open $hits: $!");

    while (<HITS>)    
    {    	
	&debug("#$.# $_");
	my @param;
	next if /^\s*$/o;
	/^>(.*)/o && do 
	{
	    $seqname = $1;
	    ($name, $start, $stop) = split /\:\[|\,|\]/, $_;
	    if ($name =~ /^>(.*?)$/)
	    {
		$name = $1;
	    }
	    next;
	};
	
	@b=();
	tr/Tt/Uu/ unless ~ /^\>/;



	@a=split /\s+/o, $_;

	s/\w/./g for @b=@a; 

	$units = @a;
	$lastunit = $units -1;
   
###
	&debug("************************\n");
	&debug("units: $units\n");
	&debug("\n@a\n");
   	&debug("\@a[$units - 4]\n");
	&debug("************************\n");
###

### Get unit numbers. May change depending on pattern
	    
	    if ($units == 12)
	    {
		$unit8 = 8;
		$unit7 = 7;
	    }
	    else # ($units == (16 || 15 || 14))
	    {
		$unit8 = $units -4;
		$unit7 = $units -5;
	    }
###

	    $b[3]='x(((('; $b[5]='xx';
	    if (length ($a[$unit8]) == 4)  { $b[$unit8]='))))' }
	    elsif (length ($a[$unit8]) == 5)  { $b[$unit8]='.))))'}
	    else { die "When not running patscan, the input file must be the output of a previous patscan session. 
Otherwise, SECISearch encounters problems with the length of the sequence.\nlength($a[$unit8])\n$units -4\n" }
	    
	$substruct = join("",@b);

### Pep's usefull debuggging options

	&debug($substruct);
	my $o=scalar( $substruct =~ s/\(/\</og );
	my $e=scalar( $substruct =~ s/\)/\>/og );
	&debug("$substruct :: $o == $e ? ".($o == $e ? "YES":"NO")."\n");
###	

	if (($upstemen = &secis_energy(3..$unit8)) <$e1 && ($fullstructen = &secis_energy(0..$lastunit))<$e2 && structure_OK() && !$R) 
	{ 
	  
	    $totalhits++;
	    open (EN ,">>$dir/$filename\.$patname\_en")|| die "can't open $dir/$filename\.$patname\_en: $!"; 
	    
	    if ($start>$stop)
	    {
		print EN ">$seqname [$start - $stop] (SECIS on complementary strand) - Stem: $upstemen   Whole structure: $fullstructen\n$currentseq\n$currentstruct\n";
	    }
	    else
	    {
		print EN ">$seqname [$start - $stop] - Stem: $upstemen   Whole structure: $fullstructen\n$currentseq\n$currentstruct\n";
	    }
	    
	    close(EN);
	    
	    push @param, $currentseq;
	    push @param, $currentstruct;
	    push @param, 15;
	    push @param, "\"$seqname".".png\"";
	    push @param, 170;
	    push @param, 250;
	    push @param, length (join "", @a[0..2])+1;
	    push @param, length (join "", @a[0..2])+4;
	    push @param, (length (join "", @a[0..2]))..((length (join "", @a[0..2]))+4);
	    push @param, (length (join "", @a[0..4]))..((length (join "", @a[0..4]))+1);
	    if (length ($a[$unit8]) == 4) 
	    {
		push @param, (length (join "", @a[0..$unit7]))..((length (join "", @a[0..$unit7]))+3);
	    }
	    else 
	    {
		push @param, (length (join "", @a[0..$unit7])+1)..((length (join "", @a[0..$unit7]))+4); 
	    }
	    &gff_output();
	    
### In the absence of I, run through the imager2 function      
	    
	    unless($I)
	    {
		imager2(@param);
	    }
	    
	    push @texta, [ @a ];
	    push @name, $name; 
	    push @start, $start;
	    push @stop, $stop;
	    push @png, "$seqname".".png";
	    push @comment, $start>$stop ? "(SECIS on complementary strand)" : "";
	    push @upstemen, $upstemen;
	    push @fullstructen, $fullstructen;
	    
	    
###display progress
	    
	    if ($v)
	    {
		print STDERR "#"; 
	    }
	    
	    
	    &debug("*");    
	}  &debug(":");
    }  &debug(".");
   
    
    ### Die if no stable structures
    
    if ($totalhits == 0)
    {
	print STDERR "\n*** Sorry, no stable structures were found. You may want to try a different energy threshold.\n";
	exit();
    }
    
}

######################################################################################################
##                           Evaluate hits thermodynamically, general mode                          ##
######################################################################################################

sub thermodynamic_eval_gen
{
    if ($v || $debug)
    {
	print STDERR "*** scan_for_matches finished, evaluating hits thermodynamically...\n";
    }

    my $seqname;
    open (HITS, "$hits");
    while (<HITS>)    
    { 
	&debug("#$.# $_");
	my @param;
	next if /^\s*$/o;
	/^>(.*):/o && do 
	{
	    $seqname = $1;
	    ($name, $start, $stop) = split /\:\[|\,|\]/, $_;
	    if ($name =~ /^>(.*?)$/)
	    {
		$name = $1;
	    }
	    next;
	};
	
	@b=();
	tr/Tt/Uu/;
	@a=split /\s+/o, $_;
	s/\w/./g for @b=@a; 

	$units = @a;
	$lastunit = $units -1;
   
###
#	&debug("************************\n");
#	&debug("units: $units\n");
#	&debug("\n@a\n");
#   	&debug("\@a[$units - 4]\n");
#	&debug("************************\n");
###

	if ($rules)
	{
	    foreach my $num (@num)
	    {
		my $length = length($a[$num]);
		$b[$num] = $hash{$num} x $length;
	    }
	}

	
	$substruct = join("",@b);

### Pep's usefull debuggging options

	&debug($substruct);
	my $o=scalar( $substruct =~ s/\(/\</og );
	my $e=scalar( $substruct =~ s/\)/\>/og );
	&debug("$substruct :: $o == $e ? ".($o == $e ? "YES":"NO")."\n");
###	
	
	if ((($energy = &general_energy()) < $e1) && !$R)
	{
	    $totalhits++;	
	    open (EN ,">>$dir/$filename\.$patname\_en")|| die "can't open $dir/$filename\.$patname\_en: $!"; 
	    if ($start>$stop)
	    {
		print EN ">$seqname [$start - $stop] (Structure on complementary strand) - Free Energy: $energy\n$currentseq\n$currentstruct\n";
	    }
	    else
	    {
		print EN ">$seqname [$start - $stop] - Free Energy: $energy\n$currentseq\n$currentstruct\n"; 
	    }
	    
	    
	    close(EN); 
	    &gff_output();
	    
	    push @texta, [ @a ];
	    push @name, $name;
	    push @start, $start;
	    push @stop, $stop;
	    push @png, "$seqname".".png";
	    push @comment, $start>$stop ? "(SECIS on complementary strand)" : "";
	    push @energy, $energy;
	    
	    ### Display Progress
	    
	    if ($v)
	    {
		print STDERR "#"; 
	    }	  
	    
	}    &debug("*");    
    }  
    
    ###display progress
    
    if ($v)
    {
	print STDERR "#"; 
    }

    &debug(":");
    
    
    
### Die if no stable structures
    
    if ($totalhits == 0)
    {
	print STDERR "\n*** Sorry, no stable structures were found. You may want to try a different energy threshold.\n";
	exit();
    }
    
}
####################################################################
##          Run input file through RNAfold   (general)            ##
####################################################################


sub general_energy
{
    if (($R) && ($v)) 
    {
	print STDERR "*** Ignoring RNAfold\n";
    }
    else
    {
	&debug("### CALC ENERGY...");    
	open(TMP, ">$dir/$$\_en"); 
	print TMP  ">$$\n", @a, "\n", @b, "\n";
	close(TMP); 
	$cmd = "$RNAfold -C < $dir/$$\_en";
	
	&debug("### RNAfold command is $cmd ###\n"); 
	open(CMD,"$cmd |");
	@_a = <CMD>; 
	close(CMD);
	
	my $line="";
	foreach my $la (@_a)
	{
	    $line = $line.$la; 
	} 
	
	(my $namee , $currentseq , $currentstruct , $energy)=split(/\s\(|\n/,$line);  
	chomp $energy; chop $energy;

	system("mv $$\_ss.ps $pngdir/$name\.ps");
	system("pstopng -margin 3 -crop -Wpnmmargin,-black $pngdir/$name\.ps") if ($I);
	unlink ("$pngdir/$name\.ps");
	unlink ("$dir/$$\_en");
	&debug("### CALC ENERGY DONE...");
	&debug ("energy is $energy\n");
	return $energy; 
    }
}



####################################################################
##              Run input file through RNAfold                    ##
####################################################################

sub secis_energy 
{
    if (($R) && ($v)) 
    {
	print STDERR "*** Ignoring RNAfold\n";
    }
    else
    {
	&debug("### CALC ENERGY...");    
	open(TMP, ">$dir/$$\_en"); 
	print TMP  ">$$\n", @a[@_], "\n", @b[@_], "\n";
	close(TMP); 
	$cmd = "$RNAfold -C < $dir/$$\_en ";
	
	
	&debug("### RNAfold command is $cmd ###\n"); 
	open(CMD,"$cmd |");
	@_a = <CMD>; 
	close(CMD);
	
	my $line="";
	foreach my $la (@_a)
	{
	    $line = $line.$la; 
	} 

	(my $namee , $currentseq , $currentstruct , $energy)=split(/\n|\s\(/,$line);  
	chomp $energy; chop $energy;


	unlink ("$dir/$$\_en");
	unlink ("$$\_ss.ps");
	
	&debug("### CALC ENERGY DONE...");
	&debug ("energy is $energy\n");
	return $energy; 

    }
}

####################################################################
##              SECIS-specific structural features                ##
####################################################################

sub structure_OK {
    &debug("### CALC STRUCTURE...");    
    my $fout=1;
    
    my $substruct=substr($currentstruct, length( join '', @a[0..3]),  length( join '', @a[4..$unit7], (length ($a[$unit8]) == 4) ? () : substr($a[$unit8], 0, 1)   )) ;
     

    my $j=$currentstruct;
    my $p=scalar( $j =~ s/\(/\</og );
    my $w=scalar( $j =~ s/\)/\>/og );

    ($p == $w) || return 0;

#Y-filter removes Y-shaped SECISes
    if ( ($substruct =~ /\)\.*\(/) && $Y_flag) { $fout=0;}
    
#O-filter removes SECISes with more than 2 unparied nucleotides in a row on any strand in the first 7nt from the quartet
    my $porog=7;
    if ( (substr($substruct, 0, $porog) =~ /\.{3,}/ ||  substr($substruct, -$porog)=~ /\.{3,}/ ) && $O_flag) { $fout=0 }
    
#B-filter removes SECISes that have 2 or more unpaired nt more on 5' side than on 3' side (visually bended to the left)
    ($substruct =~ /(.*)(\(\.*\))(.*)/);
    my $up=$1; my $down=$3;
    if (((($up =~ s/\./\./g)-($down =~ s/\./\./g)) < (-2)) && $B_flag) { $fout=0 }
    
#S-filter removes structures with less that 8 pairs
    if ((( $substruct =~ tr/\)/\)/ ) < 8 ) && $S_flag ) { $fout=0 }
    
    &debug("STRUCTURE finished: $fout");
    return $fout;

}




###############################################################################################
################                       START  imager2                          ################
###############################################################################################


sub imager2 
{ 
    &debug("### DRAWING STRUCTURE");
    
    @bold=();
    (my $seq, my $struct, $scale, $outfilename, $offx, $offy, $posbase1, $posbase2, @bold)=@_;

    @x=();
    @xnew=();
    @ynew=();
    @y=(); 
    @seq=(); 
    @xf=(); 
    @yf=();
    
    $pi=3.1415926;
    $a=$scale;
    $fontsize=$a;
    $fonttype='Courier';
    
    $seq=' '.$seq.' ';
    
    $structmod="(".$struct.")";
    
    @seq=split '', $seq;
    
    
    $x[0]=0;
    $y[0]=0;
    
    $x[$#seq]=$x[0]+$a;
    $y[$#seq]=$y[0];
    


    @pairof=pairof_from_brackets($structmod);
    
    @parts=(0);
    while (@parts) { drawpart(pop @parts) };
    
    pop @x; pop @y; pop @seq; shift @x; shift @y; shift @seq;
    @pairof=pairof_from_brackets($struct);
    position_picture($posbase1, $posbase2, $offx, $offy);
    
    
    psoutput();
##########################################################################################
    
    
    sub pairof_from_brackets 
    {
	
	&debug("#F# PAIROF FROM BRACKETS\n");
	
	my @pairof;
	my $str = shift;
	$str =~ s/\s+//og;
	my @struct=split "", $str;
		
	for my $i (0..$#struct) 
	{
	    if ($struct[$i] eq '.') { $pairof[$i]=-1 }
	    elsif ($struct[$i] eq '(') {push @left, $i}
	    elsif ($struct[$i] eq ')') {my $k=pop @left; $pairof[$k]=$i; $pairof[$i]=$k;}
	    
	    else {die "error: unrecognized symbol \"$struct[$i]\" in structure"}
	}
	
	return @pairof;
	
    }

#################################---BEGIN drawpart---################################################

sub drawpart 
{
    
    &debug("#F# DRAW PART");
    
    $i=$_[0];      
    $j=$pairof[$i];
    


    $k=1;          
    @circle=();
    @xlocal=();
    @ylocal=();
    
    $xlocal[$i]=0;
    $ylocal[$i]=0;
    
    while ($pairof[$i+$k]==($j-$k)) {

	$xlocal[$i+$k]=0;
	$ylocal[$i+$k]=$a*$k;
	($x[$i+$k], $y[$i+$k])=globalized($xlocal[$i+$k],$ylocal[$i+$k]);
	
	$xlocal[$j-$k]=$a;
	$ylocal[$j-$k]=$a*$k;
	($x[$j-$k], $y[$j-$k])=globalized($xlocal[$j-$k],$ylocal[$j-$k]);
	$k++;
    }
    
    for ($n=$i+$k; $n<=$j-$k; $n++) 
    {
	if ($pairof[$n]!=(-1)) 
	{ 
	    push @parts, $n;
	    push @circle, $n;
	    $n=$pairof[$n];
	    next if ($pairof[$n] == 0);
	}
	push @circle, $n;
	
    }
    $alpha=2*$pi/(2+@circle);
    $r=$a/(2*sin($alpha/2));
    $xcenter=$xlocal[$i+$k-1]+$r*sin($alpha/2);
    $ycenter=$ylocal[$i+$k-1]+$r*cos($alpha/2);
    
    for $ci (0..$#circle) {
	$beta=(3/2)*($pi-$alpha)-$alpha*$ci;
	$xlocal[$circle[$ci]]=$xcenter+$r*cos($beta);
	$ylocal[$circle[$ci]]=$ycenter+$r*sin($beta);
	($x[$circle[$ci]], $y[$circle[$ci]])=globalized($xlocal[$circle[$ci]], $ylocal[$circle[$ci]]);
    }
##########################
	
	sub globalized {
	    &debug("#F# GLOBALIZED");
	    my $xl=$_[0];        
	    my $yl=$_[1];
	    $rotangle=atan2( $y[$j]-$y[$i], $x[$j]-$x[$i] );
	    my $xg=sqrt( $xl**2+$yl**2 )*cos( atan2($yl, $xl)+$rotangle)+$x[$i];
	    my $yg=sqrt( $xl**2+$yl**2 )*sin( atan2($yl, $xl)+$rotangle)+$y[$i];
	    return ($xg, $yg);
	}
	
##########################
	
    }
    
###############################---END drawpart---################################################
    
    
    sub position_picture {
	&debug("#F# POSITION PICTURE");
	my ($base1, $base2, $xbase1offset, $ybase1offset) = @_;
	my $alpha=3.1415926/2-atan2 ($y[$base2]-$y[$base1], $x[$base2]-$x[$base1]) ;
	
	for $j (0..$#x) {
	    $xnew[$j]=$x[$j]-$x[$base1]; 
	    $ynew[$j]=$y[$j]-$y[$base1];
	    ($xnew[$j], $ynew[$j])= (sqrt($xnew[$j]**2+$ynew[$j]**2)*cos( atan2($ynew[$j],$xnew[$j])+$alpha),sqrt($xnew[$j]**2+$ynew[$j]**2)*sin( atan2($ynew[$j],$xnew[$j])+$alpha));
	    $xnew[$j]+=$xbase1offset;
	    $ynew[$j]+=$ybase1offset;
	}
	
	@x=@xnew;
	@y=@ynew;
    }
    
####################################################################################
    
    sub psoutput { 
	&debug("#F# PS OUTPUT");
	open(OUTF, "| $ghostvw -q -dNOPAUSE -sDEVICE=png256 -g400x600 -dBATCH -sOutputFile=$pngdir/$outfilename -");
	
	show_seq();
	show_connections();
	show_bonds();
	print OUTF "stroke\n showpage\n";
	close(OUTF); 
    }
    
####################################################################################
    
    sub show_seq {
	@nonbold=(0..$#x);
	$k=0;
	print OUTF "newpath\n /$fonttype\-Bold findfont $fontsize scalefont setfont\n";
	for $index (@bold) 
	{
	    $xf=$x[$index]-0.6*$fontsize/2;
	    $yf=$y[$index]-0.6*$fontsize/2;
	    print OUTF "$xf $yf moveto\n ($seq[$index]) show\n";
	    
	    splice @nonbold, $index-$k, 1; 
	    $k++;
	    
	}
	
	print OUTF "newpath\n /$fonttype findfont $fontsize scalefont setfont\n";
	for $index (@nonbold) 
	{
	    $xf=$x[$index]-0.6*$fontsize/2;
	    $yf=$y[$index]-0.6*$fontsize/2;
	    if (! defined $seq[$index]) { die "@seq\n@x\n"}
	    print OUTF "$xf $yf moveto\n ($seq[$index]) show\n";
	}
    }
    
####################################################################################
    
    sub show_connections {
	for my $i (0..($#x-1)) {
	    $x1=$x[$i]; 
	    $y1=$y[$i];
	    $x2=$x[$i+1];
	    $y2=$y[$i+1];
	    line_bound('black', 0.25, 0.04); 
	}
    }
    
    
####################################################################################
    
    sub show_bonds 
    {
	$gam=0.2;
	
	for my $i (0..$#x) 
	{
	    if ($i<$pairof[$i]) 
	    {
		$curpair=$seq[$i].$seq[$pairof[$i]];
		$x1=$x[$i]; $y1=$y[$i]; $x2=$x[$pairof[$i]]; $y2=$y[$pairof[$i]];
		if (($curpair eq 'GU') || ($curpair eq 'UG')) {circle_bound('lightblue', $gam)}
		elsif (($curpair eq 'AU') || ($curpair eq 'UA')) {circle_bound('mediumblue', $gam)}
		elsif  (($curpair eq 'GC') || ($curpair eq 'CG')){circle_bound('red', $gam)}
		else {diamond_bound('burgundy', $gam, 0.3)}
	    }
	}
    }

####################################################################################

sub line_bound 
{
    my @color=rgbcolor($_[0]);
    my $gamma=$_[1];
    my $teta=$a*$_[2]; 
    
    $xstart=$x1+($x2-$x1)*(1-$gamma)/2; 
    $ystart=$y1+($y2-$y1)*(1-$gamma)/2;
    $xstop=$x2-($x2-$x1)*(1-$gamma)/2;
    $ystop=$y2-($y2-$y1)*(1-$gamma)/2;
    print OUTF "newpath\n @color setrgbcolor\n $teta setlinewidth\n $xstart $ystart moveto\n $xstop $ystop lineto\n stroke\n";
}

####################################################################################

sub circle_bound 
{
    my @color=rgbcolor($_[0]);
    my $gamma=$_[1];
    my $xc=($x1+$x2)/2;
    my $yc=($y1+$y2)/2;
    my $rc=$a*$gamma/2;
    print OUTF "newpath\n @color setrgbcolor\n $xc $yc $rc 0 360 arc\n fill\n stroke\n";
}


####################################################################################

sub diamond_bound 
{
    my @color=rgbcolor($_[0]);
    my $gamma=$a*$_[1]/2;
    my $teta=$a*$_[2]/2; 
    my $rotangle=atan2 ($y2-$y1, $x2-$x1);
    my $xc=($x1+$x2)/2;
    my $yc=($y1+$y2)/2;
    my $x_1=(-1)*sin($rotangle)*$teta+$xc;
    my $y_1=cos($rotangle)*$teta+$yc;
    my $x_2=cos($rotangle)*$gamma+$xc;
    my $y_2=sin($rotangle)*$gamma+$yc;
    my $x_3=sin($rotangle)*$teta+$xc;
    my $y_3=(-1)*cos($rotangle)*$teta+$yc;
    my $x_4=(-1)*cos($rotangle)*$gamma+$xc;
    my $y_4=(-1)*sin($rotangle)*$gamma+$yc;
    print OUTF "newpath\n @color setrgbcolor\n $x_1 $y_1 moveto $x_2 $y_2 lineto $x_3 $y_3 lineto $x_4 $y_4 lineto closepath\n
fill\n stroke\n";
}

####################################################################################

sub rgbcolor 
{
    my @rgb;
    if ($_[0] eq 'lightblue') { @rgb=qw(0.7 0.7 0.7) }
    elsif ($_[0] eq 'mediumblue') { @rgb=qw(0 0 0.9) }
    elsif ($_[0] eq 'darkblue') { @rgb=qw(0 0 0.2) }
    elsif ($_[0] eq 'burgundy')  { @rgb=qw(0.5 0.1 0.1) }
    elsif ($_[0] eq 'black')  { @rgb=qw(0.5 0.5 0.5) }
    elsif ($_[0] eq 'red')  { @rgb=qw(0.9 0 0) }
    return @rgb;
}
}


###############################################################################################
################                      End of imager2                           ################
###############################################################################################




####################################################################
##    Unless -s, print SECIS elements found => $filename.secis    ##
####################################################################

sub secis_output
{

    my $mod;
    if ($gen){$mod = 'fa'}
    elsif($secis){$mod = 'secis'}
    
    if ($v)
    {
	print "\n\n";
    }

    open (HIT , ">$dir/$filename\.$patname\.$mod")|| die "can't open $dir/$filename\.$patname\.mod: $!"; 
    open (ENN,"$dir/$filename\.$patname\_en")|| die "can't open $dir/$filename\.$patname\_en: $!";
    unless ($s)
    {
	while (<ENN>)
	{
	    print HIT unless (/^\W*$/o);  ### Do not print the RNAfold constraints
	}
    }
    system("cp $dir/$filename\.$patname\_en $dir/$filename\.$patname\.en") if $en;
    close (HIT);
    close (ENN);
}

####################################################################
##               Set what will be sent to STDOUT                  ##
####################################################################
    
sub standard_out
{
    if ($output eq 'gff')
    {
	
	open (GFF, "<$dir/$filename\.$patname\_gff") || die "Cannot open $dir/$filename\.$patname\_gff: $!";
	while (<GFF>)
	{
	    print STDOUT;
	}
	close (GFF);

	if (($T || $v) && $gen)
	{
	    print STDERR "\nStructures found: $totalhits\n";
	}
	elsif (($T || $v) && $secis)
	{
	    print STDERR "\nSECIS Elements found: $totalhits\n";
	}
	
    }
    elsif ($output eq 'fs')
    {
	open (ENN,"$dir/$filename\.$patname\_en")|| die "can't open $dir/$filename\.$patname\_en: $!";
	while (<ENN>)
	{
	    unless (/^\W*$/o)
	    {
		print STDOUT;
	    }
	}
	close (ENN);
	if (($T || $v) && $gen)
	{
	    print STDERR "\nStructures found: $totalhits\n";
	}
	elsif (($T || $v) && $secis)
	{
	    print STDERR "\nSECIS Elements found: $totalhits\n";
	}
	
    }
    
}

####################################################################
##                       Produce HTML if -H                       ##
####################################################################

sub HTML_output 
{
    if ($H)
    { 
	if (@name) 
    { 
	    open (HTML, ">$dir/$filename\.$patname\.html") || die "can't open $dir/$filename\.$patname\.html: $!";
	    print HTML "<HTML><HEAD><TITLE>SECISearch Results</TITLE></HEAD><BODY BGCOLOR=lightgrey>\n";
	    print HTML "<H3>Total hits: $totalhits </H3>";
	    print HTML "<TABLE BORDER=1>";
	    for (0..$#name) {
		$hn=$_+1;
		print HTML "<TR> <TD VALIGN=top> <B> $hn </B> $name[$_]: $start[$_] - $stop[$_] $comment[$_]<BR>";
		if ($mode eq 'secis')
		{
		    print HTML "Stem: $upstemen[$_] Whole structure: $fullstructen[$_]<BR><BR>";
		}
		elsif ($mode eq 'general')
		{
		    print HTML "Free Energy: $energy[$_] <BR><BR>"; 
		}
		print HTML  <<"STOPHERE";
		@{$texta[$_]}[0..2] <FONT COLOR=red> $texta[$_][3] </FONT> $texta[$_][4] <FONT COLOR=red> 
		    $texta[$_][5] </FONT> @{$texta[$_]}[6..7]<FONT COLOR=red> $texta[$_][8] </FONT>@{$texta[$_]}[9..11]
		    <BR>
		    <TD> <IMG SRC="$pngdir/$png[$_]">
		    
STOPHERE
                            }
	    print HTML "</TABLE>";
	}
	else 
	{
	    print HTML 'Sorry, no SECIS elements were found! <br>You may want to try another pattern';    
	}
    }
}

###################################################################
##                 Create filename.hits when -t                  ## 
###################################################################

sub patscan_hits
{
    system("cp $dir/$filename\.$patname\_hits $dir/$filename\.$patname\.hits") if ($t);
}
    
   
###################################################################
##                 Create filename.gff when -g                   ## 
###################################################################


sub gff_output
{
    if (($g) || ($output eq 'gff'))   ### Necessary for the gff STDOUT output
    {
	my $frame    = '.';
	my $group    = '.';
	my $score    = '.';
	my $realname = '.';
	my $strand   = '.';
	open (GFF, ">>$dir/$filename\.$patname\_gff") || die "Cannot open $dir/$filename\.$patname\_gff: $!";
	
	$name =~ />(.*)/;
	$realname = $1;
	if ($start>$stop){$strand = '-';}
	else {$strand = '+';}
	if ($secis)
	{
	    print  GFF "$realname\tSECISearch\tSECIS\t$start\t$stop\t$score\t$strand\t$frame\t$group\n";
	}
	elsif ($gen)
	{
	    print  GFF "$realname\tSECISearch\t.\t$start\t$stop\t$score\t$strand\t$frame\t$group\n";
	}
	close (GFF);
	system("cp $dir/$filename\.$patname\_gff $dir/$filename\.$patname\.gff") if ($g);
  }   
}



#################################################################
###                     Debugging mode                        ###
#################################################################

sub debug
{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}


####################################################################
##                             Clean up                           ##
####################################################################

sub clean_up
{
    unlink("$dir/$filename\.$patname\_en");
    unlink("$dir/$filename\.$patname\_hits");
    unlink("$dir/inseq.fa");
    unlink("$dir/$filename\.$patname\_gff");
    
}

####################################################################
##                          HEEEEELP!!!!                          ##
####################################################################
sub usage()
{
    print STDOUT <<"EOF";

USAGE:
	
    anysearch [OPTIONS] -p [patscan pattern file] -f [FASTA input file]
    anysearch [OPTIONS] -p [patscan pattern file] -r [RNAfold constraints file] -f [FASTA input file]
    
    
*** If '-p' is not given, the program uses a default pattern.(SECIS)
*** If '-f' is not given, the program reads input from STDIN and expects 
    it to be in fasta format
*** For more information try anysearch -h.
EOF
}

sub help()
{
    if ($output ne 'fs' && $output ne 'gff' && $output != 0)
    {
	print STDERR "\n*** The values for the -o option are 'fs' for fasta format and 'gff' for gff format.\nTry 'SECISearch -h' for more infomation.\n"; 
	exit();
	
    }
    if ($h)
    {
    print STDOUT <<"EOF";

DESCRIPTION: 

    This program is an extension of Gregory Kryukov's SECISearch. By default it  uses patscan
 to scan nucleotide sequences for SECIS elements and then RNAfold to calculate their thermodynamic stability. 

    However, anysearch can also be used to scan nucleotide sequences for any pattern. The patscan pattern
file is specified by the -p flag. When not searching for SECIS elements, the program should be run in General 
mode (-G flag) to avoid SECIS specific features. RNAfold binding constraints can also be imposed in the General
mode (These are automatic in SECIS mode). To do so use the -r flag and a simple text file containing the rules. 
This should be in the following format:

eg.
   2       (
   5       x
   8       )

    The numbers on the left specify which of the sequence units generated by patscan this constraint applies to.
(Where the first unit is 0) The symbols on the right (separated from the number by a tab) will be multiplied by 
the length of the unit and passed to RNAfold as constraints for that unit. In the example above, the line 

    2      (

means that the third unit from the beginning should have the constraint '( x (length of unit)'. That is, if the 
unit is 5 nt long, the constraint will be '((((('. Bear in mind however, that this will not work if your pattern
has an OR with more units on one side than the other. For example with a pattern such as "
    ACCT (1...3 2..4 | 2...8) 2...3"
the positions are not standard because the matching sequence might have 3 or 4 units. This means that you cannot
really impose a constraint on positions after the OR because you will not be specifying one specific position but
two possible ones. The sames is true in patterns containing units such as 0...3 which may or may not exists
(because of the 0).
    The patscan pattern files follows the same format as in patscan itself, with the only difference that any 
rules (eg r1={au,ua,gc,cg,gu,ug}) or similar, should be on the first line of the file. The line containing the
pattern units should only contain pattern units or the program gets confused.

    Anysearch will return all the sequences which matched the pattern and passed the thermodynamic criteria as
a fasta file to STDOUT. For other output options please see below.


REQUIRMENTS:
   
    Anysearch requires RNAfold from the VIENNARNA package (http://www.tbi.univie.ac.at/RNA/) and the
scan_for_matches pattern matcher from PatScan (http://www-unix.mcs.anl.gov/compbio/PatScan/HTML/patscan.html)
Additionally, anysearch uses ghostscript (http://www.cs.wisc.edu/~ghost/) in the secis mode and pstopng 
(http://www.math.utah.edu/pub/pstopng/) in the general mode to convert the raw ps format output of RNAfold to
prettier and more easily manageable png files. The program works perfectly well however without these, you just
don't get any pretty pictures. (If you get any errors for lack of these applications, run anysearch with the -I
flag)

USAGE:
	
    anysearch [OPTIONS] -p [patscan pattern file] -f [FASTA input file]
    anysearch [OPTIONS] -p [patscan pattern file] -r [RNAfold constraints file] -f [FASTA input file]
    
    
*** If '-p' is not given, the program uses a default pattern.(SECIS)
*** If '-f' is not given, the program reads input from STDIN and expects 
    it to be in fasta format



OPTIONS:


    -c : Patscan searches complimentary strand
    -d : Debugging mode (very very very verbose)
    -e : Free energy cut-off value of the stem (or of the entire structure in general mode)
    -E : Free energy cut-off value of the entire structure (SECIS mode)
	 NOTE: If no values are given, e: -10 E: -20 (SECIS), e: -10 (General).
    -f : Fasta file of the input sequence
    -g : Create <FILENAME>.<PATNAME>.gff
    -G : Program runs in general mode, not SECIS.
    -h : Display this message and exit
    -H : Produce HTML output
    -I : Do not return images
    -n : Return <FILENAME>.en
    -o : Choice of STDOUT output format, can be 'fs' for fasta 
         or 'gff' for GFF. The default output is fasta. 
    -p : Pattern file passed to scan_for_matches. If '-p' is not given, the standard 
         (SECIS) pattern is used. Possible values are 's' for a standard pattern, 'n' 
         for a non-standard pattern, and 't' for a "twilight zone" pattern. Any bareword 
	 value is first assumed to be a file in the current directory. If no such file
         is found, the program will search for a file of that name in the
         ~/share/SECISearch directory (or wherever else the patterns happen to be stored
         on your system). Any directory address will be assumed to be the location of a 
         user defined pattern file. CAUTION: SECISearch will not warn you if you give a 
         non-existent filename as a pattern, if it cannot find the pattern it will use 
         the default SECIS pattern, so check your pattern name before running.
    -P : Do not run patscan
    -r : Specify RNAfold constraints (rules). (General mode)
    -R : Do not run RNAfold
    -s : Do not return <FILENAME>.<PATNAME>.secis
    -t : Return <FILENAME>.<PATNAME>.hits
    -T : Print total number of SECIS found (prints to STDERR, usefull when running non-verbose
         mode and/or pipes)
    -v : Verbose output
   

SECIS Structural Feature Options (ON by default, use these options to turn the filters OFF):
[These are, of course, OFF and irrelevant in the General mode]
    -B : Discards SECISes with at least 2 more unpaired nts on the 5'
         side than on the 3' side (visually bended to the left)
    -O : Discards SECISes with more than 2 consecutive,
         unpaired nts on any strand in the first 7nt after the quartet
    -Y : Discards Y-shaped SECISes 
    -S : Discards structures with less than 8 pairs
   
    
OUTPUT FILES:

    <FILENAME>.<PATNAME>.en         : Each input sequence and its RNAfold 
                                      constraints (-n flag)
    <SEQUENCE NAME>.<PATNAME>.png   : RNAfold output image beautified and in
                                      png format (SECIS mode, in general mode 
				      you still get pngs but they are not pretty... 
    <FILENAME>.<PATNAME>.html       : HTML format output (-H flag)
    <FILENAME>.<PATNAME>.hits       : Fasta file of all the sequences which
                                      matched the pattern divided into units
				      by patscan.
    <FILENAME>.<PATNAME>.secis      : Fasta file of the SECIS elements found 
                                      (<FILENAME>.<PATNAME>.fa in general mode)
    <FILENAME>.<PATNAME>.gff        : GFF file of the SECISes found (-o gff)

The default output is <FILENAME>.<PATNAME>.secis (or .fa), and the <SEQUENCE NAME>.png files. 
The latter are placed in a directory called <FILENAME>.<PATNAME>_png. All other files are optional extras.

EOF
exit();
   }
}


####################################################################################################
#######                                     END                                             ########
####################################################################################################


####################################################################################################
#                                     NOTES
# work around the pattern position problem for rule generation. Maybe count from the end? 
#
####################################################################################################
