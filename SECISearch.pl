#!/usr/bin/perl
#
# 
# $Id: SECISearch.pl,v 1.5 2004/04/16 12:37:30 cchapple Exp cchapple $
#

############################################################################
##                                                                        ##
## This script was originally developed by Gregory Kryukov and has been   ##
## modified by Charles Chapple (cchapple@imim.es).                        ##
##                                                                        ##
##                                                                        ##
##  This program is free software; you can redistribute it and/or modify  ##
##  it under the terms of the GNU General Public License as published by  ##
##  the Free Software Foundation; either version 2 of the License, or     ##
##  (at your option) any later version.                                   ##
##                                                                        ##
##  This program is distributed in the hope that it will be useful,       ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
##  GNU General Public License for more details.                          ##
##                                                                        ##
##  You should have received a copy of the GNU General Public License     ##
##  along with this program; if not, write to the Free Software           ## 
##  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             ##
############################################################################


###########################################################################

use strict;
use Getopt::Std;
use File::Copy;



#################################################################
##           Define locations of external programs             ##
#################################################################
#my $scan1GB = "/usr/local/molbio/bin/scan_for_matches"; 
my $scan1GB = "/home/terdon/bin/scan_1GB";
my $RNAfold = "/home/terdon/bin/RNAfold13"; 
my $ghostvw = "/usr/bin/gs";
my $bdir    = "/home/terdon/share/SECISearch";


#################################################################
##     Get cmd-line options and declare some variables         ##
#################################################################

my %opts;
getopts('p:f:e:o:RHITPYOBSGhvtgsndcilx',\%opts) || do { print "Invalid option, try 'SECISearch -h' for more information\n"; exit(1); };

### Input Options

my $file    = $ARGV[0];         # sequence file
my $pat     = $opts{p};         # -p specifies the pattern for scan_for_matches
my %E1      = (pat_standard_I_II => -10, pat_non_standard_I_II => -11, pat_twilight_I_II => -11);
my %E2      = (pat_standard_I_II => -11, pat_non_standard_I_II => -11, pat_twilight_I_II => -11);



### Scanning Options

my $c;              # -c scan_for_matches will scan complementary strand (ON by default)
my $Y_flag  = 1;    # S,Y,O,B are filters for specific SECIS features
my $O_flag  = 1;    #
my $B_flag  = 1;    #
my $S_flag  = 1;    #


if ($opts{Y}){$Y_flag = 0}
if ($opts{O}){$O_flag = 0}
if ($opts{B}){$B_flag = 0}
if ($opts{S}){$S_flag = 0}

### Output Options

my $I       = $opts{I} || 0;    # -I : bypass the imager2 function
my $P       = $opts{P} || 0;    # -P : bypass scan_for_matches 
my $H       = $opts{H} || 0;    # -H : produce HTML output
my $v       = $opts{v} || 0;    # Verbose
my $t       = $opts{t} || 0;    # Return <FILENAME>.hits
my $en      = $opts{n} || 0;    # Don't unlink energy param file
my $T       = $opts{T} || 0;    # Print no of SECIS found
my $g       = $opts{g} || 0;    # Produce gff output
my $secisfl = $opts{s} || 0;    # Do not return <FILENAME>.secis
my $R       = $opts{R} || 0;    # Don't run RNAfold 
my $output  = $opts{o} || 'fs';    # what is printed at the STDOUT
my $geneid  = $opts{G} || undef;

my $h       = $opts{h} || 0;    # cry for help
my $debug   = $opts{d} || 0;    # run with debugging info
my $badens  = $opts{l} || 0;    # print energies which did not pass in the logfile
my $log     = $opts{x} || 0;    # create the secislog logfile

my $cmpltry = 0;
if($opts{c}) {
    $cmpltry = 1
    }# search conmplimentary strand by default
    $log = 1 if $badens;
$g = 1 if $geneid;

#################################################################
##               Define some global variables                  ##
#################################################################

my $dir = $ENV{PWD} || './';   # get name of current directory

my $filename;
my $patname;
my $pngdir;
my $e1;
my $e2;

my $totalhits = 0;
my $units;
my $unit8;
my $lastunit;
my $unit7;
my $pngdirexists = 0;


####################################################################
##                         Define Variables                       ##
####################################################################

my ($alpha, $beta, $ci, $cmd, $curpair, $currentseq, $currentstruct, $energy, $fontsize,
    $fonttype, $fullstructen, $gam, $hn, $index, $j, $k, $n, $name, $offx, $offy, 
    $outfilename, $pi, $posbase1, $posbase2, $r, $rotangle, $scale, $start, $stop, 
    $structmod, $substruct, $upstemen, $x1, $x2, $xcenter, $xf, $xstart, $xstop, 
    $y1, $y2, $ycenter, $yf, $ystart, $ystop, $i);

my (@_a, @bold, @circle, @comment, @fullstructen, @left, @name, @nonbold, @pairof,
    @parts, @png, @seq, @start, @stop, @texta, @upstemen, @x, @xf, @xlocal, @xnew,
    @y, @yf, @ylocal, @ynew);

 
my @a = ();
my @b = ();

##############################    

    
 



###########################################################################################
###############                   Core program START                      #################
###########################################################################################


&usage();  
                      &debug("usage finished\n");
&get_params();
                      &debug("get_param finished\n");
&prepare_ground();
                      &debug("prepare_ground finished\n");
&scan_for_matches();
                      &debug("scan_for_matches finished\n");
&thermodynamic_eval();
                      &debug("thermodynamic_eval finished\n");

###########################################################################################
###############                      Core program END                     #################
###########################################################################################

                       &debug("main loop finished");
&secis_output();
                           &debug("secis_output finished\n");
&standard_out();
                       &debug("standard_out finished\n");
&unlink_en();
                       &debug("unlink_en finished\n");
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
### Really, really crude and ugly way to allow input seq from STDIN
### 
    
    unless ($ARGV[0])
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
	}
	elsif ($opts{p} !~ /\//)
	{
	    $pat = "$bdir/$opts{p}";
	}
	elsif ($opts{p} =~ /\//)
	{
	    $pat = "$opts{p}"
	}
    }
    else 
    {
	unless ($opts{P})
	{
	    print STDERR "*** no pattern specified, or specified pattern not found, using default pattern ***\n";
	    $pat = "$bdir/pat_standard_I_II";
	}
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
    if ($filename =~ /(\.fa)||(\.fasta)$/)
    {
	$filename =~ s/((\.fa)||(\.fasta))$//;
    }
    if ($pat =~ m/\/((\w*)|(\w*\.\w*))$/g)
    {
	$patname = $1;
    }
    
    else 
    {
	$patname = $pat;
    }
 
 if ($opts{P})   ### Get patname for output filenames if not running scan_for_matches
    {
	$patname = $1 if $file =~ /(pat.*?)\./o;
	$patname = 'nopat' if $file !~ /(pat.*?)\./o;
	print "\$patname is $patname\n"; 
    }

### Allow for user-defined energy cutoff values     
    
    
#   if ($opts{e}){$e1 = $opts{e};} else {$e1 = $E1{$patname} || -10}   Old energy for core
    if ($opts{e}){$e2 = $opts{e};} else {$e2 = $E2{$patname} || -11} 
    
    if ($patname =~ 'pat_standard_I_II') { $patname = 'std'}
    elsif ($patname =~ 'pat_non_standard_I_II'){$patname = 'non_std'}
    elsif ($patname =~ 'pat_twilight_I_II'){$patname = 'twil'}
    else {$patname = $patname}




###  Have scan_for_matches search complementary strand      

    if ($cmpltry){$c = '-c';}    
 
### Print pretty debugging information
    &debug("##################");
    &debug("\$dir is $dir");
    &debug( "\$filename is $filename");
    &debug( "\$file is $file");
#    &debug( "\$e1 is $e1");
    &debug( "\$e2 is $e2");
    &debug( "\$pat is $pat");
    &debug( "\$patname is $patname");
    &debug( "\$log is $log");
    &debug( "\$output is $output");
#    &debug( "\$ is $");
#    &debug( "\$ is $");
    &debug("##################\n");
    
 ### Print pretty log information
 if ($log){
     open (LOG, ">>secislog");
     my $date =  localtime;
     print LOG"##################\n";
     print LOG "$date\n\$filename : $filename\nenergy cutoff :  $e2\n\$pat : $pat\n\$patname : $patname\n";
     print LOG "##################\n";
 }


}



##################################################################
##        Delete files from a possible previous instance        ##
##################################################################


sub prepare_ground
{
    unlink("$dir/$filename\.$patname\.gff");
    unlink("$dir/$filename\.$patname\.en");
}


##################################################################
##    Running scan_for_matches with choosen pattern             ##
##################################################################


sub scan_for_matches
{

    if ($P == 0)
    {      
	if ($v || $debug)
	{
	    print STDERR "*** running scan_for_matches, please wait...\n"; 
	    print STDERR "\nNOTE :  remember that patscan will NOT scan the complementary strand\n\tunless specifically told to do so (-c)\n\n" unless $cmpltry;
	}

	&debug("### PatScan Command: $scan1GB $c $pat < $file > $dir/$filename\.$patname\.hit  ###\n");
	system("$scan1GB $c $pat < $file > $dir/$filename\.$patname\.hit") || die "cannot run scan_for_matches: $!";

	### Check that at least one hit was found
	if (-z "$dir/$filename\.$patname\.hit")
	{
	    unlink("$dir/$filename\.$patname\.hit");  ### if empty, remove patscan hits file.
	    die("*** Sorry, no matching sequences found. You may want to try another pattern ***\n");
	}
	open (HITS, "$dir/$filename\.$patname\.hit") || die "can't open $filename\.$patname\.hit: $!";
    }   
    else
    {
	if ($v)
	{
	    print STDERR "*** ignoring scan_for_matches\n";
	}
	open (HITS, "$file") || die "can't open $file: $!";
    }

###################################################################
##                 Create filename.hits when -t                  ## 
###################################################################

    if($t)
    {
	copy("$dir/$filename\.$patname\.hit","$dir/$filename\.$patname\.hits");
	print STDERR "*** Finished scanning sequence, moving on...\n" if $v;
    }    
    
}


####################################################################
##                        Ignore RNAfold                          ##
####################################################################

sub skip_rnafold
{

    open (PWP, "$file\.$patname\.hit");
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
  
    &debug("patscan_hits finished\n");
   #  if ($v)
#     {
# 	open (HITS, "$dir/$filename\.$patname\.hit")|| die "can't open $dir/$filename\.$patname\.hit: $!";
# 	while (<HITS>){ print STDOUT}
# 	close(HITS);
#     }
    
    &patscan_hits() if $t;
                      &debug("skip_rnafold finished\n");
    &thermodynamic_eval();
                      &debug("thermodynamic_eval finished\n");
    &secis_output();
                       &debug("secis_output finished\n");
    &standard_out();
                       &debug("standard_out finished\n");

    



    &unlink_en();
                     &debug("unlink_en finished\n");
    &HTML_output();
                     &debug("HTML_output finished\n");
    &clean_up();    
                     &debug("cleanup finished\n");
    exit();
}



######################################################################################################
##                                 Evaluate hits thermodynamically                                  ##
######################################################################################################

sub thermodynamic_eval
{

    if($R) ## If we don't want to run RNAFold
    {
# 	open (PWP, "$file\.$patname\.hit");
# 	while (<PWP>)
# 	{
# 	    if (/^>/)
# 	    {
# 		print; 
# 	    }
# 	    else
# 	    {
# 		s/(\w)\s(\w)/$1$2/g;
# 		print;
# 	    }
# 	}
	
# 	if ($v)
# 	{
# 	    print STDERR "\n*** WARNING: These sequences have NOT been evaluated thermodynamically\n\n";
# 	}
	
# 	&debug("patscan_hits finished\n");
# 	if ($v)
# 	{
# 	    open (HITS, "$dir/$filename\.$patname\.hit")|| die "can't open $dir/$filename\.$patname\.hit: $!";
# 	    while (<HITS>){ print STDERR}
# 	    close(HITS);
# 	}
    }
    else
    {
	if ($v || $debug)
	{
	    print STDERR "*** scan_for_matches finished, evaluating hits thermodynamically...\n";
	}
	
	my $seqname;
	my $sequence_counter = 0;
	while (<HITS>)    
	{ 
	&debug("#$.# $_");
	my @param;
	next if /^\s*$/o;
	/^>(.*)/o && do 
	{
	    $seqname = $1;
#	    ($name, $start, $stop) = split /\:\[|\,|\]/, $1;
	    $seqname =~ /^(.*):\[(\d+?),(\d+)\]/ || die("cannot match seqname $seqname : $!\n");
	    ($name, $start, $stop) = ($1,$2,$3);
	    next;
	};
	
	
	@b=();
	tr/Tt/Uu/;
	@a=split /\s+/o, $_;
	s/\w/./g for @b=@a; 

	$b[3]='x(((('; $b[5]='xx';


	$units = @a;
	$lastunit = $units -1;
	
	&debug("************************\n");
	&debug("units: $units\n");
	&debug("\n@a\n");
   	&debug("\@a[$units - 4]\n");
	&debug("************************\n");
	
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

	if (length ($a[$unit8]) == 4)  { $b[$unit8]='))))' }
	elsif (length ($a[$unit8]) == 5)  { $b[$unit8]='.))))'}
	else 
	{
	    print STDERR "When not running patscan, the input file must be the output of a previous patscan session. 
Otherwise, SECISearch encounters problems with the length of the sequence : " ;
	    print STDERR "length($a[$unit8]) : " . length($a[$unit8]) . " should be 4 or 5\n@a\n@b\n"; die(); 
	}
	

	$substruct = join("",@b);

### Pep's usefull debuggging options

	&debug($substruct);
	my $o=scalar( $substruct =~ s/\(/\</og );
	my $e=scalar( $substruct =~ s/\)/\>/og );
	&debug("$substruct :: $o == $e ? ".($o == $e ? "YES":"NO")."\n");
###	


	unless ($R)
	{ 
	    if (($fullstructen = &energy(0..$lastunit))<=$e2 && structure_OK()) 
	    { 
### Old energy stuff ($upstemen = energy(3..$unit8)) <=$E1& 		
		$totalhits++;
		open (EN ,">>$dir/$filename\.$patname\.en") || die "can't open $dir/$filename\.$patname\.en: $!"; 
		
		if ($start>$stop)
		{
		    print EN ">$name [$start - $stop] (SECIS on complementary strand) - Free Energy: $fullstructen\n$currentseq\n$currentstruct\n";
		}
		else
		{
		    print EN ">$name [$start - $stop] - Free Energy: $fullstructen\n$currentseq\n$currentstruct\n";
		}
		
		close(EN);

###
#    my $ja = join("",@b);
#    open (STR,">>$dir/structs");
#    print STR "$seq\n$struct\n$ja\n*\n";
#    close (STR);
###		
		

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
### In the absence of I, run through the imager2 function      

		unless($I)
		{
		    imager2(@param); 
		    
		    push @texta, [ @a ];
		    push @name, $name;
		    push @start, $start;
		    push @stop, $stop;
		    push @png, "$seqname".".png";
		    push @comment, $start>$stop ? "(SECIS on complementary strand)" : "";
		    # push @upstemen, $upstemen;
		    push @fullstructen, $fullstructen;
		    
		}  
###display progress
		
		if ($v)
		{
		    $sequence_counter++;
		    print STDERR "."; 
		    
		    print STDERR "[$sequence_counter]\n" if $sequence_counter % 100 ==0; 
		    
		}
		
		
		&debug("*");    
	    }
	    
	    else {print LOG ">$seqname [$start - $stop] - Free Energy $fullstructen\n$currentseq\n$currentstruct\n" if $badens && $log;}
	    
	    &debug(":");
	}  &debug(".");
	&gff_output() if (($g) || ($output eq 'gff'));
    }
	
	
### Die if no stable structures
	
	if ($totalhits == 0)
	{
	    print STDERR "\n*** Sorry, no stable structures were found. You may want to try a different energy threshold (-e).\n";
	    exit();
	}
	
    }


}


####################################################################
##              Run input file through RNAfold                    ##
####################################################################
sub energy 
{ 
    if (($R) && ($v)) 
    {
	print "Ignoring RNAfold\n";
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
    if ( (substr($substruct, 0, $porog) =~ /\.{3,}/ ||  substr($substruct, -$porog)=~ /\.{3,}/ ) && $O_flag) { $fout=0;}
    
#B-filter removes SECISes that have 2 or more unpaired nt more on 5' side than on 3' side (visually bended to the left)
    ($substruct =~ /(.*)(\(\.*\))(.*)/);
    my $up=$1; my $down=$3;
    if (((($up =~ s/\./\./g)-($down =~ s/\./\./g)) < (-2)) && $B_flag) { $fout=0;}
    
#S-filter removes structures with less that 8 pairs
    if ((( $substruct =~ tr/\)/\)/ ) < 8 ) && $S_flag ) { $fout=0;}
    
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
#    print STDERR "\$structmod is \n$structmod\n $struct\n"; die();
    
    @parts=(0);
    while (@parts) { &drawpart(pop @parts) };
    
    pop @x; pop @y; pop @seq; shift @x; shift @y; shift @seq;
    @pairof=&pairof_from_brackets($struct);
    &position_picture($posbase1, $posbase2, $offx, $offy);
    
    
    &psoutput();
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
	    my $ho = join("",@b);
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


###  Create png directory
	unless ($I || $pngdirexists) 
	{
	    mkdir("$dir/$filename\.$patname\_png", 0777) unless (-d "$dir/$filename\.$patname\_png");
	    $pngdir = "$dir/$filename\.$patname\_png";
	    $pngdirexists=1;
	}
###
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
##    if -s, print SECIS elements found => $filename.secis    ##
####################################################################

sub secis_output
{
    if($secisfl)
    {
	if($R)
	{
	    if ($v)
	    {
		print STDERR "\n\n";
	    }
	    open (HIT , ">$dir/$filename\.$patname\.secis")|| die "can't open $dir/$filename\.$patname\.secis: $!"; 
	    open (ENN,"$dir/$filename\.$patname\.hit")|| die "can't open $dir/$filename\.$patname\.hit: $!";
	    while (<ENN>)
	    {
		if (/^>/){print HIT}
		else
		{
		    my $a = uc($_);
		    $a =~ s/\s+//g;
		    my $l = length($a);
		    if ($l>60)
		    {
			my $i=0;
			while($i<$l)  
			{
			    print HIT substr($a,$i,60) . "\n"; 
			    $i=$i+60;
			}
		    }
		    else
		    {
			print HIT "$a\n";
		    }
		}
	    }
	}
	else
	{
	    if ($v)
	    {
		print "\n\n";
	    }
	    open (HIT , ">$dir/$filename\.$patname\.secis")|| die "can't open $dir/$filename\.$patname\.secis: $!"; 
	    open (ENN,"$dir/$filename\.$patname\.en")|| die "can't open $dir/$filename\.$patname\.en: $!";
	    while (<ENN>)
	    {
		unless (/^\W*$/o)
		{
		    print HIT;
		    
		}
		
	    }
	    
	    close (HIT);
	    close (ENN);
	}
    }
}

####################################################################
##               Set what will be sent to STDOUT                  ##
####################################################################

sub standard_out
{

   
    if ($output eq 'gff')
    {
	open (GFF, "<$dir/$filename\.$patname\.gff") || die "Cannot open $dir/$filename\.$patname\.gff: $!";
	while (<GFF>)
	{
	    print STDOUT;
	}
	close (GFF);
	if (($T)||($v))
	{
	    print SDTERR "\nSECIS Elements found: $totalhits\n";
	}
	unlink("$dir/$filename\.$patname\.gff") unless $g;
    }
    elsif ($output eq 'fs')
    {
	my $relevant_file;
	if($R==0){
	    $relevant_file = "$dir/$filename\.$patname\.en" 
	    }
	else{
	    $relevant_file =  "$dir/$filename\.$patname\.hit";
	}
	open (ENN,"$relevant_file")|| die "can't open $relevant_file: $!";
	while (<ENN>)
	{
	    unless (/^\W*$/o)
	    {
		if (/^>/){print STDOUT "$_"}
		else
		{
		    my $a = uc($_);
		    $a =~ s/\s+//g;
		    my $l = length($a);
		    if ($l>60)
		    {
			my $i=0;
			while($i<$l)  
			{
			    print STDOUT substr($a,$i,60) . "\n"; 
			    $i=$i+60;
			}
		    }
		    else
		    {
			print STDOUT "$a\n";
		    }
		    $totalhits++ if $R;
		}
	    }
	}
	close (ENN);
	
	if (($T)||($v))
	{
	  $R == 0 ?  print STDERR "\nSECIS Elements found: $totalhits\n" : print STDERR "\nSECIS Elements found (these sequences have NOT been thermodynamically evaluated!!):  $totalhits\n";
	}
	
    }  
    else{ die("output must be either gff or fasta, try SECISearch.pl -h for more help\n");}
    
}




#################################################################
##   Remove energy param file created by a previous            ##
##   instance, if present                                      ##
#################################################################


sub unlink_en
{
    unless ($en)
    {
	unlink("$dir/$filename\.$patname\.en");
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
		
		print HTML "Free Energy: $fullstructen[$_]<BR><BR>";
		print HTML  <<"STOPHERE";
		@{$texta[$_]}[0..2] <FONT COLOR=red> $texta[$_][3] </FONT> $texta[$_][4] <FONT COLOR=red> 
		    $texta[$_][5] </FONT> @{$texta[$_]}[6..7]<FONT COLOR=red> $texta[$_][8] </FONT>@{$texta[$_]}[9..11]
		    <BR>
		    <TD> <IMG SRC="$pngdir/$png[$_]">
		    
STOPHERE
                            }
	    print HTML "</TABLE>";
	}
	else { print HTML 'Sorry, no SECIS elements were found! <br>You may want to try another pattern';    }

	
	
    }
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
	my $score;
	if($geneid) ## geneid needs a positive score, so the thrmo score is no good
	{
	    $score = 1;
	}
	else 
	{
	    $score    = $fullstructen;
	}
	my $realname = '.';
	my $strand   = '.';
	open (GFF, ">>$dir/$filename\.$patname\.gff") || die "Cannot open $dir/$filename\.$patname\.gff: $!";
	$realname = $name;
	if ($start>$stop){$strand = '-';}
	else {$strand = '+';}
	print  GFF "$realname\tSECISearch\tSECIS\t$start\t$stop\t$score\t$strand\t$frame\n" if $strand eq '+';
	print  GFF "$realname\tSECISearch\tSECIS\t$stop\t$start\t$score\t$strand\t$frame\n" if $strand eq '-';
	close (GFF);
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
    unlink ("$dir/$filename\.$patname\.hit");
    unlink("$dir/inseq.fa");
    unless ($g)
    {
	unlink("$dir/$filename\.$patname\.gff");
    }
    if ($log) {
	print LOG "SECIS Found : $totalhits\n\n######################################################################\n";
	close(LOG);
    }
    rmdir($pngdir); ## remove png directory if empty
}


####################################################################
##                          HEEEEELP!!!!                          ##
####################################################################


sub usage()
{
    if ($output ne 'fs' && $output ne 'gff' && $output != 0)
    {
	print STDOUT "\n*** The values for the -o option are 'fs' for fasta format and 'gff' for gff format.\nTry 'SECISearch -h' for more infomation.\n"; 
	exit();
	
    }
    

    if ($h)
    {
	print STDOUT <<"EOF";

DESCRIPTION: 

    SECISearch uses patscan and RNAfold to scan nucleotide sequences for SECIS elements and evaluates
    them thermodynamically. By default, it prints the SECISes found to STDOUT in fasta format. If no 
    pattern is specified by the -p option, SECISearch uses its own default pattern. It was originally 
    developed as a web-based application by Gregory Kryukov and modified to its present form by
    Charles Chapple. 


USAGE:
	
    SECISearch [OPTIONS] -p [patscan pattern file] [FASTA input file]
    
    
*** If '-p' is not given, the program uses a default pattern.
*** If no filename is given, the program reads input from STDIN and expects 
    it to be in fasta format



OPTIONS:


    -c : Patscan will search the complimentary strand
    -d : Debugging mode (very very very verbose)
    -e : Free energy cut-off value of the entire structure
    -t : Return <FILENAME>.<PATNAME>.hits
    -g : Create <FILENAME>.<PATNAME>.gff
    -G : Create gff in geneid format (score=1), use with -g or -o gff.
    -h : Display this message and exit
    -H : Produce HTML output
    -I : Do not return images
    -l : Print structures which did not pass the thermodynamic evaluation to the logfile.
    -o : Choice of STDOUT output format, can be 'fs' for fasta 
         or 'gff' for GFF. The default output is fasta. 
    -n : Return <FILENAME>.en
    -p : Pattern file passed to scan_for_matches. If '-p' is not given, the standard 
         pattern is used. Possible values are 's' for a standard pattern, 'n' for a 
	 non-standard pattern, and 't' for a "twilight zone" pattern. Any bareword 
	 value is assumed to be one of the patterns stored in the ~/share/SECISearch
	 directory (or wherever else the patterns happen to be stored on your system).
	 Any directory address will be assumed to be the location of a user defined 
	 pattern file. CAUTION: SECISearch will not warn you if you give a non-existant
	 filename as a pattern, so check your pattern name before running.
    -P : Do not run patscan. (input file must be in the format of patscan's output) 
    -R : Do not run RNAfold
    -s : Return <FILENAME>.<PATNAME>.secis, FASTA file of SECIS elements
    -T : Print total number of SECIS found
    -v : Verbose output
    -x : Create the log file "secislog"

SECIS Structural Feature Options (ON by default, use these options to turn the filters OFF):
    
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
                                      png format
    <FILENAME>.<PATNAME>.html       : HTML format output (-H flag)
    <FILENAME>.<PATNAME>.hits       : Fasta file of all the sequences which
                                      matched the pattern divided into units
				      by patscan.
    <FILENAME>.<PATNAME>.secis      : Fasta file of the SECIS elements found 
    <FILENAME>.<PATNAME>.gff        : GFF file of the SECISes found (-o gff)
    secislog                        : Logfile(!)

The default output is <FILENAME>.<PATNAME>.secis and the <SEQUENCE NAME>.png files. The latter are placed in a directory called <FILENAME>.<PATNAME>_png. All other files are optional extras.

EOF
exit();

   }
}


####################################################################################################
#######                                     END                                             ########
####################################################################################################


####################################################################################################
#                                     NOTES
# 
#
####################################################################################################





