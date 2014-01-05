#!/usr/bin/perl  
#

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

use CGI qw(:standard);

## Variables for my laptop
# my $bdir = "/var/www";
# my $prgbin = "$bdir/cgi-bin/SECISearch";
# my $scan1GB = "$prgbin/scan_1GB";
# my $RNAfold = "$prgbin/RNAfold13";
# my $ghostvw = "/usr/bin/gs";
# my $Murl= "http://localhost/";
# my $pdir = "$bdir/htdocs/datasets/secisaln";
# #my $tmpdir = "/usr/local/Install/apache-data/genome.imim.es/cgi-bin/SECISearch/tmp"; # "$prgbin/tmp";
# my $tmpdir = "$prgbin/tmp"; # "$prgbin/tmp";
# my $pngpth = "$Murl";
# my $logfile = "$tmpdir/search$$.log";
##
##
##
## Variables for monstre1
my $type;

my $bdir = "$ENV{APACHE_ROOT}";
my $prgbin = "$bdir/cgi-bin/SECISearch";
my $scan1GB = "$prgbin/scan_1GB";
my $RNAfold = "$prgbin/RNAfold13";
my $ghostvw = "/usr/bin/gs";
$ENV{HTTP_REFERER}=~/^(.+)\//;
my $Murl= "$1";
my $pdir = "$bdir/htdocs/software/secisaln";
my $tmpdir = "$prgbin/tmp"; # "$prgbin/tmp";
my $pngpth = "$Murl";
my $logfile = "$tmpdir/search.log";
my $user_seq_name="user_seq";
my $tmp_print;
my $debug=param('debug')||0;
my $date= `date +%D" "%H:%M:%S`;
my $css;
my $added_a=0; ## set to one if sequence is too short and had top add As and Us
chomp($date);
open(EFILE,">> $logfile");
print EFILE "\n============$date============\n#<$$># STARTING NEW REQUEST #<$$>#\n";
##

my %global_pos;
my %opts;

my $spacer="left_spacer";
my $legend="legend2";
my $legend_image="type1.png";
my $lgnd_height=455;
my $lgnd_width=200;
my ($pdf_image,$png_image);
my @pred_str;
my $pat=param('pattern');



#$want_gis=param('want_gis')||undef;
my $seq=param('sequence') || param('sequence1');

my %e1      = (pat_standard_I_II => -10, pat_non_standard_I_II => -10, pat_twilight_I_II => -10);
my %E2      = (pat_standard_I_II => -17, 
	       pat_non_standard_I_II => -17, 
	       pat_twilight_I_II => -17,
	       pat_secisaln => -9,
	       pat_secisaln_loose=> -7,);

my $c = param('complstrand') ? '-c' : '';
my $e1=(param('customenergy') && param('e1')) || $e1{$pat};

#$e2=(param('customenergy') && param('e2')) || $e2{$pat};
my $e2=param('e2') || -9;


my $Y_flag=param('Y_filter');
my $O_flag=param('O_filter');
my $B_flag=param('B_filter');
my $S_flag=param('S_filter');



my $sort_by_prot=undef;
if(param('sort') eq "prot"){
    $sort_by_prot=1;
}

#### Clean files from old runs
system("ls *ps $pdir/*png $pdir/*pdf tmp/* > pslist;");
open(PS,"pslist");
while(<PS>){
    chomp;
    my $hh=(-A "$_");
    system("rm $_") if $hh > 2;
}
unlink("pslist");



#`rm /usr/local/Install/apache-data/genome.imim.es/htdocs/datasets/secisaln/*png`;


$seq = &get_sequence();
&scan_for_matches();

&SECISearch();

print EFILE "#<$$># CLOSING REQUEST #<$$>#\n";
    close(EFILE);

## Remove temp files
unlink("$tmpdir/$$\_seq");
unlink("$tmpdir/$$\_hits");





############################################################################################
################                        FUNCTIONS                         ##################
############################################################################################

sub get_sequence {
    my $c=0;
    print EFILE "Getting Sequence\n";
if ($file = param('seqfile')) { $seq=join "", <$file> };

if (! $seq) { abnormal("You must provide a query sequence") };
print EFILE "#<$$># pattern = $pat\n";

###--sequence reformatting--
$seq=~s/\r/\n/g; #changing windows curret returns to UNIX new lines
@seqarr=split /\n/, $seq;
    &abnormal("Please use a FASTA formatted sequence. Click <a href=http://www.ncbi.nlm.nih.gov/blast/fasta.shtml>here</a> for a definition of the FASTA format") unless $seqarr[0]=~/>/;
    my $hh;
    my @bb; 
    for $line (@seqarr) {
    if ($line=~/^>/) {
	$line=~/>(.{1,25})/;
	$user_seq_name=$1;
	$line=~s/ |\,/_/g;
	print EFILE "NAME : $user_seq_name\n";

    }
    else {
	$line=~s/\d| //g; #removing digits and spaces
	while( $line =~ /([^ACTGUN])/ig){push @bb, "\"$1\",";}
	$hh.=$line;
    }
}
    if($bb[0]){$bb[scalar(@bb-1)] =~ s/,/\./; &abnormal("Please submit only DNA sequences in FASTA format. Your sequence contained illegal  characters: @bb")}
    print EFILE "hh : $hh\n";

#    $seq=join "\n", @seqarr;print EFILE "ss1 : $seq\n";

    if(length($hh)<100){
	$hh= "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" . "$hh" . "uuuuuuuuuuuuuuuuuuuuuuuuuuuuuu";
	$added_a=1;
	print EFILE "saka : $hh\n";
    }                                                    
    $seq=$seqarr[0] . "\n" . $hh . "\n"; 
    	print EFILE "SEQ is now :$seq\n";
unless ($seq =~ /^>/) {$seq=">".( param('name') || "unnamed")."\n".$seq};
   
  
 return($seq);

}

sub energy {


#============================================================================================================#
# usage energy({list of numbers of PatScan pattern fields (starting from 0)})				     #
# example: energy(0..11) calculates energy for the whole structure in case of classical pattern		     #
# example: energy(3..8) calculates energy for the core SECIS ATGAN..NGAN in case of classical pattern	     #
#============================================================================================================#


    open(TMP, ">$tmpdir/$$\_en");
    print TMP  ">$$\n", @a[@_], "\n", @b[@_];
    close(TMP);
    $cmd = "$RNAfold -C < $tmpdir/$$\_en 2>>$logfile";# > $tmpdir/$$\_energy";
    print EFILE "#<$$># RUNNING RNAFOLD... $tmpdir\n";
    open(CMD,"$cmd |");
    @_a = <CMD>;
    close(CMD);
    open (OUTEN, ">$tmpdir/$$\_energy");
    print OUTEN "@_a";
    print EFILE "$cmd\n";
    print EFILE "#<$$># ...RNAFOLD DONE\n";

    ($g1, $currentseq, $currentstruct, $energy)=split(/\n|\s\(/, "@_a");
    $g1=$_a[0];
    $currentseq=$_a[1];
    $_a[2] =~ /(.+)\((.+?)\)$/;
    $currentstruct=$1;
    $energy=$2;
    chomp $energy; chop $energy;
    
    print EFILE "energy : $energy \n$_a[0],$_a[1],$_a[2]\n";
    return $energy;
}
##########################################################################################################
sub structure_OK {
    my $fout=1;
    my $substruct=substr($currentstruct, length( join '', @a[0..3]),  length( join '', @a[4..7], (length ($a[8]) == 4) ? () : substr($a[8], 0, 1)   )) ;
    # Y-filter detects Y-shaped SECISes
    if ( ($substruct =~ /\)\.*\(/) && $Y_flag) { $fout=0 }
    # O-filter detects SECISes with more than 2 umparied nucleotides in a row on any strand in the 
    # first 7nt from the quartet
    my $porog=7;
    if ( (substr($substruct, 0, $porog) =~ /\.{3,}/ ||  substr($substruct, -$porog)=~ /\.{3,}/ ) && $O_flag) { $fout=0 }
    # B-filter detects SECISes that have more than 2 unpaired nt more on 5' side that on 3' side 
    # (visually bended to the left)
    ($substruct =~ /(.*)(\(\.*\))(.*)/);
    my $up=$1; my $down=$3;
    if (((($up =~ s/\./\./g)-($down =~ s/\./\./g)) < (-2)) && $B_flag) { $fout=0 }
    
#S-filter detects structures with less that 8 pairs
    if ((( $substruct =~ tr/\)/\)/ ) < 8 ) && $S_flag ) { $fout=0 }
    
    return $fout;
}

######################Imager2 block start#####################################################################
sub imager2 {
    @bold=();
    (my $seq, my $struct, $scale, $outfilename, $offx, $offy, $posbase1, $posbase2, @bold)=@_;
    
    print EFILE "#<$$># RUNNING imager2...\n";
    print EFILE "bold : @bold\n";
    
    $struct =~ s/^\s//;
    $seq =~ s/^\s+//;
    $tmp_print="seq:$seq<br>struct:$struct-<br>, scale: $scale<br> , ourname : $outfilename<br> ofx : $offx<br> ofy: $offy<br> pos1: $posbase1<br> pos2: $posbase2<br> bold: @bold<BR>";
    
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
    &psoutput();

##########################################################################################
    sub pairof_from_brackets {
	my @pairof;
	my $str = shift;
	$str =~ s/\s+//og;
	my @struct=split "", $str;
	for my $i (0..$#struct) {
	    if ($struct[$i] eq '.') { $pairof[$i]=-1 }
	    elsif ($struct[$i] eq '(') {push @left, $i}
	    elsif ($struct[$i] eq ')') {my $k=pop @left; $pairof[$k]=$i; $pairof[$i]=$k}
	    else {die "error: unrecognized symbol \"$struct[$i]\" in structure"}
	}
	return @pairof
	}
##########################################################################################
    sub drawpart {
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
	    $k++; }
	
	for ($n=$i+$k; $n<=$j-$k; $n++) {
	    if ($pairof[$n]!=(-1)) { 
		push @parts, $n;
		push @circle, $n;
		$n=$pairof[$n];
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
	    ($x[$circle[$ci]], $y[$circle[$ci]])=globalized($xlocal[$circle[$ci]], $ylocal[$circle[$ci]])
	    }
##########################
	sub globalized {
	    my $xl=$_[0];        
	    my $yl=$_[1];
	    $rotangle=atan2( $y[$j]-$y[$i], $x[$j]-$x[$i] );
	    my $xg=sqrt( $xl**2+$yl**2 )*cos( atan2($yl, $xl)+$rotangle)+$x[$i];
	    my $yg=sqrt( $xl**2+$yl**2 )*sin( atan2($yl, $xl)+$rotangle)+$y[$i];
	    return ($xg, $yg);
	}
#########################
    }
#################################################################################################
    sub position_picture {
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
    
#################################################################################################
    sub psoutput {
	print EFILE " $ghostvw -q -dNOPAUSE -sDEVICE=png256 -g1400x2600 -r300 -dBATCH -sOutputFile=$outfilename -\n";
	open(OUTF, "| $ghostvw -q -dNOPAUSE -sDEVICE=png256 -g1400x2600 -r300 -dBATCH -sOutputFile=$outfilename -");
	show_seq();
	show_connections();
	show_bounds();
	print OUTF "stroke\n showpage\n";
    }
############################
    sub show_seq {
	
	@nonbold=(0..$#x);
	$k=0;
	
	print OUTF "newpath\n /$fonttype\-Bold findfont $fontsize scalefont setfont\n";
	for $index (@bold) {
	    $xf=$x[$index]-0.6*$fontsize/2;
	    $yf=$y[$index]-0.6*$fontsize/2;
	    print OUTF "$xf $yf moveto\n ($seq[$index]) show\n";
	    
	    splice @nonbold, $index-$k, 1; 
	    $k++;
	    
	}
	
	print OUTF "newpath\n /$fonttype findfont $fontsize scalefont setfont\n";
	for $index (@nonbold) {
	    $xf=$x[$index]-0.6*$fontsize/2;
	    $yf=$y[$index]-0.6*$fontsize/2;
	    if (! defined $seq[$index]) { die "@seq\n@x\n"}
	    print OUTF "$xf $yf moveto\n ($seq[$index]) show\n";
	}
    }
#############################
    sub show_connections {
	for my $i (0..($#x-1)) {
	    $x1=$x[$i]; 
	    $y1=$y[$i];
	    $x2=$x[$i+1];
	    $y2=$y[$i+1];
	    line_bound('black', 0.25, 0.04); 
	}
    }
#############################
    sub show_bounds {
	$gam=0.2;
	
	for my $i (0..$#x) {
	    if ($i<$pairof[$i]) {
		$curpair=$seq[$i].$seq[$pairof[$i]];
		$x1=$x[$i]; $y1=$y[$i]; $x2=$x[$pairof[$i]]; $y2=$y[$pairof[$i]];
		
		if (($curpair eq 'GU') || ($curpair eq 'UG')) {circle_bound('lightblue', $gam)}
		elsif (($curpair eq 'AU') || ($curpair eq 'UA')) {circle_bound('mediumblue', $gam)}
		elsif  (($curpair eq 'GC') || ($curpair eq 'CG')){circle_bound('red', $gam)}
		else {diamond_bound('burgundy', $gam, 0.3)}
	    }
	}
    }
###########################
    sub line_bound {
	my @color=rgbcolor($_[0]);
	my $gamma=$_[1];
	my $teta=$a*$_[2]; 
	
	$xstart=$x1+($x2-$x1)*(1-$gamma)/2; 
	$ystart=$y1+($y2-$y1)*(1-$gamma)/2;
	$xstop=$x2-($x2-$x1)*(1-$gamma)/2;
	$ystop=$y2-($y2-$y1)*(1-$gamma)/2;
	print OUTF "newpath\n @color setrgbcolor\n $teta setlinewidth\n $xstart $ystart moveto\n $xstop $ystop lineto\n stroke\n";
    }
###########################
    sub circle_bound {
	my @color=rgbcolor($_[0]);
	my $gamma=$_[1];
	my $xc=($x1+$x2)/2;
	my $yc=($y1+$y2)/2;
	my $rc=$a*$gamma/2;
	print OUTF "newpath\n @color setrgbcolor\n $xc $yc $rc 0 360 arc\n fill\n stroke\n";
    }
###########################
    sub diamond_bound {
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
###########################
    sub rgbcolor {
	my @rgb;
	if ($_[0] eq 'lightblue') { @rgb=qw(0.7 0.7 0.7) }
	elsif ($_[0] eq 'mediumblue') { @rgb=qw(0 0 0.9) }
	elsif ($_[0] eq 'darkblue') { @rgb=qw(0 0 0.2) }
	elsif ($_[0] eq 'burgundy')  { @rgb=qw(0.5 0.1 0.1) }
	elsif ($_[0] eq 'black')  { @rgb=qw(0.5 0.5 0.5) }
	elsif ($_[0] eq 'red')  { @rgb=qw(0.9 0 0) }
	return @rgb;
    }
    print EFILE "#<$$># ...imager2 DONE\n";
}
###############Imager2 block end###############################################################






######################################
sub abnormal {
    print EFILE "@_";
    print header; 
   print EFILE "\nwhowhowhowhowhwohwohwohwohw\n";    
     print  <<STOPHERE;
    <HEAD><TITLE>SECISaln</TITLE>
	
<meta name="FORMATTER" content="GNU Emacs 22.1.1">
<link rel="stylesheet" type="text/css" href="http://genome.imim.es/software/secisaln/basic.css" title="SECIS">
<link REL="SHORTCUT ICON" HREF="http://genome.imim.es/software/secisaln/favicon.ico">
</HEAD>
<BODY>
 <DIV style="width: 100%; height:100%; valign:center;">

<DIV style="position:relative; top:10em; left:20em; width:44em; border: solid #006262 ; border-width: 1 1 1 1;" >
<div id="box" style="width:99%">
<p class="box"><H2>SECI<FONT COLOR=RED>S</FONT>aln  :(</H2></p>
</div>
<div id="box" style="width:89%;">
<P class="box">$_[0]</P></div>
<img src="http://genome.imim.es/software/secisaln/images/arrowleft.gif"><A href="http://genome.imim.es/software/secisaln/">Back</a> to the SECISaln web server.
</DIV>
</DIV>



STOPHERE

    print end_html;
    exit(0);
print EFILE "whowhowhowhowhwohwohwohwohw\n";   
}

######################################


sub align{

print EFILE "sub ALIGN started...\n";    
    my $kk=0;
    my %k;
    my %str;
    my %sequence;
    my $N;
    my %seq;
    my $main;
    my $main_str;
    my @Seq;
    my @struc;
    my %opts;


my $en  = "$tmpdir/$$\_energy"; 
my $hits = "$tmpdir/$$\_hits";
     


my $apical_repetition;
my $skip=0;
print EFILE "ww $prgbin/FastaToTbl $en\n";
open(EN,"$prgbin/FastaToTbl $en |");
while(<EN>)
{
    /^(.*?)\s+([AUGCN]*)\s*([x\.\(\)]*)/;
    $k{$1} = $3;
    $sequence{$1} = $2;
}

#open(STR,"$prgbin/FastaToTbl /usr/local/Install/apache-data/genome.imim.es/cgi-bin/SECISearch/tmp/21368_hits |");

open(STR, "$prgbin/FastaToTbl $hits |");
while(<STR>)
{
    /^(.*)\t(.*)/; 
    print EFILE;
    my @structure = ();
    my @units=();
  # next unless defined($k{$1});
    $N=$$;
    $a = $2;
    $a =~ s/T/U/g;
    @units = split (/\s/, $a);
    my $l =0;
    my $skip=0;
    for(my $i=0; $i < scalar(@units); $i++ )
    {
	$skip = $skip + $l;
	$l = length($units[$i]);
	$k{$N} =~ /^.{$skip}(.{$l})/;
	$structure[$i] = $1;
    }
    $str{$N} = [@structure];
    $seq{$N} = [@units];

}
my @k = keys(%str);

# #print ">struct\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.(((((((((((((((((((.......(((((((.......)))))))..))))))))))))))))))).xxxxxxxxxxxxxxxxxxxxxxx\n";
foreach my $name (keys(%str))
{   

$skip=0;
    $kk++; # just to give unique names to the sequences
    my $aln = ""; 
    $main ="";
    $main_str = "";
  
   # next if ${$str{$name}}[6] =~ /^\.\.\.\.\.\./; 

    # fix units to be uniform
    # patscan returns varying numbers of units. However, the NGAN quartet
    # is ALWAYS the 3rd from last. Concatenate all from just asfter
    # the first AUGA quartet to the NGAN and separate later.
    # this will result in ATGA <concatenated units> NGAN being in
    # position[6]
    for (my $n=6; $n<scalar(@{$seq{$name}})-4; $n++)
    {
	$main = $main . ${$seq{$name}}[$n];
	$main_str = $main_str . ${$str{$name}}[$n];	
    }
    my @ff;
    my @str;
    my $cc = 0;
    for (my $n=0; $n<scalar(@{$seq{$name}}); $n++)
    {
	if($n >5 && $n <scalar(@{$seq{$name}})-4 && $cc == 0)
	{
	    push @ff,$main;
	    push @str,$main_str;
	    $cc=1;
	}
	elsif($n >5 && $n <scalar(@{$seq{$name}})-4 && $cc == 1)
	{next }
	else
	{
	    push @ff, ${$seq{$name}}[$n];
	    push @str, ${$str{$name}}[$n];
	}

    }
    my @Seq   = @ff;

@struc = @str;


    &debug("1 $name\n  @{$seq{$name}}\n  @Seq\n  @struc\n");   
    # check for non-standard  secis, CC in apical loop
    # or even NN
    $apical_repetition = $Seq[5];
    
    ## Sometimes a SECIS with too short a stem is predicted, skip
    $struc[4] =~ /^(.{8})/;
    #print "$struc[4]\n$Seq[4]\n";
    my @kk = split(//, $1);
    my @ll = grep(/\(/, @kk);

#    next if scalar(@ll) < 6;




# Check for type I/II. If after the apical As there are "("
# then it is type 2, else type1

    if($struc[6] =~ /^[^\(]+\)/)
    {
	$type=1;

	&align_type1($name,\@Seq,\@struc);
	
    }
    else
    {
	$type = 2;
	&align_type2($name,\@Seq,\@struc);
    }

}   

### Functions ###

sub align_type2
{


    my $name = $_[0];
    my @Seq = @{$_[1]};
    my @struc = @{$_[2]};
   
    if(length($Seq[scalar(@Seq)-4]) != 4)
    {

	$Seq[scalar(@Seq)-4] =~ s/(.*)(.GA.)(.*)/$2/;
	$Seq[scalar(@Seq)-5] =  $3 . $Seq[scalar(@Seq)-5];
	$Seq[scalar(@Seq)-3] = $Seq[scalar(@Seq)-3] . $1;

    }


# Stem 2 (after AUGA)
    if($struc[4] =~ /(\.+)$/){
	my $a = $1;
	$Seq[4] =~ s/($1)$//;
	$Seq[5] = $1 . $Seq[5];
	$struc[4] =~ s/($a)$//;
	$struc[5] = $a . $struc[5];
    }

    &debug("2 @{$seq{$name}}\n  @Seq\n  @struc");

# Struc 6 : internal loop after the as, helix 3, apical loop 
# and helix2- until core quartet. This bit separates 
# the internal loop from helix3 and joins it to struc5 which is the as
#    for my $nn(0..$#struc){print EFILE "a $nn : $struc[$nn]\n"}
    if($struc[6] =~ s/^(\.+)//){
	my $a = $1;
	$Seq[6] =~ s/^($1)//;
	$Seq[5] = $Seq[5] . $1;
	$struc[5] = $struc[5] .  $a;
    }
    &debug("3 @{$seq{$name}}\n  @Seq\n  @struc");
    my $apstem1 =0 ;
    my $apstem2=0;
    
    &debug("\@struc:@struc :: -$struc[6]-");
# special case (eg Rat_SelW,gi|56388784) where an extra bp occurs after
    # the AA and before the apical stem ==> force open and add to AAs
    if($struc[6] =~ /^\(\.\(+/)
    {
	&debug("extra bp $struc[5]:$struc[6]\n         $Seq[5]:$Seq[6]\n");
	$struc[6] =~ s/^(\(\.\(+\.+\)+\.)\)\./$1\.\./;	
	$struc[6] =~ s/^\(\.//;	
	$struc[5] = $struc[5] . "..";
	$Seq[6] =~ s/^(..)//;
	$Seq[5] = $Seq[5] . $1; 
	&debug("extra bp $struc[5]:$struc[6]\n         $Seq[5]:$Seq[6]\n");
    }
    
# Here, we separate helix3+, apical loop, and helix3-
# they are still [6] but now with space between them 
   if($struc[6] =~ /^(\(+.*?)(\.{3,7})(\)+)/){
       my $apical_loop=".";
       while($struc[6] =~ /\((\.+)\)/g){$apical_loop=$1 if length($1)>length($apical_loop)}
       $apical_loop=quotemeta($apical_loop);
       
#	if(($struc[6] =~ /^(\(+)(\.+)(\)+)/)|| ($struc[6] =~ /^(\(+\.\(+)(\.+)(\)+\.{0,1}\)+)/ ) ){
       if($struc[6] =~ /^(.+?)($apical_loop)(.+?)/){
	   $struc[6] =~ /^(.+?)($apical_loop)(.+?)$/;
	  
	   $apstem1 = length($1);
	   &debug("1,2,3,apic : $1,$2,$3,$apical_loop");
	   my $b = length($2);
	   my $qq=quotemeta($2);
	   &debug("qq  :$apical_loop\n");
	   my $kiki=$apstem1+1;
	   
	   if($struc[6] =~ /$apical_loop(\){$apstem1})/){
	       $apstem2 = length($1);
	   }
	   elsif($struc[6] =~ /$apical_loop(.{$apstem1,$kiki})/){
	       $apstem2 = length($1);
	   }
	   else{die("apstem problem\n")}
	   
	   my $f_str;

	   ## If we have unpaired residues after stem2
	   if($struc[6] =~ s/^(.{$apstem1})(.{$b})(.{$apstem2})(\.+)(.*)/$1 $2 $3 $4 $5/){
	       $f_str = $4;
	       # c,d,e,f,g = apical stem+,apical loop, apicalstem-, unpaired apical, stem-
	       ($c,$d,$e,$f,$g) = ($1,$2,$3,$4,$5);
	       &debug("c,d,e,f,g(b14) :: $c , $d , $e , $f , $g");
	       $Seq[6] =~ /^(.{$apstem1})(.{$b})(.{$apstem2})($f)(.*)/;
	       ($c,$d,$e,$f,$g) = ($1,$2,$3,$4,$5);
	   }
	   elsif($struc[6] =~ s/(.{$apstem1})(.{$b})(.{$apstem2})(.*)/$1 $2 $3 $4/){
	       $f_str = $4;
	       $Seq[6] =~ /^(.{$apstem1})(.{$b})(.{$apstem2})(.*)/;
	       # c,d,e,f = apical stem+,apical loop, apicalstem-, unpaired apical plus stem-
	       ($c,$d,$e,$f) = ($1,$2,$3,$4);
	       defined($g) ? &debug("c,d,e,f,g(b14) :: $c , $d , $e , $f , $g") :  &debug("c,d,e,f(b14) :: $c , $d , $e , $f") ;
	   }
	   else{die("apical loop parsing problem\n");}
	   &writeX(length($c),7,\$c);
	   &writeX(length($d),13,\$d);
	   
	   if(defined($g)){
	       &writeX(length($e),7,\$e,1);
	       &writeX(length($f),4,\$f);
	       &writeX(length($g),12,\$g,1); 

	   }
	   else{
	       &writeX(length($e),7,\$e,1); ## was 8
	       &writeX(length($f),12,\$f,1); 
	   }
	   #&writeX(length($f),13,\$f);
	   # $f is the 2nd helix- AND if present the rest of the
	   # apical loop, so deal with that
	   &debug("f_str : $f_str : $struc[6]\n");
	   &debug("c,d,e,f     :: $c,$d,$e,$f");
#	    if ($debug){for (my $hoho=0; $hoho<scalar(@struc); $hoho++){ print EFILE "struc[$hoho]:$struc[$hoho]\nsssec[$hoho]:$Seq[$hoho]\n"}}
	   ## If there ARE unpaired residues after stem3
	   unless(defined($g))
	   {
	       $Seq[7] = $f . $Seq[7];
	       $f = "----";
	   }
	  
	   if(defined($g)){
	       
	       my $kk=length($f);
	       $Seq[6] =~ s/^(.{$apstem1})(.{$b})(.{$apstem2})(.{$kk})(.*)/$c$d$e$f/;
	       $Seq[7] = $g . $Seq[7];
	       &debug("1c,d,e,f,g     :: $c,$d,$e,$f,$g,$apstem1,$b,$apstem2,\n6,7:$Seq[6],$Seq[7]\n");
	   }
	   else{
	       $Seq[6] =~ s/^(.{$apstem1})(.{$b})(.{$apstem2})(.*)/$c$d$e$f/;
	   }
	   print EFILE "66a66 : $Seq[6]\n"if $debug;
       
   }
    else{die("xxxxx @Seq\nxxxxx @struc\n$Seq[6]\n$struc[6]\n$name\n"); }
       
       &debug( "4 @{$seq{$name}}\n  @Seq\n  @struc");
       ## need both sides of the apical stem to have the same length
       # if($apstem1 > $apstem2)
#        {die();
# 	   my $dif = $apstem1 - $apstem2;
# 	   my @a = split(/\s+/, $Seq[6]);
# 	   die("aa : @a\n$Seq[6] : $a[2]") unless $a[2];
# 	   $struc[7] =~ s/^.{$dif}//;
# 	   $Seq[7] =~ s/^(.{$dif})//;
# 	   &debug("$a[2] : $1\n@a");
# 	   $a[2] = $a[2] . $1;
# 	   $Seq[6] = $a[0] . " " . $a[1] . " "  . $a[2]; 
	   
#        }
       
       if($struc[7] =~ /^(\.+)/){
	   &debug("7 was : -$Seq[7] $struc[7]($1)");
	   my $h = length($1); 
	   $Seq[7] =~ s/^.{$h}//;
	   &debug( "7 is  : -$Seq[7]");
	   $Seq[6] = $Seq[6] . " " . $1;
	   
	   
       }
       
       $a = $Seq[4] . " " . $Seq[5] . " " . $Seq[6] . " " . $Seq[7];
       my $b =  $struc[4] . " " . $struc[5] . " " . $struc[6] . $struc[7];
       &debug("5 @{$seq{$name}}\n  @Seq ($Seq[0]\n  @struc");
# now, align the damn things

       if(length($Seq[0]) + length($Seq[1]< 17))
       {
	   &writeX((length($Seq[0]) + length($Seq[1])),17,\$Seq[0],1);
       }




       if(length($Seq[2]) < 13)
       {
	   &writeX(length($Seq[2]),13,\$Seq[0],1);
       }
       if(length($Seq[4]) < 13)
       {
	   &writeX(length($Seq[4]),13,\$Seq[4]);
       }
       
       ## Apical loop As
       if(length($Seq[5]) < 8)
       {
	   my $beforeAs=undef;
	   if($Seq[5]=~ /^(.+?)($apical_repetition)(.*)/){
	       $beforeAs=$1; 
	   }
	   
	   if($Seq[5] =~ /^(.+?)($apical_repetition)(.*)/)
	   {
	       $Seq[5] =~ /^(.+?)($apical_repetition)(.*)/;
	       my $a = $1;
	       my $b = $3;
	       my $l;
	       length($b)>length($a) ? $l = length($b) : $l = length($a);
	       &writeX(length($a),4,\$a,1);
	       print EFILE "bbbbbb1 :$b, $Seq[5] : $apical_repetition\n" if $debug;
	       &writeX(length($b),3,\$b);
	       print EFILE "bbbbbb2 :$b\n" if $debug;
	       $a = $a . "$apical_repetition";
	       $Seq[5] =~ s/^(.+?)($apical_repetition).*/$a$b/;
	       &debug("AAs are now: $Seq[5]\n");
	   }
	   elsif(defined($beforeAs)){
	            
	       $Seq[5] = "----" . $Seq[5]; 
	       $Seq[5] =~ /^.*$apical_repetition(.*)/;
	       my $a = $1;
	       &writeX(length($a),3,\$a);
	       $Seq[5] =~ s/$apical_repetition($1)/$apical_repetition$a/;
	      
	   }
	   else
	   {
	       
 	       $Seq[5] =~ /^(.*)$apical_repetition(.*)/;
 	       my ($a,$b) = ($1,$2);
	       &writeX(length($a),4,\$a,1);
 	       &writeX(length($b),3,\$b);
	       print EFILE "a,b : $a,$b\n" if $debug;
 	       $Seq[5] = $a.$apical_repetition.$b;
#	       &writeX(length($Seq[5]),9,\$Seq[5]);
	       
	   }
	   
	   
       }
       
       # Some seqs, eg SPS2 gi|26097414 are predicted with a huge
       # apical loop and screw everything up so skip. If i do this with 
       # next, i get a warning message so i use this variable
       elsif(length($Seq[5]) > 8)
       {

	   $skip = 1; 
	   &abnormal("Weird sequence, send me a mail: charles.chapple AT crg.es")
       }
       else{
	   $Seq[5] =~ /^(.+?)($apical_repetition)(.*)/;
	   &debug("Seq5,6 : $Seq[5],$Seq[6] :: $1,$2,$3");
	   #$Seq[5] = "-" . $Seq[5]; 
	   $Seq[6] = "-" . $Seq[6]; 
	   #$Seq[6] =~ s/^-//;
	   &debug("Seq5,6 : $Seq[5],$Seq[6]");
	   
       }
       if(defined($g)){
       ## Stem, 3'
	  
	  if(length($Seq[7]) < 16)
	   {
	       &writeX(length($Seq[7]),16,\$Seq[7],1);
	   }
	
	   if(length($Seq[10]) < 21)
	   {
	       &writeX(length($Seq[10]),21,\$Seq[10]);
	   }
       }
       else{
	   #$Seq[7]="----" . $Seq[7];
	  &debug("SEQ : @Seq\n");
	   if(length($Seq[7]) < 14)
	   {
	    &writeX(length($Seq[7]),14,\$Seq[7],1);
	   }
	    
	   if(length($Seq[8]) < 4)
	   {
	       
	       #   &writeX(length($Seq[8]),4,\$Seq[8]);
	   }
	   if(length($Seq[9]) < 11)
	   {
#	    &writeX(length($Seq[9]),11,\$Seq[9],);
	   }
	   if(length($Seq[10]) < 16)
	   {
	       &writeX(length($Seq[10]),16,\$Seq[10]);
	   }
       }
       if(length($Seq[0]) < 4)
       {
	   print EFILE "WHOHAHAA : $Seq[0]\n";
	   
       }

       
      #  my $whole_seq;
# 	map{$whole_seq.=$_;}@Seq;
# 	print EFILE "aaaaaaaaaaaa", length($whole_seq), "$whole_seq\n@Seq\n"; 
# 	if(length($whole_seq) < 133){
# 	    my $a=133-length($whole_seq);
# 	    &writeX(length($Seq[10]),$a,\$Seq[10]);
	   
# 	}
       ## Sometimes there are too many gaps after end of seq
       my $whole_length;
       
       my $ll;
       map{$whole_length+=length($_); $ll.=$_;}@Seq;
       if($whole_length > 132){
	   my $extra= $whole_length - 132;
	   my $jo=scalar(@Seq)-1;
	   
	   if($Seq[scalar(@Seq)-1] =~ s/-{0,$extra}$//){}
	    else{print EFILE "PROBLEM : $Seq[scalar(@Seq)-1]: $extra : $whole_length\n@Seq\n"; die();}
	   #karo
	   
       }
       elsif($whole_length < 132){
	   my $extra= 132-$whole_length ;
	   my $jo=scalar(@Seq)-1;
	   for(my $i=0; $i<$extra;$i++){
	       $Seq[scalar(@Seq)-1] .="-";
	   }

	   # else{print "PROBLEM : $Seq[scalar(@Seq)-1]: $extra : $whole_length\n@Seq\n"; die();}
	   #karo
	   
       }
      
       
       &print_aligned_output(@Seq); 

       
     
       
   }
    else{print EFILE "SHIT struc[6] :$struc[6] :: @struc\n@Seq"; &abnormal("An error occured, please contact the webmaster charles.chapple AT crg DOT es\n")}
    
    
}


sub align_type1
{
    my $name = $_[0];
    my @Seq = @{$_[1]};
    my @struc = @{$_[2]};


    if(length($Seq[scalar(@Seq)-4]) != 4)
    {
	$Seq[scalar(@Seq)-4] =~ s/(.*)(.GA.)(.*)/$2/;
	$Seq[scalar(@Seq)-5] =  $3 . $Seq[scalar(@Seq)-5];
	$Seq[scalar(@Seq)-3] = $Seq[scalar(@Seq)-3] . $1;
    }
## Stem 2 (after AUGA)
    if($struc[4] =~ /(\.+)$/){
	my $a = $1;
	$Seq[4] =~ s/($1)$//;
	$Seq[5] = $1 . $Seq[5];
	$struc[4] =~ s/($a)$//;
	$struc[5] = $a . $struc[5];
    }
    &debug("2 @{$seq{$name}}\n  @Seq\n  @struc");

## Struc 6 : internal loop after the as, helix 3, apical loop 
## and helix2- until core quartet. This bit separates 
## the internal loop from helix3 and joins it to struc5 which is the as
    if($struc[6] =~ s/^(\.+)//){
	my $a = $1;
	$Seq[6] =~ s/^($1)//;
	$Seq[5] = $Seq[5] . $1;
	$struc[5] = $struc[5] .  $a;
    }
    &debug("3 @{$seq{$name}}\n  @Seq\n  @struc");
    my $apstem1 =0 ;
    my $apstem2=0;
    

# special case (eg ENSMUSG00000041571) where an extra bp occurs after
    # the AA and before the apical stem ==> force open and add to AAs
    if($struc[6] =~ /^\(\.\(+/)
    {
	$struc[6] =~ s/^\(\.//;	
	$struc[5] = $struc[5] . "..";
	$Seq[6] =~ s/^(..)//;
	$Seq[5] = $Seq[5] . $1; 
    }
    
# now, align the damn things



    if(length($Seq[0]) < 10)
    {
	&writeX(length($Seq[0]),10,\$Seq[0],1);
    }
     if(length($Seq[1]) < 7)
    {
	&writeX(length($Seq[1]),7,\$Seq[0],1);
    }
    if(length($Seq[2]) < 13)
    {
	&writeX(length($Seq[2]),13,\$Seq[0],1);
    }
    if(length($Seq[4]) < 13)
    {
 	&writeX(length($Seq[4]),13,\$Seq[4]);
    }
    if(length($Seq[6]) < 14)
    {
	&writeX(length($Seq[6]),14,\$Seq[6],1);
    }
  #  if(length($Seq[5]) < 19)
   # {
    &debug("xxxxxx $Seq[5] $Seq[6] Seq[7]");
    if($Seq[5] =~ /^(.+?)($apical_repetition)(.*)/)
    {
	$Seq[5] =~ /^(.+?)($apical_repetition)(.*)/;
	
	my ($l1,$l2,$l3)=(length($1),length($2),length($3));
	my $a = $1;
	my $b = $3;
	
	$struc[5] =~ /(.{$l1})(.{$l2})(.{$l3})/;
	&writeX(length($a),2,\$a,1);
	&writeX(length($b),16,\$b);
	$a = $a . "$apical_repetition";
	$Seq[5] =~ s/^(.+?)($apical_repetition).*/$a$b/;
	# must fix

    }
    elsif($Seq[5] =~ /^($apical_repetition)(.+)/)
    {	
	$Seq[5]="--" . $Seq[5];
	&writeX(length($Seq[5]),20,\$Seq[5]);

    }
    else
    {

	$Seq[5] = "--" . $Seq[5]; 
	$Seq[5] =~ /^.*$apical_repetition(.*)/;
	my $a = $1;
	&writeX(length($a),16,\$a);
	$Seq[5] =~ s/$apical_repetition($1)/$apical_repetition$a/;

    }
    
    
    #&writeX(length($Seq[5]),20,\$Seq[5]);
    #}
    
    if(length($Seq[8]) < 12)
    {
#	&writeX(length($Seq[8]),12,\$Seq[8]);
    }
    if(length($Seq[9]) < 13)
    {
	#&writeX(length($Seq[9]),8,\$Seq[9]);
    }



    my $whole_length;
       my $ll;
       map{$whole_length+=length($_); $ll.=$_;}@Seq;
       
       if($whole_length > 116){
	   my $extra= $whole_length - 116;
	   my $jo=scalar(@Seq)-1;
	   if($Seq[scalar(@Seq)-1] =~ s/-{0,$extra}$//){}
	    else{print EFILE "PROBLEM : $Seq[scalar(@Seq)-1]: $extra : $whole_length\n@Seq\n"; die();}
	   #karo
	   
       }
       elsif($whole_length < 116){
	   
	   my $extra= 116-$whole_length ;
	   my $jo=scalar(@Seq)-1;
	 
	   for(my $i=0; $i<$extra;$i++){
	       $Seq[scalar(@Seq)-1] .="-";
	   }
       }


    &print_aligned_output(@Seq); 

   
}
sub writeX
{
    
    &debug("#<$$># WritingX $_[0], $_[1], ${$_[2]}, $_[3] \n");
    my $x = $_[0];
    my $l = $_[1];
    my $res = ${$_[2]};
    my $difference = $l-$x;
    ## If there is a 4th parameter passed,
    ## place "-" to the left of the sequence, 
    ## and else to the right.
    if($_[3])
    {
	for(my $i=$x; $i<$l; $i++)
	{
	    ${$_[2]} = "-" . ${$_[2]};	    
	}
    }
    else
    {
	for(my $i=$x; $i<$l; $i++)
	{
	    ${$_[2]} = ${$_[2]} . "-";	    
	}
    }

}
sub writeX1
{
    
    &debug("#<$$># WritingX $_[0], $_[1], ${$_[2]}, $_[3] \n");
    my $x = $_[0];
    my $l = $_[1];
    my $res = ${$_[2]};
    my $difference = $l-$x;
    ## If there is a 4th parameter passed,
    ## place "-" to the left of the sequence, 
    ## and else to the right.
    if($_[3])
    {
	for(my $i=$x; $i<$l; $i++)
	{
	    ${$_[2]} = "*" . ${$_[2]};	    
	}
    }
    else
    {
	for(my $i=$x; $i<$l; $i++)
	{
	    ${$_[2]} = ${$_[2]} . "*";	    
	}
    }

}

sub debug
{
    if ($debug)
    {
	print EFILE "@_\n";
    }
}
### patscqan pattern (for reference)
    #r1={au,ua,gc,cg,gu,ug} NNNNNNNNNN p1=7...7 3...13 ATGAN p2=10...13 AA (4...12 | 0...3 p3=3...6 3...6 r1~p3 0...3) (r1~p2[2,1,1] NGAN | r1~p2[2,1,0] NNGAN) 3...10 r1~p1[1,1,1] NNNNNNNNNN  ENSRNOG00000013548    ENSPTRG00000023837


# structure
#>struct
#xxxxxxxxxx-xxxxxxx-xxxxxxxxxxxxx---.((((-(((((((((((((--.......--(((((((---.......--)))))))-..---)))))))))))))-----))))-.xxxxxxxxxxxxxxxxxxxxxxx 
print EFILE "#<$$># ...aligning DONE\n";
}

sub print_aligned_output{

    ## If this is a type1 SECIS
    if($type==1){
	$css = "http://genome.imim.es/software/secisaln/secis1.css";
	# if($want_gis){ $css = "http://genome.imim.es/datasets/secisaln/secis2.css"}
# 	else{$css = "http://genome.imim.es/datasets/secisaln/secis1.css"}
	
	    print EFILE "#<$$># printing output...\n";

	print header;
	print "<!--DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\" \n\"http://www.w3.org/TR/html4/strict.dtd\">-->\n";
print <<STOPHERE;
<HTML><HEAD><TITLE>SECISaln</TITLE>
<link href="$css" rel="stylesheet" type="text/css">
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=iso-8859-1">

<META NAME="keywords"    CONTENT="SECISaln,SECIS,selenocystein,Selenoprotein,secisearch,selenocysteine insertion sequence">
  <META NAME="description" CONTENT="This is the webserver for SECISaln, a program that produces structural alignments of eukaryotic SECIS elements">
<meta name="FORMATTER" content="GNU Emacs 22.1.1">
<style type="text/css" title="Genome" media="print">

a { text-decoration: none; color: #006262; }
 </style>
 <script language="javascript" type="text/javascript">
  
  if (window != top) { top.location.href = location.href; }
 </script>
<link REL="SHORTCUT ICON" HREF="http://genome.imim.es/software/secisaln/favicon.ico">



<SCRIPT type="text/javascript" src="http://genome.imim.es/Genome.js"></SCRIPT>

 <!--[if lte IE 7]>

<link href="http://genome.imim.es/software/secisaln/styleie.css" rel="stylesheet" type="text/css">


<![endif]-->
<script type="text/javascript">
      function OpenWindow(Url){
      WindowOpen = window.open(Url, "", "width=550,height=680")
      }

var orig;
function changeme(id,gi){
elmnt=document.getElementById(id);
orig=elmnt.innerHTML;
elmnt.innerHTML = gi;
}

function goback(id,old){
elmnt=document.getElementById(id);
elmnt.innerHTML=old;
}


</script>
<script language="Javascript">
<!--
function toggleDiv(id,flagit) {
if (flagit=="1"){
if (document.layers) document.layers[''+id+''].visibility = "show"
else if (document.all) document.all[''+id+''].style.visibility = "visible"
else if (document.getElementById) document.getElementById(''+id+'').style.visibility = "visible"
}
else
if (flagit=="0"){
if (document.layers) document.layers[''+id+''].visibility = "hide"
else if (document.all) document.all[''+id+''].style.visibility = "hidden"
else if (document.getElementById) document.getElementById(''+id+'').style.visibility = "hidden"
}
}
//-->
</script>
</HEAD><BODY>

<A NAME="TOC"></A>
<script src="http://www.google-analytics.com/urchin.js" type="text/javascript">
</script>
<script type="text/javascript">
_uacct = "UA-964675-1";
urchinTracker();

</script>
<div id="dhtmltooltip"></div>

<script type="text/javascript">

/***********************************************
* Freejavascriptkit.com
* Visit http://www.freejavascriptkit.com for more free Javascripts source code
***********************************************/
var he=1;
var offsetxpoint=0 //Customize x offset of tooltip
var offsetypoint=10 //Customize y offset of tooltip
var ie=document.all
var ns6=document.getElementById && !document.all
var enabletip=false
if (ie||ns6)
var tipobj=document.all? document.all["dhtmltooltip"] : document.getElementById? document.getElementById("dhtmltooltip") : ""

function ietruebody(){
return (document.compatMode && document.compatMode!="BackCompat")? document.documentElement : document.body
}


function ddrivetip(thetext, thecolor, thewidth){
if(he==1){	
if (ns6||ie){
if (typeof thewidth!="undefined") tipobj.style.width=thewidth+"px"
if (typeof thecolor!="undefined" && thecolor!="") tipobj.style.backgroundColor=thecolor
tipobj.innerHTML=thetext
enabletip=true
return false
}
}
    }

    function positiontip(e){
	if (enabletip){
	    var curX=(ns6)?e.pageX : event.clientX+ietruebody().scrollLeft;
	    var curY=(ns6)?e.pageY : event.clientY+ietruebody().scrollTop;
//Find out how close the mouse is to the corner of the window
var rightedge=ie&&!window.opera? ietruebody().clientWidth-event.clientX-offsetxpoint : window.innerWidth-e.clientX-offsetxpoint-20
var bottomedge=ie&&!window.opera? ietruebody().clientHeight-event.clientY-offsetypoint : window.innerHeight-e.clientY-offsetypoint-20

var leftedge=(offsetxpoint<0)? offsetxpoint*(-1) : -1000

//if the horizontal distance isn't enough to accomodate the width of the context menu
if (rightedge<tipobj.offsetWidth)
//move the horizontal position of the menu to the left by it's width
tipobj.style.left=ie? ietruebody().scrollLeft+event.clientX-tipobj.offsetWidth+"px" : window.pageXOffset+e.clientX-tipobj.offsetWidth+"px"
else if (curX<leftedge)
tipobj.style.left="5px"
else
//position the horizontal position of the menu where the mouse is positioned
tipobj.style.left=curX+offsetxpoint+"px"

//same concept with the vertical position
if (bottomedge<tipobj.offsetHeight)
tipobj.style.top=ie? ietruebody().scrollTop+event.clientY-tipobj.offsetHeight-offsetypoint+"px" : window.pageYOffset+e.clientY-tipobj.offsetHeight-offsetypoint+"px"
else
tipobj.style.top=curY+offsetypoint+"px"
tipobj.style.visibility="visible"
}
    }

    function hideddrivetip(){
	if (ns6||ie){
enabletip=false
tipobj.style.visibility="hidden"
tipobj.style.left="-1000px"
tipobj.style.backgroundColor=''
tipobj.style.width=''
}
    }

document.onmousemove=positiontip

</script>

<div align="center">
<div id="wrapper">
<!-- ############### TITLE AREA ############### -->

<div id="header">
<h1>GENOME BIOINFORMATICS RESEARCH LAB</h1>
<small>
Grup de Recerca en Inform&agrave;tica Biom&egrave;dica <br class="hhh"/>
Institut Municipal d\'Investigaci&oacute; M&egrave;dica + Universitat Pompeu Fabra + Centre de Regulaci&oacute; Gen&ograve;mica<br class="hh"  \>

</small>
<div align="center"><div id="menudiv"><table class="gblnav"><tr class="menu"><td id="Mhelp"><a href="http://genome.imim.es/software/secisaln/help.html" title="Ask for Help">Help</a></td><td id="Mnews"><a href="http://genome.imim.es/main/news.php" title="News">News</a></td><td id="Mpeople"><a href="http://genome.imim.es/main/people.php" title="Who is who in Genome Informatics Research Laboratory">People</a></td><td id="Mresearch"><a href="http://genome.imim.es/research/" title="Research at our Group">Research</a></td><td id="Msoft"><a href="http://genome.imim.es/software/" title="Software developed by our Group">Software</a></td><td id="Mpapers"><a href="http://genome.imim.es/publications/" title="Publications by our Group">Publications</a></td><td id="Mlinks"><a href="http://genome.imim.es/main/links.php" title="Interesting Links">Links</a></td></tr></tbody></table>

<table class="gblnav"><tbody><tr class="menu"><td id="Mdata"><a href="http://genome.imim.es/datasets/" title="Datasets created/curated by our Group">Resources&nbsp;&amp;&nbsp;Datasets</a></td><td id="Mpreds"><a href="http://genome.imim.es/genepredictions/" title="Gene Predictions produced by our Software">Gene&nbsp;Predictions</a></td><td id="Mseminars"><a href="http://genome.imim.es/seminars/" title="Seminars, Courses and Group Sessions">Seminars&nbsp;&amp;&nbsp;Courses</a></td></tr></tbody></table></div></div>
<div class="alt"><hr /></div>

</div></div></div>
STOPHERE

if($added_a==1){print "<DIV onMouseOver=\"ddrivetip(\'Since the submitted sequence was too short for prediction of helix 1 (less than 100nt), 30 &quot;a&quot; were added to the beginning of the sequence and the complementary 30 &quot;u&quot at the end. These are shown in lowercase in the aligned output.\',\'\',\'300\')\" OnMouseOut=\"hideddrivetip()\" id=\"warning\"><blink>IMPORTANT:</blink> Your sequence was too short for a reliable SECIS prediction. Click <a href=\"http://genome.imim.es/software/secisaln/help.html#long\" style=\"color:white; font-weight:bold;\">here</a> for more information.  </DIV>\n"}
    print <<STOPHERE;
<table class="holder"  border="0" cellpadding="0" cellspacing="0">
<tr class="nohover">
<td class="nav2">


<div id="nav2"><a href="http://www.imim.es/" title="Institut Municipal d&#39;Investigaci&oacute; M&egrave;dica">IMIM</a>&lt;<a href="http://www.upf.edu/" title="Universitat Pompeu Fabra">UPF</a>&lt;<a href="http://www.crg.es/" title="Centre de Regulaci&oacute; Gen&ograve;mica">CRG</a>&lt;<a href="http://nemo.imim.es/grib/" title="Grup de Recerca en Inform&agrave;tica Biom&egrave;dica">GRIB</a>&lt;<a href="http://genome.imim.es/" title="Genome Bioinformatics Research Laboratory HOME PAGE: You are Welcome!!!">HOME</a>&gt;<a href="http://genome.imim.es/software/" title="Software developed in our Group">Software</a>&gt;<a href="http://genome.imim.es/software/secisaln" title="SECISaln">SECISaln</a></div>
</td>
<td class="explanation"><div id ="explanation">Slide your mouse over the "S" to see the structure predicted by <a href="http://genome.unl.edu/SECISearch.html">SECISearch</a><br>Click on the "S" to see a larger image, click on the "P" for PDF format.</div></td>
<td class="explanation"><img src="http://genome.imim.es/software/secisaln/images/arrow1.gif" height="20" width="6">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>


<td class="ambiguity"><a class="ambig" href="http://genome.imim.es/software/secisaln/ambig.html">IUPAC-IUB/GCG Ambiguity Codes</a><td>
</tr><tr class="nohover"><td colspan="3">

 <table class="legend_top2" border="0"  cellpadding="0" cellspacing="0">
	<tbody><tr class="nohover"><td class="help">Help: On<input CHECKED name="help" type="radio" Onclick="he=1"> Off<input name="help" type="radio" Onclick="he=0" ></td>
STOPHERE


if(@name){## typeI
    print EFILE "**********************Got to here, type 1!!!\n" if $debug;
    unless ($skip==1){


	      $png_image = "$pngpth/$$" . ($hitnumber-1);
	      my $local_png="$pdir/$$" . ($hitnumber-1);
	      my $tmp_png=$local_png . ".png";
	      system("convert -trim /usr/local/Install/apache-data/genome.imim.es/htdocs/software/secisaln/36320.png ~/test.png 2>> ~/error");
	      my $PNG="$pdir/$$" . ($hitnumber-1);
	      $local_png=~/^.+\/(.+)/;
	      system("echo \"convert -trim $PNG.png $PNG.trimmed.png; convert $PNG.trimmed.png $PNG.trimmed.pdf  2>>$prgbin/error\" > $prgbin/command; echo \"see : $PNG.png  -trim $PNG.trimmed.png\" > $prgbin/worked; at now < $prgbin/command");
	      $png_image = "$pngpth/$$" . ($hitnumber-1) . ".trimmed.png";
	      $pdf_image = "$pngpth/$$" . ($hitnumber-1) . ".trimmed.pdf";
	     
	      open(T,"$tmpdir/$$\_energy");
	      my $structure;
	      while(<T>) {
		  next unless /\(/;
		  /^\s*([\.\(\)]+)/;
#		  @pred_str =split(//,$1);
		  $structure=$1;
	      }
	      # Using the aligned sequence units as a template
	      # add "-" to the corresponding structure units
	      # so as to print the predicted structure at the
	      # top of the web output
	      foreach my $seq_unit(@_){  ## typeI
		  $seq_unit=~/^(-*?)([^-]*)(-*)$/;
		  my $a=$1; my $b=$3;
		  my $l=length($2);
		  $structure=~s/(.{$l})//;
		  my $c=$a.$1.$b;
		  push @pred_str, $c;
	      }
	      my @tmp_array;
	      my $c=0;
	      my $b;
	      my $kk;
	      map{$kk.=$_;}@pred_str;
	      my @pred_str_split = split(//,$kk);
	      foreach my $res (@pred_str_split) {		 
		  if($res =~/\(/) {
		      $c++;
		      my $lid= "l".$c;
		      my $rid= "r".$c;
		      $b= "<td id=\"$lid\" class=\"whi\" onmouseover=\"this.className=\'red\'; document.getElementById (\'$rid\').className = \'red\'\" onmouseout=\"this.className=\'whi\'; document.getElementById (\'$rid\').className = \'whi\'\">" . $res . "</td>"; 
		      push @tmp_array, $b;
		  }
		  elsif ($res =~ /\)/){
		      
		      my $lid= "l".$c;
		      my $rid= "r".$c;
		      $b ="<td class=\"whi\" id=\"$rid\" onmouseover=\"this.className=\'red\'; document.getElementById (\'$lid\').className = \'red\'\" onmouseout=\"this.className=\'whi\'; document.getElementById (\'$lid\').className = \'whi\'\">" . $res . "</td>";
		      $c--;
		      push @tmp_array, $b;
		  }
		  else{
		      $b ="<td>" . $res . "</td>";
		      push @tmp_array, $b;
		  }
	      }
	      @pred_str=@tmp_array;

	      
	  }
	  
	  open(A, "$pdir/template_table_top1");
	  my @tt=<A>;
	  print "@tt\n";
	  close(A);
	  
	  my $a;
	  map{$a .= $_;}@_;
	  
	  ## Fix colors of input seq
	  my @tt =  &generate_table("type1_nostr.secisaln");
	  %hash = %{$tt[0]};
	  %gihash = %{$tt[1]};
	  @sequences = @{$tt[2]};
	  %col=%{&fix_colors(@sequences)};
	  print "@pred_str";
	  print "<td class=\"lastcell\">&nbsp;</td><td class=\"lastcell\">&nbsp;</td></tr><tr class=\"nohover\"><td class=\"one\" valign=\"top\" ><span class=\"name\">$user_seq_name</span></td>";
	  my @ss = split(//, $a);
	  for (my $n=0; $n<scalar(@ss); $n++){
		  print "<td class=\"user_seq\" bgcolor=\"$col{$ss[$n]}{$n}\">$ss[$n]</td>\n"
		  }
    
## generate image output page
#     my $png_page="http://genome.imim.es/software/secisaln/png_page.html";
#     open(A,">$pdir/png_page.html");
#     print A "<HTML><BODY><a href=\"\"><img src=\"$png_image\" height=\"\"></a><br><a href=\"$pdf_image\" >PDF format</BODY></HTML>\n";
#     close(A);
    
	  print <<STOPHERE;
<td class="secisimage" style="border-width: 1px 2px 0px 1px"><a class="secisimage" name="bob" href="javascript:OpenWindow('$png_image')" onmouseover="changeme('legend2','<img src=$png_image id=legend2>')" onmouseout="goback('legend2','<img src=http://genome.imim.es/software/secisaln/images/$legend_image id=legend1N height=$lgnd_height width=$lgnd_width>     ')">S</a></td>
<td class="secisimage" style="border-width: 1px 2px 0px 1px"><a class="secisimage" name="bob" href="javascript:OpenWindow('$pdf_image')">P</a></td>
</tr></tbody></table>
</td>
<td rowspan="2"  class="legend">
<div id="legend2">
<img src="http://genome.imim.es/software/secisaln/images/$legend_image" id="legend1N" height="$lgnd_height" width="$lgnd_width" alt="" style="position:relative; left:3em; top:1em; z-index:100;">
</div>

</td></tr>

<tr>

  <td class="align" colspan="3">
    <div class="scrollArea">
      <table class="align2" width="" border="0" cellpadding="0" cellspacing="0"><tbody>  
      
STOPHERE
	
	      &html_out();
	      print end_html;
	      
	  
}
	else {
	    print 'Sorry, no SECIS elements were found! <BR>';
	       if ($pat ne "pat_Sep20") {
		   print "You may want to try another pattern" 
	       } 
	}
    }
    
    ### If SECIS is type2
	else {
	# if($want_gis){ $css = "http://genome.imim.es/datasets/secisaln/secis2.css"}
# 	else{$css = "http://genome.imim.es/datasets/secisaln/secis1.css"}
	$css = "http://genome.imim.es/software/secisaln/secis1.css";
	print EFILE "#<$$># printing output...\n";
	print header;
print <<STOPHERE;


<html><head><!--DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN" "http://www.w3.org/TR/html4/strict.dtd">-->
<title>SECISaln</title>
<link href="http://genome.imim.es/software/secisaln/secis1.css" rel="stylesheet" type="text/css">
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">

<meta name="keywords" content="SECISaln,SECIS,selenocystein,Selenoprotein,secisearch,selenocysteine ins
ertion sequence">
  <meta name="description" content="This is the webserver for SECISaln, a program that produces structural
 alignments of eukaryotic SECIS elements">
<meta name="FORMATTER" content="GNU Emacs 22.1.1">
<style type="text/css" title="Genome" media="print">
body { background-color: #FFFFFF; }
	   a { text-decoration: none; color: #006262; }
 </style>
 <script language="javascript" type="text/javascript">

 if (window != top) { top.location.href = location.href; }
 </script>


<link rel="SHORTCUT ICON" href="http://genome.imim.es/software/secisaln/favicon.ico">
 <script type="text/javascript">
 function OpenWindow(Url){
      WindowOpen = window.open(Url, "", "width=450,height=620")
      }
var orig;
function changeme(id,gi){
elmnt=document.getElementById(id);
orig=elmnt.innerHTML;
elmnt.innerHTML = gi;
}
function goback(id,old){
elmnt=document.getElementById(id);
elmnt.innerHTML=old;
}
</script>
</head><body>
<script type="text/javascript" src="http://www.google-analytics.com/urchin.js" type="text/javascript">
</script>
<script type="text/javascript">
_uacct = "UA-964675-1";
	       urchinTracker();
 
</script>
<div id="dhtmltooltip"></div>

<script type="text/javascript">

/***********************************************
* Freejavascriptkit.com
* Visit http://www.freejavascriptkit.com for more free Javascripts source code
***********************************************/
var he=1;
var offsetxpoint=0 //zzCustomize x offset of tooltip
var offsetypoint=10 //zzCustomize y offset of tooltip
var ie=document.all
var ns6=document.getElementById && !document.all
var enabletip=false
if (ie||ns6)
var tipobj=document.all? document.all["dhtmltooltip"] : document.getElementById? document.getElementById("dhtmltooltip") : ""

function ietruebody(){
return (document.compatMode && document.compatMode!="BackCompat")? document.documentElement : document.body
}

function ddrivetip(thetext, thecolor, thewidth){
if(he==1){
if (ns6||ie){
if (typeof thewidth!="undefined") tipobj.style.width=thewidth+"px"
if (typeof thecolor!="undefined" && thecolor!="") tipobj.style.backgroundColor=thecolor
tipobj.innerHTML=thetext
enabletip=true
return false
}}
}

function positiontip(e){
if (enabletip){
var curX=(ns6)?e.pageX : event.clientX+ietruebody().scrollLeft;
var curY=(ns6)?e.pageY : event.clientY+ietruebody().scrollTop;
//Find out how close the mouse is to the corner of the window
var rightedge=ie&&!window.opera? ietruebody().clientWidth-event.clientX-offsetxpoint : window.innerWidth-e.clientX-offsetxpoint-20
var bottomedge=ie&&!window.opera? ietruebody().clientHeight-event.clientY-offsetypoint : window.innerHeight-e.clientY-offsetypoint-20

var leftedge=(offsetxpoint<0)? offsetxpoint*(-1) : -1000

//if the horizontal distance isn't enough to accomodate the width of the context menu
if (rightedge<tipobj.offsetWidth)
//move the horizontal position of the menu to the left by it's width
tipobj.style.left=ie? ietruebody().scrollLeft+event.clientX-tipobj.offsetWidth+"px" : window.pageXOffset+e.clientX-tipobj.offsetWidth+"px"
else if (curX<leftedge)
tipobj.style.left="5px"
else
//position the horizontal position of the menu where the mouse is positioned
tipobj.style.left=curX+offsetxpoint+"px"

//same concept with the vertical position
if (bottomedge<tipobj.offsetHeight)
tipobj.style.top=ie? ietruebody().scrollTop+event.clientY-tipobj.offsetHeight-offsetypoint+"px" : window.pageYOffset+e.clientY-tipobj.offsetHeight-offsetypoint+"px"
else
tipobj.style.top=curY+offsetypoint+"px"
tipobj.style.visibility="visible"
}
}

function hideddrivetip(){
if (ns6||ie){
enabletip=false
tipobj.style.visibility="hidden"
tipobj.style.left="-1000px"
tipobj.style.backgroundColor=''
tipobj.style.width=''
}
	       }

document.onmousemove=positiontip

</script>









<div align="center">
<div id="wrapper">
<!-- ############### TITLE AREA ############### -->

<div id="header">
<h1>GENOME BIOINFORMATICS RESEARCH LAB</h1>
<small>
Grup de Recerca en Inform&agrave;tica Biom&egrave;dica <br class="hhh">
Institut Municipal d'Investigaci&oacute; M&egrave;dica + Universitat Pompeu Fabra + Centre de Regulaci&oacute; Gen&ograve;mica<br class="hh"  />

</small>

<div align="center"><div id="menudiv"><table class="gblnav"><tbody><tr class="menu"><td id="Mhelp"><a href="http://genome.imim.es/software/secisaln/help.html" title="Ask for Help">Help</a></td>
<td id="Mnews"><a href="http://genome.imim.es/main/news.php" title="News">News</a></td>
<td id="Mpeople"><a href="http://genome.imim.es/main/people.php" title="Who is who in Genome Informatics Research Laboratory">People</a></td>
<td id="Mresearch"><a href="http://genome.imim.es/research/" title="Research at our Group">Research</a></td>
<td id="Msoft"><a href="http://genome.imim.es/software/" title="Software developed by our Group">Software</a></td>

<td id="Mpapers"><a href="http://genome.imim.es/publications/" title="Publications by our Group">Publications</a></td>

<td id="Mlinks"><a href="http://genome.imim.es/main/links.php" title="Interesting Links">Links</a></td>
</tr></tbody></table>

<table class="gblnav"><tbody><tr class="menu"><td id="Mdata"><a href="http://genome.imim.es/datasets/" title="Dat
asets created/curated by our Group">Resources&nbsp;&amp;&nbsp;Datasets</a></td>
<td id="Mpreds"><a href="http://genome.imim.es/genepredictions/" title="Gene Predictions produced by our Software">Gene&nbsp;Predicti
ons</a></td>
<td id="Mseminars"><a href="http://genome.imim.es/seminars/" title="Seminars, Courses and Grou
p Sessions">Seminars&nbsp;&amp;&nbsp;Courses</a></td>
</tr></tbody></table></div></div>
<div class="alt"><hr></div>

</div></div></div>
STOPHERE

if($added_a==1){print "<DIV onMouseOver=\"ddrivetip(\'Since the submitted sequence was too short for prediction of helix 1 (less than 100nt), 30 &quot;a&quot; were added to the beginning of the sequence and the complementary 30 &quot;u&quot at the end. These are shown in lowercase in the aligned output.\',\'\',\'300\')\" OnMouseOut=\"hideddrivetip()\" id=\"warning\"><blink>IMPORTANT:</blink> Your sequence was too short for a reliable SECIS prediction. Click <a href=\"http://genome.imim.es/software/secisaln/help.html#long\" style=\"color:white; font-weight:bold;\">here</a> for more information.  </DIV>\n"}

print <<STOPHERE;
<table class="holder" border="0" cellpadding="0" cellspacing="0">
<tr class="nohover">
<td  class="nav2">
<div id="nav2">
<a href="http://www.imim.es/" title="Institut Municipal d'Investigaci&oacute;M&egrave;dica">IMIM</a>&lt;<a href="http://www.upf.edu/" title="Universitat Pompeu Fabra">UPF</a>&lt;<a href="http://www.crg.es/" title="Centre de Regulaci&oacute; Gen&ograve;mica>CRG</a>&lt;<a href="http://nemo.imim.es/grib/" title="Grup de Recerca en Inform&agrave;tica Biom&egrave;dica">GRIB</a>&lt;<a href="http://genome.imim.es/" title="Genome Bioinformatics Research Laboratory HOME PAGE: You are Welcome!!!">HOME</a>&gt;<a href="http://genome.imim.es/software/" title="Software developed in our Group">Software</a>&gt;<a href="http://genome.imim.es/software/secisaln" title="SECISaln">SECISaln</a>

</div>
</td>
<td class="explanation"><div id ="explanation">Slide your mouse over the "S" to see the structure predicted by <a href="http://genome.unl.edu/SECISearch.html">SECISearch</a> <br>Click on the "S" to see a larger image, click on the "P" for PDF format.</div></td>
<td class="explanation"><img src="http://genome.imim.es/software/secisaln/images/arrow1.gif" height="20" width="6">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>


<td class="ambiguity"><a class="ambig" href="http://genome.imim.es/software/secisaln/ambig.html">IUPAC-IUB/GCG Ambiguity Codes</a><td>
</tr>
<tr class="nohover"><td colspan="3">
<table class="legend_top2" border="0" cellpadding="0" cellspacing="0">


<tbody><tr class="nohover"><td class="help">Help: On<input CHECKED name="help" type="radio" Onclick="he=1"> Off<input name="help" type="radio" Onclick="he=0" ></td>
STOPHERE

if(@name){## typeII
    print EFILE "**********************Got to here!!!\n" if $debug;
    unless ($skip==1){
	
	$png_image = "$pngpth/$$" . ($hitnumber-1);
	      my $local_png="$pdir/$$" . ($hitnumber-1);
	      my $tmp_png=$local_png . ".png";
	      system("convert -trim /usr/local/Install/apache-data/genome.imim.es/htdocs/software/secisaln/36320.png ~/test.png 2>> ~/error");
	      my $PNG="$pdir/$$" . ($hitnumber-1);
	      $local_png=~/^.+\/(.+)/;
	      system("echo \"convert -trim $PNG.png $PNG.trimmed.png; convert $PNG.trimmed.png $PNG.trimmed.pdf  2>>$prgbin/error\" > $prgbin/command; echo \"see : $PNG.png  -trim $PNG.trimmed.png\" > $prgbin/worked; at now < $prgbin/command");
	      $png_image = "$pngpth/$$" . ($hitnumber-1) . ".trimmed.png";
	      $pdf_image = "$pngpth/$$" . ($hitnumber-1) . ".trimmed.pdf";
	open(T,"$tmpdir/$$\_energy");
	my $structure;
	while(<T>) {
	    next unless /\(/;
	    /^\s*([\.\(\)]+)/;
	    $structure=$1;
	}

	# Using the aligned sequence units as a template
	# add "-" to the corresponding structure units
	# so as to print the predicted structure at the
	# top of the web output
	my $jaja=0;
	foreach my $seq_unit(@_){## typeII
				     
	    # if this unit has >1 set of "_"
	    if($seq_unit =~ /-[^-]+-+[^-]/){print EFILE "\$seq_unit($jaja) : $seq_unit\n" if $debug;
		my @bb=split(/-+/,$seq_unit);
		my $c=$seq_unit;
		map{
		    my $l=length($_); 
		    $structure=~s/(.{$l})//;
		    my $a=$1;
		    $c=~s/$_/$a/;
		    
		}@bb;
		push @pred_str, $c;	
	    }
	    #if unit is of the format [ACTG]---[ACTG]
	    elsif($seq_unit =~ /^([^-]+)(-+)([^-]*)/){
		print EFILE "\$seq_unit($jaja) : $seq_unit\n" if $debug;
		my $a=$1; my $b=$2; my $d=$3; 
		my $l=length($1);
		my $ll=length($3);
		$structure=~s/^(.{$l})//;
		my $c=$1.$2.$b;
		$structure=~s/^(.{$ll})//;
		$c.=$1;
		push @pred_str, $c;
		print EFILE "\CCCC     ($jaja) : $c\n" if $debug;
	    }
	    else{
		print EFILE "\$seq_unit($jaja) : $seq_unit\n" if $debug;
		$seq_unit=~/^(-*?)([^-]*)(-*)$/;
		my $a=$1; my $b=$3;
		my $l=length($2);
		print EFILE "1,2,3        : $1,$2,$3\n" if $debug;
		$structure=~s/(.{$l})//;
		my $c=$a.$1.$b;
		push @pred_str, $c;
		print EFILE "a $seq_unit\nb $c\n" if $debug;
	    }
	    $jaja++;
	}
	print EFILE "@_\n@pred_str\n" if $debug;

	my @tmp_array;
	my $c=0;
	my $C=0;
	my ($opening,$closing)=0; # number of opening/closing parentheses
	my $b;
	my $kk;
	map{$kk.=$_;}@pred_str;
	my @pred_str_split = split(//,$kk);
	print EFILE "k : $k \n @pred_str_split\n" if $debug;
	
	foreach my $res (@pred_str_split) {
	    $C++;
	    if (($opening - $closing ==0) && ($C>0) ){$c+=100; }
	    if($res =~/\(/) {
		$c++;
		$opening++;
		my $lid= "l".$c;
		my $rid= "r".$c;
		$b= "<td id=\"$lid\" class=\"whi\" onmouseover=\"this.className=\'red\'; document.getElementById (\'$rid\').className = \'red\'\" onmouseout=\"this.className=\'whi\'; document.getElementById (\'$rid\').className = \'whi\'\">" . $res . "</td>"; 
		push @tmp_array, $b;
	    }
	    elsif ($res =~ /\)/){
		$closing++;
		my $lid= "l".$c;
		my $rid= "r".$c;
		$b ="<td class=\"whi\" id=\"$rid\" onmouseover=\"this.className=\'red\'; document.getElementById (\'$lid\').className = \'red\'\" onmouseout=\"this.className=\'whi\'; document.getElementById (\'$lid\').className = \'whi\'\">" . $res . "</td>";
		$c--;
		push @tmp_array, $b;
	    }
	    else{
		$b ="<td>" . $res . "</td>";
		push @tmp_array, $b;
	    }
	
	}

	
	@pred_str=@tmp_array;

	
    }
	open(A, "$pdir/template_table_top");
	my @tt=<A>;
	print "@tt\n";
	close(A);
	
	my $a;
	map{$a .= $_;}@_;

	my @tt =  &generate_table("type2_nostr.secisaln");

	%hash = %{$tt[0]};
	%gihash = %{$tt[1]};
	@sequences = @{$tt[2]};
	%col=%{&fix_colors(@sequences)};
    print "@pred_str";
	print "<td class=\"lastcell\">&nbsp;</td><td class=\"lastcell\">&nbsp;</td></tr><tr class=\"nohover\"><td class=\"one\" style=\"border-style: solid; border-width: 0px 2px 0px 2px;\" valign=\"top\"><span class=\"name\">$user_seq_name</span></td>";
	my @ss = split(//, $a);
	
	for (my $n=0; $n<scalar(@ss); $n++){
	    print "<td class=\"user_seq\" bgcolor=\"$col{$ss[$n]}{$n}\">$ss[$n]</td>\n"
	}
print <<STOPHERE;	
<td class="secisimage" style="border-width: 1px 2px 0px 1px"><a class="secisimage" href="javascript:OpenWindow('$png_image')" onmouseover="changeme('legend2','<img src=$png_image id=legend2>')" onmouseout="goback('legend2','<img src=http://genome.imim.es/software/secisaln/images/type2.png id=legend2N height=455 width=229>     ')">S</a></td><td class="secisimage" style="border-width: 1px 2px 0px 1px"><a class="secisimage" name="bob" href="javascript:OpenWindow('$pdf_image')">P</a></td></tr></tbody></table>

</td>

<td rowspan="2"  class="legend">
<div id="legend2"><img src="http://genome.imim.es/software/secisaln/images/type2.png" alt="SECIS consensus structure"  width="229" height="455" style="position: relative; left: 0em; top: 1em;" ></div>
</td></tr>

<tr class="nohover">
<td class="align" colspan="3"><div class="scrollArea">
     <table class="align2" width="" cellpadding="0" border="0" cellspacing="0" ><tbody>
     
STOPHERE

	 &html_out();
	print end_html;
	#### ffffxxx ok so far
	  
       }
	else { 
	    &abnormal("Sorry, @name (name)no SECIS elements were found!");
	    if ($pat ne "pat_Sep20") { print "You may want to try another pattern" 
	    } 
	}
    }
# unlink("$tmpdir/$$\_hits");
# unlink("$tmpdir/$$\_en");
# unlink("$tmpdir/$$\_energy");
    } 






sub pnd{

print EFILE  "################################\nargs: @_\n";

    exit();
}




sub print_rest_of_table{
    
    if($type==1){
	open(A, "$pdir/template_table_rest1");
	my @tt=<A>;
	print "@tt\n";
	close(A);
    }
    else{print EFILE "sdafdnfñalsdjnfñalsjdnvañljdsnvñsakljdnv\n";
	open(A, "$pdir/template_table_rest");
	my @tt=<A>;
	print "@tt\n";
	close(A);
    }
}


sub generate_table{
    
    my $aln_file = $_[0];
    
    my (%hash, %gihash, %col)=();
    my $counter=0;
    open(A,"./FastaToTbl $aln_file |");
    my @sequences;
    
    while(<A>) {
	my $seq;
	my $tmp;

	if(/^(.+?):.+\t(.+)/){
	    $seq =$2;
	    $tmp = $1;  
	}
	elsif(/^(.+)\t(.+)/)	{
	    $seq =$2;
	    $tmp = $1;  
	}
	else{ print EFILE "shit, cannot match in $_\n";}
	my ($species, $prot,$gi);
	if($tmp =~ /gi/){
	    $tmp =~ /^(.+?)_(.+?)_(gi\|.+?)\|/;
	    ($species, $prot,$gi) = ($1,$2,$3);
	}
	elsif($tmp =~ /^(.+?)_(.+?)_(.+)/){
	    ($species, $prot,$gi) = ($1,$2,$3);
	}
	else{
	   
	    &abnormal("Internal server problem\n");
	     
	 }
	if (defined ($hash{$prot}{$species})){
	    
	    
	}
	push @{$hash{$prot}{$species}}, $seq;
	
#    $hash{$prot}{$species}=$seq;
	push @sequences, $seq;
	
	push @{$gihash{$prot}{$species}}, $gi;
    }
    my @tt=(\%hash, \%gihash, \@sequences);
    return  (@tt);
}



sub fix_colors{
    print EFILE "#$.# Fixing colors...\n";

    my %pos;
    my $k=0;
    foreach my $sequence (@_){
	my @seq = split(//,$sequence);
	    
	for (my $n=0; $n<scalar(@seq); $n++){
	    $pos{$seq[$n]}{$n}++;
	    $k=scalar(@seq);
	    
	}
    }
    for (my $n=$k; $n>=0; $n--){
	$pos{"A"}{$n}=0 unless defined($pos{"A"}{$n});
	$pos{"U"}{$n}=0 unless defined($pos{"U"}{$n});
	$pos{"C"}{$n}=0 unless defined($pos{"C"}{$n});
	$pos{"G"}{$n}=0 unless defined($pos{"G"}{$n});
	$pos{"-"}{$n}=0 unless defined($pos{"-"}{$n});
    }
    foreach my $nt ("A","C","U","G")
    {
	for (my $n=$k; $n>=0; $n--){
	    #print "xx \$pos{$nt}{$n} $pos{$nt}{$n}/",scalar(@_),")*100 :", ($pos{$nt}{$n}/scalar(@_))*100,"\n";
	    $col{"-"}{$n}= "";
	    if(($pos{$nt}{$n}/scalar(@_))*100>=90){
		$col{$nt}{$n}= "#6464ff";
	    }
	    elsif(($pos{$nt}{$n}/scalar(@_))*100>=60){
		$col{$nt}{$n}= "#9999ff";
	    }
	    elsif(($pos{$nt}{$n}/scalar(@_))*100>=30){
		$col{$nt}{$n}= "#ccccff";
	    }
	    else{
		$col{$nt}{$n}="";
	    }
	}
    }
for (my $n=0; $n<scalar(@ss); $n++){
	    print EFILE "<td bgcolor=\"$col{$ss[$n]}{$n}\">$ss[$n]</td>\n"
	}
    return (\%col);
    %global_pos = %pos;
    
}






sub html_out{

    
    my $logo="logo1.png";
    my $lheight="376";
    my $lwidth="680px";
    my $lcolspan=122;
    if($type==2){
	$logo="logo2.png";
	$lheight="396"; 
	$lwidth="800";
	$lcolspan=137;
	$spacer="left_spacer2";
	#$legend="legend2";
	$legend_image="type2.png";
	$lgnd_height=455;
	$lgnd_width=229;
    }
    if($sort_by_prot){
	my @proteins=sort(keys(%hash));
	    my $not_very_useful_counter=0;	
	my $last_prot="ha";
	foreach my $prot (@proteins){

	    foreach my $species (keys(%{$hash{$prot}})){

		for (my $n=0; $n<scalar(@{$hash{$prot}{$species}}); $n++)
		{
		    my $toto=$gihash{$prot}{$species}[$n];
		    $toto=~s/\|/_/g;
		    
		    my $pngname= $species . "_" .$prot."_". $toto.".png";
#		    my $pdfname= $species . "_" .$prot."_". $toto.".pdf";
		    my $image_page = $species . "_" .$prot."_". $toto.".html";		    
		    my $pdfname = $species . "_" .$prot."_". $toto.".pdf";
		    if($pngname =~/15-k/){$pngname =~s/15-k/15k/}
		    if($pdfname =~/15-k/){$pdfname =~s/15-k/15k/}

		    $counter++;
		    my $class ="one";
		    $class="border" if $prot ne $last_prot;
		    $class="one" if $not_very_useful_counter==0; 
		    $not_very_useful_counter=1;
		    $last_prot=$prot;
		    
#		    my $header ="<tr><td class=\"$spacer\">&nbsp</td><td id=\"$counter\" onMouseOver=\"changeme(\'$counter\',\'${$gihash{$prot}{$species}}[$n]\')\" onMouseOut=\"goback(\'$counter\',\'$species <span class=name>$prot</span> \')\" class=\"$class\" nowrap>$species <span class=\"name\"> $prot</span> </td>";
#my $header ="<tr><td id=\"$counter\" class=\"$class\" ><a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&cmd=search&term=${$gihash{$prot}{$species}}[$n]\" class=\"seq_link\"> $species <span class=\"name\"> $prot</span> </td>";

		    my $header ="<tr><td id=\"$counter\" onMouseOver=\"ddrivetip(\'Click through for more information on this sequence.\',\'\',\'100\')\" OnMouseOut=\"hideddrivetip()\"  class=\"$class\" ><a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&cmd=search&term=${$gihash{$prot}{$species}}[$n]\" class=\"seq_link\"> $species <span class=\"name\"> $prot</span> </td>";




		    if(${$gihash{$prot}{$species}}[$n]=~/[^gi]\|/){
			$header ="\n<tr><td id=\"$counter\" onMouseOver=\"changeme(\'$counter\',\'${$gihash{$prot}{$species}}[$n]\')\; ddrivetip(\'This sequence is from the EGO database. For more information on this sequence please copy this accession number and use the EGO search page (we provide a link to EGO on the previous page).\')\" onMouseOut=\"goback(\'$counter\',\'$species <span class=name> $prot</span>\')\; hideddrivetip()\" class=\"$class\" >$species <span class=\"name\"> $prot</span> </td>\n";
		

		    }
		    
		    print "$header";
		    
		    my @sequence = split(//,$hash{$prot}{$species}[$n] );
		    
		    for (my $i=0; $i<scalar(@sequence); $i++){
			print "<td bgcolor=\"$col{$sequence[$i]}{$i}\">$sequence[$i]</td>";
		    }
# 		open(A, ">tmp.html");
# 		print A "<html><body><a href=\"$pngpth/secis_images/$pdfname\">PDF format</a><br><img src=\"$pngpth/secis_images/backup_images/$pngname\"></body></html>";
# 		close(A);

		print "<td class=\"secisimage\"><a class=\"secisimage\" href=\"javascript:OpenWindow(\'$pngpth/secis_images/trimmed/$pngname\')\" onmouseover=\"changeme(\'$legend\',\'<img src=$pngpth/secis_images/$pngname id=legend1N>\')\" onmouseout=\"goback(\'$legend\',\'<img src=http://genome.imim.es/software/secisaln/images/$legend_image id=legend1N height=$lgnd_height width=$lgnd_width >     \')\">S</a></td><td class=\"secisimage\" style=\"border-width: 1px 2px 0px 1px\"><a class=\"secisimage\" name=\"bob\" href=\"$pngpth/secis_images/$pdfname\")\">P</a></td></tr>";
#		    print "<td class=\"secisimage\"><a class=\"secisimage\" href=\"javascript:OpenWindow(\'<html><body><a href=\"$pngpth/secis_images/$pdfname\">PDF format</a><br><img src=\"$pngpth/secis_images/backup_images/$pngname\"></body></html>\')\" onmouseover=\"changeme(\'$legend\',\'<img src=$pngpth/secis_images/$pngname id=legend1N>\')\" onmouseout=\"goback(\'$legend\',\'<img src=http://genome.imim.es/software/secisaln/images/$legend_image id=legend1N height=$lgnd_height width=$lgnd_width >     \')\">S</a></td></tr>";
#		print EFILE "AAA :  height=$lgnd_height width=$lgnd_width $type\n ";

		}
	    }
	    

    }

	print EFILE ".\n";
    }

## If we are sorting by species
    else{

	my %species_hash;
	my @proteins=sort(keys(%hash));
	
	my $last_species="ha";
	foreach my $prot (@proteins){
	    
	    foreach my $species (keys(%{$hash{$prot}})){
		$species_hash{$species}{$prot}= [@{$hash{$prot}{$species}}];
	    }
	}


	my @SPECIES= ("Human", "Chimp","M.fascicularis", "M.mulatta","Orangutan","Horse","Cow","Pig", "Sheep","C.lupus", "Dog", "Mouse", "Rat","M.auratus", "O.anatinus","Chicken","T.guttata","O.latipes", "O.mykiss",  "Salmon",   "T.nigroviridis","Zebrafish","Trout","X.laevis", "X.tropicalis","C.intestinalis","B.microplus","A.aegypti","A.gambiae","D.ananensis", "D.erecta", "D.grimshawi", "D.melanogaster", "D.mojavensis", "D.persimilis", "D.pseudoobscura", "D.sechelia", "D.simulans", "D.virilis", "D.yakuba", "Fly", "O.fasciatus","C.briggsae","C.elegans",  "E.granulosus","S.japonicum","N.vectensis","A.elegantissima","T.thermophila","C.reinhardtii");

#=keys(%species_hash);
	
	foreach my $species (@SPECIES){
	    foreach my $prot(sort(keys(%{$species_hash{$species}}))){
		for (my $n=0; $n<scalar(@{$species_hash{$species}{$prot}}); $n++)
		{
		    my $toto=$gihash{$prot}{$species}[$n];
		    $toto=~s/\|/_/g;
		    
		    my $pngname= $species . "_" .$prot."_". $toto.".png";
		    if($pngname =~/15-k/){$pngname =~s/15-k/15k/}
		    
		    $counter++;
		    my $class ="one";
		    $class="border" if $species ne $last_species;
		    $last_species=$species;
		    my $spacer;
		    if($type==2){$spacer="left_spacer2"}
		    else{$spacer="left_spacer"}
#		    my $header ="<tr><td class=\"$spacer\">&nbsp</td><td id=\"$counter\" onMouseOver=\"changeme(\'$counter\',\'${$gihash{$prot}{$species}}[$n]\')\" onMouseOut=\"goback(\'$counter\',\'$species <span class=name>$prot</span> \')\" class=\"$class\" nowrap>$species <span class=\"name\"> $prot</span> </td>";

		    my $header ="<tr><td id=\"$counter\" onMouseOver=\"ddrivetip(\'Click through for more information on this sequence.\',\'\',\'120\'))\" OnMouseOut=\"hideddrivetip()\"  class=\"$class\" ><a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=nuccore&cmd=search&term=${$gihash{$prot}{$species}}[$n]\" class=\"seq_link\"> $species <span class=\"name\"> $prot</span> </td>";
		    ## Special case for seqs with no gis (EGO)
		    if(${$gihash{$prot}{$species}}[$n]=~/[^gi]\|/){
		    $header ="\n<tr><td id=\"$counter\" onMouseOver=\"changeme(\'$counter\',\'${$gihash{$prot}{$species}}[$n]\')\; ddrivetip(\'This sequence is from the EGO database. For more information on this sequence please copy this accession number and use the EGO search page (we provide a link to EGO on the previous page).\')\" onMouseOut=\"goback(\'$counter\',\'$species <span class=name> $prot</span>\')\; hideddrivetip()\" class=\"$class\" >$species <span class=\"name\"> $prot</span> </td>\n";

		    }
		    print "$header\n";
		    my @sequence = split(//,$hash{$prot}{$species}[$n] );
		    
		    for (my $i=0; $i<scalar(@sequence); $i++){
			print "<td bgcolor=\"$col{$sequence[$i]}{$i}\">$sequence[$i]</td>\n";
		    }
		    
		    print "<td class=\"secisimage\"><a class=\"secisimage\" href=\"javascript:OpenWindow(\'$pngpth/secis_images/backup_images/$pngname\')\" onmouseover=\"changeme(\'$legend\',\'<img src=$pngpth/secis_images/$pngname id=legend2N>\')\" onmouseout=\"goback(\'$legend\',\'<img src=http://genome.imim.es/software/secisaln/images/$legend_image id=legend2N height=$lgnd_height width=$lgnd_width >     \')\">S</a></td></tr>\n";
		}
		
	    }
	    

	}

	print EFILE ".\n";

    }
print "\n<tr class=\"nohover\"><td colspan=\"$lcolspan\" class=\"spacer\">&nbsp;</td></tr>\n<tr class=\"nohover\"><td class=\"nologo\" colspan=\"$lcolspan\"><img width=\"$lwidth\" height=\"$lheight\" src=\"http://genome.imim.es/software/secisaln/images/$logo\"></td></tr><tr class=\"nohover\"><td colspan=\"$lcolspan\" class=\"spacer$type\">&nbsp;</td></tr>\n</div></td></tr></tbody></table></td></tr></tbody> </table>\n";



#    print "<img src=\"http://genome.imim.es/software/secisaln/images/logo2.png\" id=\"logo2\">";
}
	






sub print_rest_of_table_no_gis{
    if($type==1){
	open(A, "$pdir/template_table_rest_no_gis1");
	my @tt=<A>;
	print "@tt\n";
	close(A);
    }
    else{
	open(A, "$pdir/template_table_rest_no_gis");
	my @tt=<A>;
	print "@tt\n";
	close(A);
    }
}

sub color_sequence
{
    my @seq = @_;
    print EFILE "CCCCCcoloring Seq : @_ seq ----- $seq[27]\n";

    
    my @seq1;
    my $i=0;
    foreach my $res (@seq){
	$seq1[$i]="<td>$res</td>";
	$i++;
    }
    @seq = @seq1;

    if($type==2){
	if($seq[22] =~ ">G<"){$seq[22] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[23] =~ ">U<"){$seq[23] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[27] =~ ">U<"){$seq[27] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[28] =~ ">U<"){$seq[28] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[29] =~ ">A<"){$seq[29] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	## AUGA quartet
	$seq[30] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[31] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[32] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[33] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	## end quartet
	
	if($seq[35] =~ ">G<"){$seq[35] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[36] =~ ">G<"){$seq[36] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[37] =~ ">C<"){$seq[37] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[43] =~ ">C<"){$seq[43] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	
	## apical loop
	$seq[51] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[52] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	
	
	if($seq[57] =~ ">C<"){$seq[57] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[58] =~ ">C<"){$seq[58] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[74] =~ ">G<"){$seq[74] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[75] =~ ">G<"){$seq[75] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[84] =~ ">G<"){$seq[84] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[85] =~ ">G<"){$seq[85] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[90] =~ ">U<"){$seq[90] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[91] =~ ">G<"){$seq[91] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[93] =~ ">C<"){$seq[93] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[94] =~ ">U<"){$seq[94] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	
	## complement of quartet (GA)
	$seq[95] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[96] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	
	if($seq[97] =~ ">U<"){$seq[97] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[98] =~ ">G<"){$seq[98] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[100] =~">U<"){$seq[100] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[114] =~">G<"){$seq[114] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
    }
    else{

	if($seq[6] =~ ">G<"){$seq[6] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[11] =~ ">U<"){$seq[11] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[17] =~ ">U<"){$seq[17] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[19] =~ ">G<"){$seq[19] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[23] =~ ">U<"){$seq[23] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[24] =~ ">U<"){$seq[24] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[27] =~ ">U<"){$seq[27] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[28] =~ ">U<"){$seq[28] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[29] =~ ">A<"){$seq[29] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	## AUGA quartet
	$seq[30] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[31] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[32] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[33] =~ s/<td/<td bgcolor=\"\#6464ff\"/;

	if($seq[34] =~ ">U<"){$seq[34] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[35] =~ ">G<"){$seq[35] =~ s/<td/<td bgcolor=\"\#6464ff\"/}
	if($seq[36] =~ ">G<"){$seq[36] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[37] =~ ">U<"){$seq[37] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[39] =~ ">U<"){$seq[39] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[41] =~ ">U<"){$seq[41] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[42] =~ ">G<"){$seq[42] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[43] =~ ">C<"){$seq[43] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	## apical loop
	$seq[50] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[51] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	

	if($seq[52] =~ ">A<"){$seq[52] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[53] =~ ">C<"){$seq[53] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[54] =~ ">U<"){$seq[54] =~ s/<td/<td bgcolor=\"\#ccccff\"/}

	if($seq[72] =~ ">G<"){$seq[72] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[73] =~ ">G<"){$seq[73] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[74] =~ ">G<"){$seq[74] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[80] =~ ">C<"){$seq[80] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[81] =~ ">C<"){$seq[81] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[82] =~ ">U<"){$seq[82] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	
	## complement of quartet (GA)
	$seq[83] =~ s/<td/<td bgcolor=\"\#6464ff\"/;
	$seq[84] =~ s/<td/<td bgcolor=\"\#6464ff\"/;


	if($seq[85] =~ ">U<"){$seq[85] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[86] =~ ">G<"){$seq[86] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[87] =~ ">U<"){$seq[87] =~ s/<td/<td bgcolor=\"\#ccccff\"/}



	if($seq[58] =~ ">C<"){$seq[58] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[74] =~ ">G<"){$seq[74] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[75] =~ ">G<"){$seq[75] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[90] =~ ">U<"){$seq[90] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[92] =~ ">G<"){$seq[92] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[93] =~ ">A<"){$seq[93] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[94] =~ ">A<"){$seq[94] =~ s/<td/<td bgcolor=\"\#9999ff\"/}
	if($seq[97] =~ ">U<"){$seq[97] =~ s/<td/<td bgcolor=\"\#ccccff\"/}
	if($seq[99] =~ ">C<"){$seq[99] =~ s/<td/<td bgcolor=\"\#ccccff\"/}

    }
	return @seq;

}










sub scan_for_matches{
###--running scan_for_matches with choosen pattern
    open (SEQ, ">$tmpdir/$$\_seq");
    print SEQ $seq;
    close (SEQ);
    print EFILE "#<$$># RUNNING SCAN1GB...\n";
    system("$scan1GB $c $prgbin/$pat < $tmpdir/$$\_seq >$tmpdir/$$\_hits 2>>$logfile");
    print EFILE "$scan1GB $c $prgbin/$pat < $tmpdir/$$\_seq >$tmpdir/$$\_hits 2>>$logfile\n";
    print EFILE "#<$$># ...SCAN1GB DONE\n";
###--running show_hits
    print EFILE "$tmpdir/$$\_hits\n";
    
}

sub SECISearch{
    open (HITS, "$tmpdir/$$\_hits");
    
    my $CC= <HITS>;

    ## If we had no hits
    unless($CC){
	print EFILE "#<$$># No hits found\n";
	&abnormal("Sorry, no hits found. You may want to try the other pattern.\n");
    }

    $hitnumber=0;
    while (<HITS>) {
	print EFILE;
	my @param=();
	($g1, $name, $start, $stop) = split />|\:\[|\,|\]/;
	if($name=~/^(.+?)_/){
	    $name=$1; 
	}
	
	tr/Tt/Uu/;
	@a=split;
	print EFILE "aa : @a\n";
	s/\w/./g for @b=@a;
	$b[3]='x(((('; $b[5]='xx';
	my $kkk = length ($a[8]);
	
### Get unit numbers. May change depending on pattern
	$units = @a;
	$lastunit = $units -1;
	if ($units == 12)
	{
	    $unit8 = 8;
	    $unit7 = 7;
	}
	else
	{
	    for(my $n=0; $n<scalar(@a); $n++){
		if ($a[$n] =~ /^.GA.$/){
		    $unit8=$n;
		    $unit7=$n-1;
		}
		elsif($a[$n] =~ /^..GA.$/){
		    $unit8=$n;
		    $unit7=$n-1;
		}
	    }
	    
	}
###
###
	
	if (length ($a[$unit8]) == 4)  { $b[$unit8]='))))' }
	elsif (length ($a[$unit8]) == 5)  { $b[$unit8]='.))))'}
# if (length ($a[8]) == 4)  { $b[8]='))))' }
# elsif (length ($a[8]) == 5)  { $b[8]='.))))' }
	else { die(); print EFILE "acacacacac\n"; &abnormal("a8 length problem : length \$a[8]($a[$unit8])=$kkk : @a\n "); exit(); print EFILE "bbbbbbb\n";}

#if ( ($upstemen = &energy(3..8)) <$e1 && ($fullstructen = &energy(0..11))<$e2 && &structure_OK()) {
	if (($fullstructen = &energy(0..scalar(@a)))<$e2 && &structure_OK()) {
	    
	    push @param, $currentseq;
	    push @param, $currentstruct;
	    push @param, 15;
	    push @param, "$pdir/".$$.$hitnumber.".png";
	    push @param, 170;
	    push @param, 250;
	    push @param, length (join "", @a[0..2])+1;
	    push @param, length (join "", @a[0..2])+4;
	    push @param, (length (join "", @a[0..2]))..((length (join "", @a[0..2]))+4);
	    push @param, (length (join "", @a[0..4]))..((length (join "", @a[0..4]))+1);

	    ## Trick to make right residues appear in bold. Whatever the pattern,
	    ## the 3' quartet will be at position: length_of_array - 4
	    my $core_quartet=@a;
	    $core_quartet-=4;
	    &debug("\$a[\$core_quartet]=$a[$core_quartet]\n"); 
	    if (length ($a[$core_quartet]) == 4) {
		&debug("length ($a[8]) : " . length ($a[$core_quartet]) . "\n");
		push @param, (length (join "", @a[0..$core_quartet-1]))..((length (join "", @a[0..$core_quartet-1]))+3);  
	    }
	    else {
		&debug("length ($a[$core_quartet]) : " . length ($a[8]) . "\n");
		push @param, (length (join "", @a[0..$core_quartet-1])+1)..((length (join "", @a[0..$core_quartet-1]))+4) }
	    print EFILE "PARAM . @param\n";
	    &imager2(@param);
	    
	    push @texta, [ @a ];
	    push @name, $name;

	    push @start, $start;
	    push @stop, $stop;
	    push @png, "\"$pngpth/".$$.$hitnumber.".png\"";
	    
	    push @comment, $start>$stop ? "(SECIS on complementary strand)" : "";
#push @upstemen, $upstemen;
	    push @fullstructen, $fullstructen;
	    
	    $hitnumber++;
	    
	}
	else{$fullstructen = &energy(0..scalar(@a)); &abnormal("There are no predicted SECIS elements above the Free Energy threshold given ($e2). Please try a higher value"); print EFILE "full : $fullstructen\n"; }

	&align();

    }
}
