#!/usr/bin/perl -w
use strict;
use Getopt::Std;

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
getopts('dt:e:s:',\%opts);
my $debug = $opts{d} || undef;
my $en  = $opts{e} || die("need an energy file\n");
my $hits = $opts{s} || die("need a hits file\n");
my $type = $opts{t} || 2;	
my $apical_repetition;
my $skip=0;
open(EN,"FastaToTbl $en |");
while(<EN>)
{
    /^(.*?\]).*\t([AUGCN]*)([\.\(\)]*)/;
    $k{$1} = $3;
    $sequence{$1} = $2;
   
}
open(STR, "FastaToTbl $hits |");
while(<STR>)
{
   /^(.*)\t(.*)/;

    my @structure = ();
    my @units=();
   next unless defined($k{$1});
    $N=$1;
 
    $a = $2;
    $a  =~ s/T/U/g;
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

#print ">struct\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.(((((((((((((((((((.......(((((((.......)))))))..))))))))))))))))))).xxxxxxxxxxxxxxxxxxxxxxx\n";
 foreach my $name (@k)
{   
    $skip=0;
    $kk++; # just to give unique names to the sequences
    my $aln = ""; 
    $main ="";
    $main_str = "";
#next unless $name =~ /TC47246/;
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
	{next}
	else
	{
	    push @ff, ${$seq{$name}}[$n];
	    push @str, ${$str{$name}}[$n];
	}

    }
    my @Seq   = @ff;
    my @struc = @str;
    &debug("1 $name\n  @{$seq{$name}}\n  @Seq\n  @struc\n");   
    # check for non-standard  secis, CC in apical loop
    # or even NN
    $apical_repetition = $Seq[5];
    
    ## Sometimes a SECIS with too short a stem is predicted, skip
    $struc[4] =~ /^(.{8})/;
    my @kk = split(//, $1);
    my @ll = grep(/\(/, @kk);
    next if scalar(@ll) < 6;




# Check for type I/II. If after the apical As there are "("
# then it is type 2, else type1
    if($struc[6] =~ /^[^\(]+\)/)
    {
	if($type == 1){
	    &align_type1($name,\@Seq,\@struc);
	}
	else{next}
    }
    else
    {
	if($type == 2){
	    &align_type2($name,\@Seq,\@struc);
	}
	else{next}
	
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
## Here, we separate hgelix3+, apical loop, and helix3-
## they are still [6] but now with space between them 
    if($struc[6] =~ /^(\(+.*?)(\.{3,7})(\)+)/){
#    if($struc[6] =~ /^(\(+)(\.+)(\)+)/){
	&debug("1,2,3 :: $1,$2,$3");
	$apstem1 = length($1);
	my $b = length($2);
	$apstem2 = length($3);
	$struc[6] =~ s/(.{$apstem1})($2)(.{$apstem1})(.*)/$1 $2 $3 $4/;
	my $f_str = $4;
	$Seq[6] =~ /^(.{$apstem1})(.{$b})(.{$apstem1})(.*)/;
	# c,d,e,f = apical stem+,apical loop, apicalstem-, unpaired apical plus stem-
	my ($c,$d,$e,$f) = ($1,$2,$3,$4);

	&writeX(length($c),7,\$c);
	&writeX(length($d),7,\$d);
	&writeX(length($e),7,\$e,1);
#	&writeX(length($f),13,\$f,1);
	# $f is the 2nd helix- AND if present the rest of the
 	# apical loop, so deal with that
	  &debug("f_str : $f_str : $struc[6]\n");
	&debug("c,d,e,f :: $c,$d,$e,$f");
	if($f_str =~  /^(\.+)/)
	{
	    my $kk = length($1);
	    $f =~  /^(.{$kk})/;
	    my $a = $1;
	    &writeX(length($a),4,\$a);
	    $f =~ s/^.{$kk}/$a/;
	    $f =~ s/$a(.*)/$a/;
	    $Seq[7] = $1 . $Seq[7];
	    &debug("f/f_str : $f :: $f_str :: $a\n");
	}
	else
	{
	    $Seq[7] = $f . $Seq[7];
	    $f = "----";
	}
	$Seq[6] =~ s/^(.{$apstem1})(.{$b})(.{$apstem1})(.*)/$c$d$e$f/;

    }

    else{die("xxxxx @Seq\nxxxxx @struc\n$Seq[6]\n$struc[6]\n$name\n"); }
   
    &debug( "4 @{$seq{$name}}\n  @Seq\n  @struc");
    ## need both sides of the apical stem to have the same length
   #  if($apstem1 > $apstem2)
#     {
# 	my $dif = $apstem1 - $apstem2;
# 	my @a = split(/\s+/, $Seq[6]);
# 	$struc[7] =~ s/^.{$dif}//;
# 	$Seq[7] =~ s/^(.{$dif})//;
# 	&debug("$a[2] : $1\n@a");
# 	$a[2] = $a[2] . $1;
# 	$Seq[6] = $a[0] . " " . $a[1] . " "  . $a[2]; 

#     }
#     if($struc[7] =~ /^(\.+)/){
# 	&debug("7 was : -$Seq[7] $struc[7]($1)");
# 	my $h = length($1); 
# 	$Seq[7] =~ s/^.{$h}//;
# 	&debug( "7 is  : -$Seq[7]");
# 	$Seq[6] = $Seq[6] . " " . $1;
# 	print STDERR "$name : $struc[7]\n$name : $Seq[7]\n";

#     }
  
 #   $a = $Seq[4] . " " . $Seq[5] . " " . $Seq[6] . " " . $Seq[7];
 #   my $b =  $struc[4] . " " . $struc[5] . " " . $struc[6] . $struc[7];

# now, align the damn things
    if(length($Seq[2]) < 13)
    {
	&writeX(length($Seq[2]),13,\$Seq[0],1);
    }
    if(length($Seq[4]) < 12)
    {
	&writeX(length($Seq[4]),12,\$Seq[4]);
    }
   
    if(length($Seq[5]) < 8)
    {
	if($Seq[5] =~ /^(.+?)($apical_repetition)(.*)/)
	{
	   # $Seq[5] =~ /^(.+?)($apical_repetition)(.*)/;
	    my $a = $1;
	    my $b = $3;
	    my $l;
#	    length($b)>length($a) ? $l = length($b) : $l = length($a);
	    &writeX(length($a),4,\$a,1);
	    &writeX(length($b),3,\$b);
	    $a = $a . "$apical_repetition";
	    $Seq[5] =~ s/^(.+?)($apical_repetition).*/$a$b/;
	}
	else
	{
	   $Seq[5] = "----" . $Seq[5]; 
	   $Seq[5] =~ /^.*$apical_repetition(.*)/;
	   my $a = $1;
	   &writeX(length($a),3,\$a);
	   $Seq[5] =~ s/$apical_repetition($1)/$apical_repetition$a/;

	}

    }
    # Some seqs, eg SPS2 gi|26097414 are predicted with a huge
    # apical loop and screw everything up so skip. If i do this with 
    # next, i get a warning message so i use this variable
    elsif(length($Seq[5]) > 8)
    {
	$skip = 1; 
    }
    else{$Seq[5] = "-" . $Seq[5]; $Seq[6] =~ s/^-//;}
    if(length($Seq[7]) < 17)
    {
	&writeX(length($Seq[7]),17,\$Seq[7],1);
    }
 if(length($Seq[8]) < 4)
    {
#	&writeX(length($Seq[8]),4,\$Seq[8]);
    }
 if(length($Seq[9]) < 11)
    {
#	&writeX(length($Seq[9]),11,\$Seq[9]);
    }
#  if(length($Seq[10]) < 12)
#     {
# 	&writeX(length($Seq[8]),12,\$Seq[8]);
#     }

    unless ($skip==1){
	print ">$name \n";
	map{print "$_ "}@Seq;
	print "\n";
    }
}
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
    if(length($Seq[2]) < 13)
    {
	&writeX(length($Seq[2]),13,\$Seq[0],1);
    }
    if(length($Seq[4]) < 12)
    {
	&writeX(length($Seq[4]),12,\$Seq[4]);
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
	    my $a = $1;
	    my $b = $3;
	    &writeX(length($a),3,\$a,1);
	    &writeX(length($b),16,\$b);
	    $a = $a . "$apical_repetition";
	    $Seq[5] =~ s/^(.+?)($apical_repetition).*/$a$b/;
	}
	else
	{
	   $Seq[5] = "---" . $Seq[5]; 
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

    print ">$name \n";
    map{print "$_ " }@Seq;
	print "\n";
}
sub writeX
{
    my $l = $_[1];
    my $x = $_[0];
    my $res = ${$_[2]};
    my $difference = $l-$x;
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
sub debug
{
    if ($debug)
    {
	print STDERR "@_\n";
    }
}
### patscqan pattern (for reference)
    #r1={au,ua,gc,cg,gu,ug} NNNNNNNNNN p1=7...7 3...13 ATGAN p2=10...13 AA (4...12 | 0...3 p3=3...6 3...6 r1~p3 0...3) (r1~p2[2,1,1] NGAN | r1~p2[2,1,0] NNGAN) 3...10 r1~p1[1,1,1] NNNNNNNNNN  ENSRNOG00000013548    ENSPTRG00000023837


# structure
#>struct
#xxxxxxxxxx-xxxxxxx-xxxxxxxxxxxxx---.((((-(((((((((((((--.......--(((((((---.......--)))))))-..---)))))))))))))-----))))-.xxxxxxxxxxxxxxxxxxxxxxx 
