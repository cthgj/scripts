#! /usr/bin/perl -w

my $se = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC";
my $pa = "AAAAAAAC";
my @seq = split(//,$se);
my @pat = split(//,$pa);

my @netab = (-1,1,2,3,4,5,6,7);
my $iter=0;
my $i=0;

$j = 0;
$match = "";
while($j<scalar(@seq))
{
    
    $iter++;
    if($seq[$j] eq $pat[$i]){ $match=$match . $seq[$j];  $j++; $i++;}
    else{$j++; $i=$netab[$i]; }
    if($match eq $pa){ print "$match (positions " . ($j-$i+1)  . " - "  .  ($j) . ", $iter iterations)\n"; $match=""; $i=0;}
}
# print "$match (positions " . ($j-$i+1)  . " - "  .  ($j) . ", $iter iterations)\n";


print "====================\n";
$j = 0;
$iter=0;
$match = "";
my $k=0;
while($j<scalar(@seq))
{
    $iter++;
     $i=0 if $i == scalar(@pat);
    if($seq[$j+$i] eq $pat[$i]){ $match=$match . $seq[$j+$i];   $i++;}
    else{$i=0; $match="";  $j++;}
    if($match eq $pa){ print "$match (positions " . ($j+1)  . " - "  .  ($j+$i) . ", $iter iterations)\n";$match=""; $i=0; $j++;} 
}



exit();   

