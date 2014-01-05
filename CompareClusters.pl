#!/usr/bin/perl


$clusterA=$ARGV[0];
$clusterB=$ARGV[1];


open(clusterA,"$clusterA") || die "Cannot open $clusterA !!\n";

$nca=1;
while ($line=<clusterA>)
{
    chomp($line);
    @line= $line =~ /([^\s\t\s ]+)/g;
#    $nameA{$nca}=shift @line;
    @{$clusterA{$nca}}=@line;
    $nca++;
}
close clusterA;

open(clusterB,"$clusterB") || die "Cannot open $clusterB !!\n";

$ncb=1;
while ($line=<clusterB>)
{
    chomp($line);
    @line= $line =~ /([^\s\t\s ]+)/g;

#    $nameB{$ncb}=shift @line;
    @{$clusterB{$ncb}}=@line;
    $ncb++;
}
close clusterB;

$nca--;
$ncb--;


print "Number of clusters in $clusterA : $nca\n";
print "Number of clusters in $clusterB : $ncb\n";

########## DEBUT DE LA COMPARAISON ###############

for ($i=1;$i<=$nca;$i++)
{
    $dmin=2;
    $max_p = 0;
    @closecluster = ();
    @overlap = ();
    $n=0;
    $maxover=0;
    $closest=-1;
    $mind=2;

    for ($j=1;$j<=$ncb;$j++)
    {

	@result=distance(\@{$clusterA{$i}},\@{$clusterB{$j}});

	$d=shift @result ;
	$nelA = shift @result;
	$nelB = shift @result;
	$over= shift @result;
	$percent_over = $over/$nelA;

#	if ($percent_over > $max_p)

	$n+=$over;
#	if ($d < $mind)
	if ($over > $maxover)
	{
	    $maxover=$over;
	    $closest=$j;
	    $mind=$d;
	}
    }
    $maxoverlap = $maxover /$nelA;
    my $num1=scalar @{$clusterA{$i}};
    my $num2=scalar @{$clusterB{$closest}};
    print "----------------------------------------------\n";
    print "$clusterA, Class $i ($num1 proteins): @{$clusterA{$i}}\n$clusterB, Class $closest ($num2 proteins): @{$clusterB{$closest}}\nOverlap: $maxoverlap\nDistance: $mind\nShared:  $maxover of $num1\n";
#    printf("%s\t%.2f\t%.2f\n",$nameA{$i},$overlap,$maxoverlap);


}


sub distance
{
    my $refA=shift(@_);
    my $refB=shift(@_);

    my $dist;
    my $ncommon=0;
    my $nspecA=0;

    my @clusterA=@{$refA};
    my @clusterB=@{$refB};

#    print "cluster A : @clusterA\n";

    my $nelA=$#clusterA+1;
    my $nelB=$#clusterB+1;

    foreach my $elA (@clusterA)
    {
	if (grep(/^$elA$/,@clusterB))
	{
	    $ncommon++;
	}
	else
	{
	    $nspecA++;
	}
    }
    my $nspecB=$nelB-$ncommon;
    $dist = ($nspecA+$nspecB)/($nelA+$nelB-$ncommon);
    return ($dist,$nelA,$nelB,$ncommon,$nspecA,$nspecB);
}

