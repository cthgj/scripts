#!/usr/bin/perl

$db = "simct";
$user = "simct";
$pw = "simct";

#Export relevant databases

unless ($#ARGV == 1) {
    print "Usage: $0 <onto (molecular_function, biological_process, cellular_component)> <taxid>\n";
    exit;
}
$onto = $ARGV[0];
$taxid = $ARGV[1];
$outfile="out-".$onto."-".$taxid.".txt";
$ts= localtime(time);
$ts =~ s/ /-/g;

#Export annotation tables
$out1="allannot.".$ts.".txt";
print STDERR "EXPORTING RELEVANT TABLES...\n";
$query = "select distinct T1.entrezid,T2.id,T2.ancestor_id from goannotations as T1, ancestors as T2 where T1.termid=T2.id and T1.namespace=\'$onto\' and T1.taxid=\'$taxid\' " ;

print STDERR "mysql --host=10.1.1.54 -u $user --password=$pw --execute=\"$query\" $db > /tmp/$out1\n";
system("mysql --host=10.1.1.54 -u $user --password=$pw --execute=\"$query\" $db > /tmp/$out1");
print STDERR "\n\n";
$out2="ancestors.".$ts.".txt";
$query="select A.id,A.ancestor_id from ancestors as A, ontologies as O where A.oid=O.oid and O.root_name=\'$onto\' ;";
print STDERR "mysql --host=10.1.1.54 -u $user --password=$pw --execute=\"$query\" $db > /tmp/$out2\n\n";
system("mysql --host=10.1.1.54 -u $user --password=$pw --execute=\"$query\" $db > /tmp/$out2");

#========================================================
# first read list of ancestor-descendant combinations
#========================================================

print STDERR "READING ANCESTOR TABLE...\n";
open(IN,"/tmp/$out2") || die "Could not read /tmp/$out2 file !!\n";

while ($line=<IN>) {
    chomp $line;
    @line = split(/\t/,$line);
#    @a = sort ($line[0],$line[1]);
    $a = join("-",@line); # $a has the format child-parent
    $related{$a}=1 if ($line[0] ne $line[1]); # remove pairs of the form child-child
}
close IN;

#========================================================
# find how many genes have > 1 DIRECT annotations
#========================================================
print STDERR "COUNTING GENES WITH >1 DIRECT ANNOTATIONS ...\n";

# column 1 : entrezid
# column 2 : direct annotation
`cut -f1,2 /tmp/$out1 |sort -u > xxx`;
open(IN,"xxx") || die "Could not read xxx file !!\n";

$first=1;
$prec = '';
while ($line=<IN>) {
      chomp $line;
    @line = split (/\t/,$line);
    if ($line[0] eq $prec) {
	$terms{$line[1]}=1;
    }
    elsif (!$first) {

# here, we remove direct annotations that might be ancestors of other direct annotations
# (which in principle should not happen...)

	@terms = keys %terms;
	@goodterms=();
	foreach $a1 (@terms) {
	    $ok=1;
	    foreach $a2 (@terms) {
		if (exists $related{$a2."-".$a1}) {$ok=0;}
	    }
	    if ($ok) {
		push(@goodterms,$a1);
	    }
	}
	if ($#goodterms >0) {
	    $goodgene{$prec} =1;
	}
	%terms = ();
    }
    $first=0;
    $prec = $line[0];
}
@terms = keys %terms;
@goodterms=();
foreach $a1 (@terms) {
    $ok=1;
    foreach $a2 (@terms) {
	if (exists $related{$a2."-".$a1}) {$ok=0;}
    }
    if ($ok) {
	push(@goodterms,$a1);
    }
}
if ($#goodterms >0) {
    $goodgene{$prec} =1;
}
close IN;
$Ngg = scalar(keys %goodgene);
print STDERR " ==> found $Ngg genes with >1 DIRECT annotations!\n";



#========================================================
# Count combinations
#========================================================
print STDERR "COUNTING COMBINATIONS ...\n";

# column 1 : entrezid
# column 3 : direct & indirect annotations
`cut -f1,3 /tmp/$out1 |sort -u > xxx`;
open(IN,"xxx") || die "Could not read xxx file !!\n";
$nl = `cat xxx | wc -l`;
chomp $nl;

$first=1;
$prec = '';
$k=0;
open(IN,"xxx") || die "Could not read xxx !\n";
while ($line=<IN>) {
    print STDERR "$k    / $nl      \r";
    $k++;

    chomp $line;
    @line = split (/\t/,$line);

    next unless ($goodgene{$line[0]});  # we only consider genes with > 1 DIRECT annotations

    if ($line[0] eq $prec) {
	$terms{$line[1]}=1;
    }
    elsif (!$first) {
	@terms = keys %terms;

	for $i (0..$#terms-1) {
	    $nocc{$terms[$i]}++;
	    for $j ($i+1..$#terms) {
		@a = sort {$a lt $b} ($terms[$i],$terms[$j]);
		@b = sort {$b lt $a} ($terms[$i],$terms[$j]);
		$a = join("-",@a);
		$b = join("-",@b);

# we only take into account combinations of terms that are not descendant - ancestor !!
		if ((not exists $related{$a}) && (not exists $related{$b})) {
		    $n{$a}++;
		}
	    }
	}
	$nocc{$terms[$#terms]}++;
	%terms = ();
    }
    $first=0;
    $prec = $line[0];
}
for $i (0..$#terms-1) {
    for $j ($i+1..$#terms) {
	@a = sort {$a lt $b} ($terms[$i],$terms[$j]);
	@b = sort {$b lt $a} ($terms[$i],$terms[$j]);
	$a = join("-",@a);
	$b = join("-",@b);
	if ((not exists $related{$a}) && (not exists $related{$b})) {
	    $n{$a}++;
	}
    }
}

# OK, where are we now ?
# we have:
# $Ngg : number of genes with >1 DIRECT ANNOTATIONS ("good genes")
# %nocc : number of occurences of all GO terms (within the set of good genes)
# %n : number of occurences of pairs of terms;

open(OUT,">$outfile");
foreach $k (keys %n) {
    $n = $n{$k};
    @k = split(/-/,$k);
    @K = sort @k;
    $o1 = $n*$Ngg/$nocc{$K[0]}/$nocc{$K[1]};
   # $o = ($n*$Ngg)/($nocc{$K[0]}*$nocc{$K[1]});
    print OUT "$K[0]_$K[1]\t$Ngg\t$nocc{$K[0]}\t$nocc{$K[1]}\t$n\t$o1\n";#\t$o\t$o1\t$n*$Ngg/$nocc{$K[0]}/$nocc{$K[1]}\t($n*$Ngg)/($nocc{$K[0]}*$nocc{$K[1]});\n";
}
close OUT;
print STDERR "--> Written $outfile\n";
