#!/usr/bin/perl -w


# Takes a list of cluster names, extracts the sequences of all member sequences and tblastxes them
# against each other.

use strict;

my $egoseq = '/home/ug/cchapple/research/seq/4ego/ego4_060903.seq';
my %clus = ();
my %ego = ();
my $p;
my $d;



open(FILE, "/home/ug/cchapple/research/selenoproteins/4ego/erpin/ego4_060903.fa.erpin.nonsel");
my @clusters = <FILE>;
close(FILE);


open(GILE,"/home/ug/cchapple/research/seq/4ego/clusters");
while(<GILE>)
{
    foreach my $cluster (@clusters)
    {
	chomp $cluster;
	if (/^$cluster\t(.*?)$/)
	{
	    my @f = split /\s+/og, $1;
	    $clus{$cluster} = [@f];
	}
    } 
}



open(SEQ,"$egoseq") || die "0 cannot open seq: $!";
while (<SEQ>)
{
    /^(.*?)\s(.*?)$/;
    $ego{$1} = $2;     ### $ego{seq name} = sequence
}
close (SEQ);
print STDERR "*** get full seq done\n";


############################################################
############################################################



unlink("query\.tbl");  ### Delete files possibly present from previous instance 
unlink("dbase\.tbl"); 

mkdir("blast", 0777) unless (-d "blast");



foreach my $cluster (@clusters)
{
   
    my $nuum = scalar(@{$clus{$cluster}});
    my $num = $nuum / 2;
    print STDERR "********** \$num is $num **********\n";
    $p = 0;

    foreach my $seq (@{$clus{$cluster}})
    {
	open(DB, ">>dbase.tbl");
	open(TBL,">>query.tbl");
	print STDERR "\$seq : $seq\n";
	open (TBL, ">>query\.tbl")|| die "2 cannot open query: $!";
	open(DB, ">>dbase\.tbl");

	if($p < $num)
	{
	    print STDERR "\$p: $p \$num : $num\n";
	    print TBL "$seq\t$ego{$seq}\n";
	    $p++;
	}
	else
	{
	    print DB "$seq\t$ego{$seq}\n";
	    $p++;
	}
	close(TBL);
	close(DB);
    }
system ("TblToFasta query\.tbl > query\.fa") && die "shit1 : $!\n";
system ("TblToFasta dbase\.tbl > dbase\.fa") && die "shit2 : $!\n";

unlink("query\.tbl"); 
unlink("dbase\.tbl"); 
system("pressdb dbase\.fa");

system ("wu-tblastx dbase\.fa query\.fa > blast/$cluster\.out"); 

unlink("dbase\.fa\.ntb");
unlink("dbase\.fa\.csq");
unlink("dbase\.fa\.nhd");
unlink("dbase\.fa");
unlink("query\.fa");
print STDERR "****** One Down!!!!! ******\n"; 



}
