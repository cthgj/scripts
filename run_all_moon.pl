#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my %opts;
getopts('hHivo:d:',\%opts) || do { print "Invalid option, try '$0 -h' for more information\n"; exit(1); };
&usage() if $opts{h};
#&usage() unless $ARGV[0];
my $verbose=$opts{v}||undef;
my $which='';
my $only_human=$opts{H}||undef;
my $options=$opts{o}||'vc';
my $output_dir=$opts{d}||'.';
$output_dir=~s/\/$//;

my (@nets,@species);
if($only_human){
    @nets=('HSHQ.gr');
    @species=('human');
}
else{
    @nets=('CEHQ.gr','DMHQ.gr','HSHQ.gr','MMHQ.gr','SCHQ.gr');
    @species=('worm','fly','human','mouse','yeast');
}
@nets=('HSHQ.gr');
@species=('human');

my @modes=('a','b','B','c','d','i','m','p','P');
$opts{i} && do {$options .=' -i'};

for (my $i=0; $i<=$#species; $i++){
    foreach my $mode (@modes){
	my $out_file;
	$opts{i} ? ($out_file= "$output_dir/$species[$i].$mode.inter.moon") :
	    ($out_file= "$output_dir/$species[$i].$mode.moon") ;
	print STDERR "\n----------- Mode $mode on $species[$i] ----------\n" if $verbose;
	print STDERR "moonGO.pl -m $mode -$options $nets[$i] > $out_file\n" if $verbose;
	system("moonGO.pl -m $mode -$options $nets[$i] > $out_file");

   
    }

}



sub usage{
    print "Usage: run_all_moon.sh [hivod] <OUTPUT_DIR>\n";exit;

}

# moonGO.pl -vm i DMHQ.gr > results/fly.i.moon
# moonGO.pl -vm b DMHQ.gr > results/fly.b.moon
# moonGO.pl -vm d DMHQ.gr > results/fly.d.moon
# moonGO.pl -vm m DMHQ.gr > results/fly.m.moon
# moonGO.pl -vm c DMHQ.gr > results/fly.c.moon
# moonGO.pl -vm P DMHQ.gr > results/fly.P.moon
# moonGO.pl -vm p DMHQ.gr > results/fly.p.moon
# moonGO.pl -vm a DMHQ.gr > results/fly.a.moon

# moonGO.pl -vm i MMHQ.gr > results/mouse.i.moon
# moonGO.pl -vm b MMHQ.gr > results/mouse.b.moon
# moonGO.pl -vm d MMHQ.gr > results/mouse.d.moon
# moonGO.pl -vm m MMHQ.gr > results/mouse.m.moon
# moonGO.pl -vm c MMHQ.gr > results/mouse.c.moon
# moonGO.pl -vm P MMHQ.gr > results/mouse.P.moon
# moonGO.pl -vm p MMHQ.gr > results/mouse.p.moon
# moonGO.pl -vm a MMHQ.gr > results/mouse.a.moon

# moonGO.pl -vm i SCHQ.gr > results/yeast.i.moon
# moonGO.pl -vm b SCHQ.gr > results/yeast.b.moon
# moonGO.pl -vm d SCHQ.gr > results/yeast.d.moon
# moonGO.pl -vm m SCHQ.gr > results/yeast.m.moon
# moonGO.pl -vm c SCHQ.gr > results/yeast.c.moon
# moonGO.pl -vm P SCHQ.gr > results/yeast.P.moon
# moonGO.pl -vm p SCHQ.gr > results/yeast.p.moon
# moonGO.pl -vm a SCHQ.gr > results/yeast.a.moon

# moonGO.pl -vm i CEHQ.gr > results/worm.i.moon
# moonGO.pl -vm b CEHQ.gr > results/worm.b.moon
# moonGO.pl -vm d CEHQ.gr > results/worm.d.moon
# moonGO.pl -vm m CEHQ.gr > results/worm.m.moon
# moonGO.pl -vm c CEHQ.gr > results/worm.c.moon
# moonGO.pl -vm P CEHQ.gr > results/worm.P.moon
# mmoonGO.pl -vm p CEHQ.gr > results/worm.p.moon
# moonGO.pl -vm a CEHQ.gr > results/worm.a.moon

# cp results/human* results/old

# moonGO.pl -vm i HSHQ.gr > results/human.i.moon
# moonGO.pl -vm b HSHQ.gr > results/human.b.moon
# moonGO.pl -vm d HSHQ.gr > results/human.d.moon
# moonGO.pl -vm m HSHQ.gr > results/human.m.moon
# moonGO.pl -vm c HSHQ.gr > results/human.c.moon
# moonGO.pl -vm P HSHQ.gr > results/human.P.moon
# moonGO.pl -vm p HSHQ.gr > results/human.p.moon
# moonGO.pl -vm a HSHQ.gr > results/human.a.moon



# Interactome:
# moonGO.pl -vtm i DMHQ.gr > results/fly.i.inter.moon
# moonGO.pl -vtm b DMHQ.gr > results/fly.b.inter.moon
# moonGO.pl -vtm d DMHQ.gr > results/fly.d.inter.moon
# moonGO.pl -vtm m DMHQ.gr > results/fly.m.inter.moon
# moonGO.pl -vtm c DMHQ.gr > results/fly.c.inter.moon
# moonGO.pl -vtm P DMHQ.gr > results/fly.P.inter.moon
# moonGO.pl -vtm p DMHQ.gr > results/fly.p.inter.moon
# moonGO.pl -vtm a DMHQ.gr > results/fly.a.inter.moon

# moonGO.pl -vtm i MMHQ.gr > results/mouse.i.inter.moon
# moonGO.pl -vtm b MMHQ.gr > results/mouse.b.inter.moon
# moonGO.pl -vtm d MMHQ.gr > results/mouse.d.inter.moon
# moonGO.pl -vtm m MMHQ.gr > results/mouse.m.inter.moon
# moonGO.pl -vtm c MMHQ.gr > results/mouse.c.inter.moon
# moonGO.pl -vtm P MMHQ.gr > results/mouse.P.inter.moon
# moonGO.pl -vtm p MMHQ.gr > results/mouse.p.inter.moon
# moonGO.pl -vtm a MMHQ.gr > results/mouse.a.inter.moon

# moonGO.pl -vtm i SCHQ.gr > results/yeast.i.inter.moon
# moonGO.pl -vtm b SCHQ.gr > results/yeast.b.inter.moon
# moonGO.pl -vtm d SCHQ.gr > results/yeast.d.inter.moon
# moonGO.pl -vtm m SCHQ.gr > results/yeast.m.inter.moon
# moonGO.pl -vtm c SCHQ.gr > results/yeast.c.inter.moon
# moonGO.pl -vtm P SCHQ.gr > results/yeast.P.inter.moon
# moonGO.pl -vtm p SCHQ.gr > results/yeast.p.inter.moon
# moonGO.pl -vtm a SCHQ.gr > results/yeast.a.inter.moon

# moonGO.pl -vtm i CEHQ.gr > results/worm.i.inter.moon
# moonGO.pl -vtm b CEHQ.gr > results/worm.b.inter.moon
# moonGO.pl -vtm d CEHQ.gr > results/worm.d.inter.moon
# moonGO.pl -vtm m CEHQ.gr > results/worm.m.inter.moon
# moonGO.pl -vtm c CEHQ.gr > results/worm.c.inter.moon
# moonGO.pl -vtm P CEHQ.gr > results/worm.P.inter.moon
# moonGO.pl -vtm p CEHQ.gr > results/worm.p.inter.moon
# moonGO.pl -vtm a CEHQ.gr > results/worm.a.inter.moon

# moonGO.pl -vtm i HSHQ.gr > results/human.i.inter.moon
# moonGO.pl -vtm b HSHQ.gr > results/human.b.inter.moon
# moonGO.pl -vtm d HSHQ.gr > results/human.d.inter.moon
# moonGO.pl -vtm m HSHQ.gr > results/human.m.inter.moon
# moonGO.pl -vtm c HSHQ.gr > results/human.c.inter.moon
# moonGO.pl -vtm P HSHQ.gr > results/human.P.inter.moon
# moonGO.pl -vtm p HSHQ.gr > results/human.p.inter.moon
# moonGO.pl -vtm a HSHQ.gr > results/human.a.inter.moon

