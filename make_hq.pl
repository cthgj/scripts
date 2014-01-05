#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopts('sbd:',\%opts);
my $binary=$opts{b}||-1;
my $db=$opts{d}||-1;
my $print_skipped=$opts{s}||undef;
if (defined($opts{d}) && $opts{d}==0){$db=0 }
my %skipped;
my %hq=("12" => 1, "18" => 1, "30" => 1, "31" => 1, "47" => 1, "49" => 1, "53" => 1, "54" => 1, "55" => 1, "66" => 1, "71" => 1, "77" => 1, "81" => 1, "84" => 1, "89" => 1, "108" => 1, "111" => 1, "112" => 1, "114" => 1, "231" => 1, "369" => 1, "370" => 1, "397" => 1, "398" => 1, "405" => 1, "406" => 1, "412" => 1, "423" => 1, "424" => 1, "440" => 1);



my %list =("0" => "molecular interaction",
	   "4" => "affinity chromatography technology",
	   "6" => "anti bait coimmunoprecipitation",
	   "7" => "anti tag coimmunoprecipitation",
	   "13" => "biophysical",
	   "18" => "two hybrid",
	   "19" => "coimmunoprecipitation",
	   "22" => "colocalization by immunostaining",
	   "23" => "colocalization/visualisation technologies",
	   "25" => "copurification",
	   "27" => "cosedimentation",
	   "28" => "cosedimentation in solution",
	   "30" => "cross-linking study",
	   "40" => "electron microscopy",
	   "45" => "experimental interaction detection",
	   "47" => "far western blotting",
	   "49" => "filter binding",
	   "51" => "fluorescence technology",
	   "53" => "fluorescence polarization spectroscopy",
	   "59" => "gst pull down",
	   "61" => "his pull down",
	   "65" => "isothermal titration calorimetry",
	   "71" => "molecular sieving",
	   "77" => "nuclear magnetic resonance",
	   "84" => "phage display",
	   "91" => "chromatography technology",
	   "96" => "pull down",
	   "107" => "surface plasmon resonance",
	   "113" => "western blot",
	   "114" => "x-ray crystallography",
	   "226" => "ion exchange chromatography",
	   "254" => "genetic interference",
	   "339" => "undetermined sequence position",
	   "397" => "two hybrid array",
	   "398" => "two hybrid pooling approach",
	   "399" => "two hybrid fragment pooling approach",
	   "400" => "affinity technology",
	   "401" => "biochemical",
	   "403" => "colocalization",
	   "404" => "comigration in non denaturing gel electrophoresis",
	   "411" => "enzyme linked immunosorbent assay",
	   "413" => "electrophoretic mobility shift assay",
	   "415" => "enzymatic study",
	   "416" => "fluorescence microscopy",
	   "423" => "in-gel kinase assay",
	   "424" => "protein kinase assay",
	   "428" => "imaging techniques",
	   "434" => "phosphatase assay",
	   "663" => "confocal microscopy",
	   "676" => "tandem affinity purification",
	   "809" => "bimolecular fluorescence complementation",
	   "10010" => "reconstituted complex",
	   "10021" => "protein-RNA",
	   "10023" => "co-fractionation");

my %want;
while(<>){
    my @a=split(/\t/,$_);
    my $b=$a[0] . "\t" . $a[1] . "\n";
    my $c=$a[1] . "\t" . $a[0] . "\n";
    $b=$c if defined($want{$c});


    if(($binary!=-1) && ($db!=-1)){
	if(defined($hq{$a[4]})){
	    my @source=split(/\s+/,$a[6]);
	    $want{$b}++ if scalar(@source)>=$db;
	}
    }
    elsif($binary!=-1){
#	$want{$b}++ if defined($hq{$a[4]});
	$want{$b}++ if $a[4]=='493';
	$skipped{$a[4]}++;
    }
    elsif($db!=-1){
	my @source=split(/\s+/,$a[6]);
	$want{$b}++ if scalar(@source)>=$db;
    }
    else{print STDERR "Need to chose either binary (-b) interactions, or those present in multiple databases (>= -d), or both (eg,-bd 3)\n"; exit()}
    

}
if ($print_skipped){
    map{print "$_\t$list{$_}\n"}keys(%skipped);
}
else{
    print scalar(keys(%want)) . "\n";
    map{print}keys(%want);
}


## List
