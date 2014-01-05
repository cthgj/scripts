#!/usr/bin/perl
 #  $gopair e.g. : GO:0000001_GO:0050876
use Getopt::Std;
use strict;

my %opts;
getopts('mF:D:h',\%opts);
my %prob;
print STDERR "USAGE: $0 (-D gostats dir) goterms_file||go1_go2\n" if $opts{h};
exit() if $opts{h};
my $load_entire_file=$opts{F}||undef;
my @GOs;
my %pairs;
my $no_mem=$opts{m}||undef;
my %seen;
if(-e "$ARGV[0]") {
    open(F, "$ARGV[0]");
    while(<F>){
	chomp;
	if(/_/){ 
	    $pairs{$_}++;
	}
	elsif(/.+\s*GO/){
	    my @a=split(/\s+/);
	    my @b = sort {$b lt $a} (@a);
	    my $gopair=join("_",@b);
	    $pairs{$gopair}++;
	}
	else{die("Unrecognized GO pair format :$_\n");}
    }
}
else{
    if($ARGV[0]=~/_/){ 
	$pairs{$ARGV[0]}++;
    }
    elsif($ARGV[1]){
	my @b = sort {$b lt $a} ($ARGV[0],$ARGV[1]);
	my $gopair=join("_",@b);
	$pairs{$gopair}++;
    }
    else{die("Unrecognized GO pair format :$_\n");}
}

my ($high,$low);
my $stats_dir=$opts{D}||undef;
my $c=0;
my $aa= scalar(keys(%pairs));
if($load_entire_file){
    open(A,"$load_entire_file");
    while(<A>){
	chomp;
	my @a=split(/\t/);
	if (defined($pairs{$a[0]})){
	     $c++;
	     $pairs{$a[0]}='found';
	     print STDERR "$c of $aa\r";
	     print STDOUT "$a[0]\t$a[6]\n";
	}
    }
    foreach my $gopair (keys(%pairs)){
	next if $pairs{$gopair} eq 'found';
	print  STDOUT "$gopair\t1\n";
    }
}
elsif($stats_dir){die();
    my @GOPAIRS=sort {$a lt $b} keys(%pairs);
  gopair:foreach my $gopair (@GOPAIRS){#keys(%pairs)){
      $c++;
      if(defined($prob{$gopair})){print "$gopair\t$prob{$gopair}\n"; next gopair}
      print STDERR "$c of " . scalar(keys(%pairs)) . "\r";
      my ($go1, $go2)=split(/_/,$gopair);
      if ($go1 eq $go2){
	  print "$gopair\t1\n";
      }
      $gopair=~/^GO:(.....)/||die("cannot match gopair : $gopair\n");
      my $file=$1 . ".prob";
      if(-e "$stats_dir/$file.gz") {
	  open(A,"zcat $stats_dir/$file |");
      }
      elsif(-e "$stats_dir/$file") {
	  open(A,"$stats_dir/$file");
      }
      else{ ## assume overrepresented
	  die("Cannot find $file\n");
	  print $gopair, "\t1\n";
	  next;
      }
      while(<A>){
	  chomp;
	  my @a=split(/\t/);
	  my @b=sort {$a lt $b} ($gopair,$a[0]);
	  if(/$gopair/){
	      $prob{$gopair}=$a[6] unless $no_mem;
	      print STDOUT "$gopair\t$a[6]\n";
	      next gopair;
	  }
	  unless ($no_mem){$prob{$a[0]}=$a[6] if defined{$pairs{$a[0]}} ;}
	  # unless($b[0] eq $gopair){
	  #     print STDOUT "$gopair\t1\t$b[0]\n";
	  #     next gopair;
	  # }
	  
	  
      }
      print "$gopair\t1\n" ;
      close(A);
  }
    #print "$gopair\t1\n" ;
}
else{print STDERR "Read the damn source code, no help yet\n"; exit();}


