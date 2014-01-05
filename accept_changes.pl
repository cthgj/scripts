#!/usr/bin/perl -w
###########################################################
# This script will accept changes in a latex file. 	  #
#                                                         #
# \del{TEXT} is a deletion and will be deleted		  #
# \cc{TEXT} is the correction and will be kept		  #
# \com{TEXT} is a comment and will be deleted		  #
# 							  #
# Use '-i' for interactive mode, it will prompt you for	  #
# and actin for each instance. Otherwise, ALL the changes #
# will be accepted automatically.			  #
###########################################################

use strict;
use Getopt::Std;
use Term::ANSIColor; 



my %opts;
getopts('i',\%opts);


$/=undef;

open(B,"$ARGV[0]")|| die("cannot open $ARGV[0] for reading, no such file\n");
while(my $line=<B>){
    while( $line=~/(\\del\{.+?\})/s){
	$line=&replace($line,$1,"\\\\del");
    }
    while($line=~/(\\com\{.+\})/s){
	$line=&replace($line,$1,"\\\\com");
    }
    while($line=~/(\\cc\{.+?\})/s){
	$line=&replace($line,$1,"\\\\cc");
    }
    print $line;
}


sub replace{
    my $line=shift;
    my $found=shift;
    my $pat=shift;
    my @a=split(//s,$found);
    my $open=grep(/\{/s,@a);
    my $closed=grep(/\}/s,@a);
    my $i=$open-$closed;
    while($i>0){
	$found=~s/\\/\\\\/sg;
	$found=~s/\`\`/\\`\\`/sg;
	$found=~s/\'\'/\\'\\'/sg;
	$found=~s/\(/\\(/sg;
	$found=~s/\)/\\)/sg;

	$line=~/($found.*?\})/s||die("No match: $found :: $line\n");
	$found=$1;
	$found=~s/\\\\/\\/sg;
	my @aa=split(//s,$found);
	$open=grep(/\{/s,@aa);
	$closed=grep(/\}/s,@aa);
	$i=$open-$closed;
    }
    my $a=quotemeta($found);
    if($pat=~/del/s){
	if($opts{i}){
	    my $from;
	    $line=~/\n(.*?$a.*?)\n/s;
	    if($line=~/\n(.*?$a.*?)\n/s){$from=$1;}
	    elsif($line=~/^(.*?$a.*?)\n/s){$from=$1;}
	    elsif( $line=~/\n(.*?$a.*?)$/s){$from=$1;}
	    
	    my $o=0;
	    while($o==0){
		$/="\n";
		my $aa=unquotemeta($a);
		print STDERR "\nDelete '$aa'\nfrom $from? (Y/n): ";
		my $yesno=<STDIN>;
		chomp($yesno);
		print STDERR "\n\n";
		$/=undef;
		if($yesno eq 'Y' || $yesno eq 'y' || $yesno eq ''){
		    $line=~s/$a//s;
		    $o++;
		}
		elsif($yesno eq 'N' || $yesno eq 'n'){
		    $o++;
		}
		else{ print STDERR "Please type Y or N\n";}
	    }
	}
	else{
	    $line=~s/$a//s;
	}

    }
    elsif($pat=~/com/s){
		if($opts{i}){
	    my $from;
	    $line=~/\n(.*?$a.*?)\n/s;
	    if($line=~/\n(.*?$a.*?)\n/s){$from=$1;}
	    elsif($line=~/^(.*?$a.*?)\n/s){$from=$1;}
	    elsif( $line=~/\n(.*?$a.*?)$/s){$from=$1;}
	    
	    my $o=0;
	    while($o==0){
		$/="\n";
		my $aa=unquotemeta($a);
		print STDERR "\nRemove comment '$aa' (Y/n): ";
		my $yesno=<STDIN>;
		chomp($yesno);
		print STDERR "\n\n";
		$/=undef;
		if($yesno eq 'Y' || $yesno eq 'y' || $yesno eq ''){
		    $line=~s/$a//s;
		    $o++;
		}
		elsif($yesno eq 'N' || $yesno eq 'n'){
		    $o++;
		}
		else{ print STDERR "Please type Y or N\n";}
	    }
	}
	else{
	    $line=~s/$a//s;
	}
    }
    else{
	$found=~s/$pat\{(.+)\}/$1/s||die("crap: $found : $pat : $1\n");
	my $keep=$1;
	if($opts{i}){
	    my $from;
	    $line=~/\n(.*?$a.*?)\n/s;
	    if($line=~/\n(.*?$a.*?)\n/s){$from=$1;}
	    elsif($line=~/^(.*?$a.*?)\n/s){$from=$1;}
	    elsif( $line=~/\n(.*?$a.*?)$/s){$from=$1;}
	    
	    my $o=0;
	    while($o==0){
		$/="\n";
		my $aa=unquotemeta($a);
		print STDERR "\nInsert '$aa'\ninto $from? (Y/n): ";
		my $yesno=<STDIN>;
		chomp($yesno);
		print STDERR "\n\n";
		$/=undef;
		if($yesno eq 'Y' || $yesno eq 'y' || $yesno eq ''){
		    $line=~s/$a/$keep/s;
		    $o++;
		}
		elsif($yesno eq 'N' || $yesno eq 'n'){
		    $o++;
		}
		else{ print STDERR "Please type Y or N\n";}
	    }
	}
	else{
	    $line=~s/$a/$keep/s;
	}

    }
    return($line);
}


sub unquotemeta{
    my $pat=shift;
    $pat =~ s/\\\\/kikikoko123/g;
    $pat =~ s/\\//g;
    $pat =~ s/kikikoko123/\\/g;
    $pat;
}
