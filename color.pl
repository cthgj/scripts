#!/usr/bin/perl -w
use Getopt::Std;
use strict;
use Term::ANSIColor; 

my %opts;
getopts('hl:c',\%opts);
    if ($opts{h}){
	print "Use -l to specify the letter(s) to highlight. To specify more than one patern use spaces. -c makes the search case sensitive\n";
	die();
    }
my $case_sensitive=$opts{c}||undef;
my @color=("bold blue",'bold red', 'bold yellow', 'bold green', 'bold magenta','yellow on_magenta', 'white on_black', );
my @patterns;
if($opts{l}){
     @patterns=split(/\s+/,$opts{l});
}
else{
    $patterns[0]='\*';
}

# Setting $| to non-zero forces a flush right away and after 
# every write or print on the currently selected output channel. 
$|=1;

while (my $line=<>) 
{ 
    for (my $c=0; $c<=$#patterns; $c++){
	if($case_sensitive){
	    if($line=~/$patterns[$c]/){
		$line=~s/($patterns[$c])/color("$color[$c]").$1.color("reset")/ge; 
	    }
	}
	else{
	    if($line=~/$patterns[$c]/i){
		$line=~s/($patterns[$c])/color("$color[$c]").$1.color("reset")/ige; 
	    }
	}
    }
    print STDOUT $line;

}
