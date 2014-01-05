#!/usr/bin/perl
#LowestDirs02.pl
use strict;
use warnings;
my $startdir = $ARGV[0] || ".";
print_dirs_without_subdirs($startdir);

sub print_dirs_without_subdirs {
    #This subroutine calls itself for each subdirectory it finds.
    #If it finds no <strong class="highlight">directories</strong>, it prints the name of the directory in which
    # it is looking, $dir.
    my $dir = $_[0];
    opendir DH, $dir or die "Failed to open $dir: $!";
    my @d;
    while ($_ = readdir(DH)) {
        next if $_ eq "." or $_ eq "..";
        my $fn = $dir . '/' . $_;
        if (-d $fn) {
            push @d, $fn;
        }
    }
    if (scalar @d == 0) { #If no <strong class="highlight">directories</strong> found, $dir is <strong class="highlight">lowest</strong> dir in this branch
        print "$dir\n";
        return;
    }
    foreach (@d) {
        print_dirs_without_subdirs($_); #Look for <strong class="highlight">directories</strong> in directory
    }
}
