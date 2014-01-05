#!/usr/bin/perl

use Fcntl qw( SEEK_SET );

my($fh, $filename, $byte_position, $byte_value);

$filename = $ARGV[0];

open(IN, "+>", $filename);
#open IN, $filename;

seek(IN,0,SEEK_SET);
read IN, $temp, 128;

print $temp;
print "\n";

seek(IN,14,SEEK_SET);
read IN, $temp, 16;
print "Artist is :" .$temp;
print "\n";


sysseek(IN,14,SEEK_SET);
#want to replace the Artist Name with new one.
syswrite (IN,$newArtist);

print "\n";

close(IN);
