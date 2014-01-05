#!/usr/bin/perl -w

$acpi=`acpi 2>/dev/null`;
if($acpi =~/Dischar/ ){
    $acpi =~/(\d+%).\s+(.+?)\srem/;
    print "$2, $1";
}
elsif($acpi =~/(\d+%).\s+(.+?)\suntil/ ){
    print "Charging, $1 ($2)";
}
else{
    print "connected";
}
