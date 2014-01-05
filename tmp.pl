#!/usr/bin/perl -w



@names  =('Mrs Smith','Mr Jones','Ms Samuel','Dr Jansen','Sir Philip');
@medics =('Dr Black','Dr Waymour','Dr Jansen','Dr Pettle');

LBL: foreach $person (@names) {
print "$person\n";
if ($person=~/Dr /) {
foreach $doc (@medics) {
print "\t$doc\n";
last LBL if $doc eq $person;
}
}
}
