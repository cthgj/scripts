#!/usr/bin/perl -w
use Inline C;

my @a=('a','b','c','d','e','f','g','h');

our $kk=12;
our $bb='bob';
c_loop();
__END__
__C__
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
SV* c_loop(){

SV* super_i;
super_i = get_sv("kk", TRUE);
SvIV_set(super_i, 42);
SvIOK_on(super_i);
printf("Sv %d\n",SvIV(super_i));

SV* text;
text = get_sv("bb", TRUE);
 printf("Sv %s\n",SvPV(text,PL_na));




//  kkk=get_sv("kk", 0);  
  //  int bob=SvOK(*kkk);

  // printf("AA %d\n",bob);
  //  printf("Sv %d\n",SvIV(kkk));
  fflush(stdout);
}
 

