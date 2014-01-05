#!/usr/bin/perl -w

use Inline C;

my @thisone=('r','a','b','c');
my @GOs=("GO:0000001","GO:0000002","GO:0000003","GO:0000004","GO:0000005","GO:0000006","GO:0000007","GO:0000008","GO:0000009");
my $c='dabcdefg';
my @kk=(1,2,3,4);
#ff($c);
#ff("@GOs");
ff(\@GOs);
__END__
__C__

void ff(SV* array_ref){
  SV* sv_array;
  SV** sv_1;
  SV* kk;
  if (! SvROK(array_ref))
    croak("array_ref is not a reference");
  sv_array = (AV*)SvRV(array_ref);
 
  int max=sizeof(kk)/sizeof(char);
  int i;
  for(i=0;i<=max; i++){
    //    printf("AA, %d %c: %s\n",i,array_ref[i],array_ref);
    sv_1=av_fetch(sv_array,i,0); 
    kk=*sv_1;
    printf("AA, %d:%s\n",i,SvPV(kk,PL_na)); 
  }
}
