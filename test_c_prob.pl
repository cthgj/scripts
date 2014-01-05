#!/usr/bin/perl 
## Get the numbers necessary to calculate the probabilities of association of each GO pair


use Getopt::Std;
use FileHandle;
use Inline C;

my %opts;
getopts('hAaOvCo:c:s:g:S:d:',\%opts);


my @GOs=("GO:0000001","GO:0000002","GO:0000003","GO:0000004","GO:0000005","GO:0000006","GO:0000007","GO:0000008","GO:0000009");

my $tot= scalar(@GOs);
print STDERR "GOs: " . scalar(@GOs) . "\n"; 
my $golist_string=join("_",@GOs);
#while($tot>0){my $str=get_string($golist_string,scalar(@GOs)); $tot=-1;}
greet(@GOs);

__END__
__C__

void greet(SV* name1, ...) {
  //First, we need to begin our function with a "Inline_Stack_Vars" statement. This defines a few internal variables that we need to access the Stack. Now we can use "Inline_Stack_Items", which returns an integer containing the number of arguments passed to us from Perl.
  //NOTE: It is important to only use "Inline_Stack_" macros when there is an ellipsis (...) in the argument list, or the function has a return type of void.
  Inline_Stack_Vars;
  int i;
  int k;
  
  for (i = 0; i < Inline_Stack_Items; i++) {
    fprintf(stderr,"%d of %d   \r",i,Inline_Stack_Items);
       for (k = i+1; k < Inline_Stack_Items; k++) {
      char pair[100];
      char aa[10];
      strcpy(aa,"bob");

      //Get GO pair
      strcpy(pair,SvPV(Inline_Stack_Item(i),PL_na)) ; // ,"_"SvPV(Inline_Stack_Item(k),PL_na));
      strcat(pair,"_");
      strcat(pair,SvPV(Inline_Stack_Item(k),PL_na));
      //printf("aa %s\n",pair);
      printf("Hello %s %s!\n", SvPV(Inline_Stack_Item(i),PL_na),SvPV(Inline_Stack_Item(k),PL_na));
    }
  }
  
  Inline_Stack_Void;
}


/* char get_string(char* golist, long len){ */
/*     int i=0; */
/*   /\* for(i=0; i<len; i++){ *\/ */
/*   /\*   printf("HIya %d %s (%d)\n",i,golist,len); *\/ */
/*   /\* } *\/ */
/*     char normalize = 'a' ^ 'A'; */
/*     char sep =  '_'; */
/*     char *gos[100000][10000]; */
/*     int count=0; */
/*     char *c; */
/*     char *haha="haha"; */
/*     char goterm; */
/*     char test[10]; */
/*     char str='bob'; */
/*     strcat(test,"bob"); */
/*     printf("TEST is %c %c\n",test[2],str); */
/*     /\* make go list *\/ */
/*     while(c = golist[i++]) { */
/*       int aa; */
      
/*       /\* If the current char is '_' *\/ */
/*       if(aa=strcmp(sep,c) == 0){ */
/* 	count=0; */
/*       } */
/*       else{ */
/* 	//strcat(gos[i][0],c); */
/*       } */
/*       count++; */
/*       /\* printf("AA %c\n",c); *\/ */
/*       /\*  c |= normalize; *\/ */
/*       /\*  if(c == sep){printf("YES (%d) %s\n",i,golist);} *\/ */
/*       /\*  else{printf("NO (%d): %d\n",i,c);} *\/ */
/*     } */
/*     /\* int k, kk; *\/ */
/*     /\* count=0; *\/ */
/*     /\* char cc; *\/ */
/*     /\* for(k=0; k<=i; k++){ *\/ */
/*     /\*       for(kk=0; kk<=i; kk++){ *\/ */
/*     /\* 	    printf("A %s\n",gos[k][kk]); *\/ */
/*     /\* 	    //	    printf("A %d %d\n",k,kk); *\/ */
/*     /\* 	  } *\/ */
/*     /\* } *\/ */
/* } */
