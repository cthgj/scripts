#! /usr/bin/gawk -f 
BEGIN {}     

{
  
  if ( /^>/)   {
 while (sub(/	/, " " )) {} #replacing tab chars
 ind=0
 initial_title=$1
 while (arr[$1]) {$1=initial_title "_" ind; ind+=1  }  
 arr[$1]=1
}
print $0

}
