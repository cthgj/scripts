#! /usr/bin/gawk -f

BEGIN{OFS="\t"}
{if(NR==1){
    printf "%s%s",$1,OFS;
    for(i=2;i<=NF;i++){
	k[i]=$(i);
	printf "%s%s",$(i),OFS;
    }
    printf "\n"; 
    next
} 
{
    for(i=2;i<=NF;i++){s[$1][k[i]]+=$(i); names[$1]++;}
}
END{
    for(i in names){
	printf "%s%s",i,OFS; 
	for(l in s[i]){printf "%s%s", s[i][l],OFS;}
	printf "\n";
    }
}}
