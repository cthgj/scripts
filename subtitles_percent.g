#!/usr/bin/gawk -f
#define variable p to change subtitles timing in percentages . msec= msec*p          
#define variable add to move them of an number of minutes.

#e.g. subtitles_percent.g -v add=-20 -v p=-1.2    
BEGIN {FS=":" 
if (!p){p=1}  }
{
if (/[0-9]+\:/)    {

    split($0, TIMES, "-->" )
    out=""
    for (t=0; t < length (TIMES ); t++){
        $0 = TIMES [t+1]
        split($NF,  ARR , ","   )
        msec=$1*3600000+$2*60000+ARR[1]*1000+ARR[2]
	msec+=add*60000
        msec= int(msec*p)
	print "msec :",msec;
        h=int(msec/3600000)
        msec-=h*3600000
	if (h<10) {h= 0 "" h }
        m=int(msec/60000)
        msec-=m*60000
	if (m<10) {m= 0 "" m }
        s= int(msec/1000)
        msec-=s*1000
	if (s<10) {s= 0 "" s }
	while ( length(msec "" )<3 ) { msec = msec "0"} 

        out = out h ":" m ":" s "," msec 

        if (t==0){out = out " --> "}

        }
   $0= out

}

    print

 
}
