#!/bin/gawk -f

# This program extracts potential SP alignments between protein
# query and EST subject. We search all possible * in the EST 
# subject and check if they align to C in the query. 
# This shell works for multiple *.


### Define cutoff values ###

BEGIN{
    SCOREMIN=0;
    SCOREMAX=10000;
    EVALUE=1;
    IDENTITY=1;
    PAT="[\*]"; # aa in the query seq.
    bit=0;
    MATCH="\*"; ## aa in subject that PAT should match
    
}


### Extract alignment info ###

# Query ID

($1=="Query="){ # Curly brackets must go after the pattern !!!
    query_ID=$0;
    best_score=0;

}

# Identities

$1=="Identities" {
    gsub(/\,/,"",$4);
    identities=substr($4,2,length($4)-3);
}

# Length queryls -l

$2=="letters)"{

gsub(/\(/,"",$1);  
    query_len=$1;
}

# Subject ID

(substr($1,1,1)==">"){

    sbjct_ID=$0;

    while (substr($0,2,5)!="Score" && $1!="Length") { # to get all EST info

	getline;
	
	sbjct_ID=sbjct_ID"\n"$0;

    }

}

# Score and E-value



substr($0,2,5)=="Score"{
    gsub(/^ /,"",$0);
    score=$0;
    score_=$3;
    if (score_ >= best_score) {best_score =  score_; }
    #expect= substr($8,1,length($8)-1);  
    expect=$8;
}

# frame

($1=="Frame"){
  frame=$0;

}


(Identities != IDENTITY) {

# Get alignment

if ($1=="Query:"){

	query=$0;
	pos_start_query=$2;
	pos_end_query=$4;
	getline; middle_line=$0; # Read next line
	getline; sbjct_align=$0; # Read next line

### Search for Cys-Sec aln ###

# For each aln block (Query-middle-Subjct)

# Get subject * position
# Many * may exist in the EST
# We want to check all them
# and get, if any, the one
# that aligns to a C

     string=sbjct_align;
    sc=index(string,MATCH);

    out=0;
    
    while (sc != 0 && out==0)
	{
	    if (substr(query,sc,1) ~ PAT)
		out=1;

	    else
		{
		    sub(MATCH,"X",string) # mask leftmost *		    		    
		    sc=index(string,MATCH);

		    if (substr(query,sc,1) ~ PAT)
			out=1;
		}
 
	}

### Classify alignments ###

# Sec alignment (C-* or whatever we are interested)

   
  #x  print "score : " score_ "\n" "b_score : " best_score; 
    if ((identities >= IDENTITY) && (sc != 0) && (bit==0) && (score_>SCOREMIN) && (score_<SCOREMAX) && (expect+0<EVALUE) &&((substr(query,sc,1) ~ PAT)) && (sbjct_ID != query_ID_Sec) && (score_ >= best_score)) 
       { 
	     query_Sec=query; # Keep Sec alignment
	     middle_line_Sec=middle_line;
	     sbjct_align_Sec=sbjct_align;
	     query_ID_Sec=query_ID; # Identify alignment
	     sbjct_ID_Sec=sbjct_ID;# Identify alignment
	     score_Sec=score_;
	     identities_Sec=identities;
	     expect_Sec=expect;
	     pos_start_query_Sec=pos_start_query;
	     pos_end_query_Sec=pos_end_query;
	     bit=1; # last alignment was a Sec alignment, get next

          # print alignment

	     if (query_ID_last==query_ID_Sec && sbjct_ID_last==sbjct_ID_Sec && score_Sec==score_last && identities_Sec==identities_last && expect_Sec==expect_last && pos_start_query_Sec > pos_end_query_last) # check we are in the same HSP
		 {

		     print query_ID_Sec, "("query_len" aa)";
		     print "\n";
		     print sbjct_ID_Sec;
		     print frame;
		     print "\n";
		     print score, identities"%";
		     print "\n"; 
		     print query_last;		    
		     print middle_line_last;
		     print sbjct_align_last;
		     print "\n";
		     print query_Sec;
		     print middle_line_Sec;
		     print sbjct_align_Sec;
		     print "\n";

		 
		 }

	     else
		 {

		    	print query_ID_Sec, "("query_len" aa)";
			print "\n";
			print sbjct_ID_Sec;
			print "\n";
			print frame;
			print score, identities"%";
			print "\n";   
			print query_Sec;
			print middle_line_Sec;
			print sbjct_align_Sec;
			print "\n";

		 }

	 }

# last alignment was a Sec alignment, we keep next

     else if (bit==1) 
	 {
	     query_ID_next=query_ID;
	     sbjct_ID_next=sbjct_ID;
	     score_next=score_;
	     identities_next=identities;
	     expect_next=expect;
	     pos_start_query_next=pos_start_query;
	     pos_end_query_next=pos_end_query;

# Check whether next alignment belongs to the Sec alignment HSP

	     if (query_ID_Sec==query_ID_next && sbjct_ID_Sec==sbjct_ID_next && score_next==score_Sec && identities_next==identities_Sec && expect_Sec==expect_next && pos_end_query_Sec < pos_start_query_next) # check we are in the same HSP
		 {

		  	print query;
			print middle_line;
			print sbjct_align;				
			print "\n";
			print "###################################################\n"
			bit=0;

# Change last 'cos we have now printed 

	     query_last=query;
	     middle_line_last=middle_line;
	     sbjct_align_last=sbjct_align;
	     query_ID_last=query_ID;
	     sbjct_ID_last=sbjct_ID;
	     score_last=score_;
	     identities_last=identities;
	     expect_last=expect;
	     pos_start_query_last=pos_start_query;
	     pos_end_query_last=pos_end_query;
		 }

# Special case when an *-* alignment follows another *-* in different HSPs

	     else if ((identities >= IDENTITY) && (sc != 0) && (score_>SCOREMIN) && (score_<SCOREMAX) && (expect+0<EVALUE) &&((substr(query,sc,1) ~ PAT)) && (sbjct_ID!=query_ID_Sec) && (score_ >= best_score))
		 {
		     	     query_Sec=query; # Keep Sec alignment
			     middle_line_Sec=middle_line;
			     sbjct_align_Sec=sbjct_align;
			     query_ID_Sec=query_ID; # Identify alignment
			     sbjct_ID_Sec=sbjct_ID;# Identify alignment
			     score_Sec=score_;
			     identities_Sec=identities;
			     expect_Sec=expect;

		        print "###################################################\n"
		     	print query_ID_Sec, "("query_len" aa)";
			print "\n";
			print sbjct_ID_Sec;
			print frame;
			print "\n";
			print score, identities_Sec"%";
			print "\n";   
			print query_Sec;
			print middle_line_Sec;
			print sbjct_align_Sec;
			print "\n";

		 }

# keep last non Sec alignment 'cos next alignment maybe a Sec alignment

	     else 
		 {
		     query_last=query;
		     middle_line_last=middle_line;
		     sbjct_align_last=sbjct_align;
		     query_ID_last=query_ID;
		     sbjct_ID_last=sbjct_ID;
		     pos_start_query_last=pos_start_query;
		     pos_end_query_last=pos_end_query;
		     print "###################################################\n"
		     bit=0; 
		 }

	 }

# keep last non Sec alignment 'cos next alignment maybe a Sec alignment

     else 
	 {
	     query_last=query;
	     middle_line_last=middle_line;
	     sbjct_align_last=sbjct_align;
	     query_ID_last=query_ID;
	     sbjct_ID_last=sbjct_ID;
	     score_last=score_;
	     identities_last=identities;
	     expect_last=expect;
	     pos_start_query_last=pos_start_query;
	     pos_end_query_last=pos_end_query;
	     
	 }    



}

}

