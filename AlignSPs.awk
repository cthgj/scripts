#!/bin/gawk -f

# This program extracts potential SP alignments. We check for * in the query
# sequence aligned to other aa of interest (usually * or C). This shell works
# only when having one * per protein

# Define cutoff values

BEGIN{
    SCOREMIN=0;
    SCOREMAX=10000;
    EVALUE=1;
    IDENTITY=1;
    PAT="X";    ## aa in subject that PAT should match
    bit=0;
    MATCH="C";  # aa in the query seq.
}

# Extract alignment info

# Query ID

($1=="Query="){ # Curly brackets must go after the pattern !!!
    query_ID=$0;

}

# Length query

$2=="letters)"{

gsub(/\(/,"",$1);
    query_len=$1;

}

# Subject ID

(substr($1,1,1)==">"){
    sbjct_ID=$1;

}

# Length subject

$1=="Length" {
    subj_len=$3;
}

# Score and E-value

substr($0,2,5)=="Score"{
    gsub(/^ /,"",$0);
    score=$0;
    score_=$3;
    #expect= substr($8,1,length($8)-1);  
    expect=$8;
}

# Identities

$1=="Identities" {
    gsub(/\,/,"",$4);
    identities=substr($4,2,length($4)-3);
}


# Get alignment

($1=="Query:"){

    query=$0;
    getline; middle_line=$0; # Read next line
    getline; sbjct_align=$0; # Read next line

#    print query;
#    print middle_line;
#    print sbjct_align;

    sc=index(query,PAT); # Find * position in query seq

# Classify alignments

# Sec alignment (*-*, *-C or whatever we are interested)

     if ((query ~ PAT) && (bit==0) && (score_>SCOREMIN) && (score_<SCOREMAX) && (identities>=IDENTITY) && (expect+0<EVALUE) &&((substr(sbjct_align,sc,1) ~ MATCH))) 
	 {
	     query_Sec=query; # Keep Sec alignment
	     middle_line_Sec=middle_line;
	     sbjct_align_Sec=sbjct_align;
	     query_ID_Sec=query_ID; # Identify alignment
	     sbjct_ID_Sec=sbjct_ID;# Identify alignment
	     score_Sec=score_;
	     identities_Sec=identities;
	     expect_Sec=expect;
	     bit=1; # last alignment was a Sec alignment, get next

          # print alignment

	     if (query_ID_last==query_ID_Sec && sbjct_ID_last==sbjct_ID_Sec && score_Sec==score_last && identities_Sec==identities_last && expect_Sec==expect_last) # check we are in the same HSP
		 {

		     print query_ID_Sec, "("query_len" aa)";
		     print "\n";
		     print sbjct_ID_Sec, "("subj_len" nt)";
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
			print sbjct_ID_Sec, "("subj_len" nt)";
			print "\n";
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

# Check whether next alignment belongs to the Sec alignment HSP

	     if (query_ID_Sec==query_ID_next && sbjct_ID_Sec==sbjct_ID_next && score_next==score_Sec && identities_next==identities_Sec && expect_Sec==expect_next) # check we are in the same HSP
		 {

		  	print query;
			print middle_line;
			print sbjct_align;				
			print "\n";
			print "###################################################\n"
			bit=0;    
		 }

# Special case when an *-* alignment follows another *-* in different HSPs

	     else if ((query ~ PAT) && (score_>SCOREMIN) && (score_<SCOREMAX) && (identities>=IDENTITY) && (expect+0<EVALUE) &&((substr(sbjct_align,sc,1) ~ MATCH)))
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
			print sbjct_ID_Sec, "("subj_len" nt)";
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
	     
	 }    



}
