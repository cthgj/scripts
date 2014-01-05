#!/usr/bin/perl
 #  $gopair e.g. : GO:0000001_GO:0050876
use Getopt::Std;
use Math::Round;
use DBI();
use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser); 
use Time::Piece;
use 5.10.0;
use Data::Dumper;


####################################
# Generate unique output file name #
####################################
my $outfile="../tmp/$$" . "_".Time::Piece->new->strftime('%d%m%Y') . ".html";

#####################
# Various variables #
#####################
my (%prec,%prob);
my (@GOs);
my (%all_pairs, %goterm_onto);
my $have_precisions=0;
my %seen;
my $root='/var/www_8080/PoGo';
my $data_dir=$root .'/data';
my $stats_dir=$data_dir . '/gostats/human/annotations';
my $goIDs=$data_dir .'/GO.terms_alt_ids';
my $single_file=0;
my $precfile="./precision.all";


####################################
# Uncomment to run as a cgi script #
####################################
my $cgi = new CGI;
my %input;
# for my $key ( $cgi->param() ) {
#     $input{$key} = $cgi->param($key);
# }

my $sim_dissim=param('sim_dissim')||undef;
my $chosengo1=param('ChosenGO1s')||undef;
my $chosengo2=param('ChosenGO2s')||undef;
my $golist1=param('golist1')||undef;
my $golist2=param('golist2')||undef;
my $gofile1=param('gofile1')||undef;
my $gofile2=param('gofile2')||undef;
my $precision=param('precision')||0;
my $species=param('species')||undef;
my $ptype=param('ptype')||'proteome';
my $min_threshold=param('min_threshold')||1;
my $max_threshold=param('max_threshold')||1;


#########################################
# Uncomment to run from the commandline #
#########################################
# my $find_dissimilar=$opts{D}||undef;
# my $find_similar=$opts{S}||undef;
# my $outfile=$opts{o}||undef;
# my $species=$opts{s}||undef;
# my $background=$opts{b}||undef;
# my $threshold=$opts{t}||1;
# my $chosengo1=$opts{c}||undef;
# my $chosengo2=$opts{C}||undef;
# my $golist1=$opts{g}||undef;
# my $golist2=$opts{G}||undef;
# my $gofile1=$opts{f}||undef;
# my $gofile2=$opts{F}||undef;
# my $prec=$opts{P}||0;



######################
# Database variables #
######################
my $host = "localhost";
my $database = "pogo";
my $user = "root";
my $pw = "yenapas";


##############################################
# Print the placeholder page to be displayed #
# while the query is running.		     #
##############################################
print_placeholder();

##############################################
# Collect the GO terms and/or GO pairs given #
##############################################
my ($a,$b,$c)=collect_gos;
my @gos1=@{$a};
my @gos2=@{$b};
my @pairs=@{$c};

###################################
# Die if no GOs have been found   #
###################################
if($#gos2==-1 && $#gos1==-1 && $#pairs==-1){
    &usage("No valid go terms detected");
}
##################################################
# If only one set of pairs has been passed,	 #
# get all the possible combinations between them #
##################################################
@gos2=@gos1 if $#gos2==-1;
@gos1=@gos2 if $#gos1==-1;

####################################################
# Filter out pairs if any of their GOs do not pass #
# the precision threshold. Also catch unknown GOs  #
# for which we have no precision info.             #
####################################################
($a,$b,$c)=filter_by_precision(0,@gos1);
@gos1=@{$a};
my @bad_prec=@{$b};
my @no_prec=@{$c};
## Same for go2s
($a,$b,$c)=filter_by_precision(0,@gos2);
@gos2=@{$a};
push @bad_prec, @{$b};
push @no_prec, @{$c};
## Same for pairs
($a,$b,$c)=filter_by_precision(1,@pairs);
@gos2=@{$a};
push @bad_prec, @{$b};
push @no_prec, @{$c};


########################################
# Now, get all possible pairs from all #
# input fields.			       #
########################################
$a=build_pairs(\@gos1,\@gos2,\@pairs);
%all_pairs=%{$a};

#############################
# Output of the MySQL query #
############################
my @o;


############################################
# Collect those GOs that have 0 count for  #
# the current species			   #
############################################
my %has_zero;
foreach (get_counts($host, $database, $user, $pw,$species,@gos1,@gos2,@pairs)){
    my @a=split(/\t/);
    $has_zero{$a[0]}= 1 if $a[1]==0;
}

#########################################
# If we are looking for (dis)similar GOs #
#########################################

if ($sim_dissim eq "dis"){
	@o=get_mysql_output('dis',\%pairs,$host, $database, $user, $pw,$species);
}
elsif ($sim_dissim eq "sim") {
    @o=get_mysql_output('sim',\%pairs,$host, $database, $user, $pw,$species);
}
#########################################
# If we want the probabilities for the  #
# submitted GOs			        #
#########################################
else {
    
    @o=get_mysql_output('pairs',\%all_pairs,$host, $database, $user, $pw,$species);
}


#######################################################
# Create a text file with the results for downloading #
#######################################################
my $text_file=make_text_file(\@o);

###############################
# Create the output html file #
###############################
html_out(\@o,$text_file);




















############################################################################
#                             SUBROUTINES                                  #
############################################################################

##############################################
# Since we can be given GOs or pairs and by  #
# various methods, the subroutine collects   # 
# them all. It return pointers to @go1, @go2 #
# and @pairs.                                #         
##############################################
sub collect_gos{
    my (%g1,%g2,%prs);

    #########################################
    # Get the GOs entered by autocompletion #
    #########################################
    if($chosengo1){$chosengo1=~s/\r//g; $g1{$_}++ for split(/\n/,$chosengo1);}
    if($chosengo2){$chosengo2=~s/\r//g; $g2{$_}++ for split(/\n/,$chosengo2);}

    ##########################################################
    # Get the GOs pasted in as lists. These could be either  #
    # GOs or GO pairs.					     #
    ##########################################################
    if($golist1){
	## For some reason apache adds \r
	$golist1=~s/\r//g;
	foreach my $gg (split(/\n/,$golist1)){
	    ## Is this a pair or a go?
	    if (/^\s*(GO:\d+)\s*$/) {$g1{$1}++}
	    elsif (/^\s*(GO:\d+[^\w]*GO:\d+)\s*$/) {$prs{$1}++}
	    else{usage( "Invalid GO term : $_\n");}
	}
    }
    if($golist2){
	$golist2=~s/\r//g;
	foreach my $gg (split(/\n/,$golist2)){
	    ## Is this a pair or a go?
	    if (/^\s*(GO:\d+)\s*$/) {$g2{$1}++}
	    elsif (/^\s*(GO:\d+[^\w]*GO:\d+)\s*$/) {$prs{$1}++}
	    else{usage( "Invalid GO term : $_\n");}
	}
    }

    #################################
    # Get the gos uploaded as files #
    #################################
    if($gofile1){	
	open(my $fh, '<',"$gofile1")||die("Could not open $gofile1 : $!\n");
	while(<$fh>){
	    chomp;
	    s/\r//g;
	    if (/^\s*(GO:\d+)\s*$/) {$g1{$1}++}
	    elsif (/^\s*(GO:\d+[^\w]*GO:\d+)\s*$/) {$prs{$1}++}
	    else{usage( "Invalid GO term : $_\n");}
	}
    }
my @gg1=keys(%g1)
my @gg2=keys(%g2)
my @pp=keys(%prs)
return(\@gg1,\@gg2,\@pp);

}
############################################################
############################################################

#################################################
# This function will return a list of all the   #
# GO terms from the input list thaht pass the   #
# precision threshold. Input list can be either # 
# pairs or GOs.  			        #
#################################################
sub filter_by_precision{
    my (@good,@bad,@unk);
    #####################################
    # If $mode is 1, @_ contains pairs. #
    # Otherwise, it contains GOs        #
    #####################################
    my $mode=shift;
    #############################
    # Load the precision values #
    #############################
    unless ($have_precisions==1) {
	open(my $fh,'<',$precfile) or die "Could not open prec $precfile for reading : $!\n";
	while (<$fh>) {
	    chomp;
	    my @a=split(/\s+/);
	    $prec{$a[0]}=$a[1];
	}
	$have_precisions = 1;
    }
    ###############################
    # If we have been given pairs #
    ###############################
    if ($mode==1) {
	foreach my $p (@_) {
	    $p=~/^\s*(GO:\d+)[^\w]*(GO:\d+)\s*$/;
	    if ($prec{$1} >= $precision && $prec{$2} >= $precision) {
		push @good,$p;
	    }
	    elsif ($prec{$1} < $precision || $prec{$2} < $precision) {
		push @bad,$p;
	    }
	    else {
		push @unk, $p;
		print STDERR "ERROR, unknown GO term : $_\n";
	    }
	}
    }
    ###############################
    # If we have been given GOs   #
    ###############################
    else {
	foreach (map {/(GO:\d+)/g} @_){ ## remove spaces etc
	    if ($prec{$_} >= $precision) {
		push @good,$_;
	    }
	    elsif ($prec{$_} < $precision) {
		push @bad,$_;
	    }
	    else {
		push @unk, $_;
		print STDERR "ERROR, unknown GO term : $_\n";
	    }
	}
    }
    return(\@good,\@bad,\@unk);

}

####################################################################
####################################################################

################################################
# Returns counts for a list of GOs or pairs in #
# the current ontology combination	       #
################################################
sub get_counts{
    my ($host, $database, $user, $pw, $table,@gos)=@_;
    my $num_gos=scalar(keys(%gos));
    my @output;
    ## Connect to server
    my $dsn = "DBI:mysql:database=$database;host=$host;";
    ##Select database
    my $dbh=DBI->connect($dsn,$user,$pw);
    my $query="SELECT name,count from " . $table . "_counts where name=";
    for(my $i=0; $i<$#gos; $i++){
	$query.="'$gos[$i]' OR name=";
    }
    ## Add the last go
    $query.="'$gos[$#gos]' ;";
    my $result=$dbh->prepare("$query");
    $result->execute;
    while (my $ref = $result->fetchrow_hashref()) {
	push @output, "$ref->{'name'}\t$ref->{'onto'}\t$ref->{'count'}";
    }
    $result->finish();
    return(@output);
}

####################################################################
####################################################################

######################################################
# Returns a pointer to a list of pairs. It's input   #
# should be 3 pointers to  lists of GOs and/or pairs #
######################################################
sub build_pairs{
    my %pps;
    my @G1=@{shift};
    my @G2=@{shift};
    my @prs=@{shift};
    foreach my $g1 (@G1) {
	foreach my $g2 (@G2) {
	    next if $g1 eq $g2;
	    my $pair=join("_", sort  {$b lt $a} ($g1,$g2));
	    $pps{$pair}++;
	}
    }
    foreach my $p (@prs) {
	$p=~/^\s*(GO:\d+)[^\w]*(GO:\d+)\s*$/;
	my $pair=join("_", sort  {$b lt $a} ($1,$2));
	$pps{$pair}++;
	
    }

 return(\%pps)
}


####################################################################
####################################################################

#Query the server for the probabilites. Input should be:
#1. Mode, one of 'sim','dis' or 'pairs'. The first two

sub get_mysql_output{
    my $mode=shift;
    print STDERR "MMMMMODE $mode";
    my %pairs=%{shift()};
    my $num_pairs=scalar(keys(%pairs));
    my ($host, $database, $user, $pw, $table)=@_;
    my @output;
    ## Connect to server
    my $dsn = "DBI:mysql:database=$database;host=$host;";
    ##Select database
    my $dbh=DBI->connect($dsn,$user,$pw);
    my @kk=keys(%pairs);
    my $query;
    if ($mode eq 'pairs') {
	$query="SELECT gopair,tot_onto,go1_num,go2_num,overlap,P_low,P_high from $table where gopair=";
	for(my $i=0; $i<$#kk; $i++){
	    $query.="'$kk[$i]' OR gopair=";
	}
	## Add the last pair
	$query.="'$kk[$#kk]' ORDER BY P_low DESC;";
    }
    else {    print STDERR "mm : $mode @kk\n";
	$query="SELECT gopair,tot_onto,go1_num,go2_num,overlap,P_low,P_high from $table where (";
	for(my $i=0; $i<=$#kk; $i++){
	    print STDERR "$i : $kk[$i]\n";
	    if ($mode eq 'sim') {
		$query.=" (go1='$kk[$i]' OR go2='$kk[$i]')";
	    }
	    elsif ($mode eq 'dis') {
		$query.=" (go1='$kk[$i]' OR go2='$kk[$i]')";
	    }
	}
	if ($mode eq 'sim') {
	    $query.=") AND P_high > $min_threshold ORDER BY P_high DESC LIMIT 20";
	}
	elsif ($mode eq 'dis') {
	    $query.=") AND P_low < $min_threshold ORDER BY P_high DESC LIMIT 20"
	}


    }
#    print STDERR "QQ (probs) : $query\n";
    my $result=$dbh->prepare("$query");
    $result->execute;
    my $c=0;
    while (my $ref = $result->fetchrow_hashref()) {
	$c++;
	push @output, "$ref->{'gopair'}\t$ref->{'tot_onto'}\t$ref->{'go1_num'}\t$ref->{'go2_num'}\t$ref->{'overlap'}\t$ref->{'P_high'}\t$ref->{'P_low'}";
    }
    $result->finish();
    return(@output);
}
####################################################################
####################################################################


###################################################################
# This will create a tab separated text file with all the results #
###################################################################
sub make_text_file{
    $outfile=~/(.+)\./;
    my $text_file=$1 . "export.csv";
    my @out=@{shift()};
    open(my $fh,'>',"$text_file")||die "Could not open export $text_file : $!\n";
    print $fh "1st GO term\t2nd GO term\tTotal\tGO1\tGO2\tBoth\tP-value (under)  \tP-value (over)\n" if $#out>-1;
    for (my $n=0; $n<=$#out;$n++) {
	#print STDERR "aa $n : $out[$n]\n";
	my @o=split(/\t+/,$out[$n]);
	#my $prob=$out[$n+1];
	my @pairs=split(/_/,$o[0]);
	my @kk=split(/\t/,$o[1]);
	#print STDERR "aa o[0] : $o[1] : $o[2] ::: $out[$n] <= $threhold\n";
	next unless $out[7]<=$max_threshold;
	my ($onto0,$desc0)=&goterm_ontos($pairs[0]);
	my ($onto1,$desc1)=&goterm_ontos($pairs[1]);
	print $fh "$pairs[0]\t$pairs[1]\t$o[1]\t$o[2]\t$o[3]\t$o[4]\t$o[5]\t$o[6]\n";
    }
    close($fh);
    return($text_file);
}
	
####################################################################
####################################################################

#################################################
# This prints the first lines of the html file. #
# These are common for all output files.        #
# It will print until an opening <td> tag,      #
# print_html_footer() closes the table.         #
#################################################
sub print_html_header{
    my $fh=shift;
    print $fh <<EOF;
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
   "http://www.w3.org/TR/html4/loose.dtd">
<html>
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <title>PrOnto</title>
    <link rel="stylesheet" type="text/css" href="../webdev.css">
    <link rel="stylesheet" type="text/css" href="../jquery-ui-1.8.custom.cec.css">
<script type="text/javascript" src="../code/js/jquery-1.4.2.min.js"></script> 
<script type="text/javascript" src="../jquery.js"></script>
<script type="text/javascript" src="../pogo.js"></script>
<script src="../sorttable-small.js"></script>
</head>
  <body>
<br>
<div id="main-block">
  <div id="dhtmltooltipGen" class="dhtmltooltip"></div>
  <div id="dhtmltooltipBP" class="dhtmltooltip"></div>
  <div id="dhtmltooltipCC" class="dhtmltooltip"></div>
  <div id="dhtmltooltipMF" class="dhtmltooltip"></div>
<script type="text/javascript" src="../tooltip.js"></script>

<!-- added recently -->
<body>
<div id="main-block">
<table id="formTable" border="0">
  <tr class="topHeader"><td style="border-width:0px;">
<div id="pogo" style="display:inline-block; "> <img src="../images/pronto.png" alt="PrOnto logo"></div>
<div id="logos_div"><a href="http://tagc.univ-mrs.fr/welcome/"><img  id="tagc" alt="tagc logo" src="../tagc_logo.png" width="150"></a>
    <a href="http://www.univmed.fr/"><img id="univmed" alt="Université de la Mediteranée" width="140" src="../univmed.png"></a></div>
  </td></tr>
  <tr><td colspan="2" style="height:8px; border:0px solid #597e7a; border-width: 0px 0px 0px 0px;"></td></tr>
  <tr id="globalnav"><td colspan="2" class="globalnav">
    <a class="here" id="here" href="#">GO pair association Probabilities</a>
    <a class="there" href="similar.html" onMouseOver="document.getElementById('here').setAttribute('class', 'there');" onMouseOut="document.getElementById('here').setAttribute('class', 'here');" >Find (dis)similar terms</a>
  </td>
</tr>
<tr><td class="topborder">
  <table style="" class="inputTable first "><tr>
    <tr>  <td class=""  style="padding-right:20px;"></td>    
   <tr><td class="" style="background:white;"> 

EOF

};
####################################################################
####################################################################

####################################
# This prints the table of results #
####################################
sub html_out{
    my @out=@{shift()};
    my $text_file=shift;
    my $rows=$#out+2;
    open(my $fh,'>',"$outfile")||die "Could not open outfile $outfile for writing: $!\n";
    print_html_header($fh);
    ##########################################################
    # Only print the table header if results have been found #
    ##########################################################
    if ($#out>-1) {
	print $fh <<EOF;

<table id="resultsTabWrap"  border="0">
  <tr class="topHeader"><td style="border-width:0px;">
<div id="pogo" style="display:inline-block; "> <img src="../images/pronto.png" alt="PrOnto logo" height="65" ></div>
<div id="logos_div"><a href="http://tagc.univ-mrs.fr/welcome/"><img  id="tagc" alt="tagc logo" src="../tagc_logo.png" width="150"></a>
    <a href="http://www.univmed.fr/"><img id="univmed" alt="Université de la Mediteranée" width="140" src="../univmed.png"></a></div>
  </td></tr>
  <tr><td colspan="2" style="height:8px; border:0px solid #597e7a; border-width: 0px 0px 0px 0px;"></td></tr>
  <tr><td colspan="2">
 <table class="HeaderTable">
<tr  ><td  colspan="7" style="vertical-align:middle">Click on any header to sort the table by that field.</td><td style="vertical-align:middle" class="help">Help: On<input CHECKED name="help" type="radio" Onclick="he=1" style="vertical-align:middle"> Off<input name="help" type="radio" Onclick="he=0" ></td>
</tr></table>
    <table class="resultsTab sortable">
  
  <tr id="resultsTabHead"  >
    <td  onMouseover="ddrivetip('Click on a GO term to see its entry in QuickGO','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">1st GO term</td>

<td onMouseover="ddrivetip('Click on a GO term below to see its entry in QuickGO','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()" >2nd GO term</td>

<td onMouseover="ddrivetip('The number of $species proteins that have at least one direct annotation in the ontology of the first term and on in the ontology of the second. If both terms are from the same ontology, the number of $species proteins with at least two annotations in that ontology.','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">Total</td> 

<td onMouseover="ddrivetip('The number of $species proteins that are annonotated to the first GO term.','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">GO1</td> 

<td onMouseover="ddrivetip('The number of $species proteins that are annonotated to the second GO term.','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">GO2</td> 

<td onMouseover="ddrivetip('The number of $species proteins that are annotated to both terms','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">both</td>

<td onMouseover="ddrivetip('The probability of finding as many or fewer proteins annotated to both of these terms by chance.','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">P-value<br>(under)</td>

<td onMouseover="ddrivetip('The probability of finding more proteins annotated to both of these terms by chance.','', '','dhtmltooltipGen')" onMouseout="hideddrivetip()">P-value<br>(over)</td>
EOF
    };
    else {
 	print $fh <<EOF;
   No results were found
EOF
    };

    #####################
    # Print the results #
    #####################
    for (my $n=0; $n<=$#out;$n++){
	my @o=split(/\t+/,$out[$n]);
	my @pairs=split(/_/,$o[0]);
	my @kk=split(/\t/,$o[1]);
	next unless $out[7]<=$max_threshold;
	####################################################
        # Get the ontology and text description of 	   #
	# each term.					   #
        ####################################################
	my ($onto0,$desc0)=&goterm_ontos($pairs[0]);
	my ($onto1,$desc1)=&goterm_ontos($pairs[1]);
		print $fh <<EOF;
 <tr><td class="result$onto0"><span  onMouseover="ddrivetip('<b>$pairs[0]</b>\t$desc0\t$onto0','', '','dhtmltooltip$onto0')" onMouseout="hideddrivetip()"  class="result$onto0"><a class="results" href="http://www.ebi.ac.uk/QuickGO/GTerm?id=$pairs[0]" target="_blank">$pairs[0]</a></span></td>
<td  class="result$onto1"><span onMouseover="ddrivetip('<b>$pairs[1]</b>\t$desc1\t$onto1','', '','dhtmltooltip$onto1')" onMouseout="hideddrivetip()"class="result$onto1"><a class="results" href="http://www.ebi.ac.uk/QuickGO/GTerm?id=$pairs[1]"  target="_blank">$pairs[1]</a></span></td>
<td>$o[1]</td> <!-- Total -->
<td>$o[2]</td> <!--go1 -->
<td>$o[3]</td> <!-- go2-->
<td>$o[4]</td> <!--both -->
<td>$o[5]</td> <!--under -->
<td>$o[6]</td> <!--over -->
</tr>
EOF
    }
    print $fh "</table></td></tr></table>";
}
####################################################################
####################################################################

###########################################################
# This will close the table opened by print_html_header() #
# and print the links to the missing and download options #
# as well as the <body> and <html> tags. Expects an	  #
# open filehandle as $_[0]                                #
###########################################################   
 
sub print_html_footer{
    my $fh=shift;
    print $fh "</table></td></tr></table>";
    #############################################################
    # Check if anything needs to be printed to the missing file #
    # and if so, open the filehandle			        #
    #############################################################
    my $mis_fh;
    if ($no_prec[0] || $zeros[0]) {
	$outfile=~/(.+)\./;
	my $missing_file=$1 . "missing.html";
	open($mis_fh,'>',"$missing_file")||die("Could not open missing $missing_file for writing: $!\n");
    }

    #####################################################
    # If any GOs are unknown (no precision), create an  #
    # "unknown" file to list them.		        #
    #####################################################
    if ($no_prec[0]) {
	print_missing($fh,$mis_fh,"are not present in the annotations of $species and have been ignored","The following GO terms do not annotate any $species gene products:");
   }
    #################################################################
    # If any GOs have 0 count, create a "missing" file to list them #
    #################################################################
    my @zeros=keys(%has_zero);
    if ($zeros[0]) {
	print_missing($fh,$mis_fh,"are not present in the annotations of $species and have been ignored","The following GO terms do not annotate any $species gene products:");
    }
    print $fh "</div></body></html>";
    close($fh);
    close($mis_fh);
}
############################################################
############################################################


###############################################
# This will create the file with the missing  #
# GOs and will print a link to it in the file #
# handle passed.			      #
###############################################
sub print_missing{
    my $fh=shift;
    my $mis_fh=shift;
    my $msg1=shift;
    my $msg2=shift;
    $outfile=~/(.+)\./;
    my $missing_file=$1 . "missing.html";
    print $fh "<p style=\"margin-left:10px\"><a href=\"$missing_file\">", scalar(@zeros), " GO terms</a> $msg.</p>";
    open(my $mis,'>',"$missing_file")||die("Could not open missing $missing_file for writing: $!\n");
    print_html_header($mis);
    print $mis <<EOF;
<p class="missing">$msg2<br><br> 
EOF
    ########################
    # Sort GOs by ontology #
    ########################
    foreach my $go (sort_onto(@zeros)){
	my ($onto0,$desc0)=goterm_ontos($go);
	print $mis <<EOF;
<span  onMouseover="ddrivetip('<b>$go</b>\t$desc0\t$onto0','', '','dhtmltooltip$onto0')" onMouseout="hideddrivetip()"  class="result$onto0"><a class="results" href="http://www.ebi.ac.uk/QuickGO/GTerm?id=$go" target="_blank">$go</a></span>
EOF
    }
    print $mis "</tr></table></td></tr></table></body></html>";
    close($mis);
}

####################################################################
####################################################################

########################################
# Sort a list of GOs by their ontology #
########################################
sub sort_onto{
    my %sorted;
    foreach (@_) {
	my ($o,$foo)=goterm_ontos($_);
	$sorted{$o}{$_}++;
	print STDERR "\$sorted{$o}{$_}++;";
    }
    my @ret;
    foreach my $o (keys(%sorted)) {
	foreach my $go (keys(%{$sorted{$o}})) {
	    push @ret,$go;
	}
    }
return(@ret);
}

############################################################
############################################################
sub usage{
    my $msg=shift;
    print STDOUT <<EOF;
 <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
   "http://www.w3.org/TR/html4/loose.dtd">
<html>
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <title>PrOnto</title>
    <link rel="stylesheet" type="text/css" href="../pogo.min.css">
    <link rel="stylesheet" type="text/css" href="../jquery-ui-1.8.custom.cec.css">
<script type="text/javascript" src="../code/js/jquery-1.4.2.min.js"></script> 
<script type="text/javascript" src="../jquery.js"></script>
<script type="text/javascript" src="../pogo.js"></script>
</head>
  <body>
    <div id="dhtmltooltipBP" class="dhtmltooltip"></div>
    <div id="dhtmltooltipCC" class="dhtmltooltip"></div>
    <div id="dhtmltooltipMF" class="dhtmltooltip"></div>
<script type="text/javascript" src="../tooltip.js"></script>

 <div id="topHeader" >
      <div id="toplegendDIV">
	<img  src="../images/pronto.png" alt="">
      </div>
      </div>
 <table class="resultsTabWrap" ><tr><td>
    <table class="resultsTab" >
    <tr id="resultsTabHead"><td>$msg</td>

    </tr>

  </table></td><td rowspan="2" class="datadiv"></td></tr>
<tr><td class="logos" >

	<a href="http://tagc.univ-mrs.fr/welcome/"><img  alt="tagc logo" src="../tagc_tr.png" width="150"></a><span style="text-align:right"></td>
<td align="right" class="logos" >	<a href="http://www.univmed.fr/"><img alt="Université de la Mediteranée" width="170" height="40" src="../univmed.png"></a></span>

    </td></tr>
</table>
</div>


</body></html>
EOF

   

print STDERR "done done\n";
exit();
}

