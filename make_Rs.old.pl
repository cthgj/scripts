#!/usr/bin/perl 
use Getopt::Std;
require "MY_SUBS.pl";
#################################################################
# This script will create plots of my moonlighting data.        #
# The order in which the arguments are passed is essential:     #
# ARGV[0] = Candidates						#
# ARGV[1] = Hubs						#
# ARGV[2] = Multi Non Candidates				#
# ARGV[3] = Multiclassed					#
# ARGV[4] = Non Candidates					#
# ARGV[5] = Monoclassed						#
# ARGV[6] = Network						#
# ARGV[7] = Image file name					#
#################################################################
my %opts;
getopts('hpbds:t:m:',\%opts);
my $gold_file="~/research/moonlight/data/gold.names";
if ($opts{h} or not $ARGV[0]) {
    &usage();
}
## Overwrite existing files
our $force=1;

my $dist=$opts{d}||undef;
my $TITLE=$opts{t}||undef;
my $phospho=$opts{p}||undef;
## Deal with titles containing spaces
$TITLE =~ s/_/ /g;
$TITLE =~ s/"//g;
$ARGV[0]=~/(.+?)\.(.+?)\./;

my $s = $1; ## Species
my $MODE=$2; ## Mode
my @modes;
if ($opts{s}) {
    $s=$opts{s}
}
#my $pat="";
if ($opts{m}) {
    $MODE=$opts{m};
    ## Run all at once
    if ($MODE =~/^(a.*?)/i){
#	$pat=$1;
	@modes=qw(i B iB);
    }
    else {
	@modes=split(/,/,$MODE);
    }
}
my $HUBS  = $ARGV[1] || "$species.hub";
my $MONO  = $ARGV[5] || "$species.mono";
my $MULTI = $ARGV[3] || "$species.multi";
my $NETW  = $ARGV[6] || "$species.netw";


###########################
# Capitalize first letter #
###########################
my $S=$s;
$S=~/^(.)/;
my $k=uc($1);
$S=~s/^(.)/$k/;

my %names =(
	    "i"=>"Intersection",
	    "B"=>"Bridge",
	    "iB"=>"Intersection & Bridge"
);
##########################################
# These lines are the same for all modes #
##########################################
unless ($phospho) {
	print<<EEEE;
titles<-c("i"="Intersection","B"="Bridge","iB"="Intersection & Bridge")
multi<-scan("$MULTI")
mono<-scan("$MONO")
hubs<-scan("$HUBS")
netw<-scan("$NETW")
gold<-scan("gold")
EEEE
}
foreach my $m (@modes) {

    my $CANDS  =  $ARGV[0] || "$species.$m.cands"; 
    my $MNC    =  $ARGV[2] || "$species.$m.mnon"; 
    my $NON    =  $ARGV[4] || "$species.$m.non"; 
    my $stats_file="$s.$m.stats";
    ####################################################
    # Get the right file names if we are in mode 'all' #
    ####################################################
    # my $cands=$CANDS;
    # my $mnc=$MNC;
    # my $non=$NON;
    # $mnc=~s/\.$pat\./\.$m\./;
    # $non=~s/\.$pat\./\.$m\./;
    ###############
    # Graph Title #
    ###############

    my $title="$S $TITLE - $names{$m}";
    ####################################
    # Collect values for gold standard #
    ####################################
    unless (-e "gold"){
	`/bin/grep -Fwf $gold_file $s.all | gawk '{print \$2}' > gold`;
    }
    ##################################################
    # Phosphorylation data are in a different format #
    ##################################################
    if ($phospho) {
	print_phospho($stats_file);
    } 
    else {
	print<<EEEE;
cands<-scan("$CANDS")
non<-scan("$NON")
multi_non<-scan("$MNC")
EEEE
	#####################
	# Plot distribution #
	#####################
	if ($dist) {
	    print_dist($title,$stats_file);
	}
	#################
	# Make Box Plot #
	#################
	else {
	    print_box($title,$stats_file);
	}
    }
}
#########################################################
###################### SUBROUTINES ######################
#########################################################
sub print_phospho{
    my $stats_file=shift;
   # print "MM : $MNC\n";
        	print<<EEEE;
titles<-c("i"="Intersection","B"="Bridge","iB"="Intersection & Bridge")
cands<-read.table("$CANDS",header=T)
non<-read.table("$NON",header=T)
multi<-read.table("$MULTI",header=T)
mono<-read.table("$MONO",header=T)
multi_non<-read.table("$MNC",header=T)
hubs<-read.table("$HUBS",header=T)
netw<-read.table("$NETW",header=T)
labels<-c("Network","Hubs","Candidates","Multi NC","Multi","Mono","NC")
gold<-read.table("gold",header=T)
png(filename="$s.$m.events.box.png", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
boxplot(netw\$Events,hubs\$Events,cands\$Events,multi_non\$Events,multi\$Events,mono\$Events,non\$Events, names=labels, outline=F, range=1.5,notch=T,varwidth=T)
Mean <- c(mean(netw\$Events),mean(hubs\$Events),mean(cands\$Events),mean(multi_non\$Events),mean(multi\$Events),mean(mono\$Events),mean(non\$Events))
points(Mean,col="red",pch=18)
pp1<-wilcox.test(cands\$Events,hubs\$Events)
pp2<-wilcox.test(cands\$Events,multi_non\$Events)
pp3<-wilcox.test(cands\$Events,netw\$Events)
pp4<-wilcox.test(cands\$Events,multi\$Events)
pp5<-wilcox.test(cands\$Events,mono\$Events)
pp6<-wilcox.test(cands\$Events,non\$Events)
for(i in 1:length(gold\$Events)){points(2.5,gold\$Events[i],col="blue",pch=18) }
mtext(paste("(",format.pval(pp1\$p.value, eps=0.005),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2\$p.value, eps=0.005),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3\$p.value, eps=0.005),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4\$p.value, eps=0.005),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5\$p.value, eps=0.005),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6\$p.value, eps=0.005),")"),side=1,line=2,at=7)
title(paste("$S Phosphorylation Events - ",titles["$m"]))
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
dev.off()
sink("$stats_file")
cat(paste(paste("Events Cands - Hubs : ", pp1\$p.value),paste("Events Cands - Multi NC : ", pp2\$p.value),paste("Events Cands - Network : ", pp3\$p.value),paste("Events Cands - Multi : ", pp4\$p.value),paste("Events Cands - Mono : ", pp5\$p.value),paste("Events Cands - Non Cands : ", pp6\$p.value), "------------------------------------------\\n",  sep = '\\n'))




b<-wilcox.test(cands\$Sites,multi_non\$Sites)
png(filename="$s.$m.sites.png", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
boxplot(netw\$Sites,hubs\$Sites,cands\$Sites,multi_non\$Sites,multi\$Sites,mono\$Sites,non\$Sites, names=labels, outline=F, range=1.5,notch=T,varwidth=T)
Mean <- c(mean(netw\$Sites),mean(hubs\$Sites),mean(cands\$Sites),mean(multi_non\$Sites),mean(multi\$Sites),mean(mono\$Sites),mean(non\$Sites))
points(Mean,col="red",pch=18)
pp1<-wilcox.test(cands\$Sites,hubs\$Sites)
pp2<-wilcox.test(cands\$Sites,multi_non\$Sites)
pp3<-wilcox.test(cands\$Sites,netw\$Sites)
pp4<-wilcox.test(cands\$Sites,multi\$Sites)
pp5<-wilcox.test(cands\$Sites,mono\$Sites)
pp6<-wilcox.test(cands\$Sites,non\$Sites)
for(i in 1:length(gold\$Sites)){points(2.5,gold\$Sites[i],col="blue",pch=18) }
mtext(paste("(",format.pval(pp1\$p.value, eps=0.005),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2\$p.value, eps=0.005),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3\$p.value, eps=0.005),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4\$p.value, eps=0.005),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5\$p.value, eps=0.005),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6\$p.value, eps=0.005),")"),side=1,line=2,at=7)
title(paste("$S Phosphorylation Sites - ",titles["$m"]))
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
dev.off()
cat(paste(paste("Sites Cands - Hubs : ", pp1\$p.value),paste("Sites Cands - Multi NC : ", pp2\$p.value),paste("Sites Cands - Network : ", pp3\$p.value),paste("Sites Cands - Multi : ", pp4\$p.value),paste("Sites Cands - Mono : ", pp5\$p.value),paste("Sites Cands - Non Cands : ", pp6\$p.value), "------------------------------------------\\n",  sep = '\\n'))

c<-wilcox.test(cands\$Kinases,multi_non\$Kinases)
png(filename="$s.$m.kinases.png", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
boxplot(netw\$Kinases,hubs\$Kinases,cands\$Kinases,multi_non\$Kinases,multi\$Kinases,mono\$Kinases,non\$Kinases, names=labels, outline=F, range=1.5,notch=T,varwidth=T)
Mean <- c(mean(netw\$Kinases),mean(hubs\$Kinases),mean(cands\$Kinases),mean(multi_non\$Kinases),mean(multi\$Kinases),mean(mono\$Kinases),mean(non\$Kinases))
points(Mean,col="red",pch=18)
pp1<-wilcox.test(cands\$Kinases,hubs\$Kinases)
pp2<-wilcox.test(cands\$Kinases,multi_non\$Kinases)
pp3<-wilcox.test(cands\$Kinases,netw\$Kinases)
pp4<-wilcox.test(cands\$Kinases,multi\$Kinases)
pp5<-wilcox.test(cands\$Kinases,mono\$Kinases)
pp6<-wilcox.test(cands\$Kinases,non\$Kinases)
for(i in 1:length(gold\$Kinases)){points(2.5,gold\$Kinases[i],col="blue",pch=18) }
mtext(paste("(",format.pval(pp1\$p.value, eps=0.005),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2\$p.value, eps=0.005),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3\$p.value, eps=0.005),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4\$p.value, eps=0.005),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5\$p.value, eps=0.005),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6\$p.value, eps=0.005),")"),side=1,line=2,at=7)
title(paste("$S Phosphorylation Kinases - ",titles["$m"]))
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
dev.off()
cat(paste(paste("Kinases Cands - Hubs : ", pp1\$p.value),paste("Kinases Cands - Multi NC : ", pp2\$p.value),paste("Kinases Cands - Network : ", pp3\$p.value),paste("Kinases Cands - Multi : ", pp4\$p.value),paste("Kinases Cands - Mono : ", pp5\$p.value),paste("Kinases Cands - Non Cands : ", pp6\$p.value), "------------------------------------------\\n",  sep = '\\n'))
sink()
EEEE
}

sub print_dist{
    my $title=shift;
    my $stats_file=shift;
   print<<EEEE;
df<-data.frame(grp=factor(c(rep("7. Mono",each=length(mono)),rep("6. NonCands",each=length(non)),rep("5. Multi",each=length(multi)),rep("4. Network",each=length(netw)),rep("3. Multi NonCands",each=length(multi_non)),rep("2. Hubs",each=length(hubs)),rep("1. Candidates",each=length(cands)))), val=c(mono,non,multi,netw,multi_non,hubs,cands))
library(ggplot2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "magenta", "#0072B2", "#D55E00", "#CC79A7")
ggplot(df, aes(x=val,colour=grp)) + geom_density()  + scale_colour_manual(values=cbbPalette) + labs(title="$title", x="$TITLE", y="Density")
ggsave("$IMG", width=8, height=8,dpi=300)
EEEE
}



sub print_box{
    my $title=shift;
    my $stats_file=shift;
    print<<EEEE;
png(filename="$IMG", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
pp1<-wilcox.test(cands,hubs)
pp2<-wilcox.test(cands,multi_non)
pp3<-wilcox.test(cands,netw)
pp4<-wilcox.test(cands,multi)
pp5<-wilcox.test(cands,mono)
pp6<-wilcox.test(cands,non)
labels<-c("Network","Hubs","Candidates","Multi NC","Multi","Mono","NC")
boxplot(netw,hubs,cands,multi_non,multi,mono,non, names=labels, outline=F, range=1.5,notch=TRUE, varwidth=TRUE)
title("$title")
Mean <- round(c(mean(netw), mean(hubs),mean(cands),mean(multi_non),mean(multi),mean(mono),mean(non)),digits=2)
points(Mean,col="red",pch=18)
for(i in 1:length(gold)){points(2.5,gold[i],col="blue",pch=18) }
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
mtext(paste("(",format.pval(pp1\$p.value, eps=0.005),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2\$p.value, eps=0.005),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3\$p.value, eps=0.005),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4\$p.value, eps=0.005),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5\$p.value, eps=0.005),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6\$p.value, eps=0.005),")"),side=1,line=2,at=7)
dev.off()
sink("$stats_file")
cat(paste(paste("Cands - Hubs : ", pp1\$p.value),paste("Cands - Multi NC : ", pp2\$p.value),paste("Cands - Network : ", pp3\$p.value),paste("Cands - Multi : ", pp4\$p.value),paste("Cands - Mono : ", pp5\$p.value),paste("Cands - Non Cands : ", pp6\$p.value),  sep = '\\n'))
sink()

EEEE

}


sub usage{
        my $us="[options] <>";
	my $desc="";
	my %opts=(
		  "usage" => $us,
		  "desc" => $desc,
		  "p"=> "",
		  "b"=> "Make boxplots.",
		  "d"=> "Make distribution plots (density function).",
		  "s"=> "Species. Def: Parsed from \$ARGV[0]",
		  "t"=> "Graph title.",
		  "m"=> "Mode, can be one of i,B,iB or 'all'"
		 );
    print_help_exit(\%opts,0);

}
################################# R NOTES  ################################
# points(0.5,30,col="blue",pch=18)
# par(mar=c(4,3,3,10), xpd=TRUE)
# leg.txt<-c("Mean","Known MPs")
# legend("topright",leg.txt,pch =18, col=c("red","blue"))
# boxplot(cands,hubs,multi_non,netw,multi,non,mono, names=labels, outline=FALSE, varwidth=TRUE) 
# legend("topright",inset=c(-0.2,0),leg.txt,pch =18, col=c("red","blue"))
# text(1.5,Mean[1],paste("P=",round(pp1\$p.value, digits=5)))
# text(2.5,Mean[2],paste("P=",round(pp2\$p.value, digits=5)))
# text(2.5,Mean[3],paste("P=",round(pp3\$p.value, digits=5)))
#par(mar=c(4,3,3,10), xpd=TRUE)

