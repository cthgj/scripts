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
getopts('ShpaPbds:t:m:i:n:f:',\%opts);
if ($opts{h}) {
    &usage();
}
## Overwrite existing files
our $force=1;

my $dist=$opts{d}||undef;
my $box=$opts{b}||undef;
my $pie=$opts{P}||undef;
unless (defined($dist) || defined($box)) {
    $box=1;
}
my $TITLE=$opts{t}||undef;
my $phospho=$opts{p}||undef;

## Deal with titles containing spaces
$TITLE =~ s/_/ /g;
$TITLE =~ s/"//g;


my $species = $1; ## Species
my @modes;
if ($opts{s}) {
    $species=$opts{s}
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
else {
    @modes=qw(i B iB);
}


my $HUBS  = $ARGV[1] || "$species.hubs.nums";
my $MONO  = $ARGV[5] || "$species.mono.nums";
my $MULTI = $ARGV[3] || "$species.multi.nums";
my $NETW  = $opts{n} || "$species.all";

###########################
# Capitalize first letter #
###########################
my $S=$species;
$S=~/^(.)/;
my $k=uc($1);
$S=~s/^(.)/$k/;

my %names =(
	    "i"=>"Intersection",
	    "B"=>"Bridge",
	    "iB"=>"Intersection & Bridge"
);

##################
# Make Pie Chart #
##################
if ($pie) {
    foreach my $m (@modes) {
    my $IMG="$species.$m.pie.png";
    if ($opts{i}) {
	$IMG="$species.$m.$opts{i}.pie.png";
    }
    my $file="$species.$m.nums";
    $opts{f} && do {$file="$species.$m."};
    print STDERR "FF $file\n";
    print<<EndOf;
titles<-c("i"="Intersection","B"="Bridge","iB"="Intersection & Bridge")
cands<-scan("$file")
pct <- round(cands/sum(cands)*100)
lbls<-c("Neither","One", "Both")
lbls <- paste(lbls, pct)
lbls <- paste(lbls,"%",sep="")
lbls <- paste(lbls,"(", sep=" ")
lbls <- paste(lbls,cands, sep="")
lbls <- paste(lbls,")", sep="")
cbbPalette <- c("#56B4E9", "#0f4998", "#0072B2", "lightslategrey", "#CC79A7")
png(filename="$IMG", units = "px", pointsize = 17, bg = "white", width = 800, height = 800)
pie(cands,labels = lbls,col=cbbPalette)
title(paste("$TITLE\\n",sum(cands), " GO pairs", sep=""))
a<-c(paste("$S GO pairs - ",titles["$m"], " (",sum(cands), ")", sep=""),lbls)
write(a,"$species.$m.perc")


EndOf
print STDERR "made $IMG\n";
 }
}

########################
# All except pie chart #
########################
else {
    ##########################################
    # These lines are the same for all modes #
    ##########################################
    unless ($phospho) {
	print<<EEEE;
titles<-c("i"="Intersection","B"="Bridge","iB"="Intersection & Bridge")
multi<-read.table("$MULTI",sep="\\t",col.names=c("Name","val"))
mono<-read.table("$MONO",sep="\\t",col.names=c("Name","val"))
hubs<-read.table("$HUBS",sep="\\t",col.names=c("Name","val"))
netw<-read.table("$NETW",sep="\\t",col.names=c("Name","val"))
EEEE

	if ($species eq 'human' || $species eq 'vidal') {
	    print "gold<-read.table(\"gold\",sep=\"\\t\",col.names=c(\"Name\",\"val\"))\n";
	}
    }
    foreach my $m (@modes) {

	my $CANDS  =  $ARGV[0] || "$species.$m.cands.nums"; 
	my $MNC    =  $ARGV[2] || "$species.$m.mnon.nums"; 
	my $NON    =  $ARGV[4] || "$species.$m.non.nums"; 
	my $stats_file="$species.$m.stats";
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
	my $title;
	if ($opts{a}) {
	    $title=$TITLE;
	}
	else {
	    $title="$S $TITLE - $names{$m}";
	}

	##################################################
	# Phosphorylation data are in a different format #
	##################################################
	if ($phospho) {
	    print_phospho($stats_file);
	} else {
	    print<<EEEE;
cands<-read.table("$CANDS",sep="\\t",col.names=c("Name","val"))
non<-read.table("$NON",sep="\\t",col.names=c("Name","val"))
mnon<-read.table("$MNC",sep="\\t",col.names=c("Name","val"))
EEEE
	    #####################
	    # Plot distribution #
	    #####################
	    if ($dist) {
		print_dist($title,$stats_file,$m);
	    }
	    #################
	    # Make Box Plot #
	    #################
	    if ($box) {
		print_box($title,$stats_file,$m);
	    }
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
mnon<-read.table("$MNC",header=T)
hubs<-read.table("$HUBS",header=T)
netw<-read.table("$NETW",header=T)
labels<-c("Network","Hubs","Candidates","Multi NC","Multi","Mono","NC")
EEEE
    if ($species eq 'human' || $species eq 'vidal') {
	print "gold<-read.table(\"gold\",header=T)\n";
    }
    print<<EEEE;
png(filename="$species.$m.events.box.png", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
boxplot(netw\$Events,hubs\$Events,cands\$Events,mnon\$Events,multi\$Events,mono\$Events,non\$Events, names=labels, outline=F, range=1.5,notch=T,varwidth=T)
Mean <- c(mean(netw\$Events),mean(hubs\$Events),mean(cands\$Events),mean(mnon\$Events),mean(multi\$Events),mean(mono\$Events),mean(non\$Events))
points(Mean,col="red",pch=18)
pp1<-signif(wilcox.test(cands\$Events,hubs\$Events)\$p.value,digits=3)
pp2<-signif(wilcox.test(cands\$Events,mnon\$Events)\$p.value,digits=3)
pp3<-signif(wilcox.test(cands\$Events,netw\$Events)\$p.value,digits=3)
pp4<-signif(wilcox.test(cands\$Events,multi\$Events)\$p.value,digits=3)
pp5<-signif(wilcox.test(cands\$Events,mono\$Events)\$p.value,digits=3)
pp6<-signif(wilcox.test(cands\$Events,non\$Events)\$p.value,digits=3)
pp7<-signif(wilcox.test(hubs\$Events,mnon\$Events)\$p.value,digits=3)
EEEE
if ($species eq 'human' || $species eq 'vidal') {
    print "for(i in 1:length(gold\$Events)){points(2.5,gold\$Events[i],col=\"blue\",pch=18) }\n"
}
print<<EEEE;
mtext(paste("(",format.pval(pp1, eps=0.05),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2, eps=0.05),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3, eps=0.05),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4, eps=0.05),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5, eps=0.05),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6, eps=0.05),")"),side=1,line=2,at=7)
title(paste("$S Phosphorylation Events - ",titles["$m"]))
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
dev.off()
sink("$stats_file")
cat(paste(paste("Events Cands - Hubs : ", pp1),paste("Events Cands - Multi NC : ", pp2),paste("Events Cands - Network : ", pp3),paste("Events Cands - Multi : ", pp4),paste("Events Cands - Mono : ", pp5),paste("Events Cands - Non Cands : ", pp6),paste("Events Hubs - Multi NC : ", pp7)), "------------------------------------------\\n",  sep = '\\n')




b<-wilcox.test(cands\$Sites,mnon\$Sites)
png(filename="$species.$m.sites.png", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
boxplot(netw\$Sites,hubs\$Sites,cands\$Sites,mnon\$Sites,multi\$Sites,mono\$Sites,non\$Sites, names=labels, outline=F, range=1.5,notch=T,varwidth=T)
Mean <- c(mean(netw\$Sites),mean(hubs\$Sites),mean(cands\$Sites),mean(mnon\$Sites),mean(multi\$Sites),mean(mono\$Sites),mean(non\$Sites))
points(Mean,col="red",pch=18)
pp1<-signif(wilcox.test(cands\$Sites,hubs\$Sites)\$p.value,digits=3)
pp2<-signif(wilcox.test(cands\$Sites,mnon\$Sites)\$p.value,digits=3)
pp3<-signif(wilcox.test(cands\$Sites,netw\$Sites)\$p.value,digits=3)
pp4<-signif(wilcox.test(cands\$Sites,multi\$Sites)\$p.value,digits=3)
pp5<-signif(wilcox.test(cands\$Sites,mono\$Sites)\$p.value,digits=3)
pp6<-signif(wilcox.test(cands\$Sites,non\$Sites)\$p.value,digits=3)
pp7<-signif(wilcox.test(hubs\$Sites,mnon\$Sites)\$p.value,digits=3)
EEEE
if ($species eq 'human' || $species eq 'vidal') {
    print "for(i in 1:length(gold\$Sites)){points(2.5,gold\$Sites[i],col=\"blue\",pch=18) }\n";
}
print<<EEEE;
mtext(paste("(",format.pval(pp1, eps=0.05),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2, eps=0.05),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3, eps=0.05),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4, eps=0.05),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5, eps=0.05),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6, eps=0.05),")"),side=1,line=2,at=7)
title(paste("$S Phosphorylation Sites - ",titles["$m"]))
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
dev.off()
cat(paste(paste("Sites Cands - Hubs : ", pp1),paste("Sites Cands - Multi NC : ", pp2),paste("Sites Cands - Network : ", pp3),paste("Sites Cands - Multi : ", pp4),paste("Sites Cands - Mono : ", pp5),paste("Sites Cands - Non Cands : ", pp6), ,paste("Events Hubs - Multi NC : ", pp7), "------------------------------------------\\n",  sep = '\\n'))

c<-wilcox.test(cands\$Kinases,mnon\$Kinases)
png(filename="$species.$m.kinases.png", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
boxplot(netw\$Kinases,hubs\$Kinases,cands\$Kinases,mnon\$Kinases,multi\$Kinases,mono\$Kinases,non\$Kinases, names=labels, outline=F, range=1.5,notch=T,varwidth=T)
Mean <- c(mean(netw\$Kinases),mean(hubs\$Kinases),mean(cands\$Kinases),mean(mnon\$Kinases),mean(multi\$Kinases),mean(mono\$Kinases),mean(non\$Kinases))
points(Mean,col="red",pch=18)
pp1<-signif(wilcox.test(cands\$Kinases,hubs\$Kinases)\$p.value,digits=3)
pp2<-signif(wilcox.test(cands\$Kinases,mnon\$Kinases)\$p.value,digits=3)
pp3<-signif(wilcox.test(cands\$Kinases,netw\$Kinases)\$p.value,digits=3)
pp4<-signif(wilcox.test(cands\$Kinases,multi\$Kinases)\$p.value,digits=3)
pp5<-signif(wilcox.test(cands\$Kinases,mono\$Kinases)\$p.value,digits=3)
pp6<-signif(wilcox.test(cands\$Kinases,non\$Kinases)\$p.value,digits=3)
pp7<-signif(wilcox.test(hubs\$Kinases,mnon\$Kinases)\$p.value,digits=3)
EEEE
if ($species eq 'human' || $species eq 'vidal') {
    print "for(i in 1:length(gold\$Kinases)){points(2.5,gold\$Kinases[i],col=\"blue\",pch=18) }\n"
}
print<<EEEE;
mtext(paste("(",format.pval(pp1, eps=0.05),")"),side=1,line=2,at=2)
mtext(paste("(",format.pval(pp2, eps=0.05),")"),side=1,line=2,at=4)
mtext(paste("(",format.pval(pp3, eps=0.05),")"),side=1,line=2,at=1)
mtext(paste("(",format.pval(pp4, eps=0.05),")"),side=1,line=2,at=5)
mtext(paste("(",format.pval(pp5, eps=0.05),")"),side=1,line=2,at=6)
mtext(paste("(",format.pval(pp6, eps=0.05),")"),side=1,line=2,at=7)
title(paste("$S Phosphorylation Kinases - ",titles["$m"]))
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =18, col=c("red","blue"))
dev.off()
cat(paste(paste("Kinases Cands - Hubs : ", pp1),paste("Kinases Cands - Multi NC : ", pp2),paste("Kinases Cands - Network : ", pp3),paste("Kinases Cands - Multi : ", pp4),paste("Kinases Cands - Mono : ", pp5),paste("Kinases Cands - Non Cands : ", pp6), paste("Events Hubs - Multi NC : ", pp7), "------------------------------------------\\n",  sep = '\\n'))
sink()
EEEE
}

sub print_dist{
    my $title=shift;
    my $stats_file=shift;
    my $mode=shift;
    my $IMG="$species.$mode.dist.png";
    if ($opts{i}) {
	$IMG="$species.$mode.$opts{i}.dist.png";
    }
    
    
   print<<EEEE;
df<-data.frame(grp=factor(c(rep("7. Mono",each=length(mono\$val)),rep("6. NonCands",each=length(non\$val)),rep("5. Multi",each=length(multi\$val)),rep("4. Network",each=length(netw\$val)),rep("3. Multi NonCands",each=length(mnon\$val)),rep("2. Hubs",each=length(hubs\$val)),rep("1. Candidates",each=length(cands\$val)))), val=c(mono\$val,non\$val,multi\$val,netw\$val,mnon\$val,hubs\$val,cands\$val))
library(ggplot2)
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "magenta4", "#0072B2", "lightslategrey", "#CC79A7")
ggplot(df, aes(x=val,colour=grp)) + geom_density()  + scale_colour_manual(values=cbbPalette) + labs(title="$title", x="$TITLE", y="Density")
ggsave("$IMG", width=8, height=8,dpi=300)
EEEE
}

## cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ff00ff", "#0072B2", "#D55E00", "#CC79A7")



sub print_box{
    my $title=shift;
    my $stats_file=shift;
    my $mode=shift;
    my $IMG=$opts{i}||"$species.$mode.box.png";
    ## Spacers for -S
    my ($sp1,$sp2,$sp3)=("   ","    ","  ");
#    $opts{S} && do {($sp1,$sp2,$sp3)=("  ","    ","  ")};
    my $sp4="    ";

 if ($opts{i}) {
	$IMG="$species.$mode.$opts{i}.box.png";
    }
    print<<EEEE;
png(filename="$IMG", width = 800, height = 800, units = "px", pointsize = 17, bg = "white")
pp1<-signif(wilcox.test(cands\$val,netw\$val)\$p.value,digits=3)
pp2<-signif(wilcox.test(cands\$val,hubs\$val)\$p.value,digits=3)
pp3<-signif(wilcox.test(cands\$val,mnon\$val)\$p.value,digits=3)
pp4<-signif(wilcox.test(cands\$val,multi\$val)\$p.value,digits=3)
pp5<-signif(wilcox.test(cands\$val,mono\$val)\$p.value,digits=3)
pp6<-signif(wilcox.test(cands\$val,non\$val)\$p.value,digits=3)
pp7<-signif(wilcox.test(hubs\$val,mnon\$val)\$p.value,digits=3)
labels<-c("Network","Hubs","Candidates","Multi NC","Multi","Mono","NC")
colors <- c("white", "white", "#3690C0", "white", "white","white","white")
boxplot(netw\$val,hubs\$val,cands\$val,mnon\$val,multi\$val,mono\$val,non\$val, names=FALSE, outline=F, range=1.5,notch=TRUE, varwidth=TRUE,col=colors)

## Ugly hacks to try and get the spacing right for
## the box labels
y<-max(netw\$val)/-500
y
if(y>-1 && y < 0){y<- 0.7}
if(y>-6 && y < 0 ){y<- -1}
if(y>-6 && y < -1){y<- -5}
if(max(netw\$val)==0) {y <- 0.7} 

y<-min(netw\$val)-1

y<-0.7

labels=c(paste("Network$sp1",paste("",format.pval(pp1, eps=0.05),"$sp1" ),sep="\\n"),paste("Hubs$sp1",paste("",format.pval(pp2, eps=0.05),"" ),sep="\n"),paste("Candidates$sp1",sep="\\n"),paste("Multi NC $sp3\\n",paste("",format.pval(pp3, eps=0.05)," $sp3",sep="")), paste("Multi$sp4",paste("",format.pval(pp4, eps=0.05),""),sep="\\n"), paste("Mono$sp1",paste("",format.pval(pp5, eps=0.05),"" ),sep="\\n"),paste("NC$sp4",paste("",format.pval(pp6, eps=0.05),"" ),sep="\\n"))
text(1:7, par("usr")[3] - 0.25, srt = 45, adj = 1,labels = labels, xpd = NA,cex=.8)


title("$title", cex.main=1.5)
Mean <- c(mean(netw\$val), mean(hubs\$val),mean(cands\$val),mean(mnon\$val),mean(multi\$val),mean(mono\$val),mean(non\$val))
points(Mean,col="#861e23",pch=20)
EEEE


###OLD

# text(1:7, labels = paste("Network$sp1",paste("(",format.pval(pp1, eps=0.05),")$sp1" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)
# text(2,y, labels = paste("Hubs$sp1",paste("(",format.pval(pp2, eps=0.05),")" ),sep="\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)
# text(3,y, labels = paste("Candidates$sp1",sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)
# text(4,y, labels = paste("Multi NC$sp1",paste("(",format.pval(pp3, eps=0.05),")$sp3" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)
# text(5,y, labels = paste("Multi$sp4",paste("(",format.pval(pp4, eps=0.05),")" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)
# text(6,y, labels = paste("Mono$sp1",paste("(",format.pval(pp5, eps=0.05),")" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)
# text(7,y, labels = paste("NC$sp4",paste("(",format.pval(pp6, eps=0.05),")" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = NA, cex=.9)


# text(1,y, labels = paste("Network$sp1",paste("(",format.pval(pp1, eps=0.05),")$sp1" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# text(2,y, labels = paste("Hubs$sp1",paste("(",format.pval(pp2, eps=0.05),")" ),sep="\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# text(3,y, labels = paste("Candidates$sp1",sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# text(4,y, labels = paste("Multi NC$sp1",paste("(",format.pval(pp3, eps=0.05),")$sp3" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# text(5,y, labels = paste("Multi$sp4",paste("(",format.pval(pp4, eps=0.05),")" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# text(6,y, labels = paste("Mono$sp1",paste("(",format.pval(pp5, eps=0.05),")" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
# text(7,y, labels = paste("NC$sp4",paste("(",format.pval(pp6, eps=0.05),")" ),sep="\\n"), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)

if ($species eq 'human' || $species eq 'vidal') {
    print "for(i in 1:length(gold\$val)){points(2.5,gold\$val[i],col=\"#ab6e0d\",pch=20) }\n";
}
print<<EEEE;
leg.txt<-c("Mean","Known MPs")
legend("topright",leg.txt,pch =20, col=c("#861e23","#ab6e0d"))
# mtext(paste("(",format.pval(pp1, eps=0.05),")"),side=1,line=2,at=2)
# mtext(paste("(",format.pval(pp2, eps=0.05),")"),side=1,line=2,at=4)
# mtext(paste("(",format.pval(pp3, eps=0.05),")"),side=1,line=2,at=1)
# mtext(paste("(",format.pval(pp4, eps=0.05),")"),side=1,line=2,at=5)
# mtext(paste("(",format.pval(pp5, eps=0.05),")"),side=1,line=2,at=6)
# mtext(paste("(",format.pval(pp6, eps=0.05),")"),side=1,line=2,at=7)
dev.off()
sink("$stats_file")
cat(paste(paste("Cands - Hubs\\t", pp1,"\\t",mean(cands\$val), "\\t",mean(hubs\$val)),paste("Cands - Multi NC\\t", pp2,"\\t",mean(cands\$val), "\\t",mean(mnon\$val)),paste("Cands - Network\\t", pp3,"\\t",mean(cands\$val), "\\t",mean(netw\$val)),paste("Cands - Multi\\t", pp4,"\\t",mean(cands\$val), "\\t",mean(multi\$val)),paste("Cands - Mono\\t", pp5,"\\t",mean(cands\$val), "\\t",mean(mono\$val)),paste("Cands - Non Cands\\t", pp6,"\\t",mean(cands\$val), "\\t",mean(non\$val)),paste("Hubs - Multi NC\\t", pp7,"\\t",mean(hubs\$val), "\\t",mean(mnon\$val)),"",  sep = '\\n'))
sink()

EEEE

}


sub usage{
        my $us="[options] <>";
	my $desc="";
	my %opts=(
		  "usage" => $us,
		  "desc" => $desc,
		  "p"=> "process phosphorylation data",
		  "P" => "Make go term pie chart",
		  "b"=> "Make boxplots.",
                  "a" => "Make figures for an article, title is only the feature analysed.",
		  "i"=> "String to be added to the output image file name.",
		  "d"=> "Make distribution plots (density function).",
		  "s"=> "Species. Def: Parsed from \$ARGV[0]",
		  "S"=> "For some reason, the labels do not align perfectly for every dataset. For example, they work for 'Degree' but not 'Betweeness'. This is an ugly hack that just adds spaces",
		  "t"=> "Graph title.",
		  "m"=> "Mode, can be one of i,B,iB or 'all'",
		  "f"=> "This is only useful for splitting gopies according to ontologies. The file to read go pie values from will be <species>.<mode>\$opts{f}. So, for example -f F.nums",
		  "n"=> "The file containing the values for the whole network. Def: <species>.all"
		 );
    print_help_exit(\%opts,0);

}
################################# R NOTES  ################################
# points(0.5,30,col="blue",pch=18)
# par(mar=c(4,3,3,10), xpd=TRUE)
# leg.txt<-c("Mean","Known MPs")
# legend("topright",leg.txt,pch =18, col=c("red","blue"))
# boxplot(cands,hubs,mnon,netw,multi,non,mono, names=labels, outline=FALSE, varwidth=TRUE) 
# legend("topright",inset=c(-0.2,0),leg.txt,pch =18, col=c("red","blue"))
# text(1.5,Mean[1],paste("P=",round(pp1\$p.value, digits=5)))
# text(2.5,Mean[2],paste("P=",round(pp2\$p.value, digits=5)))
# text(2.5,Mean[3],paste("P=",round(pp3\$p.value, digits=5)))
#par(mar=c(4,3,3,10), xpd=TRUE)

