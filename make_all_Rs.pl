#!/usr/bin/perl -w
use Getopt::Std;
my %opts;
getopts('pd',\%opts);
my $pie;
$pie=1 if $opts{p};
my $dist=$opts{d}||undef;

my $species="human";
my $mode="i";

my %titles =(
"i"=>"Intersection",
"B"=>"Bridge",
"iB"=>"Intersection & Bridge"
);

my %EXTS=(
"degree" => "degnums",
"disorder" => "perc",
"domains" => "domnums",
"length" => "lennums",
"betweeness" => "betnums",
"shortest_paths" => "paths"
# "" => "",
# "" => "",
# "" => "",

);

my @types=qw(cands hubs multi_noncands non mono multi netw);
my %names=(
"cands" => "1. Candidates",
"hubs" => "2. Hubs",
"multi_noncands" => "3. Candidates",
"non" => "4 NonCands",
"mono" => "5. Mono",
"netw" => "6. Network"
);

foreach (keys(%EXTS)) {
    foreach my $type (@types) {
        scan($species,$mode,$_,$type)
    }
}
foreach (keys(%EXTS)) {
    make_dfs($_);
}


sub scan{
    my $species=shift;
    my $which=shift;
    my $type=shift;

    my $ext=$EXTS{$which};
    print<<EndOf;
$species.$type<-scan("$species.$type.$ext");
EndOf
}

sub make_dfs{
    my $which=shift;
    print "$which :: $names{$which} :: $which<-data.frame(grp=factor(c(rep(";
    foreach (@type) {
	print ",rep(\"$names{$which}\", each=length($which))";
    }
    print "\n";


}


exit(0);


#####################
# Plot distribution #
#####################
if ($dist) {
    my $m=$ARGV[1];
    my $S =  my $s = $ARGV[0];
    $S=~/^(.)/;
    my $k=uc($1);
    $S=~s/^(.)/$k/;
    print<<EndOf;
bet<-scan("$s.$m.cands.betnums")
deg<-scan("$s.$m.cands.degnums")
dis<-scan("$s.$m.cands.disnums")
len<-scan("$s.$m.cands.lennums")
dom<-scan("$s.$m.cands.AB.domnums")
cands<-data.frame(bet=bet,deg=deg,dis=dis,len=len,dom=dom) 

bet<-scan("$s.$m.non.betnums")
deg<-scan("$s.$m.non.degnums")
dis<-scan("$s.$m.non.disnums")
len<-scan("$s.$m.non.lennums")
dom<-scan("$s.$m.non.AB.domnums")
non<-data.frame(bet=bet,deg=deg,dis=dis,len=len,dom=dom) 

bet<-scan("$s.$m.multi_noncands.betnums")
deg<-scan("$s.$m.multi_noncands.degnums")
dis<-scan("$s.$m.multi_noncands.disnums")
len<-scan("$s.$m.multi_noncands.lennums")
dom<-scan("$s.$m.multi_noncands.AB.domnums")
multi_non<-data.frame(bet=bet,deg=deg,dis=dis,len=len,dom=dom) 

bet<-scan("$s.mono.betnums")
deg<-scan("$s.mono.degnums")
dis<-scan("$s.mono.disnums")
len<-scan("$s.mono.lennums")
dom<-scan("$s.mono.AB.domnums")
mono<-data.frame(bet=bet,deg=deg,dis=dis,len=len,dom=dom) 

bet<-scan("$s.multi.betnums")
deg<-scan("$s.multi.degnums")
dis<-scan("$s.multi.disnums")
len<-scan("$s.multi.lennums")
dom<-scan("$s.multi.AB.domnums")
multi<-data.frame(bet=bet,deg=deg,dis=dis,len=len,dom=dom) 


Mean <- round(c(mean(cands\$deg), mean(multi_non\$deg),mean(multi\$deg),mean(non\$deg),mean(mono\$deg)),digits=2)
Median <- c(median(cands\$deg), median(multi_non\$deg),median(multi\$deg),median(non\$deg),median(mono\$deg))

df_deg<-data.frame(grp=factor(c(rep("1. Candidates",each=length(cands\$deg)), rep("2. Multi NonCands",each=length(multi_non\$deg)),rep("3. Multi",each=length(multi\$deg)),rep("4. NonCands",each=length(non\$deg)),rep("5. Mono",each=length(mono\$deg)))), val=c(cands\$deg, multi_non\$deg, multi\$deg, non\$deg, mono\$deg))

df_dis<-data.frame(grp=factor(c(rep("1. Candidates",each=length(cands\$dis)), rep("2. Multi NonCands",each=length(multi_non\$dis)),rep("3. Multi",each=length(multi\$dis)),rep("4. NonCands",each=length(non\$dis)),rep("5. Mono",each=length(mono\$dis)))), val=c(cands\$dis, multi_non\$dis, multi\$dis, non\$dis, mono\$dis))

df_dom<-data.frame(grp=factor(c(rep("1. Candidates",each=length(cands\$dom)), rep("2. Multi NonCands",each=length(multi_non\$dom)),rep("3. Multi",each=length(multi\$dom)),rep("4. NonCands",each=length(non\$dom)),rep("5. Mono",each=length(mono\$dom)))), val=c(cands\$dom, multi_non\$dom, multi\$dom, non\$dom, mono\$dom))

df_len<-data.frame(grp=factor(c(rep("1. Candidates",each=length(cands\$len)), rep("2. Multi NonCands",each=length(multi_non\$len)),rep("3. Multi",each=length(multi\$len)),rep("4. NonCands",each=length(non\$len)),rep("5. Mono",each=length(mono\$len)))), val=c(cands\$len, multi_non\$len, multi\$len, non\$len, mono\$len))

df_bet<-data.frame(grp=factor(c(rep("1. Candidates",each=length(cands\$bet)), rep("2. Multi NonCands",each=length(multi_non\$bet)),rep("3. Multi",each=length(multi\$bet)),rep("4. NonCands",each=length(non\$bet)),rep("5. Mono",each=length(mono\$bet)))), val=c(cands\$bet, multi_non\$bet, multi\$bet, non\$bet, mono\$bet))

library(ggplot2)

png(filename="$s.$m.disorder.png")
ggplot(df_dis, aes(x=val)) + geom_density() + facet_grid(grp ~ .) +  labs(title="$S $names{$m} - Disorder", x="Disorder", y="Density")
dev.off()
png(filename="$s.$m.disorder.ovrl.png")
ggplot(df_dis, aes(x=val,colour=grp)) + geom_density()  +  labs(title="$S $names{$m} - Disorder", x="Disorder", y="Density") + scale_color_brewer(palette="Dark2")
dev.off()

png(filename="$s.$m.domains.png")
ggplot(df_dom, aes(x=val)) + geom_density() + facet_grid(grp ~ .) +  labs(title="$S $names{$m} - Domains", x="Domains", y="Density") + xlim(0,20)
dev.off()
png(filename="$s.$m.domains.ovrl.png")
ggplot(df_dom, aes(x=val,colour=grp)) + geom_density()  +  labs(title="$S $names{$m} - Domains", x="Domains", y="Density") + scale_color_brewer(palette="Dark2") + xlim(0,20)
dev.off()

png(filename="$s.$m.length.png")
ggplot(df_len, aes(x=val)) + geom_density() + facet_grid(grp ~ .) +  labs(title="$S $names{$m} - Length", x="Length", y="Density") + xlim(0,5000) + xlim(0,5000)
dev.off()
png(filename="$s.$m.length.ovrl.png")
ggplot(df_len, aes(x=val,colour=grp)) + geom_density()  +  labs(title="$S $names{$m} - Length", x="Length", y="Density") + scale_color_brewer(palette="Dark2") + xlim(0,5000)
dev.off()

png(filename="$s.$m.betweeness.png")
ggplot(df_bet, aes(x=val)) + geom_density() + facet_grid(grp ~ .) +  labs(title="$S $names{$m} - Betweeness", x="Betweeness", y="Density") + xlim(10000,50000)
dev.off()
png(filename="$s.$m.betweeness.ovrl.png")
ggplot(df_bet, aes(x=val,colour=grp)) + geom_density()  +  labs(title="$S $names{$m} - Betweeness", x="Betweeness", y="Density") + scale_color_brewer(palette="Dark2") + xlim(10000,50000)

dev.off()

png(filename="$s.$m.degree.png")
ggplot(df_deg, aes(x=val)) + geom_density() + facet_grid(grp ~ .) +  labs(title="$S $names{$m} - Degree", x="Degree", y="Density") + xlim(1,200)
dev.off()
png(filename="$s.$m.degree.ovrl.png")
ggplot(df_deg, aes(x=val,colour=grp)) + geom_density()  +  labs(title="$S $names{$m} - Degree", x="Degree", y="Density")  + scale_color_brewer(palette="Dark2") + xlim(1,200)
dev.off()
EndOf

exit;

}
else {
    $ARGV[0]=~/^(.+?)\.(.+?)\./;
    my $m=$2;
    my $S =  my $s = $1;
    $S=~/^(.)/;
    my $k=uc($1);
    $S=~s/^(.)/$k/;
    my $title="$S $ARGV[2] - ";


    ##################
    # Make Pie Chart #
    ##################
    if ($pie) {
	print<<EndOf;
titles<-c("i"="Intersection","B"="Bridge","iB"="Intersection & Bridge")
cands<-scan("$s.$m.nums")
pct <- round(cands/sum(cands)*100)
lbls<-c("Neither","One", "Both")
lbls <- paste(lbls, pct)
lbls <- paste(lbls,"%",sep="")
lbls <- paste(lbls,"(", sep=" ")
lbls <- paste(lbls,cands, sep="")
lbls <- paste(lbls,")", sep="")
png(filename="$ARGV[1]", units = "px", pointsize = 12, bg = "white")
pie(cands,labels = lbls)
title(paste("$title",titles["$m"], "\n",sum(cands), " GO pairs", sep=""))
a<-c(paste("$S GO pairs - ",titles["$m"], " (",sum(cands), ")", sep=""),lbls)
write(a,"$s.$m.perc")


EndOf
    }
    #####################
    # Make scatter plot #
    #####################
    # elsif ($scatter) {
	


    # }
    #################
    # Make Box Plot #
    #################
    else {
	## Special case for pfam, have A and AB
	my $stats_file="$s.$m";
	my $title="$S $ARGV[6] ";
	my $ab;
	my $means_file="$s.$m";
	if ($ARGV[0]=~/\.(AB)\./ || $ARGV[0]=~/\.(A)\./) {
	    $stats_file.=".$1";
	    $title.=" (pfam $1)";
	    $means_file.=".$1";
	}
	$stats_file.=".stats";
	$means_file.=".means";
	$title.=" - ";

	print<<EndOf;
titles<-c("i"="Intersection","B"="Bridge","iB"="Intersection & Bridge")
cands<-scan("$ARGV[0]")
non<-scan("$ARGV[1]")
multi<-scan("$ARGV[2]")
mono<-scan("$ARGV[3]")
multi_non<-scan("$ARGV[4]")
png(filename="$ARGV[5]", width = 800, height = 800, units = "px", pointsize = 12, bg = "white")
labels<-c("Candidates","Multi NonCands","Multi","NonCands","Mono")
boxplot(cands,multi_non,multi,non,mono, names=labels, outline=FALSE, varwidth=TRUE)
title(paste("$title",titles["$m"]))
Mean <- round(c(mean(cands), mean(multi_non),mean(multi),mean(non),mean(mono)),digits=2)
Median <- c(median(cands), median(multi_non),median(multi),median(non),median(mono))
points(Mean,col="red",pch=18)
names(Mean)<-labels
pp<-wilcox.test(cands,multi_non)
text(1.5,Mean[1],paste("P=",round(pp\$p.value, digits=5)))
aa<-data.frame(Mean,Median)
sink("$means_file")
aa
sink()
sink("$stats_file")
t.test(cands,non)
wilcox.test(cands,non)
t.test(cands,multi)
wilcox.test(cands,multi)
t.test(cands,mono)
wilcox.test(cands,mono)
t.test(cands,multi_non)
wilcox.test(cands,multi_non)
sink()

dev.off()
EndOf
    }
}
