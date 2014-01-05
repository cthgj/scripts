#!/usr/bin/perl

use Tk;
#use Tk::Carp qw/cluck/;
use Tk::ErrorDialog;
my %l; 
my $icon_file="/home/cchapple/mouse.gif";
$mw = MainWindow->new; 
$mw->geometry("650x70+300+380");
$mw->title("Juli Raton!");
my $icon = $mw->Photo(-file => '/home/cchapple/mouse.gif');
$mw->iconimage($icon);
my $dir='bob'; 
$mw->Label(-text=>"Use the button below to chose the directory containing the files you want to analyze",
		-font=>"helvetica 12") -> pack;
$mw->Button(-text => "Choose Directory", -command => \&choose_dir  )->pack(-side=>'left',-padx=>20 );
$mw->Button(-text => "Exit", -command => sub{ exit})->pack(-side=>'right',-padx=>20 );
 MainLoop; 

print "whohoooooo\n";
exit;

sub choose_dir{
    my @rep;
    $mw->withdraw(); ## hide main window
    ## Choose directory
    $dir = $mw->chooseDirectory(
	-initialdir => '~/',
	-title => 'Choose a directory'); 
    ## Report the dir found and ask what to do with it
    $rep_w= $mw->Dialog(
	-title => '', 
	-bitmap=>'question',
	-text=>"You have chosen \'" .  $dir . "\'. \n Is this correct?\n",
	-default_button => 'Run on this directory', 
	-buttons => [ 'Run on this directory','Choose another directory','Exit']
	)->Show();
    ## If this is the correct dir, execute program
    if ($rep_w eq 'Run on this directory'){
	## Check if another sumas.txt or todos.txt exists
	my $is_there1=`find $dir/ -name todos.txt | head -1`;
	my $is_there2=`find $dir/ -name sumas.txt | head -1`;
	chomp($is_there1,$is_there2);
	my $continue='Yes';
	if((length($is_there1)!=0)||(length($is_there2)!=0)){
	    $continue = $mw->Dialog(
		-title => 'Cuidadin!', 
		-text => "There is already a file \'$is_there1\', if you continue it (and any others) will be overwritten. \n\nContinue?\n", 
		-default_button => 'Yes', 
		-buttons => [ 'Yes', 'No'], 
		-bitmap => 'question' )->Show( ); 
	}
	## If everything is OK, continue.
	if ($continue eq 'Yes') {	  
	    ## Find lowest level dirs
	    my @low=@{&lowest_dir($dir)};
	    foreach my $dd (@low){   
		my ($res,$res1)=&juli_raton($dd);
		push @rep,$res;
		push @rep,$res1;
	    }
	    $mw->deiconify( ); 
	    $mw->raise( ); 
	}
	else{
	    $mw->deiconify( ); 
	    $mw->raise( ); 
	}
    }
    ## Else, chose another dir
    elsif ($rep_w eq 'Choose another directory'){
	&choose_dir;
	
    }
    else{exit};
    &print_report(\@rep);
}



sub lowest_dir{
    my $dir=shift;
    %l=();
    &find_lowest_dirs($dir);

    sub find_lowest_dirs{
	my $dir=shift;
	opendir DH, $dir or die "Failed to open $dir: $!";
	my @d;
	while ($_ = readdir(DH)) {
	    next if $_ eq "." or $_ eq "..";
	    my $fn = $dir . '/' . $_;
	    if (-d $fn) {
		push @d, $fn;
	    }
	}
	if (scalar @d == 0) { #If no directories found, $dir is lowest dir in this branch
	    $l{$dir}++;
	}
	foreach (@d) {
	    &find_lowest_dirs($_); #Look for subdirectories in directory
	}
    }
    my @lowest=keys(%l);
    return(\@lowest);
}

sub juli_raton{
    my $name;
    my (%mm,%data,%totals,%data2);
    my $dir=shift;

    my $outfile1="$dir/todos.txt";
    my $outfile2="$dir/sumas.txt";
    open(OUT,">$outfile1");
    open(OUT2,">$outfile2");
    my $is_there=`find $dir -type f -name "*.E*"`;
    return unless $is_there;
    open(IN, "cat $dir/*.E*|");
    while(<IN>){
	next if /^Point/;
	if(/^..(.{6}).+\.E/){
	    $name=$1;
	    $mm{$name}++;	
	}
	elsif(/^totals\s*(\d+)/i){
	    $totals{$name}=$1;
	}
	else{
	    chomp;
	    /(\d+)\s+(\d+)/;
	    my @a;
	    $data{$1}{$name}=$2;
	}
    }
    my @mice=sort(keys(%mm));
    my @oo=keys(%data);
    foreach my $time (sort{ $a <=> $b }keys(%data)){
	foreach my $mouse (@mice){
	    $data{$time}{$mouse}='n/a' unless defined($data{$time}{$mouse});
	}
    }
    foreach my $mouse (@mice){
	print OUT "\t$mouse";
    }
    print OUT "\n";
    my @times=sort{ $a <=> $b } keys(%data);
    foreach my $time (@times){
	print OUT "$time";
	foreach my $mouse (@mice){
	    print OUT "\t$data{$time}{$mouse}";
	}
	print OUT "\n";

    }
    print OUT "Tot\t";
    foreach my $mouse (@mice){
	print OUT "$totals{$mouse}\t";
    }
    print OUT "\n";
    my $rep="$outfile1\n";
    foreach my $mouse (@mice){
	print  OUT2 "\t$mouse";
    }
    print OUT2 "\n";
    my $rep1="$outfile2\n";

#    my @times2=sort{ $a <=> $b }keys(%{$data2{$mouse}});
    for (my $n=1; $n<scalar(@times); $n=$n+2){
	print OUT2 "$times[$n-1]+$times[$n]";
	foreach my $mouse (@mice){
	    $data{$times[$n]}{$mouse}=0 if $data{$times[$n]}{$mouse}=~'n/a';
	    $data{$times[$n-1]}{$mouse}=0 if $data{$times[$n-1]}{$mouse}=~'n/a';
	    my $sum=$data{$times[$n]}{$mouse} + $data{$times[$n-1]}{$mouse};
	    print OUT2 "\t$sum "  ;
	}
	print OUT2 "\n";
    }

    return($rep,$rep1);
 }

sub print_report{
    $mw->withdraw(); ## hide main window

    my $rreport=shift;
    my @report=@{$rreport};
    print "rr : @report\n";
    my $rep= MainWindow->new; 
    $rep->geometry("650x370+300+380");
    my $icon = $rep->Photo(-file => $icon_file);
    $rep->iconimage($icon);
    $rep->title("Juli Raton: Report");
    my $icon = $rep->Photo(-file => $icon_file);
    $rep->iconimage($icon);

    $rep->Label(-text=>"The following files were created:\n\n@report\n",
		-font=>"helvetica 12") -> pack;
    $rep->Button(
	-text=>"OK",
	-command => sub{
	    $rep->withdraw(); ## hide main window
	    $mw->deiconify( ); 
	    $mw->raise( ); 
	}
	)->pack(-side=>'left',-padx=>120 );
    $rep->Button(
	-text=>"Exit",
	-command => sub{exit}  
	)->pack(-side=>'right',-padx=>120 );
}
