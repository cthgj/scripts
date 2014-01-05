#!/usr/bin/perl

use strict;
use Gtk2 '-init';
use Glib qw/TRUE FALSE/; 
#standard window creation, placement, and signal connecting
my $window = Gtk2::Window->new('toplevel');
$window->signal_connect('delete_event' => sub { Gtk2->main_quit; });
$window->set_border_width(5);
$window->set_position('center_always');
$window->set_title (my $title='Juli Raton!');
my $label_feedback;

#add and show the vbox
open(AA,">/home/cchapple/hahahaha");
print AA "WUWUWUWUWUWUW\n";

$window->set_icon_from_file ("/home/cchapple/mouse.gif"); 
$window->add(&ret_vbox);
$window->show();

#our main event-loop
Gtk2->main();

sub ret_vbox {
    
    my $vbox = Gtk2::VBox->new(FALSE,5);
    	$label_feedback = Gtk2::Label->new();
	my $image=Gtk2::Image->new_from_file("/home/cchapple/mouse_tr.gif");
    $vbox->pack_start($image,FALSE,FALSE,6);

    #***************************************
    #Show the filechooserbuttons (open and select-folder action types)
    #Select Folder---->
    my $hbox_fl_chooser_dialog = Gtk2::HBox->new(FALSE,5);
    $hbox_fl_chooser_dialog->set_border_width(10);
    my $frm_fl_chooser_dialog = Gtk2::Frame->new('Select the folder containing your data:');
    my $folder;
    my $btn_select_folder   = Gtk2::Button->new('Select a folder');
    $btn_select_folder->signal_connect('clicked' => 
				       sub{ $folder=&show_chooser('File Chooser type select-folder','select-folder'); });
    
    # my $btn_select_folder =Gtk2::FileChooserButton->new ('select a folder' , 'select-folder');
    my $e_button=Gtk2::Button->new('Exit');
    $e_button->signal_connect("clicked" => sub{exit}, "Exit");    
    $hbox_fl_chooser_dialog->pack_start($btn_select_folder,TRUE,TRUE,6);
    $hbox_fl_chooser_dialog->pack_start($e_button,TRUE,TRUE,6);
    $frm_fl_chooser_dialog->add($hbox_fl_chooser_dialog);
    $vbox->pack_start($frm_fl_chooser_dialog,FALSE,FALSE,6);

    

    $vbox->show_all();
    return $vbox;
}

sub show_chooser {
    #---------------------------------------------------
    #Pops up a standard file chooser--------------------
    #Specify a header to be displayed-------------------
    #Specify a type depending on your needs-------------
    #Optionally add a filter to show only certain files-
    #will return a path, if valid----------------------
    #---------------------------------------------------
    $window->hide();
    my($heading,$type,$filter) =@_;
    #$type can be:
    #* 'open' 
    #* 'save' 
    #* 'select-folder'
    #* 'create-folder' 
    my $file_chooser =  Gtk2::FileChooserDialog->new ( 
	$heading,
	undef,
	$type,
	'gtk-cancel' => 'cancel',
	'gtk-ok' => 'ok'
	);
    (defined $filter)&&($file_chooser->add_filter($filter));
    
    my $filename;
    
    if ('ok' eq $file_chooser->run){    
        $filename = $file_chooser->get_filename;
    }
     $file_chooser->destroy;
    my (@rep,@res);
    my $response=&show_message_dialog($window,'question',"You have chosen <b>$filename</b>, is this correct?\n\nClick 'OK' to run on this directory, or 'Cancel' to choose another.\n\n",'ok-cancel');
    if($response eq 'ok'){
	my $is_there=`find $filename -type f -name "*.E??"`;
	## If no appropriate files are found, warn
	unless ($is_there){
	    my $response=&show_message_dialog($window,'warning',"No appropriate file were found in folder $filename\n\n The files should have an extension like .E12\n\n",'ok');
	    $window->show();
	    return;
	}
	## Check if another sumas.txt or todos.txt exists
	my $is_there1=`find $filename/ -name todos.txt | head -1`;
	my $is_there2=`find $filename/ -name sumas.txt | head -1`;
	chomp($is_there1,$is_there2);
	if((length($is_there1)!=0)||(length($is_there2)!=0)){
	    my $response=&show_message_dialog($window,'question',"There is already a file \'$is_there1\', if you continue it (and any others) will be overwritten. \n\nContinue?\n",'ok-cancel');
	    if($response eq 'cancel'){ 
		$window->show();
		return;
	    }
	}
	
	my @low=@{&lowest_dir($filename)};
	foreach my $dd (@low){
	    
	    print "dd : $dd\n";
	    my ($res,$res1)=&juli_raton($dd);
	    push @rep,$res;
	    push @rep,$res1;
	}
	print STDERR "Done...\n";
	&show_message_dialog($window,'info',"The following files were created:\n\n@rep\n",'ok');
	$window->show();
    }
    else{$window->show();}
    return;
}


sub show_message_dialog {
    #---------------------------------------------------
    #you tell it what to display, and how to display it
    #$parent is the parent window, or "undef"
    #$icon can be one of the following: a) 'info'
    #                   b) 'warning'
    #                   c) 'error'
    #                   d) 'question'
    #$text can be pango markup text, or just plain text, IE the message
    #$button_type can be one of the following:  a) 'none'
    #                       b) 'ok'
    #                       c) 'close'
    #                       d) 'cancel'
    #                       e) 'yes-no'
    #                       f) 'ok-cancel'
    #---------------------------------------------------
    
    my ($parent,$icon,$text,$button_type) = @_;
    
    my $dialog = Gtk2::MessageDialog->new_with_markup ($parent,
						       [qw/modal destroy-with-parent/],
						       $icon,
						       $button_type,
						       sprintf "$text");
  
    my $retval = $dialog->run;
    $dialog->destroy;
    return $retval;
}


sub lowest_dir{
    my $dir=shift;
    my %l=();
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

