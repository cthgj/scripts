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
$window->set_icon_from_file ("/home/cchapple/mouse.gif"); 
$window->add(&ret_vbox);
$window->show();

#our main event-loop
Gtk2->main();

sub ret_vbox {
    
    my $vbox = Gtk2::VBox->new(FALSE,5);
    	$label_feedback = Gtk2::Label->new();
	
    #***************************************
    #Show the filechooserbuttons (open and select-folder action types)
    #Select Folder---->
    my $hbox_fl_chooser_dialog = Gtk2::HBox->new(FALSE,5);
    $hbox_fl_chooser_dialog->set_border_width(10);
    my $frm_fl_chooser_dialog = Gtk2::Frame->new('Select the folder containing your data:');
    my $folder;
    #my $btn_select_folder   = Gtk2::Button->new();
    #$btn_select_folder->signal_connect('clicked' => 
#				       sub{ $folder=&show_chooser('File Chooser type select-folder','select-folder'); });
    
    my $btn_select_folder =Gtk2::FileChooserButton->new ('select a folder' , 'select-folder');
    my $filename = $btn_select_folder->get_filename;    print STDERR "ff $filename\n";die();

    my $e_button=Gtk2::Button->new('Exit');
    $e_button->signal_connect("clicked" => sub{exit}, "Exit");
    
    my $run_button=Gtk2::Button->new('Run');
    $run_button->signal_connect('clicked' => sub{&run_juli_raton($filename)}, 'Run');
# this calls our box creating function
   # my $box = xpm_label_box("/home/cchapple/mouse.xpm", 'Select Folder');

# pack and show all our widgets
   # $box->show();
    #$btn_select_folder ->add($box);
    #$button->show();
    #$vbox->add($button);
    
    $hbox_fl_chooser_dialog->pack_start($btn_select_folder,TRUE,TRUE,6);
    $hbox_fl_chooser_dialog->pack_start($e_button,TRUE,TRUE,6);
    $hbox_fl_chooser_dialog->pack_start($run_button,TRUE,TRUE,6);
    
    #Create Folder---->
#    $vbox->pack_end(Gtk2::HSeparator->new(),FALSE,FALSE,4);

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
        print "filename $filename\n";
    }
     $file_chooser->destroy;
    my (@rep,@res);
    my $response=&show_message_dialog($window,'question',"You have chosen <b>$filename</b>, is this correct?\n\nClick 'OK' to run on this directory, or 'Cancel' to choose another.\n\n",'ok-cancel');
    print STDERR "RR : $response\n";
    if($response eq 'ok'){
	my @low=@{&lowest_dir($filename)};
	foreach my $dd (@low){   
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


sub ret_png_filter {
    #----------------------------------------
    #Returns a filter, filtering only png files
    #----------------------------------------
    
    my $filter = Gtk2::FileFilter->new();
    $filter->set_name("Images");
    $filter->add_mime_type("image/png");
    
    return $filter;
}


sub lowest_dir{
    my $dir=shift;
    my %l=();
    &find_lowest_dirs($dir);

    sub find_lowest_dirs{
     print STDERR "..1 $dir..\n";

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





# Create a new hbox with an image and a label packed into it and return
# the box.
sub xpm_label_box {
	my ($xpm_filename, $label_text) = @_;

	# create box for image and label 
	my $box = Gtk2::HBox->new(FALSE, 0);
	$box->set_border_width(2);

	# now on to the image stuff
	my $image = Gtk2::Image->new_from_file($xpm_filename);
	
	# Create a label for the button
	my $label = Gtk2::Label->new($label_text);
	
	# pack the image and label into the box
	$box->pack_start($image, FALSE, FALSE, 3);
	$box->pack_start($label, FALSE, FALSE, 3);

	$image->show;
	$label->show;

	return $box
}


sub run_juli_raton{
    my $filename=shift;
    print STDERR "ff $filename\n";die();
    my (@rep,@res);
    # my $response=&show_message_dialog($window,'question',"You have chosen <b>$filename</b>, is this correct?\n\nClick 'OK' to run on this directory, or 'Cancel' to choose another.\n\n",'ok-cancel');
    # if($response eq 'ok'){
	my @low=@{&lowest_dir($filename)};
     print STDERR "....\n";
	foreach my $dd (@low){  
	    my ($res,$res1)=&juli_raton($dd);
	    push @rep,$res;
	    push @rep,$res1;
	}
	print STDERR "Done...\n";
	
	&show_message_dialog($window,'info',"The following files were created:\n\n@rep\n",'ok');
		$window->show();
	
    # }
    # else{$window->show();}
}
