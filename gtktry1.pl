#!/usr/bin/perl

use Glib qw/TRUE FALSE/;
use Gtk2 '-init';

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

# our usual callback function
sub callback {
	my $widget = shift;
	my $data   = shift;
	printf "Hello again - %s was pressed\n", $data;
}

##################
# main starts here

# create a new window
my $window = Gtk2::Window->new('toplevel');

$window->set_title("Pixmap'd Buttons!");

# # it's a good idea to do this for all windows
# $window-&gt;signal_connect("destroy" => sub { Gtk2->main_quit; });
# $window-&gt;signal_connect("delete-event" => sub { Gtk2->main_quit; });

# sets the border width of the window
$window->set_border_width(10);

# create a new button
my $button = Gtk2::Button->new();

# connect the 'clicked' signal of the button to our callback
#$button->signal_connect("clicked" => \&callback, "cool button");

# this calls our box creating function
my $box = xpm_label_box("mouse.gif", "cool button");

# pack and show all our widgets
$box->show();

$button->add($box);

$button->show();

$window->add($button);

$window->show();

# rest in the GTK main loop and wait for the fun to begin!
Gtk2->main;

0;
	
