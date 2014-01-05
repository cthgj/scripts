#!/usr/bin/perl
## simple daemon to allow my user access to the Xserver

my @hosts;
my %allowed;
my $user=`echo \$USER`;
chomp($user);
my $who="bob";
my $old_who;
&get_hosts();

while(1){
    $old_who=$who;
    $who=`who | grep $user | gawk '{print \$NF}'`;
    if($who ne $old_who){
	$old_who=$who;
	@a=split(/\n/,$who);
	my @WHO;
	map{/\(([\w-\.]+)/; push @WHO,$1;}@a;
	map{
	    unless (defined($allowed{$_})){
		`xhost $_ 2>>/home/juka/xerr` ;
		&get_hosts();
	    }
	}@WHO; 
    }
    sleep(3);
}

sub get_hosts{
    $a=`xhost | grep -v access`;
    @hosts=split(/\n/,$a);
    foreach my $host (@hosts){
	$allowed{$host}++;
    }

}
