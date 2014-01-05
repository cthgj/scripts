#!/usr/bin/perl

use feature "switch";
use Encode;
use Text::Wrap;
use Getopt::Std;

my $weatherline;

# This script was written by lvleph and inspired by the original conky weather script written by azhag (azhag@bsd.miki.eu.org)
## It has been significantly modified by terdon, (singollo@yahoo.com)

# Modyfied by LazarusHC to list more details


#getopts('l:,h',\%opts) || do { print "Invalid option\n"; usage(); };
my $offset=1;
my $night;
my $location=$ARGV[0] || &usage(); #zipcode or weather.com city code
my ($htemp,$ltemp);
# ($code,$code2)=split(/-/,$ARGV[0]); # weather.com and weather.noaa.gov codes
my $system =$ARGV[1]|| &usage(); #f for imperial c for metric
$what=$ARGV[2]|| &usage();	 #what are we looking for?
$file="/tmp/weather.html";	 #temp holding weather
my @lines;
$update=3600; #time in seconds to update $file if set to 0 don't use $file
my ($LOC,$code1,$code2);
$hleadspace="    ";		#spacing before each high
$htrailspace="    ";		#spacing after each high.

$lleadspace="       ";		#spacing before each  low
$ltrailspace=" ";		#spacing after each low.
$leadspace="   ";		#spacing before each high low
$trailspace="    ";		#spacing after each high low.
$fspaces="  ";			#spacing between condition symbols.
$dspaces="";			#spacing between each day
$lines="\n\n\n\n"; #each \n represents one line between the days and temps

$Text::Wrap::columns = 58;
$initial_tab="";      #tab before first line in weather output
$subsequent_tab="\t"; #tab before each subsequet line in weather output

$degree= encode_utf8( "\x{00B0}" ); #give me the degree symbol, not everyone has same locale


#ensure user inputs proper system
if ($system !=~ "c" || $system !=~ "f") {
    $what=0;
}				#this will give usage error 
## Get location

if ($location =~ /bcn/i) {
    $LOC="Barcelona";
    $code1="spain/catalonia/barcelona-753692";
    $code2="LEBL";
} elsif ($location=~ /ath/i) {
    $LOC="Athens";
    $code1="greece/attiki/athens-946738";
    $code2="LGAV";
} elsif ($location=~ /mar/i) {
    $LOC="Marseille";
    #    $code1="france/provence-alpes-cote-d%27azur/marseille-610264";
    #    $code1="france/provence-alpes-c%c3%b4te-dazur/marseille-610264/";
    $code1="062/c01055.htm";
    $code2="LFML";
} else {
    $LOC="Athens";
    $code1="greece/attiki/athens-946738";
    $code2="LGAV";
    
}

given($what) {
    when (/"c"/) 
    {
	&file_op;		#save weather to $file
	while (<FILE>) {	#cycle through file
	    if (/<dt>Feels Like:<\/dt><dd>([\+-\d]+)/) { #found feels like temp
		($tmf) = $1;
	    }
	    if (/<dt>Humidity:<\/dt>/) { #found current humidity
		($hmt) = /<dd>(\d+\s*\%)/; #save current humidity					
	     
	    }
	    if (/<dt>Wind:<\/dt><dd>(.+?)<\/dd>/) { #found wind conditions
		($wnd) = $1;	#save wind conditions
		#do we have current conditions?
	    }
	    if (/>Sunrise:<.+?(\d+:\d+)/i) {
		($sunrise) = $1;
	    }
	    if (/<dt>Sunset:<.+?(\d+:\d+)/) {
		($sunset) = $1; 
	    }
	    
	    #	    if( $tmf && $hmt && $wnd && $sunset && $sunrise ){
	    if ( /id=.yw-temp.>([-\d]+)/) {
		$realtemp=$1;
		print "       Temperature: $realtemp $degree ($tmf)\n";
		print "       Humidity:        $hmt\n"; 
		print "       Wind:        $wnd\n"; 
		print "       Sunrise:     $sunrise\n"; 
		print "       Sunset:      $sunset\n"; 	       
		
		exit;
	    }
	}
    }
    when (/[1-5]dp$/)
    {		  #display the conditions from today through day $days
	&file_op; #save weather to $file
	my $day=(split "dp", $what)[0]; #how many days are we looking for
	my $flag=0;	#set flag for when we find start of conditions
	my $count=0; 
	while (<FILE>) {
	    if ($_=~/id=.yw-fivedayforecast.+?class=\"(night)\"/) {
		$night=1; 
	    }
	    if (/class=\"fiveday/) {
		$flag=1;$weatherline .=$_;
	    }			#found the start of conditions
	    elsif ($flag && /^\s*<td/) {
		$weatherline .=$_;
		$weatherline=~s/(\w+)-\n/$1/;
	    } elsif ($flag && $weatherline=~ /<\/div><\/td>/ && ++$count<=$day) { #look for the conditions upto specified day
		$weatherline .=$_;
		$weatherline=~s/\s*(\w+)-\n/$1/;
		@cnd=($weatherline=~/<br\/>(.+?)<\/div><\/td>/g);
	    } elsif ($flag && /SpaceID=0/) { #don't keep looking if everything has been found 
		map{&cond_symb($_);$night=0}@cnd;
		exit;
	    } elsif (/fiveday-temps/) {
		$flag=0;
	    }
	}
    }

    when (/[1-5]d$/)
    {
	&file_op;		#save weather to $file
	my $count1 = my $count2=0;
	my $day=(split "dt", $what)[0]; #how many days are we looking for
	my $flag=0;			#print days once
	while (<FILE>) { 
	    #get the high temp
	    (my $high) = /High:\s*(\d+)/;
	    #get the low temp
	    (my $low) = />Low:\s*(\d+)/;
	    #print the high and low temp for the specified day
	    if (/yw-fivedayforecast/ && ++$count1<=$day) { #look for the conditions upto specified daya
		$_=~/id=.yw-fivedayforecast.+?>(.+nobg.>)/; 
		$kk=$1;
		$kk=~s/<.+?>/ /g;
		$kk=~s/^\s+//;
		@days=split(/\s+/,$kk);
		#print "dd : @days\n";
		map{&day_space($_); $count++; $count2++}@days;
	    } elsif ($high=~/\d+/ && $low=~/\d+/ && ++$count2<=$day) {
		print "c : $count2\n";$ttext.=$leadspace.$high.$degree."-".$low.$degree.$trailspace;
	    } elsif ($count2>=$day) {
		print "  $dtext\n"; exit;
	    }	      #don't keep looking if everything has been found
	}
    }
    when (/[1-5]dt/)
    {
	&file_op;		#save weather to $file
	my $count1 = my $count2=0;
	my $day=(split "dt", $what)[0]; #how many days are we looking for
	my $flag=1;			#print days once
			       
	my $found=0;
	my $hl=0;
	my @forecast_days;
	while (<FILE>) {
	    if (/forecast_header_new/) {
		$found++;
	    } 
	    elsif (/FDAY/) {
		$found=0;
	    }
	    next unless $found>0;
	    ## Check if we are at high or low
	    if (/#0000ff/) {
		$hl=0;
	    }
	    if (/#ff3300/) {
		$hl=1;
	    }
	    # Collect days
	    if (/\(\s*(\w{3})\s*/) {
		push @forecast_days, $1;
	    }
	    # Get high and low
	    if (/^\s*(\d+)\s*$/) {
		push @temps, $1
	    }
	}
	print "@forecast_days\n@temps\n";
	print "   \n";
	for (my $nn=0;$nn<scalar(@temps);$nn+=2) {
	    if (length($temps[$nn])<2) {
		$temps[$nn]=" " . $temps[$nn];
	    }
	    print "     $temps[$nn] ";
	}
	print "   \n";
	for (my $nn=1;$nn<scalar(@temps);$nn+=2) {
	    if (length($temps[$nn])<2) {
		$temps[$nn]=" " . $temps[$nn];
	    }
	    print "     $temps[$nn] ";
	}
	# # Get lows
	# #get the high temp
	# (my $high) = /High:\s*(\d+)/;
	# #get the low temp
	# (my $low) = />Low: \-?(\d+)/;
	
	# ## deal with 1 or 2 digit temps
	# if (length($high)<2) {
	#     $high=" ".$high;
	# }
	# if (length($low)<2) {
	#     $low=" ".$low;
	# }
	# #print the high and low temp for the specified day
	#     if (/class=\"fiveday-temps.>(.+)<div\sid=/) { 
	# 	$kk=$1;
	# 	$kk=~s/<.+?>/ /g;
	# 	$kk=~s/^\s+//;
	# 	$kk=~s/\&.+?\;//g;
	# 	#print "kk : -$kk-\n";
	# 	@temps=split(/\s+/,$kk);
	# 	@temps=grep(/\d+/,@temps);
	# 	print "   ";
	# 	print "\n";
	#     } elsif ($high=~/\d+/ && $low=~/\d+/ && ++$count2<=$day) {
	# 	$htemp.=$hleadspace.$high.$degree.$htrailspace;
	# 	$ltemp.=$lleadspace.$low.$degree.$ltrailspace; 
		
	#     }
	#     #			elsif($high=~/\d+/ && $low=~/\d+/ && ++$count2<=$day){$ttext.=$leadspace.$high.$degree."\n".$low.$degree.$trailspace ."\n";  }
	#     elsif ($count2>=$day) {
	# 	print "$htemp\n  $ltemp\n$dtext"; exit;
	#     }	      #don't keep lopking if everything has been found
	    
	#}
    }
    when (/l/)
    {
	print "$LOC\n";
    }
    when (/C/)
    {
	&file_op;		#save weather to $file
	while (<FILE>) {
	    if (/forecast-icon(.*)/) {
		$1=~/src=\'(.+\/(.+?))\'/ || die("shit\n");
		print STDERR  "$1 :::::::::: $2\n";
		print STDERR "wget $1; mv $2 ~/conkyimages/current.png ;convert ~/conkyimages/current.png ~/conkyimages/current.gif\n";
		system("wget $1; mv $2 ~/conkyimages/current.png ;convert ~/conkyimages/current.png ~/conkyimages/current.gif");

	    }
	}
    }

}


sub file_op {		     #do file operations
    if (-e $file ) {	     #does the file exist and is it not empty?
	my $size=`stat -c %s $file`;
	if ($size >= 1000) {
	    my $date=`date -u +%s`; #get current date in seconds
	    my $created=`stat -c %Y $file`; #get creation date of file in seconds
	    $age=$date - $created;	    #determine age of file
	} else {
	    $age=$update+1;
	}
    } else { #if file doesn't exist make it and set to update the file
	`touch $file`;
	$age=$update+1;
    }
    $a=`/sbin/ifconfig | grep Bcast |wc`;
    #only get a new file every hour and only if connected
    if (($age>=$update) && ($a != 0)) { 
	#obtain the weather forecast and store it in $file
	#	`wget -O - http://weather.yahoo.com/$code1/?unit=c > $file`;
	`wget -O - http://www.worldweather.org/$code1 > $file`;
    }
    open(FILE, $file) or die "Could not open file $file: $!\n";
    @lines=<FILE>;
    close(FILE);
    open(FILE, $file) or die "Could not open file $file: $!\n";

}

sub usage {	   #if correct options haven't been passed usage error
    print "Usage error weather.pl <citycode> <system> <option>\n";
    print "weather.pl [lc] <option>\n";
    print "\t<citycode> : Location (bcn|ath|mar)\n";
    print "\t<system>   : c for metric or f for imperial\n";
    print "\t<option>   : Only one option can be entered at a time\n";
    print "\t\tc displays current conditions\n";
    print "\t\tw displays list of current conditions\n";
    print "\t\tcp displays current conditions symbol\n";
    print "\t\tt displays current temp in chosen system\n";
    print "\t\t[1-5]d displays the days up to specified day\n"; 
    print "\t\t[1-5]dp displays condition symbol for days up to specified day\n";
    print "\t\t[1-5]t displays high/low temp in chosen system up to specified day\n";
    print "\t\t[1-5]dt displays days and then high/low temp in chosen system up to specified day\n";
    print "\t\t[1-7]w displays the weather in words up number specified\n";
    print "\t\t[1-5]p displays conditions for specified day\n";
    print "\t\t[1-5] displays high/low temp in chosen system for specified day\n";
}



sub cond_symb {	    #translates conditions into symbol in weather font
    # if ($_ =~ "Partly Cloudy"){$_="c";}
    # elsif($_ =~ /Mostly\sSunny/ || $_ =~ /Mostly\sClear/ ){$_="B"}
    # elsif ($_ =~ /Mostly\sCloudy/){$_="d";}
    # elsif ($_ =~ "Cloud" || $_ =~ "Fog"){$_="e";}
    # elsif ($_ =~ "Storm" || $_ =~ "Thunder" || $_ =~ "T-"){$_="i";}
    # elsif ($_ =~ "Snow" || $_ =~ "Flurries" || $_ =~ "Wintry"){$_="k";}
    # #elsif(/Shower/ && /Wind/){$_="f"}
    # elsif ($_ =~ "Fair"  || $_ =~ "Sunny"|| $_ =~ "Sun" || $_ =~ "Clear"){$_="A";}
    # elsif(/Scattered\sShowers/){$_="g"}
    # elsif ($_ =~ /Rain/ || /Drizzle/){$_="i";}
    # elsif ($_ =~ /Shower/){$_="h"; }
    if ($night && $_ =~/Clear/) {
	$_="~/conkyimages/moon.gif";
    } elsif ($night && $_ =~/Fair/) {
	$_="~/conkyimages/mclouds1.gif";
    } elsif ($night && $_ =~/Cloudy/) {
	$_="~/conkyimages/mclouds2.gif";
    } elsif ($_ =~ "Partly Cloudy") {
	$_="~/conkyimages/partly_cloudy.gif";
    } elsif ($_ =~ /Wind/ ) {
	$_="~/conkyimages/wind.gif";
    } elsif ($_ =~ /Mostly\sSunny/ || $_ =~ /Mostly\sClear/ ) {
	$_="~/conkyimages/mostly_sunny.gif";
    } elsif ($_ =~ /Mostly\sCloudy/) {
	$_="~/conkyimages/mostly_cloudy.gif";
    } elsif ($_ =~ "Cloud" || $_ =~ "Fog") {
	$_="~/conkyimages/cloudy.gif";
    } elsif ($_ =~ "Storm" || $_ =~ "Thunder" || $_ =~ "T-") {
	$_="~/conkyimages/storm.gif";
    } elsif ($_ =~ "Snow" || $_ =~ "Flurries" || $_ =~ "Wintry") {
	$_="~/conkyimages/snow.gif";
    }
    #elsif(/Shower/ && /Wind/){$_="f"}
    elsif ($_ =~ "Fair"  || $_ =~ "Sunny"|| $_ =~ "Sun" || $_ =~ "Clear") {
	$_="~/conkyimages/sunny.gif";
    } elsif (/AM\sShowers/) {
	$_="~/conkyimages/am_showers.gif";
    } elsif (/Scattered\sShowers/) {
	$_="~/conkyimages/showers.gif";
    } elsif ($_ =~ /Rain/ || /Drizzle/) {
	$_="~/conkyimages/rain.gif";
    } elsif ($_ =~ /Showers/) {
	$_="~/conkyimages/showers4.gif";
    }
    $ctext.="\${image " . $_ . " -p $offset,755}\n";
    my $imgname="image" . $offset . ".gif";
    $offset++;
    system("rm ~/conkyimages/$imgname; ln -s $_ ~/conkyimages/$imgname");
    print STDERR "rm ~/conkyimages/$imgname; ln -s $_ ~/conkyimages/$imgname\n";
    $night=1;
}

sub day_space {			#Adds spaces for aligment
    if ($_ =~ "Today") {
	$_="  Tdy  ";
    }
    # 	elsif ($_ =~ "Tonight"){$="Tonight  ";}
    # 	elsif ($_ =~ "Tomorrow"){$_=" Tomorrow";}
    elsif ($_ =~ "Tonight") {
	$_="  Tng  ";
    } elsif ($_ =~ "Tomorrow") {
	$_="  Tmr  ";
    } elsif ($_ =~ "Thu") {
	$_="   Thu  ";
    } elsif ($_ =~ "Fri") {
	$_="   Fri  ";
    } elsif ($_ =~ "Sat") {
	$_="   Sat  ";
    } elsif ($_ =~ "Sun") {
	$_="   Sun  ";
    } elsif ($_ =~ "Mon") {
	$_="   Mon  ";
    } elsif ($_ =~ "Tue") {
	$_="   Tue  ";
    } elsif ($_ =~ "Wed") {
	$_="   Wed  ";
    }
    $dtext.=$_.$dspaces;

}


