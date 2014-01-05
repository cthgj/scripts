#!/usr/bin/perl 

use Date::Calc qw(Add_Delta_DHMS);

my $pepito=0;
print  "pepito tiene valor : $pepito\n";
my ($stime,$start_time,$start_date,$etime, $mouse_number);
my %values;
my @medidas=("Activ.","Stere.","Locom.","P.Time","PT","Res.T.","RT","Mov.S.","MS","Mov.F.","MF","N.Rea.","M.Rea.");

while(<>){

    if(/Subject\s+Track\s+Number.+(\d+)/){
	$mouse_number=$1;
	print "===========================================\n";
	}

    if(/^Track\s+Date.+?:\s+(.+)/){
	($start_date,$start_time)=split(/\s+/,$1);
	print "time date : $start_time, $start_date\n";
    }

    if (/From/){
	$pepito=$pepito+1 ;
	/From\s+(.+?)\s+to\s+(.+)/;
	$stime=$1; 
	$etime=$2;	
    }

    if($pepito>1){
	if (/TOTAL/){
	    my @valores=split;
	    for(my $n=0; $n<scalar(@medidas); $n++){
		$values{$mouse_number}{$medidas[$n]}=$valores[$n+1];

	    }
	    print "MMM : $mouse_number\n";
	  #  map{print "$_ : $values{$mouse_number}{$_}\n"}keys(%{$values{$mouse_number}});
	    
	    &day_night();
	    # for(my $n=0; $n<scalar(@a); $n++){
# 		print "$n : $a[$n]\n";

# 	    }
# 	    die();
	    
	}


    }



}



sub day_night{
    my $day=0;
 
 my ($hour, $minute, $second)=split(/:/,$start_time);
    my $days_offset=0;
    ($day, $month, $year) = split(/\//,$start_date);
    print "aaa : $year, $month, $day\n";    
    my ($hour_offset, $minute_offset, $second_offset)=split(/:/,$stime);
    $stime=~/^(\d+):/;
    my $start=$1+$hour;
        
    my ($year2, $month2, $day2, $STIME, $m2, $s2) = Add_Delta_DHMS( $year, $month, $day, $hour, $minute, $second,
                $days_offset, $hour, $minute, $second );

    ($year2, $month2, $day2, $h2, $m2, $s2) = 
    Add_Delta_DHMS( $year, $month, $day, $hour, $minute, $second,
                $days_offset, $hour_offset, $minute_offset, $second_offset );
    
#    print STDOUT "aa $hour:$minute:$second\n";
    my @aa= ($hour,$minute,$second,$STIME);
 #   print STDOUT "a $hour, $minute, $second $year2, $month2, $day2, $h2:$m2:$s2\n";

    
    my $realtime=$h2+$stime;
    print "RE : $realtime\n";
    if($realtime>8

    return($day);

    die() if $pepito>7;
    
}
