#!/usr/bin/perl -w

use strict;
#my @temps=split(/\s+/,`sensors| grep Core | gawk '{i=i+substr(\$3,2,length(\$3)-3); }END{print i/2,substr(\$6,2,length(\$6)-4),substr(\$9,2,length(\$9)-4)}'`);
my @temps=split(/\s+/,`sensors| grep Core|sensors| grep Core | gawk '{print substr(\$3,2,length(\$3)-3) }END{print substr(\$6,2,length(\$6)-4),substr(\$9,2,length(\$9)-4)}'`);

#$temps[0]:cpu1 
#$temps[1]:cpu2
#$temps[2]:max 
#$temps[3]:crit
my $temp;
$temps[0] > $temps[1] ? ($temp=$temps[0]) : ($temp=$temps[1]);
if ( $temp >= 2*($temps[2]/3)+5 ){
    print  "\${color darkorange}$temp\${color}";

}
elsif ( $temp >= $temps[2] ){
    print  "\${font Arial:bold:size=8}\${color red}$temp\${color}\${font}";
}
elsif ( $temp >= $temps[3]-5 ) {
    print  "\${font Arial:bold:size=18}\${color firebrick}$temp\${color}\${font}";
}
else{
    print  "\${color white}$temp\${color}";
#    print  "\${font Arial:bold:size=8}\${color white}$temp\${color}\${font}";

}

#my @a=($tt=~/\+(\d+)/g); print "@a\n";' | gawk '{k=k+$1}rint k/2"\n",$2"\n",$3"\n"}'`

# ;





############### Old Version ###########################
# temps[0]="THM0";
# temps[1]="THM1";
# a=0;
# for (( i=0;i<2;i++)); do
#     if (( $i == 1 )); then
# 	echo -ne "/"
#     fi
#     name=${temps[${i}]}
#     temp=`cat /proc/acpi/thermal_zone/"$name"/temperature|gawk '{print $2}'`
#     if (( $temp > 79 )); then
# 	a=1
# 	echo -ne "\${color red}$temp\${color}";
#     else
# 	echo -ne "$temp";
#     fi
# done 

# exit
