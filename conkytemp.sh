#!/bin/bash

temps=(`sensors| grep Core |perl -ne '@a=(/\+(\d+)/g); print "@a\n";' | gawk '{k=k+$1}END{print k/2"\n",$2"\n",$3"\n"}'`);
echo "${temps[0]} ${temps[1]} ${temps[2]}"
if [[ ${temps[0]} > ${temps[1]} ]]; then
    echo -e "\${color red}${temps[0]}\${color}";
elif [[ ${temps[0]} > ${temps[2]} ]]; then
    echo -e "\${color red}${temps[0]} (Critical)\${color}";
else
    echo -e "${temps[0]}";
fi





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
