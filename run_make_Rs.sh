#!/bin/bash
## cands, nc
# declare -a dirs=(betweeness  degree  disopred domains length shortest_paths classes isoforms);
# declare -a titles=(Betweeness Degree Disorder Domains Length Shortest_Paths Classes Isoforms)
# declare -a exts=(betnums degnums perc domnums lennums paths clascount count)
 declare -a dirs=(betweeness  degree  disorder  domains  length  classes);
# declare -a titles=(Betweeness Degree Disorder Domains Length Classes)
 declare -A titles;
titles[betweeness]="Betweeness"
titles[degree]="Degree"
titles[disopred]="Disorder"
titles[domains]="Domains"
titles[length]="Length"
titles[classes]="Classes"
titles[shortest_paths]="Shortest_Paths"


declare -A dirnames;
dirnames[betweeness]="betweeness"
dirnames[degree]="degree"
dirnames[disopred]="disopred"
dirnames[domains]="domains"
dirnames[length]="length"
dirnames[classes]="classes"
dirnames[shortest_paths]="paths"

declare -A exts;
exts[betweeness]="betnums"
exts[degree]="degnums"
exts[disopred]="perc"
exts[domains]="domnums"
exts[length]="lennums"
exts[classes]="clascount"
exts[shortest_paths]="count"

species=$1;
mode=$2;
declare -a files=($species.$mode.cands $species.hub $species.$mode.multi_noncands $species.multi $species.$mode.non $species.mono $species $species.$mode)

#echo ${Unix[@]:3:2}
aa=(${@});
## If more than 2 arguments are given, take the rest as dir names
if (( $# > 2 )); then
    dirs=();
    dirs=${aa[@]:2:$#}
fi
echo "Dirs : ${dirs[@]}" 1>&2

echo "Building boxplots for $species $mode" 1>&2

want_phospho=0;
for name in ${dirs[@]}
do
    ## Check if we want to run the phosphorylation data
    if [ $name -eq "phospho" ]
    then
	want_phospho=1;
    fi
    echo "-------------------- $name  $species $mode ----------------------------";
    cd $name;
    com="make_Rs.pl -m $mode -t \"${titles[$name]}\" ";
    echo "$com" 1>&2;
    for f in ${files[@]}
    do
	com+=" $f.${exts[$name]}"
    done
    com+=".box.png"
    echo "$com";
    $com > $species.$mode.box.R && Rscript $species.$mode.box.R
    cd - >/dev/null

done

echo "Building density plots for $species $mode" 1>&2
for name in ${dirs[@]}
do
    cd $name;
    com="make_Rs.pl -m $mode -dt \"${titles[$name]}\" ";
    for f in ${files[@]}
    do
	com+=" $f.${exts[$name]}"
    done
    com+=".distr.png"
    echo $com;
    $com > $species.$mode.distr.R && Rscript $species.$mode.distr.R
    cd - >/dev/null
done

### Special case for the phosphorylation data
if [[ $want_phospho == 1 ]]
then
    echo "-------------------- Phospho $species $mode----------------------------" 1>&2;
    cd phospho
    echo "cd phospho; make_Rs.pl -m $mode -pt \"Phosphorylation\"  $species.$mode.cands.ph $species.hub.ph $species.$mode.multi_noncands.ph $species.multi.ph $species.$mode.non.ph $species.mono.ph $species.all.stats ";
 
    make_Rs.pl -m $mode -pt "Phosphorylation" $species.$mode.cands.ph $species.hub.ph $species.$mode.multi_noncands.ph $species.multi.ph $species.$mode.non.ph $species.mono.ph $species.all.stats  > $species.$mode.R; Rscript $species.$mode.R; 
    cd -
fi
