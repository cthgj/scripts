#!/bin/bash
## This script should be run in a folder containing the various results of moonlight candidate
## analyses, it expects to find a file named <species>.<mode>.<stats>, e.g.: human.i.stats in 
## the various subfolders. It will print a latex table.



## LaTeX:
#\usepackage{booktabs}
#\usepackage{color}
# \definecolor{green}      {cmyk}{0.74,0.00,0.76,0.28}
# \definecolor{red}      {cmyk}{0,0.93,0.93,0.22}
# \def\grAr{\large{\textcolor{green}{$\Uparrow$}}}
# \def\reAr{\large{\textcolor{red}{$\Downarrow$}}}
# \def\grar{\large{-}}
# \def\rear{\large{-}}

cat<<EOF
\begin{table}[htb]
\begin{centering}
\begin{tabular}{lccc}
\toprule
 & \multicolumn{1}{c}{\textbf{Multi NC}} & \multicolumn{1}{c}{\textbf{Hubs}} & \multicolumn{1}{c}{\textbf{Network}} \\\\
EOF
for i in `find  . -name "$1*.$2.*stats" | grep -v old`; do 
    foo=`echo "$(basename $(dirname $i))" | sed 's/_/ /g'`; 
    echo -n "${foo^}"
    for k in NC Hubs Net; do 
	cands=`grep $k $i | gawk '{print $(NF-1)}'` 
	group=`grep $k $i | gawk '{print $(NF)}'`; 
	pval=`grep $k $i | gawk '{print $(NF-2)}'`; 

	## If the p-val is significant
	test=`echo $pval | gawk '{if($1<=0.05){print "1\n"}else{print "-1"}}'`
	if [ $test == 1 ]; then
	    ## If the cands' mean is greater than that of the current group
	    test=`echo "$cands $group" | gawk '{if($1>$2){print "1\n"} else if($1<$2){print "-1"} else{print "0"}}'`
#	    echo -e "\nTEST ($k) is $test for  $cands > $group";
	    if [ $test == 1 ]; then
		ar="\grAr";
	    elif [ $test == -1 ] ; then  
		ar="\reAr";
	    else
		ar="CHECK";
	    fi

        ## if the p-value is not significant
	else
            ## If the cands' mean is greater than that of the current group
	    test=`echo "$cands $group" | gawk '{if($1>$2){print "1\n"} else if($1<$2){print "-1"} else{print "0"}}'`
	    if [ $test == 1 ]; then
		ar="\upar";
	    elif [ $test == -1 ] ; then  
		ar="\dnar";
	    else
		ar="CHECK";
	    fi

	fi
	echo -n " & $ar "
    done
    echo "\\\\"
done

cat<<EOF
\bottomrule
\end{tabular}
\caption{}
\end{centering}
\end{table}

EOF
