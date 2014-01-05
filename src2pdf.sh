#!/usr/bin/env bash


cat<<EOF >/tmp/all.tex   ## Print the tex file header
\documentclass{article}
\usepackage{listings}
\usepackage[usenames,dvipsnames]{color}
\lstdefinestyle{customasm}{
  belowcaptionskip=1\baselineskip,
  xleftmargin=\parindent,
  language=C++,
  breaklines=true, %%%%%  This is useful if you have very long lines in your code
  basicstyle=\footnotesize\ttfamily,
  commentstyle=\itshape\color{Gray},
  stringstyle=\color{Black},
  keywordstyle=\bfseries\color{OliveGreen},
  identifierstyle=\color{blue},
  xleftmargin=-8em,
}

\usepackage[colorlinks=true,linkcolor=blue]{hyperref} 
\begin{document}
\tableofcontents

EOF


find . -type f ! -regex ".*/\..*" ! -name ".*" ! -name "*~" | 
sed 's/^\..//' |                 ## Change ./foo/bar.src to foo/bar.src
while read  i; do                ## Loop through each file
    echo "\newpage" >> /tmp/all.tex   ## start each section on a new page
    echo "\section{$i}" >> /tmp/all.tex  ## Create a section for each file
    echo "\lstinputlisting[style=customasm]{$i}" >> /tmp/all.tex ## Include the file

done &&
echo "\end{document}" >> /tmp/all.tex &&
pdflatex /tmp/all.tex -output-directory . && pdflatex /tmp/all.tex -output-directory .

