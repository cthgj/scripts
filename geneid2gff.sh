#!/bin/sh


gawk '
BEGIN{ OFS="\t" }

(/\# S/){ k=$3 } 
($5 == "-" || $5 == "+"){
    print k, "geneid", $1, $2, $3, $4, $5, $6, "."
} ' "$@"
