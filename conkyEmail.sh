#!/bin/bash

`$HOME/scripts/conkyEmail.py --servertype=IMAP --servername=imap.googlemail.com -u karolos.chapple -p asemenampw --ssl >$HOME/.conky_emails 2>>$HOME/.conky/.conky_emails`
emails=`wc $HOME/.conky_emails | gawk '{print $1}' `
if (( $emails == 1 )); then
    emails=`cat $HOME/.conky_emails`
    ## If there is no internet connection, conkyEmail.py
    ## will return '?'. At the moment, I have no internet 
    ## connection. Therefore, I cannot search online for how
    ## the FUCK I can match '?' in bash. So, grep to the rescue.
    aa=`grep "?" .conky_emails |wc -l`;
    if (( $aa == 1)); then
	exit;
    fi
    if (( $emails == 1 )); then
	echo "\${color orange}$emails\${color white} email"
    
    elif (( $emails > 1 )); then
	echo "\${color orange}$emails\${color white} emails"
    fi

fi
