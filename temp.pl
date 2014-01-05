#!/usr/bin/perl -w

%ucnt=(); 
for (system("who")){ 
    s/\s.*\n//; 
    $ucnt{$_}++;
} 
@users = sort keys %ucnt; 
print "users: @users\n"
    
