#!/usr/bin/perl -w

my $pid=`ps aux | egrep 'vuze|azureus'  | gawk '\$11~/java/{print \$2}'` || 0;
 chomp($pid);

$ARGV[0]="-h" unless defined $ARGV[0];
my $file=$ARGV[1] || "";
#$pid == 0 && do{print "Vuze is not running\n"; exit();};
if($ARGV[0] eq "status" || $ARGV[0] eq "-s" || $ARGV[0] eq "--status"){
    if($pid == 0){
	print STDOUT "Vuze is not running\n";
	exit(0);
    }
    else{
	 print STDOUT "Vuze is running with PID: $pid\n";
	 exit(0);
    }
}
elsif($ARGV[0] eq "kill" || $ARGV[0] eq "-k" || $ARGV[0] eq "--kill"){
    if($pid == 0){
	print STDERR "Vuze is not running.\n" 
    }
    else{
	print STDERR "Killing vuze...\n"; 
	`kill -9 $pid`;
    }
    exit();
}
elsif($ARGV[0] eq "force" || $ARGV[0] eq "-f" || $ARGV[0] eq "--force"){
    if($pid == 0){
	print STDERR "Vuze is not running.\n" 
    }
    else{
	print STDERR "Killing vuze...\n"; 
	`kill -9 $pid`;
    }
    exit();
}
elsif($ARGV[0] eq "restart" || $ARGV[0] eq "-r" || $ARGV[0] eq "--restart"){
  $pid == 0 ?  print STDERR "Vuze is not running, launching now...\n" : print STDERR "Killing vuze...\n"; 
  if($pid != 0) { `kill $pid`;   sleep 10;};
    system("vuze $file&");
}
elsif($ARGV[0] eq "help" || $ARGV[0] eq "-h" || $ARGV[0] eq "--help"){
    usage();
}
elsif($ARGV[0] eq "remote" || $ARGV[0] eq "-R" || $ARGV[0] eq "--remote"){
    $pid == 0 ?  print STDERR "Vuze is not running, launching now...\n" : print STDERR "Killing vuze...\n"; 
  if($pid != 0) { `kill $pid`;   sleep 10;};
     system("export DISPLAY=:0.0; vuze $file&");
}
else{
    usage();
}

sub usage{
 print <<EndPrint;
  manage_vuze.pl status|kill|force|restart|remote|help

 Options can be given as -h, -help or --help

 -f : Force vuze to exit (kill -9).
 -k : Kill currently running instance of vuze (if any).
 -r : Kill currently running instance of vuze (if any) and launch a new one displayed on the local computer.
 -R : As above but will kill currently running instance of vuze (if any) and launch a new one displayed on the remote computer.
 -s : Check if vuze is running.


EndPrint
    exit();
}
