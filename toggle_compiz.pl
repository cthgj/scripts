#!/usr/bin/perl -w




my $pid=`ps aux | grep metacity  | gawk '\$11~/metaci/{print \$2}'` || `ps aux | grep -w enlightenment  | gawk '\$11~/enlightenment/{print \$2}'` || `ps aux | grep xfce4-session  | gawk '\$11~/xfce4/{print \$2}'`  || 0;
chomp($pid);

if($pid != "0"){  
    print STDERR "Killing window manager, starting compiz ...\n";
#    system("kill $pid; compiz --replace&");
    system("compiz --replace&");
    system("gtk-window-decorator --replace&");
    system("sleep 5; killall gnome-panel");
    exit();

}
open(A,"ps aux | grep compiz |");
$pid=0;
while(<A>){$pid=1 if /\scompiz\s/ || /--replace/;}

if($pid != "0"){
    print STDERR "Killing compiz, starting metacity...\n";
    system("killall compiz.real; metacity --replace &");
    system("gtk-window-decorator --replace&");
    exit();
}
#print "pp : -$pid-\n";
