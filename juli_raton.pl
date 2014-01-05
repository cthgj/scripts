#!/usr/bin/perl -w

unless($ARGV[0]){
    print STDERR<<ENDofHELP;

Notas: '> todos.txt' quiere decir que ponga la salida del juli_raton en el archivo todos.txt

Caso 1: Tienes UNA carpeta con muchos FICHEROS de datos. La carpeta es To_Exp_juli/BL6/ADQ

 pasos: cd hoy/NOMBRE_DE_CARPETA
        juli_raton.pl * > todos.txt


ENDofHELP
exit;
}


my $name;
my (%mm,%data,%totals,%data2);

my $outfile1="todos.txt";
my $outfile2="sumas.txt";

open(OUT,">$outfile1");
open(OUT2,">$outfile2");
while(<>){
    next if /^Point/;

    if(/^..(.{6}).+\.E/){
	$name=$1;
	$mm{$name}++;	
    }
    elsif(/^totals\s*(\d+)/i){
	$totals{$name}=$1;

    }
    else{
	chomp;
	/(\d+)\s+(\d+)/;
	my @a;
	$data{$1}{$name}=$2;
#	$data2{$name}{$1}=$2;

    }
}

my @mice=sort(keys(%mm));
my @oo=keys(%data);
#print "x kk :  @oo  \n"; die();
foreach my $time (sort{ $a <=> $b }keys(%data)){
    foreach my $mouse (@mice){
	$data{$time}{$mouse}='n/a' unless defined($data{$time}{$mouse});
    }
}
foreach my $mouse (@mice){
	print OUT "\t$mouse";
}
print OUT "\n";
my @times=sort{ $a <=> $b } keys(%data);
foreach my $time (@times){
    print OUT "$time";
    foreach my $mouse (@mice){
	print OUT "\t$data{$time}{$mouse}";
    }
    print OUT "\n";

}
    print OUT "Tot\t";
foreach my $mouse (@mice){
    print OUT "$totals{$mouse}\t";
}
print OUT "\n";

foreach my $mouse (@mice){
	print  OUT2 "\t$mouse";
}
print OUT2 "\n";

#    my @times2=sort{ $a <=> $b }keys(%{$data2{$mouse}});
for (my $n=1; $n<scalar(@times); $n=$n+2){
    print OUT2 "$times[$n-1]+$times[$n]";
    foreach my $mouse (@mice){
	$data{$times[$n]}{$mouse}=0 if $data{$times[$n]}{$mouse}=~'n/a';
	$data{$times[$n-1]}{$mouse}=0 if $data{$times[$n-1]}{$mouse}=~'n/a';
	my $sum=$data{$times[$n]}{$mouse} + $data{$times[$n-1]}{$mouse};
	print OUT2 "\t$sum "  ;
    }
    print OUT2 "\n";
}
