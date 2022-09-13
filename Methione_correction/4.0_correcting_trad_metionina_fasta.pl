#!usr/bin/perl

###2-12-2020
### vinicius parreira
## This script replace in the FASTA file the sequences that have the first methione 
#  corrected based on SwissProt

$arquivo1 = $ARGV[0]; #corrected file
$arquivo2 = $ARGV[1]; #FASTA file

open (IN, "$arquivo1") || die "nao abriu $arquivo1\n"; # abre o arquivo
open (IN2, "$arquivo2") || die "nao abriu $arquivo2\n"; # abre o arquivo

$/=">"; # muda o input separator para > por ser fasta

$count = 0;
%hash_correcao = ();

while (<IN>){ # le linha por linha 
    $count++;
    next if $count <=1;    
    chomp $_;
    $line = $_;
    @array1 = split (/\n/, $line);
    $array1[0] =~ s/_novaseq//;
    $hash_correcao{$array1[0]} = $array1[1];
#    print "$array1[1]\n"; # debug ---ok
}
$count2 = 0;
while (<IN2>){ # le linha por linha
    $count2++;
    next if $count2 <=1;
    chomp $_;
    $line2 = $_;
    @array2 = split (/\n/, $line2);
    @ENST = split (/\]/, $array2[0]);
    $teste =0;
    foreach $a (@ENST){
	$a =~ s/\[//g;
	if ($hash_correcao{$a}){
	    print ">$array2[0]\*\n";
	    print "$hash_correcao{$a}\n";
	    $teste = 1;
	    last;
	}
    }
    if ($teste == 0){
	print ">$array2[0]\n";
	print "$array2[1]\n";
    }
}
