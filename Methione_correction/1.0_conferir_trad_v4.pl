#!usr/bin/perl

###21-12-2020
### vinicius parreira
## version 04 : generate a multifasta file
# This Script compares Spliceprot and Swissprot translation.
$arquivo1 = $ARGV[0]; #SpliceProt translation
$arquivo2 = $ARGV[1]; #SwissProt sequence file
$arquivo3 = $ARGV[2]; #Parser file tha make the link between ENST and SwissProt ID 

open (IN, "$arquivo1") || die "nao abriu $arquivo1";
open (IN2, "$arquivo2") || die "nao abriu $arquivo2";
open (IN3, "$arquivo3") || die  "nao abriu $arquivo3";


%hash_parser = ();

#print"====== Lendo $arquivo3 ======\n";

while (<IN3>) {
    chomp $_;
    $line_parser = $_;
    next if $line_parser =~ /^Entry name/;
    @array_parser = split (/\t/, $line_parser);
#    print "$array_parser[6]\n";
    $ac_parser = $array_parser[0];
    @isoformas = split (/\;/, $array_parser[6]);
    $num_isoformas = scalar @isoformas;
    for ($a=0; $a<$num_isoformas; $a++){
	if ($isoformas[$a] !~ /\[/g){
	   $hash_parser{$isoformas[$a]} = $ac_parser."-1";
	}else{
	    $isoformas[$a] =~ /(.*?) \[(.*?)\]/;
	    $enst_parser = $1;
	    $idd_parser  = $2;
	    $hash_parser{$enst_parser} = $idd_parser;
	    
	}
    }
}

=c # imprime o arquivo parser de dicionario ENST  vs Swissprot ID
foreach $chave (keys %hash_parser){
    print "$chave;$hash_parser{$chave}\n";
}
=cut


%hashdef = (); # hash definitivo recebe trasncrito e sequencia
%hashseqSw = (); # hash temporario, avalia cada entrada do uniprot 

#print"====== Lendo $arquivo2 ======\n";

$seq_sw = "";
while (<IN2>){
    chomp($_);
    $lineSw = $_;
    if ($lineSw =~ />/) {
	$lineSw =~ />sp\|(.*?)\|(.*?)/;
	$ac_sw = $1;
	if ($ac_sw !~ /\-/){
	    $ac_sw =  $ac_sw."-1";
	}
	    
#	print"$lineSw\t$ac_sw\n"; # ---debug ok
	unless ($before_ac eq ""){
	$hashseqSw{$before_ac} = $seq_sw;
#	print "$before_ac\n"; # --- debug ok
	$seq_sw = "";
	}
    }else{
	$before_ac = $ac_sw;
	$seq_sw = $seq_sw.$lineSw;
    }

}


=c # imprime o AC e a sequencia da proteina do uniprot, canonica e alternativa
foreach $chave2 (keys %hashseqSw){
    print "$chave2;$hashseqSw{$chave2}\n";
}
=cut


%hash_my=();

#print"====== Lendo $arquivo1 ======\n";
while (<IN>){
    chomp($_);
    $line_my = $_;
    if ($line_my =~ /^>/){
	@arraymy = split(/\[/);
	$tamanho = scalar @arraymy;
#	print"$tamanho\n";
    }else{
	$seqmy = $line_my;
       	for ($i = 1; $i<=$tamanho ;$i++) {

	    $arraymy[$i] =~ s/\]//;
	    $ENSTmy = $arraymy[$i];
	    $hash_my{$ENSTmy} = $seqmy 
	}
    }
}
#print "acabou o segundo arquivo\n";

foreach $keys (keys %hash_my){
    $keys2 = $hash_parser{$keys}; # repsectiva Ac do Swissprot para o ENST correspondente

    if ($hash_my{$keys} eq $hashseqSw{$keys2}){
	print"[$keys][$keys2]:My_dataset vs SwissProt = TRU3\n>[spliceprot]\n$hash_my{$keys}\n>[swissprot]\n$hashseqSw{$keys2}\n";
	print "----------\n";
    }elsif ($hash_my{$keys} =~ $hashseqSw{$keys2} and $hashseqSw{$keys2} =~ /\w/g){
	print "[$keys][$keys2]:My_dataset vs SwissProt = Swiss_E_my\n>[spliceprot]\n$hash_my{$keys}\n>[swissprot]\n$hashseqSw{$keys2}\n";
    print "----------\n";
    }elsif ($hashseqSw{$keys2} =~ $hash_my{$keys} and $hashseqSw{$keys2} =~ /\w/g){
        print "[$keys][$keys2]:My_dataset vs SwissProt = My_E_swiss\n>[spliceprot]\n$hash_my{$keys}\n>[swissprot]\n$hashseqSw{$keys2}\n";
    print "----------\n";
    }elsif ($hashseqSw{$keys2}){
	print"[$keys][$keys2]:My_dataset vs SwissProt = FALS3\n>[spliceprot]\n$hash_my{$keys}\n>[swissprot]\n$hashseqSw{$keys2}\n";
    print "----------\n";
    }else{
	print"[$keys][$keys2]:My_dataset vs SwissProt = NON3\n>[spliceprot]\n$hash_my{$keys}\n>[swissprot]\n$hashseqSw{$keys2}\n"; 
    }
    print "----------\n";
}
