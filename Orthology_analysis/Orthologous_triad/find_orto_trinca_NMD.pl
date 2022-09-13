#!usr/bin/perl
#vinicius parreira
#24-05-22
# encontra trinca de ortologos que foram encontrados pelo pattern lab pelo fasta de NMD


$arq1 = $ARGV[0]; # fasta NMD
=c
>Hs.100836_vn.7_::_vn.2  [ENST00000553960]  [ENST00000556809]
MAAAAAAAAAAGAAGGRGSGPGRRRHLVPGAGGEAGEGAPGGAGDYGNGLESEELEPEELLLEPEPEPEPEEEPPRPRAPPGAPGPGPGSGAPGSQEEEEEPGLVEGDPGDGAIEDPVREGGRASRPAAGRPGHWRPRARASGWQAGGGVGRGITWLGPGRAGDGSAITTRGPTGLIRASRLPSVVLERVAFLLSRPSHGAREMAEHG
>Hs.168090_vn.7  [ENST00000468499]
MAAAAAAAAATNGTGGSSGMEVDAAVVPSVMACGVTGSVSVALHPLVILNISDHWIRMRSQEGRPVQVIGALIGKQEGRNIEVMNSFELLSHTVEEKIIIDKEYYYTKEEQCV
=cut


$arq2 = $ARGV[1]; # lista de trincas ortologas
=c
";"human_id";"mouse_id";"rat_id"
"1";"ENSG00000188976-vn.1";"ENSMUSG00000095567-vn.3";"ENSRNOG00000021392-vn.1"
"2";"ENSG00000187961-vn.4";"ENSMUSG00000078485-vn.1_vn.4";"ENSRNOG00000020302-vn.1"
"3";"ENSG00000187642-vn.2";"ENSMUSG00000078486-vn.1";"ENSRNOG00000020244-vn.1"
"4";"ENSG00000187608-vn.2";"ENSMUSG00000035692-vn.1";"ENSRNOG00000021802-vn.1"
"5";"ENSG00000188157-vn.4";"ENSMUSG00000041936-vn.7";"ENSRNOG00000020205-vn.2"
"6";"ENSG00000131591-vn.9_vn.2";"ENSMUSG00000059939-vn.5_vn.6";"ENSRNOG00000020199-vn.1"
=cut

$sp = $ARGV[2]; # Mm Rn Hs

open (IN, "$arq1") || die "nao abre $arq1";

%hash_parser =();
%hash_parser_ensts= ();
    
while (<IN>){
    chomp $_;
    $line = $_;
    $line =~ s/\_\:\:\_/\_/g;
    $line =~ s/\"//g;
    $line =~ s/\]/ /g;
    $line =~ s/\[/ /g;
    $line =~ s/\>//g;
    @array1 = split (/\s+/, $line);
    

    foreach $a (@array1){
	if (($a =~ /^Hs/) || ($a =~ /^Mm/) || ($a =~ /^Rn/)){
	   
	    $hash_parser{$a} = 1;
#	    print "$a $hash_parser{$a}\n";
	    
	}else {next;}
    }
  
}


open (IN2, "$arq2") || die "nao abre $arq2\n";

while(<IN2>){
    chomp $_;
    chop $_;
    
    $line2 = $_;
    $tmp = $line2;
    $tmp =~ s/"//g;
#    $tmp =~ s/0{1,6}//; ##### resolver isso 
    $tmp =~ s/\_\:\:\_/\_/g;
    $tmp =~ s/\-/_/g;
    @array3 = split (/\;/, $tmp);
    #print "$tmp\n";

    if ($sp =~ /Hs/){
	$array3[1] =~ s/ENSG0{1,6}/Hs./g;
	$array3[1] =~ s/\s+//g;
	print "$array3[1]\t";
	
	if ($hash_parser{$array3[1]}) {
	    print "match\t$line2;$hash_parser{$array3[1]}";
	}
    }elsif ($sp =~ /Mm/){
	$array3[2] =~ s/ENSMUSG0{1,6}/Mm./g;
	$array3[2] =~ s/\s+//g;
	print "$array3[2]\t";
	if ($hash_parser{$array3[2]}) {
	    print "match\t$line2;$hash_parser{$array3[2]}";
	}
    }elsif ($sp =~ /Rn/){
	$array3[3] =~ s/ENSRNOG0{1,6}/Rn./g;
	$array3[3] =~ s/\s+//g;
	print "$array3[3]\t";
	if ($hash_parser{$array3[3]}) {
	    print "match\t$line2;$hash_parser{$array3[3]}";
	   
	}
    }else { die "digite o simbolo da especie\n"
    }
    print "\n";    
}



