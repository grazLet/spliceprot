#!usr/bin/perl
#vinicius parreira
#24-05-22
# encontra trinca de ortologos que foram encontrados pelo pattern lab


$arq1 = $ARGV[0]; # parser
=c
","TheProteinLoci","Sequence","SpectrumCount","HeaderOriginal"
"1","49564","QQSHFAMMHGGTGFTRIDSSSPEV","1","Hs.251484_vn.1 [ENST00000429752]"
"2","27761","KPGGFDLSLFYR","7","Hs.237417_vn.1 [ENST00000420392]"
"3","6009","EDAQEIADTPSGDKTSLETR","1","Hs.197043_vn.5 [ENST00000521512]"
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
    $line =~ s/"//g;
    $line =~ s/IDENTICA//g;
    $line =~ s/SPLICEPROT_SUB_UNIPROT//g;
    $line =~ s/METHIONINE_SPLICEPROT_SUB_UNIPROT//g;
    $line =~ s/UNIPROT_SUB_SPLICEPROT//g;
    $line =~ s/METHIONINE_UNIPROT_SUB_SPLICEPROT//g;
    $line =~ s/\*//g;
    $line =~ s/\_\:\:\_/\_/g;
    $line =~ s/\"//g;
    $line =~ s/\]/ /g;
    $line =~ s/\[/ /g;
    $line =~ s/\>//g;
    @array1 = split (/\,/, $line);
    $protein_loci = $array1[1];
    $header =  $array1[4];
    
    @array_to_hdSplice = split (/\;/, $header); # para encontrar somentes os id do splice

    foreach $a (@array_to_hdSplice){
	if (($a =~ /^Hs/) || ($a =~ /^Mm/) || ($a =~ /^Rn/)){
	    @array2 = split (/\s+/, $a);
	    $id = $array2[0];
	    shift @array2;
	    
	    $ensts = join( ',', @array2);

	    $hash_parser{$id} = $protein_loci;
	    $hash_parser_ensts{$id} = $ensts;
#	    print "$id $hash_parser{$id} $hash_parser_ensts{$id}\n";
	}else {next;}
    }
  
    #    $id = $array2[0];
    #   $hash_parser{$id} = $protein_loci;
#    print "$id $hash_parser{$id}\n";
    
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
	    print "$hash_parser_ensts{$array3[1]}\tmatch\t$line2;$hash_parser{$array3[1]}";
	}
    }elsif ($sp =~ /Mm/){
	$array3[2] =~ s/ENSMUSG0{1,6}/Mm./g;
	$array3[2] =~ s/\s+//g;
	print "$array3[2]\t";
	if ($hash_parser{$array3[2]}) {
	    print "$hash_parser_ensts{$array3[2]}\tmatch\t$line2;$hash_parser{$array3[2]}";
	}
    }elsif ($sp =~ /Rn/){
	$array3[3] =~ s/ENSRNOG0{1,6}/Rn./g;
	$array3[3] =~ s/\s+//g;
	print "$array3[3]\t";
	if ($hash_parser{$array3[3]}) {
	    print "$hash_parser_ensts{$array3[3]}\tmatch\t$line2;$hash_parser{$array3[3]}";
	   
	}
    }else { die "digite o simbolo da especie\n"
    }
    print "\n";    
}



