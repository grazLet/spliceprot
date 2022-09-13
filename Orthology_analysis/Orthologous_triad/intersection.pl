#!usr/bin/perl
#vinicius parreira
#05-06-22
# encontra intesecao de ortologos que foram encontrados pelo pattern lab


$arq1 = $ARGV[0]; # humano
$arq2 = $ARGV[1]; # camundongo
$arq3 =  $ARGV[2]; #rato

open (IN, "$arq1") || die "nao abre $arq1";
open (IN2, "$arq2") || die "nao abre $arq2";
open (IN3, "$arq3") || die "nao abre $arq3";


%hash_human =();
%hash_mouse =();
%hash_rat =();

%HH_loci=();
%HM_loci=();
%HR_loci=();

%HH_trans=();
%HM_trans=();
%HR_trans=();

while (<IN>){
    chomp $_;
    $line_hum = $_;

    $line_hum =~ s/-/_/g;
    $line_hum =~ s/"//g;
    @array_hum = split (/\t/, $line_hum);
    @array1 = split (/\;/, $line_hum);
    $enst_human = $array_hum[1];
    $H_hum = $array1[1];
    $H_mou = $array1[2];
    $H_rat = $array1[3];
    $hash_human{$H_hum} = 1;
    $HH_loci{$H_hum} = $array1[4];
    $HH_trans{$H_hum} = $enst_human;
#     print "$H_hum $hash_human{$H_hum}\n";
}



while (<IN2>){
    chomp $_;
    $line_cam = $_;

    $line_cam =~ s/-/_/g;
    $line_cam =~ s/"//g;
    @array_mouse = split (/\t/, $line_cam);
    @array2 = split (/\;/, $line_cam);
    $enst_mouse = $array_mouse[1];
    $C_hum = $array2[1];
    $C_mou = $array2[2];
    $C_rat = $array2[3];
    $hash_mouse{$C_mou} = 1;
    $HM_loci{$C_mou} = $array2[4];
    $HM_trans{$C_mou} = $enst_mouse;
#    print "$C_mou $hash_mouse{$C_mou}\n";
}



while (<IN3>){
    chomp $_;
    $line_rat = $_;

    $line_rat =~ s/-/_/g;
    $line_rat =~ s/"//g;
    @array_rat = split (/\t/, $line_rat);
    @array3 = split (/\;/, $line_rat);
    $enst_rat = $array_rat[1];
    $R_hum = $array3[1];
    $R_mou = $array3[2];
    $R_rat = $array3[3];
    $loci_rat= $array3[4];
    $hash_rat{$R_rat} = 1;
    
    #    print "$R_rat $hash_rat{$R_rat}\n";
    if (($hash_human{$R_hum} == 1) && ($hash_mouse{$R_mou})) {
	print "$R_hum\t($HH_loci{$R_hum})\t[$HH_trans{$R_hum}]\t$R_mou\t($HM_loci{$R_mou})\t[$HM_trans{$R_mou}]\t$R_rat\t($loci_rat)\t[$enst_rat]\n";


    }
}
