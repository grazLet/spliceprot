#!usr/bin/perl

###21-12-2020
### vinicius parreira
## This Script are based on Clustaw alignment file for each sequence, the output is a unique file 
#  with sequences that have been corrected.

#$diretorio = "/home/users/vinicius.parreira/projeto_vinicius/results/spliceProt/clutalw_comparacao_enst/"; #humano
$diretorio = "/home/users/vinicius.parreira/projeto_vinicius/results/spliceProt/camundongo/clutalw_comparacao_enst_camundongo/"; #camundongo
#print "$diretorio\n"; # debug --ok
opendir (DIR, "$diretorio") || die "nao abriu $diretorio\n";
$separator = "/^$/";


while (readdir DIR){
 #   chomp $_;
    $file = $_;
 #   print "$file\n"; #debug --ok
    next if ($file eq "\." || $file eq "\.\.");
    if ($file =~".aln"){
#	print "$file\n"; # debug --ok
	$caminho_in = $diretorio.$file;
#	print "$caminho_in";
	open (IN, "$caminho_in") || die "nao abriu $file\n"; # abre o arquivo
	$enst = $file;
	$enst =~ s/.aln//;
	$seq_splice = "";
	$seq_swiss = "";
	while (<IN>){ # le linha por linha 
	    chomp $_;
	    $line = $_;
	    if ($line =~ /splice/){
		$line =~ s/\s+/;/;
		@aln_splice = split (/\;/, $line);
		$seq_splice = $seq_splice.$aln_splice[1]
		
	    }elsif ($line =~ /swiss/){
                $line =~ s/\s+/;/;
                @aln_swiss = split (/\;/, $line);
		$seq_swiss = $seq_swiss.$aln_swiss[1];
	    } #tenho as sequencias separadas
	    
	}#leu todo o arquivo
	#print "$seq_splice\n$seq_swiss\n";#debug --ok
	@AA_splice = split (//, $seq_splice); #quebra a string por caracter
	@AA_swiss = split (//, $seq_swiss);
	$count = 0;
	$pula= 0;	
	if ($seq_swiss =~ /^-/){
	    for ($i = 0; $i<= length $seq_splice; $i++){
#		print "$AA_splice[$i]\t$AA_swiss[$i]\n"; #debug --ok
		if ($AA_swiss[$i] eq "-" or $AA_swiss[$i] eq "M"){
		    if ($AA_splice[$i] ne $AA_swiss[$i]){
		    #print "$AA_splice[$i]\t$AA_swiss[$i]\n"; #debug --ok
			$count++;
			$pula = 0;		    
			next;
			
		    }elsif ($AA_splice[$i] eq $AA_swiss[$i] and $AA_splice[$i] =~ /M/){
#		    print "$AA_splice[$i]\t$AA_swiss[$i]\t"; #debug --ok
			#print"$AA_swiss[$i]\n";
			$pula = 0;
			last;
		    }
		}elsif($AA_swiss[$i] ne "-" and $AA_swiss[$i] ne "M"){
		    $pula = 1;
		    last;
		}
	    }	
	    #print "$pula\n";
	    if ($pula == 0){
	    $Tam = scalar @AA_splice;
	    @new_seq_splice = @AA_splice[$count..$Tam];
	    print ">$enst"."_novaseq\n";
	    $newseq = "";
	    foreach $a (@new_seq_splice){
		$newseq = $newseq.$a;
	    }
	    $newseq =~ s/-//g;
	    print  $newseq;
	    print "\n";
	    }elsif($pula==1){
		next;
	    }
	    
	}
	
#	exit;
    }

    
}



