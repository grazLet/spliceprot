
#!usr/bin/perl
#vinicius parreira
#05.12.2019
##the output is a file to populate the clust_seq table in the mapping schema of mouse 
#24.08.2020
#script create the  select_consenso file of mouse 
#id_tsl  seq_acc sequence        tsl_info        gene_id
#    1       ENST00000496771 GTCCCCGTCTCGGACAGTTATAAAGGTGATGACCATTA       3       ENSG00000000003


$in = $ARGV[0];
$oi = $ARGV[1];
 
open (IN, "$in") || die "did not open $in\n";
print "id_select_consenso\tnm_transcrito\tcluster_id\ttsl_true\tlength\ttsl\n";
$count=0;
%hashseq=();
%hashtsl=();
%hashenst=();
%hashlista=();
%hashgene=();

while (<IN>){
    chomp ($_);
    $count++;
    next if ($count<=1);
    $line = $_;
    @array = split (/\t/, $line);
    $enst3 = $array[1];
#    $enst3 =~ m/([A-Z]*)(\d*)/ig;
 #   $order = $2;
    $seq = $array[2];
    $tsl = $array[3];
    if ($tsl eq "NULL" || $tsl eq "NA"){
	$tsl = 7;
    }
    $numero = length ($seq);
    $enst = $tsl."00".$numero.$enst3;
    $gene = $array[4];
    $gene =~ m/([A-Z]*)([0]*)(\d*)/ig;
    $gene = $3;
#    print "$gene\t$enst\t$tsl\n"; #ok
    $hashlista{$gene}{$enst}=$seq;
#    print "$hashlista{$gene}{$enst}\n"; #ok
    $hashenst{$gene}{$seq}= $enst;
    $hashseq{$enst}=$seq;
    $hashtsl{$enst}=$tsl;
    $hashgene{$enst}=$gene;
    
}
@consenso=();
    
foreach $gene1 (sort {$a <=> $b }keys %hashlista){
#    print "$gene1\n"; #ok
    @transcritos_1=();
    @transcritos_2=();
    @transcritos_3=();
    @transcritos_4=();
    @transcritos_5=();
    @transcritos_0=();
    
    foreach $enst1 (sort {$a <=> $b} keys %{$hashlista{$gene1}}){
#	print "$hashtsl{$enst1}\n"; #ok 
	
	if ($hashtsl{$enst1} == 1){
	    push @transcritos_1, $enst1;
	}elsif ($hashtsl{$enst1} == 2){
	    push @transcritos_2, $enst1;
	}elsif ($hashtsl{$enst1} == 3){
	    push @transcritos_3, $enst1;
	}elsif ($hashtsl{$enst1} == 4){
            push @transcritos_4, $enst1;
	}elsif ($hashtsl{$enst1} == 5){
            push @transcritos_5, $enst1;
	}else{
	    push @transcritos_0, $enst1;
	}
    }
    $num_trans_1 = scalar @transcritos_1;
    $num_trans_2 = scalar @transcritos_2;
    $num_trans_3 = scalar @transcritos_3;
    $num_trans_4 = scalar @transcritos_4;
    $num_trans_5 = scalar @transcritos_5;
    $num_trans_0 = scalar @transcritos_0;
    $tamanho_avaliar_transcrito_old;
        
    #print "$num_trans_1\n$num_trans_2\n$num_trans_3\n$num_trans_4\n$num_trans_5\n$num_trans_0\n"; #ok
    
    if ($num_trans_1 > 0){
	#	print "$num_trans_1\n";
	if ($num_trans_1 < 2){
	    foreach $escolha1 (@transcritos_1){
		push @consenso, $escolha1;
	    }
	}else{
	    $tamanho_avaliar_transcrito_old=0;
	    foreach $escolher1 (@transcritos_1){
#		print "$escolher1\t"; #ok
		$tamanho_avaliar_transcrito = length($hashseq{$escolher1});
#		print "$tamanho_avaliar_transcrito\n"; #ok
		if ($tamanho_avaliar_transcrito > $tamanho_avaliar_transcrito_old){
		    $meu_transcrito1  =  $escolher1;
		}
		$tamanho_avaliar_transcrito_old =  $tamanho_avaliar_transcrito;
	    }
	    push @consenso, $meu_transcrito1;
	}next;
    }elsif($num_trans_2 > 0){
	if ($num_trans_2 < 2){
	    foreach $escolha2 (@transcritos_2){
		push @consenso, $escolha2;
	    }
	}else{
	    $tamanho_avaliar_transcrito_old=0;
	    foreach $escolher2 (@transcritos_2){
		$tamanho_avaliar_transcrito = length($hashseq{$escolher2});
		if ($tamanho_avaliar_transcrito > $tamanho_avaliar_transcrito_old){
		    $meu_transcrito2  =  $escolher2;
		}
		$tamanho_avaliar_transcrito_old =  $tamanho_avaliar_transcrito;
	    }
	    push @consenso, $meu_transcrito2;
	}next;
    }elsif($num_trans_3 > 0){
	if ($num_trans_3 < 2){
	    foreach $escolha3 (@transcritos_3){
		push @consenso, $escolha3;
	    }
	}else{
	    $tamanho_avaliar_transcrito_old=0;
	    foreach $escolher3 (@transcritos_3){
		$tamanho_avaliar_transcrito = length($hashseq{$escolher3});
		if ($tamanho_avaliar_transcrito > $tamanho_avaliar_transcrito_old){
		    $meu_transcrito3 =  $escolher3;
		}
		$tamanho_avaliar_transcrito_old =  $tamanho_avaliar_transcrito;
	    }
	    push @consenso, $meu_transcrito3;
	}next;
    }elsif($num_trans_4 > 0){
	if ($num_trans_4 < 2){
	    foreach $escolha4 (@transcritos_4){
		push @consenso, $escolha4;
	    }
	}else{
	    $tamanho_avaliar_transcrito_old=0;
	    foreach $escolher4 (@transcritos_4){
		$tamanho_avaliar_transcrito = length($hashseq{$escolher4});
		if ($tamanho_avaliar_transcrito > $tamanho_avaliar_transcrito_old){
		    $meu_transcrito4 =  $escolher4;
		}
		$tamanho_avaliar_transcrito_old =  $tamanho_avaliar_transcrito;
	    }
	    push @consenso, $meu_transcrito4;
	}next;
    }elsif($num_trans_5 > 0){
	if ($num_trans_5 < 2){
	    foreach $escolha5 (@transcritos_5){
		push @consenso, $escolha5;
	    }
	}else{
	    $tamanho_avaliar_transcrito_old=0;
	    foreach $escolher5 (@transcritos_5){
		$tamanho_avaliar_transcrito = length($hashseq{$escolher5});
		if ($tamanho_avaliar_transcrito > $tamanho_avaliar_transcrito_old){
		    $meu_transcrito5 =  $escolher5;
		}
		$tamanho_avaliar_transcrito_old =  $tamanho_avaliar_transcrito;
	    }
	    push @consenso, $meu_transcrito5;
	}next;
    }elsif($num_trans_0 > 0){
	if ($num_trans_0 < 2){
	    foreach $escolha0 (@transcritos_0){
		push @consenso, $escolha0;
	    }
	}else{
	    $tamanho_avaliar_transcrito_old=0;
	    foreach $escolher0 (@transcritos_0){
		$tamanho_avaliar_transcrito = length($hashseq{$escolher0});
		if ($tamanho_avaliar_transcrito > $tamanho_avaliar_transcrito_old){
		    $meu_transcrito0 =  $escolher0;
		}
		$tamanho_avaliar_transcrito_old =  $tamanho_avaliar_transcrito;
	    }
	    push @consenso, $meu_transcrito0;
	}next;
    }
}                	

$chave=0;
%hashcheck=();
foreach $i (@consenso){
    $chave++;
    $i =~ m/(.*?)ENSMUST(\d*)/ig;
    $i2="ENSMUST".$2;
    $tamanho= length $hashseq{$i};
    print "$chave\t$i2\t$hashgene{$i}\t$hashtsl{$i}\t$tamanho\ttrue\n";

}

