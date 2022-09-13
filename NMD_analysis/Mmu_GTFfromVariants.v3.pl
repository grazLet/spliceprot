#usr/bin/perl

##parreira.vsc@gmail.com
#24.09.2019
#programa que utiliza as informacoes das variantes geradas pelo spliceprot e pega o consenso para gerar um gtf

print "==========================\nexecute perl matrizes.csv consenso.csv cluster.csv anotacao.gtf\n===========================\n";

$variant_file = $ARGV[0];
$consenso_file = $ARGV[1];
$cluster_file = $ARGV[2];
$anotacao_file = $ARGV[3];

print "LENDO: ... $variant_file\t$consenso_file\t$cluster_file\t$anotacao_file\n";

###vinicius parrreira versao2 - informacoes adicionais

open (GTF, "<$anotacao_file") || die "did not open $anotacao_file";

=c
1       havana  gene    11869   14409   .       +       .
       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1       havana  transcript      11869   14409   .       +
       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";


=cut

while ($anotacao= <GTF>){
    chomp($anotacao);
    @infoanotacao=split (/\t/, $anotacao); #separa linha do arquivo anotacao
    if ($infoanotacao[2]=~"transcript"){ #executa somente se for uma anotacao para o transcrito
    #$infoanotacao[-1] =~ /(.*?); transcript_id \"(.*?)\"; (.*?); transcript_biotype \"(.*?)\"; (.*?); transcript_support_level \"(\w)\"/;
	$infoanotacao[-1] =~ m/transcript_id \"(.*?)\"\;/; #busca o ENST
	$trans_id = $1;
	$infoanotacao[-1] =~ m/transcript_biotype \"(.*?)\";/; #busca o tipo do ENST
	$gene_tipo = $1;	
	unless($infoanotacao[-1] =~/transcript_support_level/){ #se nao tiver um TSL usa NULL
	    $TSL = "NULL";
	}else{ #se tiver um TSL captura a indormacao da anotacao
	    $infoanotacao[-1] =~ /transcript_support_level \"(\w{1,2})(.*?)\"/;
	    $TSL = $1;
	}
	$variant_name = $trans_id; #nome da variante ou ENST
	$hashtsl{$variant_name} = $TSL; #essa variante vai estar associada a um TSL
	$hashfunction{$variant_name} = $gene_tipo; #essa variante esta associada a algo na anotação

#	print "$trans_id\t$gene_tipo\t$TSL\n"; # - debug
    }
}   




open (VARIANT,"<$variant_file") || die "did not open $variant_file\n";# arquivos de variantes para pegar as matrizes preditas
open (CONSENSO,"<$consenso_file") || die "did not open $consenso_file\n";# arquivo consenso para pegar as coordenadas
open (CLUSTER, "<$cluster_file") || die "did not open $cluster_file\n"; #arquivo cluster para pegar a orientacao do gene
open (OUT, ">Variantes_Ensembl.gtf") || die "did not create Variantes_SpliceProt.gtf\n";


############# CLUSTER ###################
#header
#504505  X       ENST00000612152 100633931       100634029       675     577     -       100     3
##########################################

%hash=();


print "\n==========\ntrabalhando arquivo de cluster\n===========\n";
 
while (<CLUSTER>){ #recupera informacao da fita em que o gene é transcrito
    chomp($_);
    if ($_ =~ /^\d/){ #nao le a header
	@cluster = split /\t/, $_; #quebra a linha
	$fita=$cluster[7]; #orinetacao
	$ENSG_idd=$cluster[-1];#ENST do transcrito
	$hash{$ENSG_idd}=$fita;
#	print "$_\n"; #--debug
    }
}


############# CONSENSO ###################
#id_consenso     cluster_id      vl_inicio       vl_fim
#6758719 3       100627108       100629986
#6758720 3       100630759       100630866
#6758721 3       100632063       100632068 
#6758722 3       100632485       100632568
#6758723 3       100633405       100633539
#6758724 3       100633931       100634029
#6758725 3       100635178       100635252
#6758726 3       100635558       100635746
#6758727 3       100636191       100636607 
#6758728 3       100636608       100636792
#6758729 3       100636793       100636806
#6758730 3       100636807       100637104
#6758731 3       100639945       100639991
##########################################

$count_var_name=1; #variant name
%hashchr=();
%hashconsenso_inicio=();
%hashconsenso_fim=();
%exonposition_in=();
%exonposition_end=();

# antes usava consenso, agora usa cluster, pq preciso de todas as coordenadas, não somente as coordenadas de TSL1
#OS ARQUIVOS TEM QUE ESTAR OOOOOORRRRRRDDDDEEEENNNAAAAAAAADOOOOOOOOOOOS!!!!!!!!!!!!!!!!!!! variantes ordenados pelo cluster_id

print "\n==========\ntrabalhando arquivo de consenso\n===========\n";
open (CONSENSO,"<$consenso_file") || die "did not open $cluster_file\n"; #arquivo cluster para pegar a 

#open (CLUSTER2, "<$cluster_file") || die "did not open $cluster_file\n"; #arquivo cluster para pegar a orientacao do gene
#while ($consenso_line = <CLUSTER2>){
while ($consenso_line = <CONSENSO>){    
    chomp($consenso_line);
#    print "$consenso_line\n";
    if ($consenso_line =~ /^\d/){ #nao le a header
#    $count_con++; next if ($count_con<=1); # not header
#      print "$consenso_line\n"; #debug
    @consenso_array = split(/\t/, $consenso_line); #quebra a linha 
    $ENSG_consenso = $consenso_array[1]; #ENSG do transcrito
    $chr_in_consenso =$consenso_array[2]; #posicao do comeco do transcrito
    $chr_end_consenso =$consenso_array[3]; #posicao do fim do trasncrito
# print "$ENSG_consenso\t$chr_in_consenso\t$chr_end_consenso\n";   
    
    if ($ENSG_consenso != $new_ENSG_consenso){
	$count_consenso_exon=0;   #quando comecar a avaliar outro ENSG comeca a contar o exons novamente
    }
    
    $count_consenso_exon++;    # contar o numero do exon no consenso, esse numero de exon eh teorico pq se as coord forem consecutivas eh um exon so

  
    $hashconsenso_inicio{$ENSG_consenso}{$count_consenso_exon}=$chr_in_consenso;
    #                        o gene        o exon teorico         in. loc. chr
    
#    print "$ENSG_consenso\t$count_consenso_exon\t====>$chr_in_consenso\n";   #-debug
   

    $hashconsenso_fim{$ENSG_consenso}{$count_consenso_exon}=$chr_end_consenso;
#                        o gene        o exon teorico         fim. loc. chr
    
 #   print "$ENSG_consenso\t$count_consenso_exon\t====>$chr_end_consenso\n\n";  #  -degub
    $new_ENSG_consenso = $ENSG_consenso; #avaliar depois se ainda esta vendo um mesmo gene ou se mudou
#print "$hashconsenso_inicio{$ENSG_consenso}{$count_consenso_exon}=$chr_in_consenso\t$hashconsenso_fim{$ENSG_consenso}{$count_consenso_exon}=$chr_end_consenso";    
    }    
} #para o gene eu tenho o numero do exon e onde ele comeca e termina

    
#print "\n==========\ntrabalhando arquivo de variantes\n===========\n";
    
    

    

$count_var=0; #header kill
    
############# VARIANT ####################
#id_matrizes     chr     seq_acc matriz  cluster_id      variante
#1   X    ENST00000373020 | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 0110 | 0 |    3     | 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 0111 | 1 |
#2   X    ENST00000494424 | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 0011 | 1 |    3     | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0011 | 1 |
#3   X    ENST00000496771 | 0 | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1100 | 0 |    3     | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1111 | 1 |
##########################################


while ($variant_line = <VARIANT>){
    chomp($variant_line);
    $count_var++; next if ($count_var<=1); # not header
    #    print "$variant_line\n"; #debug
    @variant_array = split(/\t/, $variant_line); #quebra a linha do arquivo
    $variant_matrix = $variant_array[3]; # coluna matriz (v3)
#    $variant_matrix = $variant_array[5]; # coluna variante
    $ENSG_variant = $variant_array[4]; #gene_id
    $hashchr{$ENSG_variant}=$variant_array[1]; #em que cromossomo está o gene
    
#    if ($NEW_ENSG_variant != $ENSG_variant){ #avaliar a cada novo ENSG(gene_id)
#	$count_var_name=1; # variavel que serve para criar um nome numerico para cada variante (ENST) (v3 - não usa mais pq esta usando o ENST sem o 'ENST')
 #   } 
    
  #  $variant_name = $ENSG_variant."var".$count_var_name; (v3 - não usa mais pq esta usando o ENST sem o 'ENST')

    #03.10 parreira.vsc v2 - inserir ESNT ao invés de nome teorico
    $variant_name = $variant_array[2]; #captura o seq_acc (ENST)
    $variant_name =~ s/ENSMUST//;

 #   $NEW_ENSG_variant=$ENSG_variant;  #variavel para avaliar depois se ainda esta vendo um mesmo gene ou se mudou
        
#    $count_var_name++;
 #       print "$variant_matrix\t$ENSG_variant\t$variant_name\n"; #debug

    @each_column_variant = split(//, $variant_matrix); #pega cada elemento da matriz de variantes
    @each_exon_variant=(); #zera o array que pega as informacoes para cada exon (1 ou 0)
    foreach (@each_column_variant){
#	print "$_\n"; #debug
	if ($_ =~ /\d/){ # elimina o que é "|" do array, so vai usar 1 e 0 para saber se a coordenada esta presente ou nao
	    push (@each_exon_variant, $_);
	}
    }
    #foreach (@each_exon_variant){ #debug
    #	    print "$_\n"; #debug
    #} #debug 

    #agora neste ponto eu tenho um array que contem todos os digitos presentes na matriz de variantes, removi o pipe, os espaços, logo, quem esta no mesmo exon, ou seja, nao era separado por pipe esta junto com todos os outros, eu vou ter que deduzir que eles eram do mesmo exon vendo que eles sao consecutivos, caso fique ruim eu volto nesta etapa e marco esses que estarao dentro de um mesmo par de pipe.
    

    $count_exon_hash=0; # um numero para identificar cada exon pq 0 e 1 vai repetir

    foreach $exon_teorico (@each_exon_variant){
	#$exon_teorico - informacao de ausencia ou presenca do exon 1 ou 0
	$count_exon_hash++; #numero do éxon teórico
	
	$exonposition_in{$ENSG_variant}{$exon_teorico}{$count_exon_hash} =  $hashconsenso_inicio{$ENSG_variant}{$count_exon_hash};
		
	$exonposition_end{$ENSG_variant}{$exon_teorico}{$count_exon_hash} = $hashconsenso_fim{$ENSG_variant}{$count_exon_hash};
	
	# neste hash eu tenho cada numero (exon teorico) da matriz de variantes associado a uma coordenada de inicio e uma coordenada de fim.

=c - debug
	print "$ENSG_variant\t$exon_teorico\t$count_exon_hash\t=======>$hashconsenso_inicio{$ENSG_variant}{$count_exon_hash}\n";
	print "$ENSG_variant\t$exon_teorico\t$count_exon_hash\t=======>$hashconsenso_fim{$ENSG_variant}{$count_exon_hash}\n";


=cut - debug
	

    }



### 03.10 vinicius parreira v2 insercao da informacao de TSL e anotacao
    #criados
    # hashtsl{$variant_name}
    #"$hashfunction{$variant_name}
    
    
    $count_exon_validate=0; # vai contar para pegar os exons teoricos 
    $posicaofinal=0; #auxiliar para trabalhar com as consecutivas
    $exon_number=0; # vai entrar o verdadeiro numero do exon no transcrito

    
    $i = scalar @each_exon_variant; #tamanho do array (numero de 0  e 1) (o mesmo que o numero do ultimo exon)
    @each_exon_variant2 = @each_exon_variant; 

   
    
    while (@each_exon_variant2){
	if ($each_exon_variant2[-1] eq "1"){
	    $j = scalar @each_exon_variant2;
	    last;
	}
	$lasts =  pop @each_exon_variant2;
#	print "$lasts\n";
	if ($lasts eq "0"){next;}
	else{
	    $j = scalar @each_exon_variant2;
	    last;
	}
    }#esse while serve para identificar o real tamanho do array, ou seja, tem exon ate onde tem o ultimo 1
   
    @each_exon_variant3 = @each_exon_variant;
    $q=1;
    foreach (@each_exon_variant3){
	if ($_ eq "1"){
	    last;
	}else{ $q++;}
    }
	
#    print "$i\t$j\n"; #- debug
    $variant_name2 = "ENSMUST".$variant_name;
    #    print "\n\n$i\n"; - debug
    #    print OUT "$hashchr{$ENSG_variant}\tSpliceProt\ttranscript\t$exonposition_in{$ENSG_variant}{1}{1}\t$exonposition_end{$ENSG_variant}{1}{$i}\t1000\t$hash{$ENSG_variant}\t.\tgene_id \"$ENSG_variant\"; transcript_id \"$variant_name\"; transcript_suport_level \"$hashtsl{$variant_name}\"; annotation \"$hashfunction{$variant_name}\";\n";
        print OUT "$hashchr{$ENSG_variant}\tSpliceProt\ttranscript\t$exonposition_in{$ENSG_variant}{1}{$q}\t$exonposition_end{$ENSG_variant}{1}{$j}\t1000\t$hash{$ENSG_variant}\t.\tgene_id \"$ENSG_variant\"; transcript_id \"$variant_name\"; transcript_suport_level \"$hashtsl{$variant_name2}\"; annotation \"$hashfunction{$variant_name2}\";\n";
    
    # $exonposition_in  tem a key2 e key 3 igual a um pq como eh referente ao trasncrito inteiro ele vai do primeiro ao ultimo exon ja $exonposition_end tem a key 3 igual a $i pq seu valor eh o numero do ultimo exon.
    
    foreach $leitura (@each_exon_variant){
	chomp ($leitura);
	$count_exon_validate++;
	if ($leitura == 1){ #tem que comecar com '1' pq eh quando tem o exon
	   
	    if ($posicaofinal != 0){
		
		if ($posicaofinal == $exonposition_in{$ENSG_variant}{$leitura}{$count_exon_validate}){  #esse if so vai poder ocorrer depois de o else abaixo "$posicaofinal != 0"
		    $posicaofinal = $exonposition_end{$ENSG_variant}{$leitura}{$count_exon_validate};
		    $posicaofinal++;
		    $posicaofinal_Eanterior = $exonposition_end{$ENSG_variant}{$leitura}{$count_exon_validate};
		    if ($count_exon_validate == $j){ #resolve o caso de ter exons teoricos consecutivos em coordenadas no final da variante
			print OUT "$posicaofinal_Eanterior\t1000\t$hash{$ENSG_variant}\t.\tgene_id \"$ENSG_variant\"; transcript_id \"$variant_name\"; exon_number \"$exon_number\";\n";
		    }
		    next;
		}else{
		    print OUT "$posicaofinal_Eanterior\t1000\t$hash{$ENSG_variant}\t.\tgene_id \"$ENSG_variant\"; transcript_id \"$variant_name\"; exon_number \"$exon_number\";\n";
		}
	    }
	    $exon_number++;
	    print OUT "$hashchr{$ENSG_variant}\tSpliceProt\texon\t$exonposition_in{$ENSG_variant}{$leitura}{$count_exon_validate}\t";
	    #a posicao final
	    $posicaofinal = $exonposition_end{$ENSG_variant}{$leitura}{$count_exon_validate};
	    $posicaofinal_Eanterior = $exonposition_end{$ENSG_variant}{$leitura}{$count_exon_validate}; #criado para avaliar se as posicoes sao consecutivas (mesmo exon)
	    $posicaofinal++;
#	    print  "#$exonposition_end{$ENSG_variant}{$leitura}{$count_exon_validate}\n";
	    #	    print "#$posicaofinal\n";
	        if ($count_exon_validate == $j){ #resolve a ultima coordenada que nao vai ter como comparar pq nao vem nada depois
		    print OUT "$posicaofinal_Eanterior\t1000\t$hash{$ENSG_variant}\t.\tgene_id \"$ENSG_variant\"; transcript_id \"$variant_name\"; exon_number \"$exon_number\";\n";
	    }
	   	    
	} # leitura = 1
	#	else ($count_exon_validate == $i){  print OUT "$posicaofinal_Eanterior\t1000\t$hash{$ENSG_variant}\t.\tgene_id \"$ENSG_variant\"; transcript_id \"$variant_name\"; exon_number \"$exon_\number\";\n";
	
    } #foreach
}

print "\n==========\ngerou GTf...   Variantes_SpliceProt.gtf   \n===========\n";

close(VARIANT);
close(CONSENSO);
close(CLUSTER);
close(OUT);
close(GTFFILE);
exit;
