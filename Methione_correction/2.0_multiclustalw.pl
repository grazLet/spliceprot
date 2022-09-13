#!usr/bin/perl

###21-12-2020
### vinicius parreira
## programa criado para rodar o clustalw para as comparacoes entre swissprot e spliceprot

$arquivo1 = $ARGV[0]; #arquivo de comparacao

open (IN, "$arquivo1") || die "nao abriu $arquivo1";

$/ = "----------\n"; #mudanca do input record sep. para adapatar o meu arquivo 
#$path = "/home/users/vinicius.parreira/projeto_vinicius/results/spliceProt/clutalw_comparacao_enst/"; #resultados humano
$path = "/home/users/vinicius.parreira/projeto_vinicius/results/spliceProt/camundongo/clutalw_comparacao_enst_camundongo/"; #resultados camundongo
while (<IN>){
    chop($_);
    #print"$_\n";
    @entradas = split (/\n/, $_);
    $resultado = $entradas[0];
    $splice_seq = $entradas[2];
    $swiss_seq = $entradas[4];
    $resultado =~ /\[(.*?)\]\[(.*?)\]\:(.*?)/;
    $enst = $1;
    $uniprot = $2;
#    print "$resultado , $splice_seq , $swiss_seq  , $enst , $uniprot \n"; #debug ok
    $splice = $enst."_splice";
    $swiss = $enst."_swiss";
    $output = $path.$enst.".fasta";
    print" $output\n";
    open (OUT1, ">$output") || die "nao criou fasta do $enst\n"; #arquivo contendo o fasta das duas sequencias
    print OUT1 ">$splice\n$splice_seq\n";
    print OUT1 ">$swiss \n$swiss_seq";
    $log = $enst."_log.txt";
    system "clustalw -infile=$output > $path$log";

}
