 #!usr/bin/perl
use strict;
use warnings;

####
# vinicius parreira
# 9-04-2019
#this scripts uses the PSL file to search the entries in the blast 9 alignment file and get the % of align identity. The output is a file to populate a SQL Database
####


print "\n=================\nperl FormatSQLtable.BLAST9ver.pl dicio.tsv in.psl in.blast9 out.sql\n=================\n";
my $dicionario = $ARGV[0];
my $infile = $ARGV[1];
my $blast9file = $ARGV[2];
my $outfile = $ARGV[3];
my @info;
my $count=0;
my $primkey=0;
my @tid;
my @transcripts;
my %hash;


###########################
#dicionario de identidades do alinhamento
=c
# BLAT 36x2 [2009/02/26]
# Query: ENST00000390583.1
# Database: /data1/users/viniciusparreira/Chr_seq_fa/hg38.fa
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
ENST00000390583.1	chr14	100.00	31	0	0	1	31	105904527	105904497	2.9e-09	60.0
=cut

open (IDENT, "$blast9file") || die "Did not open $blast9file\n";
my $blastline;
my $string;
my $idtrans;
my @blast;
my $enst;
my %hashidentid;
while ($blastline =  <IDENT>){
    chomp ($blastline);
#    print "$blastline\n";
    next if ($blastline=~/^\#/);
    @blast = split(/\s+/, $blastline);
    $blast[0] =~ s/\.\d{1,2}//;
#    print "$blast[0]\n";
    $string = "$blast[1], $blast[0], $blast[8], $blast[9]";
    $hashidentid{$string} = $blast[2];
}
@blast=();
close (IDENT);
###########################

###########################
#dicionario ENST vs ENSG

open (DICIO,"<$dicionario") || die "did not open $dicionario\n";

my %hashdicio;
my $ensembl;
my @splithash;


while ($ensembl = <DICIO>){
    chomp($ensembl);
    @splithash = split(/\s+/, $ensembl);
#    print "$splithash[0]    $splithash[1]\n";                                                                                                                                                              
    $hashdicio{$splithash[0]} = $splithash[1];
}
close(DICIO);

###########################

open (IN, "$infile") || die "did not open $infile\n"; ##psl
open (OUT, ">$outfile") || die "did not open $outfile\n";


print OUT "id_clusters\tchr\tseq_acc\tstart_chr\tend_chr\tstart_est\tend_est\tstrand\tidentity\tcluster_id\n";


while(<IN>){
    chomp $_;
    $count++;
    next if ($count <= 5);
    my @split = split(/\t/,$_);
    my $chrid = $split[13];
    my $seq_acc = $split[9];
    $seq_acc =~ (s/\.\d{1,2}//);
    my $cluster_id = $hashdicio{$seq_acc};
    $cluster_id =~ m/([A-Z]*)([0]*)(\d*)/ig; # marca 3 subvariaveis $1, $2, S3, aqui, segue (sequencias de letras tamanho var.)(sequencia de 0s tamanho var.)(sequencia de numero tamanho var.)
    my $newnameCluster = $3;
    print "$cluster_id\n";
    my $strand = $split[8];
    my @blocks = split(/\,/,$split[18]);	
    my @transcriptstarts = split(/\,/,$split[19]);
    my @starts = split(/\,/,$split[20]);
    my $blocknum = scalar @blocks;
    my $start_chr;
    my $end_chr;
    my $start_est;
    my $end_est;
    my $identity;
    my %hash=();
    my $bloco;
    if ($strand eq "+") {

	for(my $i = 0; $i < $blocknum; $i++){ 
	    $bloco=($blocks[$i]-1);
	    $start_chr = ($starts[$i]+1);# soma mais 1 no que tá no psl pq ele considera o primeiro nucleotideo como poiscao 0
	    $end_chr = ($start_chr+$bloco);# mas para ter o final do trasncrito eu somo o tamanho com a posicao inicial no cromossomo
	    $start_est = ($transcriptstarts[$i]+1); # posicao inicial do exon do transcrito
	    $end_est = ($start_est+$bloco); #somo o tamanho do bloco(exon) para ter a posicao final do exon no transcrito
	    $identity = "$chrid, $seq_acc, $start_chr, $end_chr";
	    $primkey++;
	    if ($hashidentid{$identity}){
		print OUT"$primkey\t$chrid\t$seq_acc\t$start_chr\t$end_chr\t$start_est\t$end_est\t$split[8]\t$hashidentid{$identity}\t$newnameCluster\n";	    
	    }else{
		print OUT"$primkey\t$chrid\t$seq_acc\t$start_chr\t$end_chr\t$start_est\t$end_est\t$split[8]\t99\t$newnameCluster\n";	    
	    }
	    
	}    
	
    }else{ ## negativa

	my $total =0;#	my @reverseblocks = reverse @blocks; #for antisense
	my @reversetranscriptstarts = reverse @transcriptstarts; #inversion for antisense
	for(my $i = 0; $i < $blocknum; $i++){ 
	    $bloco=($blocks[$i]-1);
	    #$bloco=($blocks[$i]);
	    $total=0;
	    $primkey++;
	    unless ($i >= 1){
		foreach my $soma (@blocks){
		    chomp ($soma);
		    $total=$total+$soma;# psl add +1 positions due to initialization in position 0
		}
		$start_est = ($total);
		$end_est = ($total - $bloco);
#		if ($end_est == 0){ $start_est++; $end_est++;}
		$start_chr = ($starts[$i]+1);# psl add +1 positions due to initialization in position 0		
		$end_chr = ($bloco+$start_chr);# exon end is the length plus initial position
		$identity = "$chrid, $seq_acc, $end_chr, $start_chr";
#		print "(-) first $identity\n\n";
		if ($hashidentid{$identity}){
		    print OUT"$primkey\t$chrid\t$seq_acc\t$start_chr\t$end_chr\t$start_est\t$end_est\t$split[8]\t$hashidentid{$identity}\t$newnameCluster\n";	    
		}else{
		    print OUT"$primkey\t$chrid\t$seq_acc\t$start_chr\t$end_chr\t$start_est\t$end_est\t$split[8]\t99\t$newnameCluster\n";
		}
	    }else{
		$start_est = ($end_est-1);
		$end_est = ($start_est - $bloco);
		$start_chr = ($starts[$i]+1);# psl add +1 positions due to initialization in position 0		
		$end_chr = ($bloco+$start_chr);#  exon end is the length plus initial position
		$identity = "$chrid, $seq_acc, $end_chr, $start_chr";
#		print "(-) others $identity\n\n";
		if ($hashidentid{$identity}){
		    print OUT"$primkey\t$chrid\t$seq_acc\t$start_chr\t$end_chr\t$start_est\t$end_est\t$split[8]\t$hashidentid{$identity}\t$newnameCluster\n";	    
		}else{
		    print OUT"$primkey\t$chrid\t$seq_acc\t$start_chr\t$end_chr\t$start_est\t$end_est\t$split[8]\t99\t$newnameCluster\n";
		}
		
	    }	

	}

    }
}

print "============END============\n\nCriou $outfile :)\n\n";
close (IN);
close (OUT);

exit;
