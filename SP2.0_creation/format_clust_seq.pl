#!usr/bin/perl

###parreira.vcs@gmail.com
#12/07/2019

#format clust_seq table to SpliceProt database postgre update data


#####
#cluster_id    sequence    seq_acc   tissue    TSL    cds_start    cds_end
#(ENSG)       (fasta_seq)   (ENST)   (NULL)    (TSL)     (start)            (end)
#####


#input arquives:
#ENST_seq_TSL.csv>>
#id_tsl  seq_acc        sequence        tsl_info        gene_id

#input_clusters_table.csv>>
#id_clusters     chr     seq_acc start_chr       end_chr start_est       end_est strand  identity       cluster_id



#input_clusters_table.csv - info about  ENSG, ENST, START, END
#ENST_seq_TSL.csv -info about TSL and sequence


$arquive1=$ARGV[0];
$arquive2 =$ARGV[1];

#print "$arquive1\t$arquive2\n"; 

open (IN1, "<$arquive1") || die "did not open input 1\n"; 
open (IN2, "<$arquive2") || die "did not opne input 2\n"; 
open (OUT1, ">cluster_seq_ENSEMBL.csv") || die "did not create output 1\n"; 
open (OUT2, ">seqfile_ENSEMBL.csv") || die "did not create output 2\n"; 


%hash=();
$count=0;

print OUT1 "cluster_id\tsequence\tseq_acc\ttissue\ttsl\n";
while ($a=<IN1>){
    $count++;
    next if $count<=1;
    chomp($a);
    @array1 = split (/\t/, $a);
    $cluster_id=$array1[4];
    $cluster_id=~ m/([A-Z]*)([0]*)(\d*)/ig;
    $newnameCluster = $3;
    $ENST = $array1[1];
    $sequence = $array1[2];
    $tsl = $array1[3];
    $tissue = "NULL";
    $cds_start="";
    $cds_end="";

    print OUT1  "$newnameCluster\t$sequence\t$ENST\t$tissue\t$tsl\n";
#    $hash{$ENST}= "$sequence\t$ENST\t$tissue\t$tsl"; 

#print "$ENST\t$sequence\t$tsl\t$tissue\n"; 
}



$count=0;

while ($b=<IN2>){
    $count++;
    next if $count<=1;
    chomp($b);
    @array2 = split(/\t/, $b);
    $seq_acc = $array2[2];
    $transcript_start = $array2[5];
    $transcript_end = $array2[6];
    $cluster_id = $array2[9];
    if ($hash{$seq_acc}){
	print OUT1 "$cluster_id\t$hash{$seq_acc}\t$cds_start\t$cds_end\n";
    }else{
	push @teste, $seq_acc;
    }
}
foreach $c (@teste){print OUT2 "$c\n";}


close (IN1);
close (IN2);
close (OUT1);
close (OUT2);
