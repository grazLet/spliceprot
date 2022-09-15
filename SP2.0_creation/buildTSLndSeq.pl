#!usr/bin/perl

#use strict;
#use warnings;

#parreira.vsc@gmail.com
#03/06/2019 - ICC - Fiocruz
# Search in the genome the TSL information and make a tabel with ENST and primary key to upload in a database.

###### version 2.0, sequence associated ######

my $cdnafile = $ARGV[0];

open (CDNA, "<$cdnafile") || die "did not open $cdnafile (first file). Call perl:\n===================\nperl buildTSLndSeq.pl cdna.fa genome.gtf\n===================\n";


$/="\>"; # changed the input record separator to read a fasta archive. The identificator and the sequence will be a same element (line) when reading the arquive (together in the $getseq variable).

my %seqhash;
my $count2=0;
my @onlyENST;


while (my $getseq = <CDNA>){
    $countNs=0;
    chomp($getseq);
    $count2++;
    
    next if $count2<=1; 
    
    @onlyENST = ();
    
    if ($getseq =~m/\n/){ #this option make a first '\n' in the $getseq eq to '\s'. The fasta sequence will be in the same line of informations separeted by '\s'.
	$getseq =~s/\n/ /;
    }
    $getseq =~s/\n//g; #the first '\n' was eliminated above, now this exclude the '\n's that breaks the fasta sequence.

#    print "$getseq\n";
 
    @onlyENST = split /\s+/, $getseq; #separe all the information of fasta header and sequence that now is the last element 
    
#    print "$onlyENST[0]\n";

    $onlyENST[0] =~ s/\.\d{1,2}//; # remove the last digits of an ENST

#   print "$onlyENST[0]\t";
#    print "$onlyENST[-1]\n";

    $seqhash{$onlyENST[0]}=$onlyENST[-1];
}

close (CDNA);

###### ======== ###### End of sequence seeking 

###### ======== ###### Start of TSL seeking

$/ = "\n"; # return the input record separator to '\n' to read each line.

my$genome = $ARGV[1];
open (DATA, "<$genome") || die "did not open $genome (second file). Call perl:\nperl \
buildTSLndSeq.pl cdna.fa genome.gtf\n";

#1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";

#1       havana  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";                                               

#1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";

my @linharray; my @infoarray; my @tslarray; my $genomeinfo; my $tsl; my @onlytsl; my %exontsl; my @arrayinfo; my @arrayENSG;
my $ENST; my $ENSG;

open (OUT1, ">ENST_seq_TSL.csv") || die "did not create tabelabanco.csv";
open (OUT2, ">semsequencias.csv") || die "did not create semsequencias.csv";


print OUT1 "id_tsl\tseq_acc\tsequence\ttsl_info\tgene_id\n";

while ($genomeinfo = <DATA>){ #open all line of genome archive 
    chomp($genomeinfo);
    next if ($genomeinfo =~ /^\#/); #ignore comments in the genome archive
    @linharray=(); #split the line
    @arrayinfo = (); #split nineth element of the line (info)
    @linharray = split(/\t/, $genomeinfo);
    my $type = $linharray[2]; # gene || transcript || exon
    my $comeco = $linharray[3]; # chr start position
    my $final = $linharray[4]; #  chr end position
    my $fita = $linharray[6]; # + || -
    my $info = $linharray[8]; # aditional information
    @arrayinfo = split (/;\s+/, $info);
    $ENSG = ""; # ENSG identificator
    $ENST = ""; # ENST identificator
    if ($type eq "transcript") {
	@onlytsl=(); # TSL information
	$arrayinfo[0] =~ /gene_id "(.*?)"/; #select only ENSG
	$ENSG = $1;
	$arrayinfo[2] =~ /transcript_id "(.*?)"/; #select only ENST
	$ENST = $1;
	
	if  ($arrayinfo[-1] =~ m/transcript_support_level/g){ #some info did not contain TLS information, so if it has execute...
	    
	    $arrayinfo[-1] =~ /transcript_support_level "(.*?)"/; #select only TSL number, some have aditional information
	    @onlytsl = split (/\s/, $1); #separe the information about TSL to catch only the number
	    $tsl = $onlytsl[0];
	}
	else {
	    $tsl="NULL";# para transcritos q nao tem TSL.
	}
	
	$exontsl{$ENSG}{$ENST}=$tsl; # for each ENSG and ENST I have a TSL number based on transcript information
    }
}

my $gene; # primary key
my $transc; # secondary key
my @listagene; #total ENSG 
my $count="0";
foreach $gene (sort keys %exontsl){ #keys = primary key
    unshift @listagene, $gene; 
    foreach $transc (keys %{$exontsl{$gene}}){ # keys = secondary key
	#print "\t$transc=$exontsl{$gene}{$transc}";

        #dark side of programming!!!

#	@$gene = "$transc=$exontsl{$gene}{$transc}"; # the $gene-named array contain every transcript and the TSL information about it.
	if ($seqhash{$transc}) {
	    $count++;
	    print OUT1 "$count\t$transc\t$seqhash{$transc}\t$exontsl{$gene}{$transc}\t$gene\n"; # ENST and it own TSL information
	}else{
	    print OUT2 "$transc\t$gene\t$exontsl{$gene}{$transc}\n";

	}
    }
}

=c
foreach my $line (@listagene) { #for each gene that exist in the genome 
    foreach my $line2 (@$line) { #Now $line is some ENSG and @$line ir a $gene-named array
	if ($line2 =~ m/'1'/g){ #Look only for ENSG that has at least one TLS=1
	    print "$line\n";  #ENSG
	} 
    }
=cut  

close(DATA);
close (OUT1);
close(OUT2);
exit;
