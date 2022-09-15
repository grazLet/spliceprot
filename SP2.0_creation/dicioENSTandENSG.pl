#! usr/bin/perl
use strict;
use warnings;
####Vinicius Parreira 22.02.19####
###Dictionary ###
###Leticia Mattos 13.09.19###
###Mouse atualization###
#print "===== Download of cDNA list from  GRCh38 release 100 Ensembl=====\n";

print "==========\nperl dicioENSTandENSG.pl filein.fa fileout.tsv\n==========\n";

# you need download the fasta file for each sepecie to generate each dictionary.
#system "wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";(versao0.1)

#system "wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"; 

#system "wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.fa.gz";

my $filein = $ARGV[0];
open (IN, "<$filein") || die "did not open $filein";

#open (IN, "Homo_sapiens.GRCh38.ncdna.all.fa") || die "Não abriu Homo_sapiens.GRCh38.cdna.fa";(versão0.1)
my $fileout = $ARGV[1];

open (OUT, ">$fileout.tsv") || die "did not create $fileout.tsv";

my @head;

print "=====ENST vs ENSG list=====\n";
while (my $line = <IN>){
    chomp($line);
    if ($line =~ /^\>/){
       	@head = split(/\s+/, $line);
	$head[3]=~ (s/gene://);
	$head[3]=~ (s/\.\d{1,2}//);
	$head[0]=~ (s/>//);
	$head[0]=~ (s/\.\d{1,2}//);
	print OUT "$head[0] $head[3]\n";
    }
}
print "===Create $fileout.tsv===\n";
close (IN);
close (OUT);
