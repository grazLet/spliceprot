#! usr/bin/perl

=c
Vinicius Parreira - Fiocruz PR
26-04-2019
This Script aims to remove some duplicate entry in the psl file.
All the transcripts that have more than one entry in that file was removed.
=cut
#perl removeduplicates.pl input.BEST.psl

use strict;
use warnings;

my $file1= $ARGV[0];
my $conferir="a"; # variable initialization 
my @repetidos;
my $count = 0;
my %hash = ();
my $count2 = 0;

open (IN, "$file1") ||  die "did not open $file1\n";
open (OUT, ">repeticoes.txt") || die "did not create.txt\n";
open (OUT2, ">$file1.out.semrep.psl") ||die "did not create $file1.out.semrep.psl\n";
open (OUT3, ">tempfile.psl") || die "did not create tempfile.psl\n";

#PSL file in ascending order
while (my $input = <IN>){
    chomp ($input);
    if ($count2 < 5){
	$count2++;
	print OUT2 "$input\n";
	next; }
    print OUT3 "$input\n";
}

close (OUT3);

system "sort -k10 tempfile.psl > tempfile.sorted.psl";
open (IN2, "tempfile.sorted.psl") || die "did not open tempfile.psl\n"; 

while (my $line1 = <IN2>){
    chomp ($line1);
    #if ($count < 5){ $count++; next; }
    my @array1 = split (/\s+/, $line1);
    my $name = $array1[9];
#    print "$name\n";
    if ($conferir eq $name) {
	push @repetidos, $name;
    }
    $conferir = $array1[9];
}
my $line3="";
foreach my $line2 (@repetidos){
    chomp ($line2);
    #remove the duplicates in the redundance file
    if ($line2 eq $line3){next;}else{
	print OUT "$line2\n";
	$hash{$line2}=1;
	$line3 = $line2;
    }
}

#foreach my $key (keys %hash){
#    print "$key == $hash{$key}\n";
#}
$count =0;
close (OUT);

open (IN, "$file1") ||  die "did not open $file1\n"; 
while (my $line3 = <IN>){
    chomp ($line3);
    if ($count < 5){ $count++; next;  }
#    print "$line3\n";
    my @array2 = split (/\s+/, $line3);
    my $name2 = $array2[9];
    #  print "$name2\n";
    #   print "$hash{'ENST00000631869.1'}\n";
    
    if ($hash{$name2}){
	next;#print "=====hash\n";
    }else{
	print OUT2 "$line3\n";
    }
    
}
close (IN);
close (IN2);
close (OUT2);
system "rm tempfile.psl";
system "rm tempfile.sorted.psl";
exit;
