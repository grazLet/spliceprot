#!usr/bin/perl
#vinicius parreira
#13-10-2021#
#removing NMD-target 

#you have to change the IDS when change your analysed species 

$nmd_list = $ARGV[0]; # from NMD classifier
$my_fasta = $ARGV[1];

open (IN, $nmd_list) || die "did not open $nmd_list\n";
open (IN2, $my_fasta) || die " did not open $my_fasta\n";

while (<IN>){
    chomp ($_);
    $line= $_;
    @targets = split (/\t/, $line);
    $ids = $targets[1]; ## human or mouse
   # $ids = "ENSRNOT". $targets[1]; ## rat
    $hash_nmd{$ids} = "1";
#    print "$ids\n"; # debug 
}



while (<IN2>){
    chomp $_;
    $line2 = $_;
    
    if ($line2 =~ /^>/){
	$header = $line2;
	$header =~ s/\]/\] /g;
	$header =~ s/\[/ \[/g;
     }else{
	$sequencia = $line2;
	@cabecalho = split (/\s+/, $header);
	@ids=();
	foreach $a (@cabecalho){
	    if ($a =~ /ENS/){
		$a =~ s/\[//g;
		$a =~ s/\]//g;
		#$a =~ s/\s+//g;
	#	print "$a\t"; # ===debug ok
		push @ids, $a;
	    }
	}
	# into the else conditioning - we capture the right header and sequence
 
	#	print "\n"; # === debug ok
	foreach $b (@ids){
	    $test = 0;
#	    print "$b ===> $hash_nmd{$b}\n"; ##debug
	    if ($hash_nmd{$b} eq "1") {
		$header = "NMDT".$header;

		next;
	    }else{
		next;
	    }next;
	}
	#out of the @ids foreach loops
	
	if ($header =~ /NMDT/){
	    $header =~ s/NMDT//g;
#	    print "$header\n$sequencia\n";#debug
	    $header = "";
	    $sequencia= "";
	    next;
	}
	else{
	    print "$header\n$sequencia\n";
	    $header = "";
	    $sequencia= "";
	}next;
    }# closing else fo $sequencia and return to while loop to get another header 
   
    next;

}


exit;
