#! perl -w    
        open LOGP, ">$ARGV[0].log" || die "$!";
		#$startalignment=$ARGV[2];#1 is min.
		#$endalignment=$ARGV[3]; #The last character which is aligned.
		#$endnucl = '';
		#$test = $endnucl . $startalignment . $endalignment;
		#print "$test\n";
		#print "Please make sure a single fasta only has two lines\n";
#usage: perl hmmreformat.pl HMMSTO hmmin.align hmmout.fasta
#reformat hmm out to fit MIPE		
#>SRR949316.363770432  
#..atgc..---...-AGAG...GTACCGG-----cag...........
#$i $startposition is 15, $endposition is 28
#012345678901234567890123456789012345678901234567
#alignemnt $startalignment is 9, $endalignment 34
#123456789012345678901234567890123456789012345678
#	 system("hmmbuild --dna $pname:tmp_primer_template.'hmm' $pname:tmp_primer_template");
#	system("hmmalign -o $pname:tmp_pbsite_candidate.sto $pname:tmp_primer_template.hmm $pname:tmp_pbsite_candidate.fa");

open HMMSTO,"$ARGV[0]" || die "$!";		
$stogcrf = '';
while($sto=<HMMSTO>){
if($sto =~ /^#=GC RF(.+)/){
$stogcrf = $stogcrf . $1;
chomp($stogcrf);
$stogcrf =~ s/[\n\r\s]//g;
}
}
close(HMMSTO);

@gcrf=split(/x/,$stogcrf);
	$startalignment=length($gcrf[0])+1;#1 is min.
	$endalignment=length($stogcrf)-length($gcrf[-1]);#The last character which is aligned.
   print LOGP "$sto";
   print LOGP "startalignment: $startalignment endalignment: $endalignment\n";
   print "startalignment: $startalignment endalignment: $endalignment\n";

   #=GC RF                                                                                          ..............xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.............
  
#--------------------------------------------------------------------------------

		
system("/home/server/hmmer-3.1b2-linux-intel-x86_64/easel/miniapps/esl-reformat afa $ARGV[0] > $ARGV[1]"); #change PATH based on your system.
#system("perl quhuanhang.zb.pl $pname:tmp_pbsite_candidate_normal.align $pname:tmp_pbsite_candidate_one.align");

#--------------------------------------------------------------------------------

	
	
	
	
	
$file=$ARGV[1];
$fileout=$ARGV[1] . '_one';
open(HMMALIGN, "$file") || die "$!";
open(OUT, ">$fileout");
$fasta = <HMMALIGN>;
chomp($fasta);
$fasta =~ s/[\n\r]//;
print OUT "$fasta\n";
while($fasta=<HMMALIGN>){
   chomp($fasta);  
  if($fasta =~ /^>/){
  print OUT "\n";
  print OUT "$fasta\n";
  }else{
    print OUT "$fasta";
  }  
}
close(HMMALIGN);
close(OUT);
#--------------------------------------------------------------------------------

open HMMOUT, ">$ARGV[2]" || die "$!";
open(HMMIN, "$fileout");
while ($in=<HMMIN>){
    #print title
	if($in =~ /^>/){ 
    print HMMOUT "$in";
	print LOGP "$in";
    $in = '';
  }else{
   #get nucl and find err.
    chomp($in);
    @line=split //,$in;
	$lengthline=length($in);
	if($in =~ /[^\.\-ATCGUNatcgun]/){print "ERR: illegal nucleotide:$in\n ";}
	#print $lengthline;
	
	#step1: extract aligned part.
	#1.1 find the position of first nucl.
for(my $i=$startalignment-1;$i<$endalignment;$i++)
	{
	   if($line[$i] =~ /[ATCGUNatcgun]/){
	   $startposition =$i;
	   print LOGP "startposition: $startposition+1\n";
	   last;
       }
    }
	#1.2 find the position of last nucl.
	for(my $i=$endalignment-1;$i>=$startalignment-1;$i--)
	{
	   if($line[$i] =~ /[ATCGUNatcgun]/){
	   $endposition =$i;
	   print LOGP "endposition: $endposition+1\n";
	   last;
       }
    }
	#1.3 tr/\./-/
	$alignment ='';
	for(my $i=$startposition;$i<=$endposition;$i++)
	{
	   $alignment=$alignment . $line[$i];
     }
	 $alignment=~tr/\./-/; 
     print LOGP "alignment: $alignment\n";
     #step2: extract nucls before aligned part.
	$startnucl = '';
    for(my $i=0;$i<$startalignment;$i++)
	{
	   if($line[$i] =~ /[ATCGUNatcgun]/){
	   $startnucl=$startnucl . $line[$i];
       }
    }
     print LOGP "startnucl: $startnucl\n";
	 @startnucl_rev=split //, reverse($startnucl);
	 $length_start=length($startnucl);
	 $beforealignment_rev ='';
	 $j=$startposition;
	 for(my $i=0;$i<$length_start;$i++)
	{
	   $j--;
	   if($j>=$startalignment-1){
	   if($line[$j] =~ /\./){$beforealignment_rev=$beforealignment_rev . '-';$i--;}else{$beforealignment_rev=$beforealignment_rev . $startnucl_rev[$i];}
	   }else{
	   $beforealignment_rev=$beforealignment_rev . $startnucl_rev[$i];
	   }
    }
	    $beforealignment_rev=$beforealignment_rev . ('.' x ($startposition-length($beforealignment_rev)));
		print LOGP "beforealignment: " . reverse($beforealignment_rev) . "\n";
	
	 
	 
	 
	 
    #step3: extract nucls after aligned part.
	$endnucl = '';
    for(my $i=$endalignment;$i<$lengthline;$i++)
	{
	   if($line[$i] =~ /[ATCGUNatcgun]/){
	   $endnucl=$endnucl . $line[$i];
       }
    }
     print LOGP "endnucl: $endnucl\n";
	 @endnucl=split //, $endnucl;
	 $length_end=length($endnucl);
	 $afteralignment ='';
	 $j=$endposition;
	 for(my $i=0;$i<$length_end;$i++)
	{
	   $j++;
	   if($j<=$endalignment-1){
	   if($line[$j] =~ /\./){$afteralignment=$afteralignment . '-';$i--;}else{$afteralignment=$afteralignment . $endnucl[$i];}
	   }else{
	   $afteralignment=$afteralignment . $endnucl[$i];
	   }
    }
	    $afteralignment=$afteralignment . ('.' x ($lengthline-$endposition-1-length($afteralignment)));
		print LOGP "endalignment: $afteralignment\n";
	#step4: cat 3 parts
	$alignmentcat=reverse($beforealignment_rev) . $alignment . $afteralignment;
	$alignmentcat =~ tr/atcgun/ATCGUN/; #If u comment out this line, you can find whether unaligned nucls are well arranged. 
   print HMMOUT "$alignmentcat\n";
  }
}
close HMMIN;
close HMMOUT;
close LOGP;

#>SRR949316.363770432  
#..............................................................................................atgc..................................------------------------------------...---------------------------...---------------------...---...---------------...---------------...-----------------------------------------...-TGGATCGATTGGAAGGATCGGCAGTTTTGGGTCACGGTCACACCAATCGTAGAG...GTGATGTATCCCGGGGCAATCATGTACTACTTCTGGACGTTCTACCGG---------------------------......---------------------------------------------...---------------------------------------------------------------------------------------------------------------------------------------------------------------...------...------------------------------------------------------------------------------------------------...---------------------...------------...---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------cag...............................................................................................................................................................................................................................................................................................................
#>SRR949316.364902962  
#....................................................................................................................................------------------------------------...---------------------------...---------------------...---...---------------...---------------...-----------------------------------------...-------------------------------------------------------...------------------------------------GTGAACTATCGCCAGCCCTTTGGCGCGACGATCACCATT......CTGGCGTTGCTGGCCGGGAAATGGGTCACGATCGTGGCGGCCTGG...TGGTGGTGGTCGAACTACCCCTATAACTTC---------------------------------------------------------------------------------------------------------------------------------...------...------------------------------------------------------------------------------------------------...---------------------...------------...---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------..................................................................................................................................................................................................................................................................................................................
#>NSFW
#....................................................................................................................................ATGTTTAGAACCGATGAGATTATTAAGGCCTCGAAG...TTGCCGCCGGAAGGCGTGGCGATGTCC...AGGCACATCGATTACATTTAT...TTCattCCGATCTTGTTTGTG...ACCATCGTGGGCACG...TTCCACATGCACTTCGACCTGCTGGCGGGCGACTGGGATTT...CTGGATTGATTGGAAGGATCGGCAGTGGTGGCCGATCGTGACGCCCGTCACGGCG...ATTACCTTTTGTGCGGCTCTCCAATATTACAATTGGGTGAACTATCGCCAGCCCTTTGGGGCGACGATCACCATC......CTGGCGCTGCTGGCCGGGAAGTGGGTGACGATCGTGGCCGCCTGG...TGGTGGTGGTCGAACTATCCGTATAACTTCGTCATGCCGTCCACTCTGCTGCCGAGCGCGCTGGTGCTGGACATCGTCTTGTTGTTGACCCGGAACTGGACGTTGACGGCGGTGATCGGCGCGTGGCTGTTTGCCGCGTTGTTCTATCCGACCAACTGG...GCCATC...TTCGCCTATAGCCATACGCCGCTCGTCATCGACGGCACCTTGCTGTCGTGGGCGGACTACATGGGCTTTGCGTATGTTCGGACCGGGACGCCGGAG...TACATTCGGATGATCGAAGTG...GGCTCGCTACGC...ACGTTCGGCGGGCACAGCACCATGATCTCCTCGTTCTTCGCGGCGTTCGCGTCGTCGTTGATGTACATCCTGTGGTGGCAGTTCGGGAAGTTCTTCTGCACATCCTACTTCTACCTGACGGATGACCGGCAGCGGACGACCAAAGTGTACGATGTGTTCGCGTACGCC---------------------acgttggcgaaccaggacaaggccaaagtcggagggaaagca........................................................................................................................................................................................................................................................................
#>AMOA_CXWL01000529.1cladeb
#....................................................................................................................................---------------------------GCCAGCAAG...CTGCCGCCGGAAGGCGTCAAACTGTCA...AGGATGGTAGACGCGGTGTAC...TTCCTTCCAATTGCATACCTG...GGAGTATTCGCGACC...TTTCATATGCACTTCGACTTGCTTGCGGGAGACTGGGATTT...CTGGATCGATTGGAAAGACCGGCAGTTCTGGGTAACGGTAACTCCTATTGTGGAG...GTCATGTACCCAGGTGCGATCATGTACTATTTCTGGACATTCTACCGGCAACCCTTTGGTGCAACGTTAAGCATC......ACAGGGCTCCTTGTCGGAAAGTGGATTACGGTCGTGTTGGCCTGG...TACTGGTGGGCAAATTTCCCGGTAAACTTCGTCATGCCTGCCACCATGGTCTCCAGCGCATTAATCCTGGACTGCACCCTGTTGCTGACCAGGAGCTGGATGATGACAGCGATATTCGGAGTATGGGCATTTGCGATGATGTTCTATCCTACGCAATAC...GCCCTG...TTCGGATACAGTAAGCAACCGATGGTCGTGGATGGTCAGCTGCTTTCCTTGGCCGACTACATGGGCTTCACCTATGTCCGGACGGGAACGCCTGAG...TATATCCGCATCATCGAAGTC...GGTTCGTTGCGT...ACATTCGGCGGACACACCGTTTGGATCTCCGCATTTTTCTCAGCGTTTGTTTCCATGTTGGTCTTCACCATTTGGTGGCAAGTCGGGCACTTCGTCGGCCAATCATTCTTTTGGGTCAAGGGCAAACGTGGCGTGCTGTCACGGAAAAATGACATCATGATCTGCGGC---------------------ACGGATGAGTTTGCCGCGAGACTTGGTACACCGGATGGGGTTGCACCAAGACTGACAGCAGCAGGCACACTTGCGGAAG...................................................................................................................................................................................................................................
