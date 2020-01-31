#! perl

###################################################################################################################################
#MIPE: Metagenome based community structure explorer and primer evaluation tool
#version 2.0.0
#    This program is free software.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    If you use MIPE in your work, please cite these manuscripts, as appropriate.
#    Mao DP, Zhou Q, Chen CY, Quan ZX (2012). Coverage evaluation of universal bacterial primers using the metagenomic datasets. BMC Microbiol.3;12:66.
#    Zou, B., Li, J., Zhou, Q., & Quan, Z. X. (2017). MIPE: A metagenome-based community structure explorer and SSU primer evaluation tool. PloS one, 12(3), e0174609.
#
###################################################################################################################################

###################################################################################################################################
##Dependence: Linux, mothur v1.33.3 
# Copy mothur into the working directory, then input: chmod 755 mothur
#
#
##Database in need: v119silva.SSU.fasta, v119silva.LSU.fasta, v119fangsilva.SSU.silva.tax, v119silva.SSU_PICK.fasta. They should be put in the directory of ./database. 
# If you want to change them, plasea change their pathways in MIPE_alignment_taxonomy.pl.
# Put input file, primer file, reference file into the working directory.
#
#
##Usage example: perl MIPE_matchtype_stat.pl -i input.fasta -r Ref_16S.fasta -p primer1.txt -w ssu_merge_88888.wang.taxonomy -n 88888
# -i is the input fasta file.
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
# -w is the taxonomy file which should be named as ssu_merge_xxxxx.silva.wang.taxonomy
#
##Output file: 
# 88888 is an example. Please change it.
####################################################################################################################################


# The modules used: 
# Getopt::std   parameters input
# Thread   multiple-thread computation in the alignment




use Getopt::Std;
use Thread;
use POSIX;

#---Subroutine1: Do alignment via Mothur----------------------------------------------------- ------------------------#

sub mothur_align
{
     my ($candidate,$template)=@_;
     system("./mothur \"#align.seqs(candidate=$candidate,template=$template)\"");
}

#---Subroutine2: Get reverse complement of the dna sequence--------------------------------------------------------#

sub get_rev_complement
{
     my $sequence=$_[0];
     $sequence=~tr/ACGTacgt/TGCATGCA/;      # tr operator
     my @nucleotides=split //,$sequence;
     join '',reverse @nucleotides;     
}

#---Subroutine3: Get the primer suite from the degenerate primer----------------------------------------------------#

sub get_primer_suite 
{
     my %degenerate=
     (
      "R" => "A,G",
      "Y" => "C,T",
      "M" => "A,C",
      "K" => "G,T",
      "S" => "G,C",
      "W" => "A,T",
      "H" => "A,T,C",
      "B" => "G,T,C",
      "V" => "G,A,C",
      "D" => "G,A,T",
      "N" => "A,G,C,T",
    );
  
   # $_[0] is like:  '8F:AGACK....'

   my $input_primer_seq = (split /:/,$_[0])[1];
   my @bases= split //,$input_primer_seq;

   # Create @suite to save all the nondegenerate primer formulations

   my @suite;
   $suite[0]="";
 
   # Add a base to the primer seqs each time. If the base is a degenerate one, then multiply the number of primer seqs by its degenerate degree,
   # and add its degenerate compositions respectively. 

   foreach my $base(@bases) 
   {
    if (exists $degenerate{$base})
    {
         my $num_of_primers=$#suite+1;
         my @degenerate_bases=split /,/,$degenerate{$base};
         my $degree=$#degenerate_bases; 
         
    	 for (my $i=0;$i<=$num_of_primers -1 ;$i++ ) 
	 {
             my $part_of_seq=$suite[$i];
             for (my $j=0;$j<=$degree;$j++)
	     {
               $suite[$i+$num_of_primers*$j]=$part_of_seq.$degenerate_bases[$j];
             }           
         }
    }
    else
    {
      for(my $i=0;$i<=$#suite;$i++ ) {$suite[$i] .= $base;}	
    }
  } 
  
  return @suite;
}

#---Stage 1: Get the primer-binding-sites and taxonomic information of each sequence-----------------------------------#


#---Step 1: Declare my variables and define Getopt parameters---------------------------------------------------------#

# HELP
=head1 Description

MIPE: Metagenome based community structure explorer and primer evaluation tool. Version 2.0.0

This script is used to analyze match types and statistics.

=head1 Usage

    $0 -i <input> -p <primer_file> -n <sign> -r <reference_file> -w <silva_wang_taxonomy>
    
=head1 Parameters

    -i  [str]   Input multiple nucleotide fasta file.
    -p  [str]   The file of primers. Primer format:"primer_name:NNNNNNNN". Forward primer name: end with F and no R or r is allowed in name. Reverse primer name: end with R and no F or f is allowed in name.
    -r  [str]   Reference fasta file which contains only one reference sequence whose name cannot appear in <input>. 
    -w  [str]   The taxonomy file which should be named as ssu_merge_xxxxx.silva.wang.taxonomy.
    -n  [int]   The sign of this sample.it should be a positive integer. 10000-99999 is recommanded.

=cut


getopt ('ipwnr'); 

# -w is the taxonomy file which should be named as ssu_merge_xxxxx.silva.wang.taxonomy
# -n is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.


my $input_fasta=$opt_i;
my $input_primer=$opt_p;
my $silva_wang_taxonomy=$opt_w;
my $sign=$opt_n;
my $input_ref_seq=$opt_r;

die `pod2text $0` if ((!$input_fasta) or (!$input_primer) or (!$sign) or (!$input_ref_seq) or (!$silva_wang_taxonomy));


$log_file="$sign"."$input_fasta";
open LOG, ">>$log_file.log";

#---Step 2: According to the information from the checkbox, append primers to the input primer file-------------------------#

#---Step 3: Read the input primer file and translate the degenerate primers-----------------------------------------------------------------#

# Create three hashes to save the information of all the primers
# The keys (primer_name) are the names of each non-degenerate primer
# e.g.  8F0, 8F1, 8F2, 8F3
# %primer_length: primer_name => primer_length
# %primer_sequence: primer_name => primer_sequence
# %degenerate_primer: primer_name=> the number of all its compositions, that is, the number of all the non-degenerate primers in its primer suite
# e.g.  primer1:AATTK  $degenerate_primer{'primer1'} equals 2

my %primer_length;      
my %primer_sequence;         
my %degenerate_primer;  

# Create a new primer file (SUITE) which contain all the non-degenerate primers

open PRIMER, "$input_primer" || die "$!";
open SUITE, ">complete_$input_primer" || die "$!";
my @primer_names;

while (<PRIMER>)
{
    print LOG "Get the primer suite of the degenerate primers:\n";
    chomp;
    s/\s+//g;
    $_=uc($_);
    next unless /[a-zA-Z]/;
    my @name_and_seq = split /:/,$_;
    my @all_primers = get_primer_suite ($_);
    $degenerate_primer{$name_and_seq[0]} = $#all_primers;   
    for (my $i=0;$i<=$#all_primers;$i++) 
    {
       $name = $name_and_seq[0] . $i; 
       $primer_length{$name} = length ($all_primers[$i]);
       print SUITE "$name:$all_primers[$i]\n";
       push @primer_names, $name;

	   print LOG "$name\t";
    }
   print LOG "\n";
   push @primer_names, @degenerate_primer_names;

}

close SUITE;
close PRIMER;


#---Step 4: Align the fasta file with LSU/SSU in both original and reverse directions (multiple threads may be used)----------------------------------#


#---Step 5: Select bacteria sequences according to the bootstrap value----------------------------------------------#

my %taxo_info;  

open TAXO,"$silva_wang_taxonomy" || die "$!"; #To cover the bug of mothur v1.33, ".silva.wang.taxonomy" is canned format.
while(my $line=<TAXO>)
{
	my @item=split(/\s+/,$line);
	my @taxo=split(/;/,$item[1]);
	$taxo_info{"$item[0]"}="$taxo[0];$taxo[1];$taxo[2];$taxo[3];$taxo[4];$taxo[5]"; ##$taxo_info{"$item[0]"}="$taxo[0];$taxo[1];$taxo[2];$taxo[3];$taxo[4];$taxo[5]"; is also OK.
}
close TAXO;




#---Step 6: Add reference sequence to the picked bacteria/archaea fasta---------------------------------------------------------#
my $ref_seq_name;
open REF,"$input_ref_seq" || die "$!";
while(my $line=<REF>)
{

	if ($line =~ />/)
	{
	    chomp ($line);
	    $line=~ s/\s+/_/g;
	    $line=~/>(.*)/;
	    $ref_seq_name=$1;
	}
}
close REF;


#---Step 7: Align the primer with the reference sequence via Mothur-------------------------------------------------------------# 

#---Step 8: Find primer-binding-sites according to the align result-------------------------------#

#---Step 9: Write primer binding site positions into log file-----------------------------------------------------------------------#

#---Step 10: Extract primer binding sites from the filter.fasta file and save in the .fa files--------------------------------------------------------------------#

#---Stage 2: output match_type of each sequence for each primer-----------------------------------------------------------------------#

#---Step 1: Make primer into template files and align the primer-binding-sites via Mothur-------------------------------------------------------------------#	

#---Step 2: Compare the primer with each primer-binding-site in aligned format and output matching type-----------------------------------#			


my @in_path=split(/\//,$input_fasta);
(my $fasta_name="$sign"."_$in_path[$#in_path]")=~s/\.fasta$//;

my $new_dir_name= 'result_'. $fasta_name;
my %primer_sequence;
open SUITE, "complete_$input_primer" || die "$!";


while (<SUITE>)
{
    chomp;
    s/\s+//g;
    $_=uc($_);
    next unless /[a-zA-Z]/;
    my @name_and_seq = split /:/,$_;
    if ($name_and_seq[0]=~/R/i)
    {
      my $rev_seq = get_rev_complement ($name_and_seq[1]);
      $primer_sequence{$name_and_seq[0]}=$rev_seq;
    }
	else
    {
     $primer_sequence{$name_and_seq[0]}=$name_and_seq[1];
   }
   my $pname=$name_and_seq[0];

   my $pbsite_file_name=$new_dir_name . "/$pname" . '_primer_binding_sites.fa';
   
      
	my %primer_binding_sites;

   open PBSITE, "$pbsite_file_name" || die "$!";
   print "$pbsite_file_name\n";
        my @lines=<PBSITE>;
        for (my $i=0; $i<$#lines; $i+=2)
	{
	    chomp ($lines[$i]);
	    chomp ($lines[$i+1]);
	    (my $pbsite_name=$lines[$i])=~s/>//;
            $primer_binding_sites{$pbsite_name}=$lines[$i+1];
   }
   close PBSITE;
   

        open FILT,"$pname:tmp_pbsite_candidate.filter.fasta" || die "$!";
        open LOGP, ">$pname.log" || die $!;

        my @names_and_pbsites=<FILT>;
	close FILT;
        print LOG "Create:$pname,match_type";
	open RES,">$new_dir_name/${pname}_match_type" || die $!;

	print RES "Seq_name\tForward or reverse\tMatch_type\tPrimer-binding site\tMisMatches in last 4 nt\t".
	"Mismatches\tIf matching\tCompleteness\tTaxonomic information\n";
	
	my $aligned_primer_template=$names_and_pbsites[$#names_and_pbsites];
	chomp($aligned_primer_template);
	
	my @primer_components=split //, $aligned_primer_template;
        print LOGP "template:$aligned_primer_template\n";
	print LOGP "number of fa lines: $#names_and_pbsites\n";
        
	my $if_reverse;
	if ($pname=~/R/){$if_reverse='Reverse';}else{$if_reverse='Forward';}

	for(my $i=1;$i<$#names_and_pbsites;$i+=2)
	{
	        print LOGP "$i\n";
		chomp($names_and_pbsites[$i-1]);		
		print RES "$names_and_pbsites[$i-1]\t$if_reverse\t";			
		chomp($names_and_pbsites[$i]);

# $mismatch is the number of mismatches of each primer-binding-site.
# $mismatch_in_last_4 is the number of mismatches in last 4 nt. 
# The "last" is relative because it is related with the primer direction (forward or reverse).
# In the output match types:
#   = the same nt as the primer
#   A(TCG) subsititution
#   a(tcg) insertion
#   d deletion
#   . missing
                
		print LOGP "$names_and_pbsites[$i-1]\n$names_and_pbsites[$i]\n";
        
		my $aligned_length=length($names_and_pbsites[$i]);
		my $aligned_length_4=$aligned_length-4;
		my $mismatch_in_last_4=0;
		my $mismatch=0;
		my $match_type;
		my $valid_bases=$primer_length{$pname};

                print LOGP "Aligned_length:$aligned_length\tvalid_bases:$valid_bases\n";

		if($aligned_length > $#primer_components+1) {
		  print LOG "$pname:aligned_length=$aligned_length\n#primer_components=$#primer_components\n";
		  $aligned_length=$#primer_components+1;
		}  #necessary?

		for(my $j=0;$j<$aligned_length;$j++)
		{
			my $single_character=substr($names_and_pbsites[$i],$j,1);
			if($single_character eq $primer_components[$j])
			{
				next if($primer_components[$j] eq '-');
				next if($primer_components[$j] eq '.');#add this line to ignore . 20180110
				$match_type .= '=';
			}
			else                      
			{
			  if (($primer_components[$j] eq "-") && ($single_character =~ /[ATCGU]/))
			  {   
				 my $insertion=lc($single_character);
				 $match_type .= "$insertion";
                                 $mismatch ++;
			  }elsif(($primer_components[$j] =~ /[ATCGU]/) && ($single_character eq "-")) #elsif($single_character eq "-")  change it
			  {
                                 $mismatch ++;
			         $match_type .= "d";
			  }elsif ( ($primer_components[$j] =~ /[ATCGU]/) && ($single_character=~/[ATCGU.]/)  )
			  {
			         $match_type .= "$single_character";
                                 $mismatch ++;
			  }
			  if (($primer_components[$j] =~ /[ATCGU]/) && ($single_character eq ".")){
			         $valid_bases --;
			  }
			}
		}
		

    
		my @match_type=split //, $match_type;
		my $type_length=length($match_type);

# $type_length does not necessarily equal to the length of the primer, because there may be insertions. 
# Count the mismatches in the last 4 nucleotides:		
		
		for ( my $k=0; $k<=$#match_type; $k++) 
	        { 
	  	    if ($match_type[$k] eq '=') {next;}	  	
	            elsif ($if_reverse eq "Forward" && $k>=$type_length -4)
		    {
	    	        $mismatch_in_last_4 ++;
	    	        print LOGP "Mismatch:$names_and_pbsites[$i-1]\t$k\t$type_length\t$match_type[$k]\n";
	    	    }elsif($if_reverse eq "Reverse" && $k<=3)
		    {
		        $mismatch_in_last_4 ++;
		    }
		}
		
		
		my $if_matching='N';
		if($mismatch<=1 && $mismatch_in_last_4 == 0)
		{
		    $if_matching='Y';
		}		
		
		my $percent_of_valid_bases=$valid_bases / $primer_length{$pname};
		(my $pbsite_name = $names_and_pbsites[$i-1])=~s/>//;
		
		if ($pbsite_name=~/$ref_seq_name/)
		{
		print RES "$match_type\t$primer_binding_sites{$pbsite_name}\t$mismatch_in_last_4\t$mismatch\t$if_matching\t".
		"$percent_of_valid_bases\n";		
	        }
		else
		{
		print RES "$match_type\t$primer_binding_sites{$pbsite_name}\t$mismatch_in_last_4\t$mismatch\t$if_matching\t".
		"$percent_of_valid_bases\t$taxo_info{$pbsite_name}\n";
                }


	}
	close LOGP;
	close RES;

}
close SUITE;
#---Step 3: Combine the match types of non-degenerate primers for the same position-----------------------------------------------------------#


my @degenerate_primer_names; 
@degenerate_primer_names = keys %degenerate_primer;

# degenerate_primer_names are like '8F' '338F'
# primer_names are like '8F0' '8F1' '8F2' '8F3' '338F0'


foreach my $pname (@degenerate_primer_names)
{	

	my $if_reverse;
        if ($pname=~/R/){$if_reverse='Reverse';}else{$if_reverse='Forward';}

        my @reference_bases;
	
# Create degenerate reference primer bases

	for (my $i=0;$i<=$degenerate_primer{$pname};$i++)
	{	
		my $non_degenerate_name = "$pname" . "$i";
		print LOG "727:$pname\t$i\t$non_degenerate_name\t$primer_sequence{$non_degenerate_name}\n";
		my @primer_bases = split //, $primer_sequence{$non_degenerate_name};
                print LOG "729:@primer_bases\n"; 
                for (my $n = 0; $n<=$#primer_bases; $n++){     	
    	             if (! (defined $reference_bases[$n]))
		     {
    	                 $reference_bases[$n]=$primer_bases[$n];
                     }
		     else
		     {
      	                 my $single_base = $reference_bases[$n];
      	                 if(! ($primer_bases[$n] =~ /[$single_base]/) )
			 {
                            $reference_bases[$n] .= $primer_bases[$n];
                         }          
                     }
                }   
        }  

# @Reference_bases are like ATCGG[AT]GC...
# [AT] means a degenerate base
# Open each non-degenerate primer match type file and combine them into one.
# The match type is changed due to the degenerate bases.
# The number of mismatches and the number of mismatches in last 4 are recounted.
# Whether the primer-binding-site is matching is rejudged.

         for (my $i=0;$i<=$degenerate_primer{$pname};$i++)
	 {
	     my $match_type_file=$pname . $i . '_match_type';
             open TYPE, "$new_dir_name/$match_type_file";
             open TMPTYPE, ">$new_dir_name/${match_type_file}_tmp";
   	     print TMPTYPE "Seq_name\tForward or reverse\tMatch_type\tPrimer-binding site\tMisMatches in last 4 nt\t".
	     "Mismatches\tWhether matching\tCompleteness\tTaxonomic information\tScore\n";	       
             while (<TYPE>)
	     {
    	         next if (/^Seq_name/);
    	         chomp;
		 my @items = split /\t/,$_;
    	         my @pbsite_bases = split //, $items[2];
    	         my $insertion_count = 0;
    	         my $mismatch_in_last_4 = $items [4];
    	         my $mismatch = $items[5];
    	         my $if_matching = $items[6];
		 print LOG "Invalid pbsite:$items[2]\n" unless (@pbsite_bases);
    	
    	         for (my $n = 0; $n<=$#pbsite_bases; $n++)
	         {
    		     if ($pbsite_bases[$n]=~ /[atcg]/)
		     {
    	                 $insertion_count ++;
    	             }		
    		     elsif ($pbsite_bases[$n] =~ /[ATCG]/)
		     {
    		         $degenerate=$reference_bases[$n-$insertion_count];
    		         print LOG "$pname\t$i\t$n\t$degenerate\n";    		  
    		         if ($pbsite_bases[$n] =~ /[$degenerate]/)
			 {
    		             $pbsite_bases[$n] = '=';
               		     $mismatch --;  #mismatches change thus
    		             if ($if_reverse eq "Forward" && $n > $#pbsite_bases - 4){$mismatch_in_last_4 --;}
    		             if ($if_reverse eq "Rerverse" && $n< 4 ){$mismatch_in_last_4 -- ;}
    		         }    		
    	             }
    	 	}

                if ( $mismatch_in_last_4 == 0 && $mismatch<=1 ) {$if_matching='Y';}
        	$items[4]=$mismatch_in_last_4;
    	        $items[5]=$mismatch;
    	        $items[6]=$if_matching;
         	$items[2] = join '', @pbsite_bases;
        	$item = join "\t", @items; 	
             
	        my $match_score=( 1 - ( $items[4]+$items[5] + 10 - ($items[7] * 10) ) / ( $#reference_bases + 14) ) * 100;    
                printf TMPTYPE "$item\t%6.2f\n", $match_score;    	
             }  

             close TYPE;
             close TMPTYPE;
        }
}


#---Step 4: Compare the match scores of each primer-binding site and pick out the best matching one------------------------------------------------------#	
foreach my $pname (@degenerate_primer_names)
{
	open TYPE, ">$new_dir_name/${pname}_match_type" || die "$!";	
	print TYPE "Seq_name\tForward or reverse\tMatch_type\tPrimer-binding site\tMisMatches in last 4 nt\t".
	"Mismatches\tWhether matching\tCompleteness\tTaxonomic Information\tScore\n";  
        my $primer_name_0="$pname" . '0' ; 
        my $word_count=`wc -l $new_dir_name/${primer_name_0}_match_type`;
        $word_count =~ /^([0-9]+)\s/;
        my $number_of_lines= $1;
##########
        my @lines, @match_score;
		$lines[0]=0;
		$lines[1]=0;
		$match_score[0]=0;
		$match_score[1]=0;

            my $j=0;
            my $i=2;
	        my $match_type_file=$pname . $j . '_match_type';
                open TMPTYPE, "$new_dir_name/${match_type_file}_tmp";	
                <TMPTYPE>;
				while (<TMPTYPE>){

				$lines[$i] = $_ ;
	  	      /(-?[0-9]+\.[0-9]+)$/m;##/\t([0-9]+.[0-9]+)\n$/ changed as /(-?[0-9]+\.[0-9]+)$/m
	  		    $match_score[$i] = $1;
	      	      $i=$i+1;
				  }  	  	  		
                close TMPTYPE;  	
                	
		###################  
# mark which non-degenerate primer

#	    my $primer_number=0;
#	    my $match_score=0;
#	    my @lines;  	
      	    for (my $j=1;$j<=$degenerate_primer{$pname};$j++){

  	        my $i=2;
            
	        my $match_type_file=$pname . $j . '_match_type';
                open TMPTYPE, "$new_dir_name/${match_type_file}_tmp";	
                       <TMPTYPE>;
				while (<TMPTYPE>){
	    	      #next unless ($. == $i);
				  my $linetmp=$_;
	  	      /(-?[0-9]+\.[0-9]+)$/m;##/\t([0-9]+.[0-9]+)\n$/ changed as /(-?[0-9]+\.[0-9]+)$/m
              
			  if ($match_score[$i] < $1)
		      {
	  		    $match_score[$i] = $1;
	  	        $lines[$i] =~ /^(.+?)\t/;
				my $name1=$1;
				$linetmp =~ /^(.+?)\t/;
				my $name2=$1;
                unless($name1 eq $name2){
				print TYPE "ERROR: match type files are not in the same order. See $match_type_file\n"; 
				print LOG "ERROR: match type files are not in the same order. See $match_type_file\n"; 
				die "ERROR: match type files are not in the same order. See $match_type_file\n";
				}
				$lines[$i] = $linetmp ;
				  }  	  	  		
	      	      $i=$i+1;
	        
			}
                close TMPTYPE;  	
                	
}         

        for (my $i=2; $i<=$number_of_lines; $i++){print TYPE "$lines[$i]";}

	 close TYPE;
}


push @primer_names, @degenerate_primer_names;

foreach my $pname (@degenerate_primer_names)
{
       for (my $i = 0; $i<=$degenerate_primer{$pname};$i++)
       {
           my $non_degenerate_name = $pname . $i;
      	   $primer_sequence{$pname} .= "$primer_sequence{$non_degenerate_name}\n";
  	   $primer_length{$pname} = $primer_length{$non_degenerate_name}; 	
       }  
}

#---Step 5: Output stats--------------------------------------------------------------------------------------------------------------#

foreach my $pname (@primer_names)
{
   open STAT, ">$new_dir_name/$pname.stat";
   open RES, "$new_dir_name/${pname}_match_type";

# @mismatch saves the numbers of pbsites with certain mismatches
# @complete saves the numbers of pbsites with completeness over certain value
# $matching saves the number of matching primer-binding-sites

   my @mismatch; 
   my $matching = 0;
   my @completeness; 
 
   print STAT "$pname:\n$primer_sequence{$pname}\nLength:$primer_length{$pname}\n\n";
 
   while (<RES>)
   {
      next if (/Seq_name\tForward/);
      my @line=split /\t/,$_;
   
      my $mismatch = $line[5];
      my $if_matching = $line[6];
      my $completeness = $line[7];
      $mismatch[$mismatch] ++;
      $matching ++ if $if_matching eq "Y";
      $completeness *= 10;
 
      for (my $i=0;$i<=10;$i++)
      {
           $completeness[$i]++ if $completeness>=$i;   
      }  
   }
 
   $matching=0 unless defined $matching; 

#---Table 1: Primer-binding-site numbers with given mismatches--------------------------------------------------#
 
    print STAT "Table 1\nMatching\t$matching\n";
 
    for (my $i=0;$i<=$#mismatch;$i++)
    {
 	$mismatch[$i] = 0 unless defined $mismatch[$i];
        print STAT "$i"." mismatch(es)\t$mismatch[$i]\n";  
    }

#---Table 2: Primer_binding-site numbers with completeness over given value--------------------------------------------------------------#
    print STAT "\nTable 2\n".
    "Completeness cutoff\tNumber of seq\tCoverage rates\n";
 
    for (my $j=1;$j<=10;$j++)
    {
        my $j_10 = $j/10;
        eval 
	{ 
	   my $coverage_rate = $matching / $completeness[$j];
           print STAT "$j_10\t$completeness[$j]\t$coverage_rate\n";
	};  
    } 
 
    eval {my $coverage = $matching / $completeness[0];}; 
    print STAT "Total\t$completeness[0]\t$coverage\n\n\n";
 
    close RES;
 
#---Table 3: Most frequent mismathes----------------------------------------------------------------------------# 

    open RES, "$new_dir_name/${pname}_match_type";
 
    my @match_type;  
    my %match_type;
    my $word_count=`wc -l $new_dir_name/${pname}_match_type`;
    $word_count =~ /^([0-9]+)\s/;
    my $number_of_pbsites= $1-1;
    
    while (<RES>)
    {
	  	next if /^Seq_name\t/;
	  	my @items = split /\t/,$_;
	  	push @match_type, $items[2];
	  	$match_type{$items[2]} ++ ;
    }
    close RES;


    my (@subA, @subT, @subG, @subC, 
       @insA, @insT, @insG, @insC, @del, @pts);
  	 	 
    foreach my $type (@match_type)
    {
        my  @type_bases = split //, $type;
        my $dot_count=0;
        foreach $type_base(@type_bases){$dot_count ++ if $type_base eq ".";}	
        my $completeness_tmp=1 - ($dot_count / $primer_length{$pname});
        my $insertion_count =0;
        for (my $v = 0; $v <= $#type_bases; $v++)
	{
   	    next if  $completeness_tmp<0.8;
   	    if ($type_bases[$v] =~ /[atcg]/){$insertion_count ++ ;}
   	    $position= $v -$insertion_count;
   	    if ($type_bases[$v] eq "A"){$subA[$position] ++;}
   	    if ($type_bases[$v] eq "T"){$subT[$position] ++;}   	
   	    if ($type_bases[$v] eq "G"){$subG[$position] ++;}   	
   	    if ($type_bases[$v] eq "C"){$subC[$position] ++;}
   	    if ($type_bases[$v] eq "a"){$insA[$position] ++;}   	
   	    if ($type_bases[$v] eq "t"){$insT[$position] ++;}
   	    if ($type_bases[$v] eq "c"){$insC[$position] ++;}   	   	   	
   	    if ($type_bases[$v] eq "g"){$insG[$position] ++;}
   	    if ($type_bases[$v] eq "d"){$del[$position] ++;}
   	    if ($type_bases[$v] eq "."){$pts[$position] ++;}   	   	   	
        }	
     }
  
     print STAT "Table 3: Most frequent mismatch\n";
 
  
     for (0 .. ($primer_length{$pname}-1))
     {
 	    my $position = $_ + 1;
 	    if ($subA[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tA-substitution\t$subA[$_]\n";}
 	    if ($subT[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tT-substitution\t$subT[$_]\n";}
 	    if ($subG[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tG-substitution\t$subG[$_]\n";}
 	    if ($subC[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tC-substitution\t$subC[$_]\n";}
 	    if ($insA[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tA-insertion\t$insA[$_]\n";}
 	    if ($insT[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tT-insertion\t$insT[$_]\n";}
 	    if ($insG[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tG-insertion\t$insG[$_]\n";}
 	    if ($insC[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tC-insertion\t$insC[$_]\n";} 	 	 	
 	    if ($del[$_] >= 0.05 * $number_of_pbsites){print STAT "$position\tDeletion\t$del[$_]\n";}	 	 	 	
     }     
 
#---Table 4: Most frequent unmatching primer-binding-sites/match_types--------------------------------------------------------------------------# 

    print STAT "\nTable 4: Most frequent unmatching primer-binding sites\n";
  
    foreach my $match_type ( keys %match_type )
    {
 	  if ( $match_type{$match_type} >= 0.02 * $number_of_pbsites)
	  {
 	  	print STAT "$match_type\t$match_type{$match_type}\n"; 	
 	  }
    }   
#---Table 5: Number of mismtaches at each position-------------------------------------------------------------------------------# 

   print STAT "\nTable 5: Number of mismatches at each position\n";

   for (1 .. $primer_length{$pname}){print STAT "\t$_";}
   print STAT "\nA-substitution";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$subA[$_]";}
   print STAT "\nT-substitution";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$subT[$_]";} 
   print STAT "\nC-substitution";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$subC[$_]";} 
   print STAT "\nG-substitution";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$subG[$_]";} 
   print STAT "\na-insertion";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$insA[$_]";} 
   print STAT "\nt-insertion";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$insT[$_]";} 
   print STAT "\nc-insertion";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$insC[$_]";}
   print STAT "\ng-insertion";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$insG[$_]";}
   print STAT "\ndeletion";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$del[$_]";}
   print STAT "\nmissing";
   for (0 .. ($primer_length{$pname}-1)){print STAT "\t$pts[$_]";}
   print STAT "\n\n";


#----Table 7: coverage rates in each phylogroup; (Completeness >= 100%)-----------------------------------------------#

   my @matching_taxo;
   my @unmatching_taxo;
   my %matching_count;
   my %unmatching_count;

   open RES, "$new_dir_name/${pname}_match_type";
   while (<RES>)
   {
       next if /^Seq_name\t/; 
       chomp;
       my @items = split /\t/,$_;
       next if $items[7]<1;##Completeness can be changed here. Now it is 100%.
       if ($items[6] eq 'Y')
       {
          push @matching_taxo, $items[8];
          print LOG "$pname:m:$items[0]\n";
          
       }
       else
       {
          push @unmatching_taxo, $items[8];
	  print LOG "$pname:n:$items[0]\n";
       }
   }
   close RES;
   
   my @total_taxo = ([@matching_taxo],[@unmatching_taxo]);
   my $which_group=0;

   foreach my $taxo_ref (@total_taxo)
   {
       foreach my $taxo (@{$taxo_ref})
       {
       	   print LOG "$taxo\t";
           my @taxo_split=split /;/, $taxo;
           my $taxo_order=1;
           foreach my $single_taxo (@taxo_split)
           {
              print LOG "$single_taxo\t";
              $single_taxo=~ s/\([0-9]{1,3}\)//;
	      $single_taxo .= "-$taxo_order";
	            print LOG "$single_taxo\n";
              $matching_count{$single_taxo}++ if $which_group == 0;
   	      $unmatching_count{$single_taxo}++ if $which_group ==1;
	      $taxo_order ++; 
	   }
      }
      $which_group ++;
   }

   my @count_record;
   foreach my $taxo (keys %matching_count)
   {
      my $match=$matching_count{$taxo};
      my $unmatch;
      
      if (exists $unmatching_count{$taxo})
      {
      $unmatch = $unmatching_count{$taxo};
      }
      else
      {
      $unmatch=0; 
      }
      
      my $percent=$match/($match+$unmatch);
      my $count_record=join "\t", $taxo, $match, $unmatch, $percent;
      push @count_record, $count_record;
   }
   foreach my $taxo (keys %unmatching_count)
   { 
      next if exists $matching_count{$taxo};
      my $unmatch=$unmatching_count{$taxo};
      my $count_record=join "\t", $taxo, '0', $unmatch, '0';
      push @count_record,$count_record;

   }
   
   for (my $i=1;$i<=6;$i++) ##for (my $i=1;$i<=6;$i++) ##if you change line 354 to get deeper taxonomy, you have to change this line.
   {      
      print STAT "Taxonomy\tmatch\tunmatch\tmatch_rates\n";
	  foreach my $line (@count_record)
      {
          next unless $line=~/[a-zA-Z]/;
          next unless $line=~/\-$i\t/;
          print STAT "$line\n";

      }
      print STAT "\n";
   } 
   close STAT;
}

system ("rm *.report *.logfile .*.swo .*.swp *8mer* *train* *sum* "); ##If you want to get more files, comment out the line.
#system("rm *tmp_*.* *tmp_*");

close LOG;


