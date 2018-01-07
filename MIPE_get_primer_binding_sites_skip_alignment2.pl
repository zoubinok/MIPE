#! perl -w

###################################################################################################################################
#MIPE: Metagenome based community structure explorer and primer evaluation tool
#version 1.0.0
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
##Usage example: perl MIPE_get_primer_binding_sites_skip_alignment2.pl -i input.fasta -r Ref_16S.fasta -p primer1.txt -n 88888  -l 4 
# -i is the input fasta file.
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
#
##Output: manual_complete_primer1.txt
# 88888 is an example. Please change it.
####################################################################################################################################


# The modules used: 
# Getopt::std   parameters input
# Thread   multiple-thread computation in the alignment




use Getopt::Std;
use Thread;

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

getopt ('iprnl'); 

# -l -k are used in the extraction of primer-binding-sites. -k is the min percent of pbsite length to filt sequences not aligned in this part (now - is also regarded as aligned, but . is not).
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -n is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.

my $input_fasta=$opt_i;
my $input_primer=$opt_p;
my $sign=$opt_n;
my $input_ref_seq=$opt_r;
my $length_of_elongation=$opt_l;   




# Declare the template and taxonomy files used in Mothur jobs
# These names should change when analyzing functional gene sequences.
# Create a log file to record information and for debugging
# Log file name is like:  3423_test.fasta.log



system ("rm *tmp_* *.log*"); ##If you want to get more files, comment out the line.

#---Step 2: According to the information from the checkbox, append primers to the input primer file-------------------------#


#---Step 3: Read the input primer file and translate the degenerate primers-----------------------------------------------------------------#
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





my @refalign=split(/\.fa/,$input_ref_seq);
my $candidate_align = 'candidate_' . $refalign[0] . '.align';


my @in_path=split(/\//,$input_fasta);
(my $fasta_name="$sign"."_$in_path[$#in_path]")=~s/\.fasta$//;

my $new_dir_name= 'result_'. $fasta_name;
my %primer_sequence;
open SUITE, "$candidate_align" || die "$!"; #candidate_Ref_16S2.align
open MPRIMER, ">manual_complete_$input_primer" || die "$!";
my $n_elongation = 'N'x$length_of_elongation;
my @name_and_seq;
<SUITE>;
<SUITE>;
while ($name_and_seq[0]=<SUITE>)
{
    
	chomp($name_and_seq[0]);
	$name_and_seq[0] =~ s/[\n\r]//g;
    $name_and_seq[0] =~ s/\s+//g;
    $name_and_seq[0]=uc($name_and_seq[0]);
	
	$name_and_seq[1]=<SUITE>;
	chomp($name_and_seq[1]);
	$name_and_seq[1] =~ s/[\n\r]//g;
    $name_and_seq[1] =~ s/\s+//g;
	$name_and_seq[1] =~ s/\.+//g;

    $name_and_seq[1]=uc($name_and_seq[1]);
	
    if ($name_and_seq[0]=~ />(.+)/){$name_and_seq[0]=$1;}else {print "ERR:$name_and_seq[0]. A sigle sequence in $candidate_align should only has two lines and there should be no null string.\n"; die "$!";}
	 
	 $primer_sequence{$name_and_seq[0]}= $n_elongation . $name_and_seq[1] . $n_elongation;
   
   my $pname=$name_and_seq[0];

   my $pbsite_file_name=$new_dir_name . "/$pname" . '_primer_binding_sites.fa';
   open PBSITE, "$pbsite_file_name" || die "$!";
   print "$pbsite_file_name\n";
   
   my $sequence_name='';
   while (my $line=<PBSITE>)
{
  
  chomp ($line);
  if ($sequence_name =~ /$ref_seq_name/)
  {	
      print MPRIMER "Reference$pname:$line\n";
	  print MPRIMER "$pname:$primer_sequence{$pname}\n";

	  $line='';
  }
  if ($line =~ s/>//)
  {
      $sequence_name=$line;
  }
   }
   $pbsite_file_name='';
   @name_and_seq='';
   close PBSITE;
}

close SUITE;
close MPRIMER;
print "Please edit manual_complete_$input_primer\n";


#---Step 4: Align the fasta file with LSU/SSU in both original and reverse directions (multiple threads may be used)----------------------------------#


# Get the fasta file name and add a random number to make the job unique and easy to recognize:


# Create a result directory to save the files.


#----Step 5: According to the searchscores in the report files, pick seqs-----------------------------------------------------------#


#---Step 6: Add reference sequence to the picked bacteria/archaea fasta---------------------------------------------------------#


#---Step 7: Align the primer with the reference sequence via Mothur-------------------------------------------------------------# 

#---Step 8: Find primer-binding-sites according to the align result-------------------------------#

#---Stage 2: output match_type of each sequence for each primer-----------------------------------------------------------------------#


#---Step 1: Make primer into template files and align the primer-binding-sites via Mothur-------------------------------------------------------------------#	



system ("rm  *.report *.logfile .*.swo .*.swp *8mer* *train* *sum*"); ##If you want to get more files, comment out the line.



