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
##Usage example: perl MIPE_get_primer_binding_sites.pl -i input.fasta -r Ref_16S.fasta -p primer1.txt -n 88888 -l 4 -k 01 
# -i is the input fasta file.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
# -l -k are used in the extraction of primer-binding-sites. -l should be 0-5. -k should be (0,1].
#
##Output directory: ./result_88888_input
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

getopt ('iprlkn'); 

# -l -k are used in the extraction of primer-binding-sites. -k is the min percent of pbsite length to filt sequences not aligned in this part (now - is also regarded as aligned, but . is not).
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -c is used to determine whether the mode of multiple threads is on.
# -n is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.

my $input_fasta=$opt_i;
my $input_primer=$opt_p;
my $input_ref_seq=$opt_r;
my $length_of_elongation=$opt_l;   
my $min_percent_of_pbsite_length=$opt_k;  
my $sign=$opt_n;


my $candidate_primer_to_ref = 'candidate_' . $input_ref_seq;


# Declare the template and taxonomy files used in Mothur jobs
# These names should change when analyzing functional gene sequences.
# Create a log file to record information and for debugging
# Log file name is like:  3423_test.fasta.log


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
       print LOG "$name\t";
    }
   print LOG "\n";
}

close SUITE;
close PRIMER;

#---Step 4: Align the fasta file with LSU/SSU in both original and reverse directions (multiple threads may be used)----------------------------------#


# Get the fasta file name and add a random number to make the job unique and easy to recognize:

my @in_path=split(/\//,$input_fasta);
(my $fasta_name="$sign"."_$in_path[$#in_path]")=~s/\.fasta$//;

# Create a result directory to save the files.

my $new_dir_name= 'result_'. $fasta_name;
if (-e "$new_dir_name"){die "ERR: Please make sure that $new_dir_name does not exist.$!\n";}
print LOG "ERR: Please make sure that $new_dir_name does not exist.$!\n";
mkdir ("$new_dir_name", 0755);
system ("rm *tmp_* *.log*"); ##If you want to get more files, comment out the line.

#----Step 5: According to the searchscores in the report files, pick seqs-----------------------------------------------------------#


#---Step 6: Add reference sequence to the picked bacteria/archaea fasta---------------------------------------------------------#


#---Step 7: Align the primer with the reference sequence via Mothur-------------------------------------------------------------# 

open PRIMER, "complete_$input_primer" || die "$!";
while (<PRIMER>)
{
  chomp;
  s/\s+//g;
  my @name_and_seq= split /:/,$_;
	
  if ($name_and_seq[0]=~/R/i)
  {
    my $rev_seq = get_rev_complement ($name_and_seq[1]);
    $primer_sequence{$name_and_seq[0]}=$rev_seq;
  }
  else
  {
    $primer_sequence{$name_and_seq[0]}=$name_and_seq[1];
  }
}

close PRIMER;

#---Step 8: Find primer-binding-sites according to the align result-------------------------------#

# Create hashes to save the positions of primer-binding-sites.
# first_base_on_ref & last_base_on_ref indicate the starting and ending base of the primer on the reference sequence.
# primer_name=>base_number
# start_position & end_position indicate the starting and ending position of primer binding site of each sequence in the 
# filter.fasta file (in the aligned format). 
# frist_base_on_ref & last_base_on_ref are not used in later steps, but may be used in the future.


my %first_base_on_ref;
my %last_base_on_ref;
my %start_position;   
my %end_position;   

(my $primer_align_result = $candidate_primer_to_ref) =~ s/\.fasta/.align/;

open ALIGN, "$primer_align_result" || die "$!"; 

my $sequence_name = "";
my @ref_bases;
my $ref_seq;

while (my $line=<ALIGN>)
{
  
  print LOG "Enter ALIGN candidate_align file\n";
  chomp ($line);
  if ($sequence_name =~ /Reference/)
  {	
      $ref_seq=$line;
      print LOG "ref_seq:\t$sequence_name\t$ref_seq\n";
      @ref_bases = split //, $line;
      $sequence_name="";
  }
  elsif ($sequence_name ne "")
  {
      @primer_bases = split //, $line;
      my $total_base=$#ref_bases+1;
      my $first_base=0;  
      my $start_position=0;
      my $end_position=$total_base;
      my $base_for_count=$ref_seq;
      my $last_base=($base_for_count=~tr/ATCGUatcgu//);
      print LOG "ref_base:$#ref_bases\nlast_base:$last_base\n"; 
     
      for (my $i=0;$i<=$#ref_bases;$i++)
      {
          if ($ref_bases[$i]=~ /[ATCGU]/i)
	  {
              $first_base ++ ; 
          }
          if ($primer_bases[$i]=~ /[ATCGU]/i)
	  {        
             $start_position=$i +1;
             last;
          }
      }
    
      for (my $i=$#ref_bases;$i>=0;$i--)
      {
          if ($primer_bases[$i]=~ /[ATCGU]/i)
	  {
             $end_position=$i +1;
             last;
          }
         if ($ref_bases[$i]=~ /[ATCGU]/i)
	 {
             $last_base -- ;
         }        
      } 
    
      $first_base_on_ref{$sequence_name}=$first_base;
      $last_base_on_ref{$sequence_name}=$last_base;
      $start_position{$sequence_name} = $start_position;
      $end_position{$sequence_name} = $end_position;
      print LOG "$sequence_name:\tstart:$start_position\tend:$end_position\n";

      if ($start_position > $length_of_elongation ) 
      {
         $start_position{$sequence_name} = $start_position - $length_of_elongation ;   
      }
      else
      {
         $start_position{$sequence_name} = 1;
      }
    	
      if ($end_position <= $total_base -$length_of_elongation)
      {	
         $end_position{$sequence_name} = $end_position + $length_of_elongation ;   
      }
      else
      {
         $end_position{$sequence_name} = $total_base;
      }

      $sequence_name ="";   
  }
  
  if ($line =~ s/>//)
  {
      $sequence_name=$line;
      print LOG "$sequence_name\n";
  }
}

print LOG "Leave ALIGN\n";
close ALIGN;

#---Step 9: Write primer binding site positions into log file-----------------------------------------------------------------------#

# Create an array to save the name of each nondegenerate primer.
# This array will be frequently used in later steps.

my @primer_names = keys %start_position;
foreach my $pname (@primer_names) 
{
    print LOG "$pname\tPosition in the aligned format\tBase number on the reference sequence\n".
    "Start:$start_position{$pname}\t$first_base_on_ref{$pname}\tEnd:$end_position{$pname}\t$last_base_on_ref{$pname}\n";
}

#---Step 10: Extract primer binding sites from the filter.fasta file and save in the .fa files--------------------------------------------------------------------#

open FILT, "ssu_merge_$sign.pick.filter.fasta" || die "$!";
{
local $/= ">";	
while (my $line=<FILT>)
{
    print LOG "At  FILT $.\n";
    chomp($line);
    next unless $line;
    my @name_and_seq= split /\n/, $line; 
    my $seq_name=$name_and_seq[0];
    unless ($seq_name)
    {
      print LOG "Invalid:$line\n";
      next;
    }
    $seq_name =~ s/\s+/_/g;
    $seq_name =~ s/_$//;
    (my $aligned_seq=$name_and_seq[1])=~ s/\s//g;
    my $length_of_aligned_seq = length ($aligned_seq);	 
    unless ( $seq_name || $aligned_seq )
    {
      print LOG "Invalid sequence:$seq_name\n$aligned_seq\n";
      next;
    }
   
   print LOG "primer_name:@primer_names\n";   
   foreach my $pname (@primer_names) 
   {
      my $pbsite_file = "$pname".'_primer_binding_sites.fa';
      my $no_pbsite_file = "$pname" . '_no_pbsite_sequence.list';

      open PBSITE, ">>$new_dir_name/$pbsite_file" || die "$!";
      open NOPB, ">>$new_dir_name/$no_pbsite_file" || die "$!";
      
      print LOG "Enter: $pbsite_file\n";

      if ($start_position{$pname}<$length_of_aligned_seq)
      {
          my $seq_extracted=substr($aligned_seq,$start_position{$pname}-1,$end_position{$pname}-$start_position{$pname}+1);
          $seq_extracted= uc ($seq_extracted);

          my $seq_for_count=$seq_extracted;
          my $total_nt= ($seq_for_count=~tr/ATCGUatcgu//); 
    
          if ($seq_extracted =~ s/[^\.\-ATCGU]//g)   #delete N and other characters
	  {
                print LOG "Unrecognized character found in seq:$seq_name\n";
          }
    
          if ($total_nt >=  $min_percent_of_pbsite_length * $primer_length{$pname})

	  {
               print PBSITE ">$seq_name\n$seq_extracted\n";
          }
	  else
	  {   
               print NOPB ">$seq_name\n";
          }
      }else{print PBSITE "\n";}
      close PBSITE;
      close NOPB;
  } 
    
}
}


close FILT;



#---Stage 2: output match_type of each sequence for each primer-----------------------------------------------------------------------#


#---Step 1: Make primer into template files and align the primer-binding-sites via Mothur-------------------------------------------------------------------#	


print LOG "primer_names:@primer_names\nStart to output match types\n";

foreach my $pname(@primer_names)
{
	
	print LOG "Enter\n";
	my $if_reverse;
	if ($pname=~/R/){$if_reverse='Reverse';}else{$if_reverse='Forward';}

	open TEMPLATE,">$pname:tmp_primer_template";
        
        # Insert gaps between nucleotides of primers. Because there may be insertions in some sequences' primer-binding-sites.

        my @primer_bases=split //,$primer_sequence{$pname};
	my $template_primer=join '----', @primer_bases;
        $template_primer= 'NNNNN'.$template_primer.'NNNNN';

        print TEMPLATE ">$pname\n$template_primer\n";
	close TEMPLATE;

        my $pbsite_file_name=$new_dir_name . "/$pname" . '_primer_binding_sites.fa';
        print LOG "Open:$pbsite_file_name\n";

	system("cp $pbsite_file_name $pname:tmp_pbsite_candidate.fa");
	
	open CANDI, ">>$pname:tmp_pbsite_candidate.fa" || die "$!";
	print CANDI ">primer_sequence\n$template_primer\n";
	close CANDI;

	
	system("mothur \"#align.seqs(candidate=$pname:tmp_pbsite_candidate.fa,template=$pname:tmp_primer_template)\"");
	system("mothur \"#filter.seqs(fasta=$pname:tmp_pbsite_candidate.align)\"");
}



system ("rm  *.report *.logfile .*.swo .*.swp *8mer* *train* *sum*"); ##If you want to get more files, comment out the line.

close LOG;


