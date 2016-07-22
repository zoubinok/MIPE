#! perl 

###################################################################################################################################
#MIPE: Metagenome based community structure explorer and primer evaluation tool
#Copyright (C) 2014-2015 Bin Zou, Jiefu Li  (14210700126@fudan.edu.cn)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
#    If you use MIPE in your work, please cite these manuscripts, as appropriate.
#    Mao DP, Zhou Q, Chen CY, Quan ZX (2012). Coverage evaluation of universal bacterial primers using the metagenomic datasets. BMC Microbiol.3;12:66.
#
###################################################################################################################################


###################################################################################################################################
##Dependence: Linux, mothur v1.33.3 
# Copy mothur into the working directory, then input: chmod 755 mothur
#
#
##Database in need: v119silva.SSU.fasta, v119silva.LSU.fasta, v119fangsilva.SSU.silva.tax, v119silva.LSU.silva.tax, v119silva.SSU_PICK.fasta.
# Put them and input file, primer file, reference file into the working directory.
#
#
##Usage: perl MIPE_main_program_metatranscriptome.pl -i input.fasta -p primer1.txt -r Ref_16S.fasta -l 4 -k 01 -s 10 -t 50 -c 0 -F 0000000 -b b
# -i is the input fasta file.
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -l -k are used in the extraction of primer-binding-sites.
# -s is used in determining the direction and source of each sequence in th input fasta file.
# -t is used in determining the reliability of the taxonomic classification of each sequence.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -c is used to determine whether the mode of multiple threads is on.
# -F is the return value of the checkbox: seven digits; 1 means checked; 0 means unchecked
# -b is the domain you want to evaluate primers. b is Bacteria, e is Eukaryota, a is Archaea.
#
#
##Output file: ssu_merge_XXXX.fasta, ssu_merge_XXXX.silva.wang.taxonomy, lsu_merge_XXXX.fasta, lsu_merge_XXXX.silva.wang.taxonomy and directory of result_$$_input
#
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

getopt ('iprlkstcbF'); 

# -l -k are used in the extraction of primer-binding-sites.
# -s is used in determining the direction and source of each sequence in th input fasta file.
# -t is used in determining the reliability of the taxonomic classification of each sequence.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -c is used to determine whether the mode of multiple threads is on.
# -F is the return value of the checkbox: 
# -F seven digits; 1 means checked; 0 means unchecked 

my $input_fasta=$opt_i;
my $input_primer=$opt_p;
my $input_ref_seq=$opt_r;
my $length_of_elongation=$opt_l;   
my $min_percent_of_pbsite_length=$opt_k;  
my $cutoff_of_search_score=$opt_s;  
my $cutoff_of_taxo_bootstrap=$opt_t;
my $if_multiple_threads=$opt_c;
my $checkbox_input=$opt_F;
my $domain_selected=$opt_b;


# Declare the template and taxonomy files used in Mothur jobs
# These names should change when analyzing functional gene sequences.

my $ssu_template='v119silva.SSU.fasta';    ##If you want to change database, change these four file names.
my $lsu_template='v119silva.LSU.fasta';   ##If you do not want to do LSU alignment, do not add this file into the working directory.
my $taxo_ref='v119fangsilva.SSU.silva.tax'; ##".silva.tax" is canned format, do not change it.
my $greengene='v119silva.SSU_PICK.fasta';   ##It is a custom SILVA database to evaluate 27F and 1492R, greengene is also available.
my $taxo_lsu_ref='v119silva.LSU.silva.tax';  ##2015.09.01 for LSU tax.

# Create a log file to record information and for debugging
# Log file name is like:  3423_test.fasta.log

$log_file="$$"."$input_fasta";
open LOG, ">$log_file.log";

# Check the -F value

$checkbox_input=~m/[01]{7}/ or die "option F error.\n";
print LOG "Checkbox:$checkbox_input\n";


#---Step 2: According to the information from the checkbox, append primers to the input primer file-------------------------#

# The default seven primers that are shown on the webpage

my @seven_primers = qw/8F 338F 519F 907F 1100R 1390R 1492R/;
my @seven_primer_seqs = qw/
AGAGTTTGATCATGGCTCAG 
ACTCCTACGGGAGGCAGC 
CAGCAGCCGCGGTAATAC 
AAACTCAAATGAATTGACGG 
CAACGAGCGCAACCCT 
TTGTACACACCGCCCGT 
AAGTCGTAACAAGGTA/;

# If the checkbox is checked, that is, if the value equals 1, then append the primer's name and sequence to the input primer file.

open PRIMER, ">>$input_primer" || die "$!";
print PRIMER "\n";
for(my $i=0;$i<7;$i++)
{
	if(substr($checkbox_input,$i,1) ==1)
	{
          print PRIMER "$seven_primers[$i]:$seven_primer_seqs[$i]\n";
	  print LOG "Primer checked and added:$seven_primers[$i]\n";
	}
}
close PRIMER;


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
(my $fasta_name="$$"."_$in_path[$#in_path]")=~s/\.fasta$//;

# Create a result directory to save the files.

my $new_dir_name= 'result_'. $fasta_name;
mkdir ("$new_dir_name", 0755);

# Create four fasta files for each alignment

$ssu_forward=$fasta_name . '.ssu.fasta';
$lsu_forward=$fasta_name . '.lsu.fasta';
$ssu_reverse=$fasta_name . '.ssu.rc.fasta';
$lsu_reverse=$fasta_name . '.lsu.rc.fasta';

my @candidate=($ssu_forward, $ssu_reverse, $lsu_forward, $lsu_reverse);
my @template=($ssu_template, $lsu_template);

system("cp $input_fasta $ssu_forward");
system("./mothur \"#reverse.seqs(fasta=$ssu_forward)\"");
system("cp $ssu_forward $lsu_forward");
system("cp $ssu_reverse $lsu_reverse");

# Whether multiple threads:

if($if_multiple_threads)
{
	my @threads;
	for (my $i=0;$i<=3;$i++)
	{
          $threads[$i]=Thread->new(\&mothur_align,"$candidate[$i]","$template[$i/2]");
	}
	foreach my $thread(@threads){$thread->join();}
}
else
{
	for (my $i=0;$i<=3;$i++)
	{
          system("./mothur \"#align.seqs(candidate=$candidate[$i],template=$template[$i/2])\"");
        }
}

#----Step 5: According to the searchscores in the report files, pick seqs-----------------------------------------------------------#

# The array @report contains four elements.
# Each element is a reference to another array that saves all the lines in an .align.report file

my @report;
my @candi_name;

for(my $i=0;$i<=3;$i++)
{
        ($candi_name[$i]=$candidate[$i])=~s/\.fasta$//;
	open REPORT,"$candi_name[$i].align.report" || die "$!";
	my @sequence_records=<REPORT>;
	push(@report,[@sequence_records]);   #create anonymous array references and push them into @report
        close REPORT;
}

# Find sequences with searchscores bigger than the cutoff and that belong to SSU rDNA
# Open two files to write down these sequences' names
# One for SSU-forward alignment, one for SSU-reverse alignment 

open(SSUF,">ssu_forward_$$") or die $!;
open(SSUR,">ssu_reverse_$$") or die $!;
open(LSUF,">lsu_forward_$$") or die $!;##2015.09.01print lsu
open(LSUR,">lsu_reverse_$$") or die $!;##2015.09.01print lsu

for( my $i=1;$i<=$#{$report[0]};$i++)
{
	my @item_ssu_f=split(/\s+/,$report[0]->[$i]);  
	my $seq_name=$item_ssu_f[0];
	my $search_score=$item_ssu_f[5];

        # $order means which file the max of 4 searchscores comes from.

	my $order=0;
	
	for(my $j=1;$j<=3;$j++)
	{
		my @item_other=split(/\s/,$report[$j]->[$i]);
		if($item_other[5] > $search_score)
		{
			$search_score=$item_other[5];
			$order=$j;
		}
	}
	if($search_score>$cutoff_of_search_score)
	{
		if($order == 0) {print SSUF "$seq_name\n";}
		elsif($order == 1) {print SSUR "$seq_name\n";}
		elsif($order == 2) {print LSUF "$seq_name\n";}##2015.09.01print lsu
		elsif($order == 3) {print LSUR "$seq_name\n";}##2015.09.01print lsu
	}
}
close SSUF;
close SSUR;
close LSUF;##2015.09.01print lsu
close LSUR;##2015.09.01print lsu


# The merge file below is used in the next step.

system("./mothur \"#get.seqs(accnos=ssu_forward_$$,fasta=$candi_name[0].align)\"");
system("./mothur \"#get.seqs(accnos=ssu_reverse_$$,fasta=$candi_name[1].align)\"");
my $ssu_merge='ssu_merge_'."$$".'.align';
if(-z "ssu_forward_$$"){system ("cp ssu_forward_$$ $candi_name[0].pick.align");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.
if(-z "ssu_reverse_$$"){system ("cp ssu_reverse_$$ $candi_name[1].pick.align");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.


system("./mothur \"#merge.files(input=$candi_name[0].pick.align-$candi_name[1].pick.align,output=$ssu_merge)\"");


# This merge file below is for later use in Step 6.


system("./mothur \"#get.seqs(accnos=ssu_forward_$$,fasta=$candidate[0])\"");
system("./mothur \"#get.seqs(accnos=ssu_reverse_$$,fasta=$candidate[1])\"");
my $ssu_merge_fasta='ssu_merge_'."$$".'.fasta';
if(-z "ssu_forward_$$"){system ("cp ssu_forward_$$ $candi_name[0].pick.fasta");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.
if(-z "ssu_reverse_$$"){system ("cp ssu_reverse_$$ $candi_name[1].pick.fasta");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.

system("./mothur \"#merge.files(input=$candi_name[0].pick.fasta-$candi_name[1].pick.fasta,output=$ssu_merge_fasta)\"");


system("./mothur \"#classify.seqs(fasta=$ssu_merge_fasta,reference=$ssu_template,taxonomy=$taxo_ref,cutoff=$cutoff_of_taxo_bootstrap)\"");

##merge--classify--lsu
system("./mothur \"#get.seqs(accnos=lsu_forward_$$,fasta=$candidate[2])\"");
system("./mothur \"#get.seqs(accnos=lsu_reverse_$$,fasta=$candidate[3])\"");
my $lsu_merge_fasta='lsu_merge_'."$$".'.fasta';
if(-z "lsu_forward_$$"){system ("cp lsu_forward_$$ $candi_name[2].pick.fasta");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.
if(-z "lsu_reverse_$$"){system ("cp lsu_reverse_$$ $candi_name[3].pick.fasta");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.

system("./mothur \"#merge.files(input=$candi_name[2].pick.fasta-$candi_name[3].pick.fasta,output=$lsu_merge_fasta)\"");


system("./mothur \"#classify.seqs(fasta=$lsu_merge_fasta,reference=$lsu_template,taxonomy=$taxo_lsu_ref,cutoff=$cutoff_of_taxo_bootstrap)\"");

system ("rm *.align ssu_*_$$ *.logfile *.flip.accnos *$$*.ssu.fasta *$$*.lsu.fasta *$$*.rc.fasta"); ##If you want to get more files, comment out the line.



#---Step 5: Select bacteria sequences according to the bootstrap value----------------------------------------------#

# Create a hash to save the taxonomy information of each sequence.
# The domain_name can be Bacteria or Archaea, according to the user's purpose.

my %taxo_info;  

(my $ssu_merge_name=$ssu_merge)=~s/\.align$//;
my $domain_name='Bacteria' if $domain_selected eq 'b';
$domain_name='Archaea' if $domain_selected eq 'a';
$domain_name='Eukaryota' if $domain_selected eq 'e';

open TAXO,"$ssu_merge_name.silva.wang.taxonomy" || die "$!"; #To cover the bug of mothur v1.33, ".silva.wang.taxonomy" is canned format.
open DOMAIN, ">pick_${domain_name}_$$" || die "$!";
while(my $line=<TAXO>)
{
	my @item=split(/\s+/,$line);
	my @taxo=split(/;/,$item[1]);
	if($taxo[0]=~m/^$domain_name\(([0-9]{1,4})\)/){
		if($1 >= $cutoff_of_taxo_bootstrap){
			print DOMAIN "$item[0]\n";
			$taxo_info{"$item[0]"}="$taxo[0];$taxo[1];$taxo[2];$taxo[3];$taxo[4];$taxo[5]"; ##$taxo_info{"$item[0]"}="$taxo[0];$taxo[1];$taxo[2];$taxo[3];$taxo[4];$taxo[5]"; is also OK.
		}
	}
}
close TAXO;
close DOMAIN;


system("./mothur \"#get.seqs(accnos=pick_${domain_name}_$$,fasta=$ssu_merge_fasta)\"");
unlink("pick_${domain_name}_$$");


#---Step 6: Add reference sequence to the picked bacteria/archaea fasta---------------------------------------------------------#

# Add the reference sequence into the pick.fasta file. The ref seq is used in finding the primer binding sites.

my $ref_seq_name;
my $domain_picked='ssu_merge_'."$$".'.pick';

open PICK,">>$domain_picked.fasta" || die "$!";
open REF,"$input_ref_seq" || die "$!";
while(my $line=<REF>)
{

	if ($line =~ />/)
	{
	    chomp ($line);
	    $line=~ s/\s+/_/g;
	    $line=~/>(.*)/;
	    $ref_seq_name=$1;
	    print PICK "$line\n";  
	}
	else 
	{	
	    print PICK $line;  
	}
}
	print LOG "Ref_Name:$ref_seq_name";	 
close REF;
close PICK;

# Align again with the greengene template file, 
# because silva template sequences do not contain the primer-binding-sites of 8F and 1492R.

system("./mothur \"#align.seqs(candidate=$domain_picked.fasta,template=$greengene)\"");
system("./mothur \"#filter.seqs(fasta=$domain_picked.align)\"");



#---Step 7: Align the primer with the reference sequence via Mothur-------------------------------------------------------------# 

# Candidate: the reference sequence and the primers
# Template: the ref seq
# This alignment is run to find the primer-binding-sites on the reference sequence
# The ref seq here comes from the filter.fasta file

my $candidate_primer_to_ref = 'candidate_' . $input_ref_seq;
my $template_primer_to_ref='reference_' . "$input_fasta";
open FILT, "$domain_picked.filter.fasta" || die "$!";
open CANDI, ">$candidate_primer_to_ref" || die "$!"; 
open PRIMER, "complete_$input_primer" || die "$!";
open REF, ">$template_primer_to_ref" || die "$!";

print CANDI ">Reference\n";
print REF ">Reference\n";

my $mark_if_ref=0;

while (<FILT>)
{
  if ($mark_if_ref==1)
  {
  	print CANDI $_;
  	print REF $_;  	
  	$mark_if_ref=0;
  }
  if (/>($ref_seq_name)/)
  {
        $mark_if_ref=1;
  }	      
}

while (<PRIMER>)
{
  chomp;
  s/\s+//g;
  my @name_and_seq= split /:/,$_;
	
  if ($name_and_seq[0]=~/R/i)
  {
    print CANDI ">$name_and_seq[0]\n";
    my $rev_seq = get_rev_complement ($name_and_seq[1]);
    $primer_sequence{$name_and_seq[0]}=$rev_seq;
    print CANDI "$rev_seq\n";
  }
  else
  {
    print CANDI ">$name_and_seq[0]\n$name_and_seq[1]\n";
    $primer_sequence{$name_and_seq[0]}=$name_and_seq[1];
  }
}

close FILT;
close PRIMER;
close CANDI;
close REF;

system("mothur \"#align.seqs(candidate=$candidate_primer_to_ref,template=$template_primer_to_ref)\"");


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

open FILT, "$domain_picked.filter.fasta" || die "$!";
open SUPERV, ">supervision.txt";
open SV2, ">superv2.txt";

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
    
          if ($seq_extracted =~ s/[^\.\-ATCGU]//g)
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

close SUPERV;
close FILT;
close SV2;




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

        # Create a hash to save all the primer_binding_sites of the primer.

        my %primer_binding_sites;
        open CANDI, "$pname:tmp_pbsite_candidate.fa" || die "$!";
        my @lines=<CANDI>;
        for (my $i=0; $i<$#lines; $i+=2)
	{
	    chomp ($lines[$i]);
	    chomp ($lines[$i+1]);
	    (my $pbsite_name=$lines[$i])=~s/>//;
            $primer_binding_sites{$pbsite_name}=$lines[$i+1];
	}
	close CANDI;
	
	system("mothur \"#align.seqs(candidate=$pname:tmp_pbsite_candidate.fa,template=$pname:tmp_primer_template)\"");
	system("mothur \"#filter.seqs(fasta=$pname:tmp_pbsite_candidate.align)\"");

#---Step 2: Compare the primer with each primer-binding-site in aligned format and output matching type-----------------------------------#			

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
				$match_type .= '=';
			}
			else                      
			{
			  if (($primer_components[$j] eq "-") && ($single_character =~ /[ATCGU]/))
			  {   
				 my $insertion=lc($single_character);
				 $match_type .= "$insertion";
                                 $mismatch ++;
			  }elsif($single_character eq "-")
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
        system("rm *tmp_*.* *tmp_*");

}

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
  
        for (my $i=2; $i<=$number_of_lines; $i++)
	{
# mark which non-degenerate primer

	    my $primer_number=0;
	    my $match_score=0;
	    my @lines;  	
  	
      	    for (my $j=0;$j<=$degenerate_primer{$pname};$j++)
	    {
	        my $match_type_file=$pname . $j . '_match_type';
                open TMPTYPE, "$new_dir_name/${match_type_file}_tmp";	
                while (<TMPTYPE>){
	    	      next unless ($. == $i);
	  	      $lines[$j] = $_ ;
	  	      /(-?[0-9]+\.[0-9]+)$/m;##/\t([0-9]+.[0-9]+)\n$/ changed as /(-?[0-9]+\.[0-9]+)$/m
                      if ($match_score < $1)
		      {
	  		    $match_score = $1;
	  		    $primer_number = $j;
	      	      }  	  	  		
	        }
                close TMPTYPE;  	
            }    	
            print TYPE $lines[$primer_number];	
         }
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
   
   for (my $i=1;$i<=6;$i++) ##for (my $i=1;$i<=6;$i++) ##if you change line 426 to get deeper taxonomy, you have to change this line.
   {      
      print STAT "Taxonomy\tmatch\tunmatch\tmatch_percent\n";
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

system ("rm *.pick.fasta  *.report *.logfile .*.swo .*.swp *8mer* *train* *sum* candi*"); ##If you want to get more files, comment out the line.

close LOG;


