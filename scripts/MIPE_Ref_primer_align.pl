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
##Usage example: perl MIPE_Ref_primer_align.pl -i input.fasta -r Ref_16S.fasta -p primer1.txt  -n 88888
# -i is the input fasta file.
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
#
##Output file: candidate_Ref_16S.align
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

This script is used to align reference fasta (from aligned ssu_merge_xxxxx.pick.filter.fasta) and non-degenerate primers. You need to prepare an <input> multiple fasta file, an aligned ssu_merge_xxxxx.pick.filter.fasta, a formatted <primer_file> and a proper <reference_file>.

=head1 Usage

    $0 -i <input> -p <primer_file> -n <sign> -r <reference_file>
    
=head1 Parameters

    -i  [str]   Input multiple nucleotide fasta file.
    -p  [str]   The file of primers. Primer format:"primer_name:NNNNNNNN". Forward primer name: end with F and no R or r is allowed in name. Reverse primer name: end with R and no F or f is allowed in name.
    -r  [str]   Reference fasta file which contains only one reference sequence whose name cannot appear in <input>. 
    -n  [int]   The sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
=cut

getopt ('iprn'); 

# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.

my $input_fasta=$opt_i;
my $input_primer=$opt_p;
my $input_ref_seq=$opt_r;
my $sign=$opt_n;

die `pod2text $0` if ((!$input_fasta) or (!$input_primer) or (!$sign) or (!$input_ref_seq));

print "the aligned reference and primer file is candidate_Ref_name.align, you can make it in several ways\n";

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
my @primer_names;
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
	   push @primer_names, $name;

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


#----Step 5: According to the searchscores in the report files, pick seqs-----------------------------------------------------------#

#---Step 6: Add reference sequence to the picked bacteria/archaea fasta---------------------------------------------------------#

my $ref_seq_name;
my $domain_picked='ssu_merge_'."$sign".'.pick';

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
	print LOG "Ref_Name:$ref_seq_name";	 
close REF;

#---Step 7: Align the primer with the reference sequence via Mothur-------------------------------------------------------------# 

# Candidate: the reference sequence and the primers
# Template: the ref seq
# This alignment is run to find the primer-binding-sites on the reference sequence
# The ref seq here comes from the filter.fasta file

my $candidate_primer_to_ref = 'candidate_' . $input_ref_seq;
my $template_primer_to_ref='reference_' . "$input_fasta"; #unly put the aligned reference into it.
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

my @refalign=split(/\.fa/,$input_ref_seq);
my $candidate_align = 'candidate_' . $refalign[0] . '.align';
open CANDIALI, "$candidate_align" || die "$!"; 
my %alignrefprimer;
while(my $line=<CANDIALI>)
{

	if ($line =~ />/)
	{
	    chomp ($line);
	    $line=~/>(.*)/;
	    $ref_seq_name=$1;
	}
	else 
	{	
        my @primeralignstart=split(/[ATCGatcg\-]/,$line);
		$alignrefprimer{$ref_seq_name}=length($primeralignstart[0]);
		$ref_seq_name='';
	}
}
close CANDIALI;

open PRIMER, "$input_primer" || die "$!";

while (<PRIMER>)
{
    chomp;
    s/\s+//g;
    $_=uc($_);
    next unless /[a-zA-Z]/;
    my @name_and_seq = split /:/,$_;
    my @all_primers = get_primer_suite ($_);
	my $name = $name_and_seq[0] . '0'; 
    for (my $i=0;$i<=$#all_primers;$i++) 
    {
       my $alignrefprimer_old=$alignrefprimer{$name};
	   $name = $name_and_seq[0] . $i; 
       unless ($alignrefprimer{$name}==$alignrefprimer_old){print "ERR: $name_and_seq[0] does not align well. \n";}
    }
}

close PRIMER;





print "If it did not aligned well, please correct the .algn file. But do not change >Reference as it is extracted from ssu_merge_$sign.pick.filter.fasta.\n";

system ("rm *.report *.logfile .*.swo .*.swp *8mer* *train* *sum* "); ##If you want to get more files, comment out the line.

close LOG;


