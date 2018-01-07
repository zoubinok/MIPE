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
##Usage example: perl MIPE_primer_format.pl -p primer1.txt -n 88888
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
#
##Output file: complete_primer1.txt
# 88888 is an example. Please change it.
####################################################################################################################################



# The modules used: 
# Getopt::std   parameters input
# Thread   multiple-thread computation in the alignment




use Getopt::Std;
use Thread;

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

getopt ('pn'); 


my $input_primer=$opt_p;
my $sign=$opt_n;


#$log_file="$sign"."$input_fasta";
#open LOG, ">>$log_file.log";

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
    #print LOG "Get the primer suite of the degenerate primers:\n";
    chomp;
    s/\s+//g;
    $_=uc($_);
    next unless /[a-zA-Z]/;
    my @name_and_seq = split /:/,$_;
    my @all_primers = get_primer_suite ($_);
	#print @all_primers;
    $degenerate_primer{$name_and_seq[0]} = $#all_primers; 
    print '$degenerate_primer{$name_and_seq[0]}:';
	print "$#all_primers\n";
	
    for (my $i=0;$i<=$#all_primers;$i++) 
    {
       $name = $name_and_seq[0] . $i; 
       $primer_length{$name} = length ($all_primers[$i]);
       print SUITE "$name:$all_primers[$i]\n";
       #print LOG "$name\t";
	   push @primer_names, $name;

    }
#   print LOG "\n";
}
   print @primer_names;

close SUITE;
close PRIMER;
