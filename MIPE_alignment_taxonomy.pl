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
##Usage example: perl MIPE_MIPE_alignment_taxonomy.pl -i input.fasta -r Ref_16S.fasta -s 10 -t 50 -c 4 -n 88888
# -i is the input fasta file.
# -s is used in determining the direction and source of each sequence in th input fasta file.
# -t is used in determining the reliability of the taxonomic classification of each sequence.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -c is used to determine the number of multiple threads in Mothur.
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
#
##Output file: ssu_merge_88888.fasta, ssu_merge_88888.silva.wang.taxonomy and ssu_merge_88888.filter.fasta
# 88888 is an example. Please change it.
####################################################################################################################################

# The modules used: 
# Getopt::std   parameters input
# Thread   multiple-thread computation in the alignment




use Getopt::Std;
use Thread;


#---Stage 1: Get the primer-binding-sites and taxonomic information of each sequence-----------------------------------#


#---Step 1: Declare my variables and define Getopt parameters---------------------------------------------------------#

getopt ('irsctn'); 

# -s is used in determining the direction and source of each sequence in th input fasta file.
# -t is used in determining the reliability of the taxonomic classification of each sequence.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -c is multiple threadsis used in Mothur.
# -n is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.

my $input_fasta=$opt_i;
my $input_ref_seq=$opt_r;
my $cutoff_of_search_score=$opt_s;  
my $multiple_threads=$opt_c;
my $cutoff_of_taxo_bootstrap=$opt_t;
my $sign=$opt_n;

# Declare the template and taxonomy files used in Mothur jobs
# These names should change when analyzing functional gene sequences.

my $ssu_template='./database/v119silva.SSU.fasta';    ##If you want to change database, change these four file names.
my $lsu_template='./database/v119silva.LSU.fasta';   ##If you do not want to do LSU alignment, do not add this file into the working directory.
my $taxo_ref='./database/v119fangsilva.SSU.silva.tax'; ##".silva.tax" is canned format, do not change it.
my $greengene='./database/v119silva.SSU_PICK.fasta';   ##It is a custom SILVA database to evaluate 27F and 1492R, greengene is also available.

# Create a log file to record information and for debugging
# Log file name is like:  3423_test.fasta.log

$log_file="$sign"."$input_fasta";
open LOG, ">>$log_file.log";


#---Step 2: According to the information from the checkbox, append primers to the input primer file-------------------------#


#---Step 3: Read the input primer file and translate the degenerate primers-----------------------------------------------------------------#


#---Step 4: Align the fasta file with LSU/SSU in both original and reverse directions (multiple threads may be used)----------------------------------#
print "Align the fasta file with LSU/SSU in both original and reverse directions (multiple threads may be used)\n";
# Get the fasta file name and add a random number to make the job unique and easy to recognize:

my @in_path=split(/\//,$input_fasta);
(my $fasta_name="$sign"."_$in_path[$#in_path]")=~s/\.fasta$//;

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


	for (my $i=0;$i<=3;$i++)
	{
          system("./mothur \"#align.seqs(candidate=$candidate[$i],template=$template[$i/2],processors=$multiple_threads)\"");
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
        #print "@report\n";
		close REPORT;
}

# Find sequences with searchscores bigger than the cutoff and that belong to SSU rDNA
# Open two files to write down these sequences' names
# One for SSU-forward alignment, one for SSU-reverse alignment 

open(SSUF,">ssu_forward_$sign") or die $!;
open(SSUR,">ssu_reverse_$sign") or die $!;

for( my $i=1;$i<=$#{$report[0]};$i++)
{
	#print "$report[0]          $#{$report[0]}                    $report[0]->[$i]\n";
	my @item_ssu_f=split(/\s+/,$report[0]->[$i]);  #$repoet[0]的第$i个元素
	#print @item_ssu_f;
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
	}
}
close SSUF;
close SSUR;

# The merge file below is used in the next step.

system("./mothur \"#get.seqs(accnos=ssu_forward_$sign,fasta=$candi_name[0].align)\"");
system("./mothur \"#get.seqs(accnos=ssu_reverse_$sign,fasta=$candi_name[1].align)\"");
my $ssu_merge='ssu_merge_'."$sign".'.align';
if(-z "ssu_forward_$sign"){system ("cp ssu_forward_$sign $candi_name[0].pick.align");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.
if(-z "ssu_reverse_$sign"){system ("cp ssu_reverse_$sign $candi_name[1].pick.align");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.


system("./mothur \"#merge.files(input=$candi_name[0].pick.align-$candi_name[1].pick.align,output=$ssu_merge)\"");


# This merge file below is for later use in Step 6.


system("./mothur \"#get.seqs(accnos=ssu_forward_$sign,fasta=$candidate[0])\"");
system("./mothur \"#get.seqs(accnos=ssu_reverse_$sign,fasta=$candidate[1])\"");
my $ssu_merge_fasta='ssu_merge_'."$sign".'.fasta';
if(-z "ssu_forward_$sign"){system ("cp ssu_forward_$sign $candi_name[0].pick.fasta");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.
if(-z "ssu_reverse_$sign"){system ("cp ssu_reverse_$sign $candi_name[1].pick.fasta");}#To cover the bug of mothur v1.33, add a empty file if mothur reports error.
print "Make ssu_merge file\n";
system("./mothur \"#merge.files(input=$candi_name[0].pick.fasta-$candi_name[1].pick.fasta,output=$ssu_merge_fasta)\"");


system("./mothur \"#classify.seqs(fasta=$ssu_merge_fasta,reference=$ssu_template,taxonomy=$taxo_ref,cutoff=$cutoff_of_taxo_bootstrap,processors=$multiple_threads)\"");

system ("rm *.align ssu_*_$sign *.logfile *.flip.accnos *$sign*.ssu.fasta *$sign*.lsu.fasta *$sign*.rc.fasta"); ##If you want to get more files, comment out the line.



#---Step 5: Select bacteria sequences according to the bootstrap value----------------------------------------------#

# Create a hash to save the taxonomy information of each sequence.
# The domain_name can be Bacteria or Archaea, according to the user's purpose.


#---Step 6: Add reference sequence to the ssu_merge fasta---------------------------------------------------------#
print "Add reference sequence to the ssu_merge fasta\n";
# Add the reference sequence into the pick.fasta file. The ref seq is used in finding the primer binding sites.

my $ref_seq_name;
my $domain_picked='ssu_merge_'."$sign".'.pick';
system("cp $ssu_merge_fasta $domain_picked.fasta");

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

system("./mothur \"#align.seqs(candidate=$domain_picked.fasta,template=$greengene,processors=$multiple_threads)\"");#keep this file so that you can add new reference sequence.
system("./mothur \"#filter.seqs(fasta=$domain_picked.align)\"");
close LOG;

