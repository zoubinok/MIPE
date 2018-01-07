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
##Usage example: perl MIPE_check_import_format.pl -i input.fasta -r Ref_16S.fasta -p primer1.txt -n 88888
# -i is the input fasta file.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -p is the file of primers. The primer should be wrote like this:"name:NNNNNNNN", the name of forward primer name cannot contain R or r,
#    and reverse primer name should contain R or r. 
# -n is is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
#
##Output file: check88888_input.fasta checkRef_16S.fasta checkinput88888_input.txt checkprimer1.txt
# 88888 is an example. Please change it.
# These files can be used as input of downstream analysis. but I suggest you to adjust your own input so that they can be used directly.
####################################################################################################################################


# The modules used: 
# Getopt::std   parameters input
# Thread   multiple-thread computation in the alignment




use Getopt::Std;
use Thread;


getopt ('iprn'); 

# -l -k are used in the extraction of primer-binding-sites. -k is the min percent of pbsite length to filt sequences not aligned in this part (now - is also regarded as aligned, but . is not).
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.
# -c is used to determine whether the mode of multiple threads is on.
# -n is the sign of this sample.it should be a positive integer. 10000-99999 is recommanded.
# -r is the reference fasta file which contains only one reference sequence. Default is the 16s rDNA sequence of E.coli.

my $input_fasta=$opt_i;
my $input_primer=$opt_p;
my $sign=$opt_n;
my $input_ref_seq=$opt_r;



my @in_path=split(/\//,$input_fasta);
(my $fasta_name="$sign"."_$in_path[$#in_path]")=~s/\.fasta$//;
my $seq_name;
open CHK, ">checkinput$fasta_name.txt";


my $ref_seq_name;
open REF,"$input_ref_seq" || die "$!";
open OUTREF,">check$input_ref_seq";

while(my $line=<REF>)
{
	    chomp ($line);
        $line =~ s/[\r\n]//g;
	next unless $line =~/[a-zA-Z0-9]/;
		

	if ($line =~ /^>(.+)/)
	{
	     $ref_seq_name=$1;
		if($ref_seq_name=~ s/[^a-zA-Z0-9]//g){print CHK "The name of $ref_seq_name has the special character(s), I delete it. I SUGGEST YOU ONLY USE LETTERS AND NUMBERS.\n";}
	    print OUTREF ">$ref_seq_name\n";
	
	}else{
		if($line=~ s/[^ATCGatcg]//g){print CHK "The sequence of $ref_seq_name has the illegal character(s), I delete it. I SUGGEST YOU CHANGE THE REFERENCE.\n";}
	    print OUTREF "$line\n";	
	}
}
close REF;
close OUTREF;

open PRIMER, "$input_primer" || die "$!";
open OUTPRIMER, ">check$input_primer" || die "$!";
print CHK "\n\n\nI SUGGEST YOU DO NOT USE LOWERCASES IN THIS FILE.\n";
while (<PRIMER>)
{
    chomp;
	s/[\r\n]//g;
	next unless /[a-zA-Z]/;

    my @name_and_seq = split /:/,$_;
    my $pname=$name_and_seq[0];
    	if(($pname =~ /[Ff]/) && ($pname =~ /[Rr]/)){
		print CHK "$pname has F&R. As F&R is used to recognize forward primers or reverse primers, they cannot be in a single name. PLEASE CHANGE IT AND I SUGGEST YOU DO NOT USE LOWERCASES IN THIS FILE.\n";
		print OUTPRIMER "ERR plases check checkinput$fasta_name.txt\n";
		}elsif($pname =~ /[0-9]$/){
		print CHK "The last character of $pname has a number. You'd better use F or R as the end. I SUGGEST YOU DO NOT END WITH A NUMBER. PLEASE CHANGE IT AND I SUGGEST YOU DO NOT USE LOWERCASES IN THIS FILE.\n";
		print OUTPRIMER "ERR plases check checkinput$fasta_name.txt\n";
		}elsif($name_and_seq[1] =~ s/[^ATCGRYMKSWHBVDN]//gi){
		print CHK "$pname has illegal characters. PLEASE CHANGE IT AND I SUGGEST YOU DO NOT USE LOWERCASES IN THIS FILE.\n";
		print OUTPRIMER "ERR plases check checkinput$fasta_name.txt\n";
		}
		else{
		$_=uc($_);
		s/\s+//g;

		print OUTPRIMER "$_\n";
		}

}

close PRIMER;
close OUTPRIMER;

open IN,"$input_fasta" || die "$!";
open OUTIN,">check$fasta_name.fasta";
while(my $line=<IN>)
{
	chomp ($line);
	$line =~ s/[\r\n]//g;
	next unless $line =~/[a-zA-Z0-9]/;

	if ($line =~ /^>/)
	{
	    $seq_name=$line;
		if($line=~ s/\s+/_/g){print CHK "$seq_name has the blank space(s), I replaced it with \"_\". I SUGGEST YOU DO NOT USE THE BLACK SPACE.\n";}
		if($line=~ s/:+/_/g){print CHK "$seq_name has \":\", I replaced it with \"_\". I SUGGEST YOU DO NOT USE \":\".\n";}
		print OUTIN "$line\n";
		}else{
		if($line=~ s/[^\.\-ATCGUatcgu]//){print CHK "$seq_name has illegal characters. I delete it.\n";}
		print OUTIN "$line\n";
		}
}
close IN;
close OUTIN;
close CHK;
print "Please check checkinput$fasta_name.txt to solve problems\n";


