#!/usr/bin/perl -w

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
##Dependence: Linux, mothur v1.33, ncbi blast tools. 
# Copy mothur into the working directory, then input: chmod 755 mothur
# NCBI BLAST tools should be in environment variable.
#
##Put "silva.SSU75seeds.fasta" and input file into the working directory.
#
##Usage: perl MIPE_pre_program.pl input.fasta 0 10
# Here input.fasta is the inputfile, and should be end in .fasta.
# 0 is Percentage of alignment, the minimun is 0.
# 10 is the cutoff of Evalue in blast.
#
#
##Output file: input_picked.fasta. It is used as input of MIPE_main_program.pl
#
####################################################################################################################################







my @fastaname=(split /\.f/, $ARGV[0]);
system("formatdb -i $ARGV[0] -p F -o T");
system("blastall -p blastn -d $ARGV[0] -i silva.SSU75seeds.fasta -o $ARGV[0].out -v 65535 -b 65535 -F F");
die "Usage:	$0	\$BlastResultFile	\$DescriptionOutputFile	\$WantAlignPercentage \$eValue[?]\n" if(@ARGV != 3);
use strict;
my $WantPercentage = $ARGV[1];                          
my $wantEvalue = $ARGV[2] || 1e-15;
my ($MasterGene,$LastMasterGene,$MasterGeneLength,$LastDescription,$Description,$LastHitLength,$HitLength);
my ($EvalueFlag,$PrintFlag,$FirstHsp,$Score,$EValue,$Identity);
my (%QueryStoreInfo,%SbjctStoreInfo,%TotalAlignHash,@TempArray);
my ($QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd)=("","",1000000000,0,1000000000,0);
my $strandflag=0;
my $count;
open(BLAST,"$ARGV[0].out") or die "Can't open file $ARGV[0]: $!";     # input
open(OUT,">$ARGV[0].out1") or die "Can't create file $ARGV[1]: $!";    # output
open(ACCONS,">$fastaname[0].accnos");
print OUT "Query Name\tQuery Length\tSbjct Length\tQuery Alignment\tSbjct Alignment\tAnnotation\tScore\tE Value\tIdentity\tIdentity_Rate\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tQueryStrand\tSubjectStrand\n";
Abstract();
close(ACCONS);
close(BLAST);
close(OUT);
print "$ARGV[0]\t$count\n";	
#-------------------------------------------------------------------------------------------------------------------
sub Abstract{
		while(<BLAST>){
		chomp;
		if(/Query= (.+)/){
			my $QueryAlignLength = keys %QueryStoreInfo;        
			my $SbjctAlignLength = keys %SbjctStoreInfo;
			my $TotalAlignLength = keys %TotalAlignHash;
			if($QueryAlignLength > 0){
				#$PrintFlag = PrintAnno($QueryAlignLength,$SbjctAlignLength,$TotalAlignLength,$MasterGeneLength,$LastHitLength,$PrintFlag);
				$PrintFlag = PrintAnno($QueryAlignLength,$SbjctAlignLength,$TotalAlignLength,$MasterGeneLength,$LastHitLength,$PrintFlag,$QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd);
			}		
			%QueryStoreInfo = ();
			%SbjctStoreInfo = ();
			%TotalAlignHash = ();	
			($QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd)=("","",100000000,0,100000000,0);
			$strandflag=0;			
			$MasterGene = $1;
			$PrintFlag = 1;
			my $TempLine;
			while($TempLine = <BLAST>){
				chomp($TempLine);
				if($TempLine =~ /\((\S+)\s+letters\)/){
					$MasterGeneLength = $1;              
					$MasterGeneLength =~ s/,//g if($MasterGeneLength =~ /,/);
					last;
				}	
				$TempLine =~ s/^\s+/ /;
				$MasterGene .=" $TempLine";
			}
			$MasterGene =~s/\s+/ /g;			
		}elsif(/^>(.+)/){
			$Description = $1;
			$FirstHsp = 1;
			my $TempLine;
			while($TempLine = <BLAST>){
				chomp($TempLine);
				if($TempLine =~ /Length = (\S+)/){
					$HitLength = $1;
					($HitLength =~ s/,//g) if($HitLength =~ /,/);
					last;
				}	
				$TempLine =~ s/^\s+/ /;
				$Description .= $TempLine;
			}	
			
			my $QueryAlignLength = keys %QueryStoreInfo;       
			my $SbjctAlignLength = keys %SbjctStoreInfo;
			my $TotalAlignLength = keys %TotalAlignHash;
			if($QueryAlignLength > 0){
				#$PrintFlag = PrintAnno($QueryAlignLength,$SbjctAlignLength,$TotalAlignLength,$MasterGeneLength,$LastHitLength,$PrintFlag);
				$PrintFlag = PrintAnno($QueryAlignLength,$SbjctAlignLength,$TotalAlignLength,$MasterGeneLength,$LastHitLength,$PrintFlag,$QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd);
			}		
			$LastDescription = $Description;
			$LastHitLength = $HitLength;
			%QueryStoreInfo = ();
			%SbjctStoreInfo = ();
			%TotalAlignHash = ();	
			($QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd)=("","",100000000,0,100000000,0);
			$strandflag=0;		
		}
			elsif(/^ Strand = (\S+)\s+\S+\s+(\S+)/){
			#print $_;
			if($QueryStrand && $SubjectStrand)
			{
				if($QueryStrand ne $1 ||$SubjectStrand ne $2)
				{$strandflag=0;}
			}
			else
			{
				$QueryStrand = $1;
				$SubjectStrand = $2;
				$strandflag=1;
			}
		}
		# Score = 77.8 bits (190), Expect = 1e-16,   Method: Composition-based stats.
		elsif(/Score\s*=\s*(\S+) bits.+, Expect\S* = (\S+)/ && $PrintFlag == 1){                #control e_value
			$EvalueFlag = 0;
			my $TempScore = $1;
			my $e_value = $2;
			$e_value =~s/\W+$//;
			$e_value = '1'.$e_value if($e_value =~ /^e/);
			$EvalueFlag = 1 if($e_value <= $wantEvalue);
				
			if($FirstHsp == 1){
				$FirstHsp = 0;
				$Score = $TempScore;
				$EValue = $e_value;
			}	
		}	
		elsif(/Query:/&& $EvalueFlag == 1 && $PrintFlag == 1){
			my $Query = $_;
			%QueryStoreInfo = QueryStat($Query,%QueryStoreInfo);
			#Query: 1    aacgaacgctggcggcaggcttaacacatgcaagtcgagcgccccgcaaggggagcggca 60
			if($strandflag==1 && $Query=~/^Query:\s+(\S+)\s+\S+\s+(\S+)/)	
			{
				if($QueryStrand eq "Plus")
				{
					$QueryStart=$1 if($QueryStart>$1);
					$QueryEnd=$2 if($QueryEnd<$2);
					#print "$QueryStart\t$QueryEnd\n";
				}
				if($QueryStrand eq "Minus")
				{
					$QueryStart=$2 if($QueryStart>$2);
					$QueryEnd=$1 if($QueryEnd<$1);				
				}
			}	
		}
		elsif(/Sbjct:/&& $EvalueFlag == 1 && $PrintFlag == 1){
			my $Sbjct = $_;
			%SbjctStoreInfo = SbjctStat($Sbjct,%SbjctStoreInfo);
			if($strandflag==1 && $Sbjct=~/^Sbjct:\s+(\S+)\s+\S+\s+(\S+)/)	
			{
				if($SubjectStrand eq "Plus")
				{	
					$SubjectStart=$1 if($SubjectStart>$1);
					$subjectEnd=$2 if($subjectEnd<$2);
				}
				if($SubjectStrand eq "Minus")
				{
					$SubjectStart=$2 if($SubjectStart>$2);
					$subjectEnd=$1 if($subjectEnd<$1);				
				}									
			}		
		}
		elsif(/^  Database:/){	
			
			my $QueryAlignLength = keys %QueryStoreInfo;        
			my $SbjctAlignLength = keys %SbjctStoreInfo;
			my $TotalAlignLength = keys %TotalAlignHash;
			#print "$MasterGene\t$QueryAlignLength\n";
			if($QueryAlignLength > 0){
				#$PrintFlag = PrintAnno($QueryAlignLength,$SbjctAlignLength,$TotalAlignLength,$MasterGeneLength,$LastHitLength,$PrintFlag);
				$PrintFlag = PrintAnno($QueryAlignLength,$SbjctAlignLength,$TotalAlignLength,$MasterGeneLength,$LastHitLength,$PrintFlag,$QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd);
			}	
			%QueryStoreInfo = ();
			%SbjctStoreInfo = ();
			%TotalAlignHash = ();	
			($QueryStrand,$SubjectStrand,$QueryStart,$QueryEnd,$SubjectStart,$subjectEnd)=("","",100000000,0,100000000,0);
			$strandflag=0;
		}	
	}
}
#--------------------------------------------------------------------------------------
sub QueryStat{
	my ($QueryLine,%TempQueryStoreInfo) = @_;
	
	if($QueryLine =~ /Query: (\d+)(\s*)(.+)\s(\d+)$/){
		my $BlankLength = length($1) + length($2) + 7;
		my $QueryBegin = $1;
		my $QueryEnd = $4;
		my $QueryContent = $3;
		my %QueryStore;
		
		%QueryStore = StorePlus($QueryContent,$QueryBegin);
		
		for(my $i = $QueryBegin; $i <= $QueryEnd; $i ++){
			$TotalAlignHash{$i} = 1;		
		}	
		
		my $MiddleLine = <BLAST>;
		chomp($MiddleLine);
		$MiddleLine =~ s/^\s{$BlankLength}//;
		@TempArray = split(//,$MiddleLine);
		
		my $StartPos = $QueryBegin;
		my $StatPos = $StartPos;
		foreach my $ele (@TempArray){
			($StartPos ++,next) if(exists $QueryStore{$StartPos});
			$TempQueryStoreInfo{$StatPos} = 1 if($ele eq '|');
			$StatPos ++;
			$StartPos ++;
		}
	}else{
		print "Wrong! Can't Match Query Line.\n";
	}
	return %TempQueryStoreInfo;		
}
#--------------------------------------------------------------------------------------
sub SbjctStat{
	my($SbjctLine,%TempSbjctStoreInfo) = @_;
	my $Strand;
	if($SbjctLine =~ /Sbjct: (\d+)(\s*)(.+)\s(\d+)$/){
		my $SbjctBegin = $1;
		my $SbjctEnd = $4;
		my $SbjctContent = $3;
		my %SbjctStore;
		$Strand = ($SbjctBegin < $SbjctEnd ? 1 : 0);	
		if($Strand == 1){
			%SbjctStore = StorePlus($SbjctContent,$SbjctBegin);
			
			my $StartPos = $SbjctBegin;
			my $StatPos = $StartPos;
			foreach my $ele (@TempArray){
				($StartPos ++,next) if(exists $SbjctStore{$StartPos});
				$TempSbjctStoreInfo{$StatPos} = 1 if($ele eq '|');
				$StatPos ++;
				$StartPos ++;
			}
		}else{
			%SbjctStore = StoreMinus($SbjctContent,$SbjctBegin);
			
			my $StartPos = $SbjctBegin;
			my $StatPos = $StartPos;
			foreach my $ele (@TempArray){
				($StartPos --,next) if(exists $SbjctStore{$StartPos});
				$TempSbjctStoreInfo{$StatPos} = 1 if($ele eq '|');
				$StatPos --;
				$StartPos --;
			}
		}		
	}else{
		print "Wrong! Can't Match Sbjct Line.\n";
	}
	return %TempSbjctStoreInfo;	
}	
#--------------------------------------------------------------------------------------
sub StorePlus{
	my ($content,$Begin) = @_;
	my %Store;
	while($content =~ /-+/g){
		my $Position = pos($content);
		my $Length = length($&);
		my $PositionStart = $Position - $Length + 1;
		for(my $i = $PositionStart; $i <= $Position; $i ++){
			$Store{$i + $Begin - 1} = 1;
		}	
	}
	return(%Store);
}
#--------------------------------------------------------------------------------------
sub StoreMinus{
	my ($content,$Begin) = @_;
	my %Store;
	while($content =~ /-+/g){
		my $Position = pos($content);
		my $Length = length($&);
		my $PositionStart = $Position - $Length + 1;
		for(my $i = $PositionStart; $i <= $Position; $i ++){
			$Store{$Begin - $i + 1} = 1;
		}	
	}
	return(%Store);
}
#--------------------------------------------------------------------------------------
sub PrintAnno{
	my($qL,$sL,$tL,$mL,$hL,$flag,$QueryStrandtmp,$SubjectStrandtmp,$QueryStarttmp,$QueryEndtmp,$SubjectStarttmp,$SubjectEndtmp) = @_;
	my $QueryAlignPercentage = sprintf("%.2f",$qL / $mL * 100);
	my $SbjctAlignPercentage = sprintf("%.2f",$sL / $hL * 100);			
	$Identity = sprintf("%.2f",$qL / $tL * 100);
	if(($QueryAlignPercentage >= $WantPercentage || $SbjctAlignPercentage >= $WantPercentage) && $flag == 1){
		$count++;
		print OUT "$MasterGene\t$mL\t$hL\t$QueryAlignPercentage\t$SbjctAlignPercentage\t$LastDescription\t$Score\t$EValue\t$qL/$tL\t$Identity\t$QueryStarttmp\t$QueryEndtmp\t$SubjectStarttmp\t$SubjectEndtmp\t$QueryStrandtmp\t$SubjectStrandtmp\n";
		print ACCONS "$LastDescription\n";   ##change "$LastDescription\t" to "$LastDescription\n"#
	}
	return $flag;
}	
system("./mothur \"#get.seqs(accnos=$fastaname[0].accnos,fasta=$ARGV[0])\"");
system("mv $fastaname[0].pick.fasta  $fastaname[0]_picked.fasta");
