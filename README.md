# MIPE
MIPE: Metagenome based community structure explorer and primer evaluation tool

MIPE (MIcrobiota Prior knowledge Explorer) is a metagenome based community structure explorer and primer evaluation tool. With a metagenome dataset input, MIPE demonstrates the detailed composition taxonomy information, and the coverage of the primer of interest.

MIPE contains three perl scripts.

MIPE_pre_program.pl does Step0: Pre-Propcessing. It uses BLAST tools to abstract target gene from metagenome data set. Metatranscriptome data can skip this step.

MIPE_main_program.pl does SSU taxonomy by mothur v1.33.3 and evaluates SSU primer binding sites and calculate primer coverage.

MIPE_main_program_metatranscriptome.pl does SSU and LSU taxonomy with mothur v1.33.3 and evaluates SSU primer binding sites and calculate primer coverage.

########################################################
#Platform: Linux
#Installation: MIPE contains somes perl scripts. No need to install MIPE but you need to download MIPE dependences.

#MIPE dependences downloading
#If you use the SILVA reference files you should be aware of their dual-use license https://www.arb-silva.de/silva-license-information

#All MIPE files and dependences must be in the same directory.
#Copy mothur v1.33.3 into the working directory, then input: 
$chmod 755 mothur

#Download Silva.seed_v119.tgz from mothur website. I use mothue modified Silva reference files release 119 to build v119silva.SSU.fasta, v119silva.SSU_PICK.fasta and v119fangsilva.SSU.silva.tax .
$wget http://www.mothur.org/w/images/5/56/Silva.seed_v119.tgz
$tar xvzf Silva.seed_v119.tgz
$mv silva.seed_v119.align v119silva.SSU.fasta
$mv silva.seed_v119.tax v119fangsilva.SSU.silva.tax #I recommend you to change all parentheses for brackets.
$perl SILVA_ENLONG.pl #To build v119silva.SSU_PICK.fasta

#Download SILVA_119_LSURef_tax_silva_full_align_trunc.fasta.gz from SILVA website. I clustered SILVA lSU database (v119) with Usearch (v5.2.32) at a 99% sequence identity level to decrease LSU database and made a taxonomy file (v119silva.LSU.silva.tax) for mothur.
$wget https://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/SILVA_119_LSURef_tax_silva_full_align_trunc.fasta.gz
$gunzip -c SILVA_119_LSURef_tax_silva_full_align_trunc.fasta.gz > SILVA_119_LSURef_tax_silva_full_align_trunc.fasta
$./mothur "#get.seqs(fasta=SILVA_119_LSURef_tax_silva_full_align_trunc.fasta,accnos=v119silva.LSU.silva.tax)" 
$mv SILVA_119_LSURef_tax_silva_full_align_trunc.pick.fasta v119silva.LSU.fasta
