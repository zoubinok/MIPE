# MIPE
MIPE: Metagenome based community structure explorer and primer evaluation tool

MIPE (MIcrobiota Prior knowledge Explorer) is a metagenome based community structure explorer and primer evaluation tool.With a metagenome data set input, MIPE demonstrate the detailed composition taxonomy information, into class level, and the coverage of the primer of interest.

MIPE contains three perl scripts.

MIPE_pre_program.pl does Step0: Pre-Propcessing.it use BLAST tools to abstract target gene from metagenome data set. Metatranscriptome data can skip this step.

MIPE_main_program.pl does SSU taxonomy by mothur v1.33, find SSU primer binding sites and calculate primer coverage.

MIPE_main_program_metatranscriptome.pl does SSU and LSU taxonomy by mothur v1.33, find SSU primer binding sites and calculate primer coverage.
