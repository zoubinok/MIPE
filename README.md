# MIPE
MIPE: Metagenome based community structure explorer and primer evaluation tool

MIPE (MIcrobiota Prior knowledge Explorer) is a metagenome based community structure explorer and primer evaluation tool.With a metagenome data set input, MIPE demonstrate the detailed composition taxonomy information, into class level, and the coverage of the primer of interest.

MIPE contains two perl scripts.

MIPE_pre_program.pl does Step0: Pre-Propcessing.it use BLAST tools to abstract target gene from metagenome data set.

MIPE_main_program.pl do taxonomy by mothur v1.33, find primer binding sites and calculate primer coverage.
