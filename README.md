MIPE: Metagenome based community structure explorer and primer evaluation tool  
Version 2.0.0  
  
MIPE (MIcrobiota Prior knowledge Explorer) is a metagenome based community structure explorer and primer evaluation tool. With a metagenome dataset input, MIPE demonstrates the detailed composition taxonomy information, and the coverage of the primer of interest.  

Based on previous experiences, I divided it into several parts and focused on primer evaluation so that all kinds of genes can use MIPE now. You can also find old version in MIPE_old.  

#######################################################################################################################################  

1.Installation  
Platform: Linux  
Language: Perl  
Installation: Clone or download this repository. MIPE contains some perl scripts. No need to install MIPE but you need to download MIPE dependences.  
  
$chmod 755 mothur #Mothur version 1.33.3  
  
You can also use other versions of mothur but I haven't tested them. Please check the names of mothur output.  
Download HMMER 3.1b2 (http://eddylab.org/software/hmmer/hmmer-3.1b2.tar.gz) and install it into the environment variable.  

MIPE dependences downloading If you use the SILVA reference files you should be aware of their dual-use license https://www.arb-silva.de/silva-license-information  

Download Silva.seed_v119.tgz from mothur website. I use mothur modified Silva reference files release 119 to build v119silva.SSU.fasta, v119silva.SSU_PICK.fasta and v119fangsilva.SSU.silva.tax .  

$cd database  
$wget http://www.mothur.org/w/images/5/56/Silva.seed_v119.tgz  
$tar xvzf Silva.seed_v119.tgz  
$mv silva.seed_v119.align v119silva.SSU.fasta  
$mv silva.seed_v119.tax v119fangsilva.SSU.silva.tax #I recommend you to change all parentheses for brackets.  
$perl ../scripts/SILVA_ENLONG.pl #To build v119silva.SSU_PICK.fasta  

Download SILVA_119_LSURef_tax_silva_full_align_trunc.fasta.gz from SILVA website. I clustered SILVA LSU database (v119) with Usearch (v5.2.32) at a 99% sequence identity level to decrease LSU database and made a taxonomy file (v119silva.LSU.silva.tax) for mothur.  
  
$wget https://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/SILVA_119_LSURef_tax_silva_full_align_trunc.fasta.gz  
$gunzip -c SILVA_119_LSURef_tax_silva_full_align_trunc.fasta.gz > SILVA_119_LSURef_tax_silva_full_align_trunc.fasta  
$../mothur "#get.seqs(fasta=SILVA_119_LSURef_tax_silva_full_align_trunc.fasta,accnos=v119silva.LSU.silva.tax)"  
$mv SILVA_119_LSURef_tax_silva_full_align_trunc.pick.fasta v119silva.LSU.fasta   
$cd ..  

You can also use other versions of database but I haven't tested them. Please check whether the primer site is in the alignment.  
  
  
2.Usage  

You can just type the script without any options to get help.  
  
2.1 Pipeline: SSU primer evaluation  

Prepare an <input> multiple fasta file, a formatted <primer_file>, a proper <reference_file> and put them in MIPE directory.  
To choose a proper reference file is the most important. Please make sure the reference fasta and primer can get well aligned as this position is the standard to find primer binding site.     
  
The file of primers should be formatted. Primer format:"8F:AGAGTTTGATYMTGGCTCAG". Forward primer name: end with F and no R or r is allowed in name. Reverse primer name: end with R and no F or f is allowed in name.    
  
Example:  
$cp ./mock_dataset/mockinput.fasta ./ #end with .fasta  
$cp ./mock_dataset/Ref_16S2.fasta ./ #end with .fasta   
$cp ./mock_dataset/mockprimer8.txt ./  
$perl ./scripts/MIPE_check_import_format.pl -i mockinput.fasta -p mockprimer8.txt -r Ref_16S2.fasta -n 88888 #Change your input file if any ERR occurs.  
$perl ./scripts/MIPE_alignment_taxonomy_hmm.pl -i mockinput.fasta -r Ref_16S2.fasta -s 10 -c 4 -t 50 -n 88888  
$perl ./scripts/MIPE_Ref_primer_align.pl -i mockinput.fasta -p mockprimer8.txt -r Ref_16S2.fasta -n 88888  
$perl ./scripts/MIPE_get_primer_binding_sites_hmm_bigdata.pl -i mockinput.fasta -p mockprimer8.txt -r Ref_16S2.fasta -l 4 -k 01 -n 100  
$perl ./scripts/MIPE_matchtype_stat_bigdata.pl -i mockinput.fasta -p mockprimer8.txt -w ssu_merge_88888.silva.wang.taxonomy -n 88888 -r Ref_16S2.fasta  
$rm -f *88888* *tmp* 88888* *.log   

You can delete any file you want but you'd better keep ssu_merge_88888.silva.wang.taxonomy and ssu_merge_88888.pick.filter.fasta as they are the alignment and taxonomy files. The output files are in result_88888_mockinput. All *F_match_type *R_match_type *F.stat *R.stat are results of degenerate primers.  
  
2.2 Pipeline: any gene primer evaluation  
  
Prepare an <input> multiple fasta file, a formatted <primer_file>, a proper <reference_file>, aligned file by yourself (please add the reference fasta at the end before aligning and change the name of aligned file as ssu_merge_xxxxx_.pick.filter.fasta, every fasta-format gene should only has two lines: the first line is title,like '>XXX', the second line is nucleotide alignment), a taxonomy file in mothur format(you can edit it yourself and change the name as ssu_merge_xxxxx_.silva.wang.taxonomy) and put them in MIPE directory.  
  
To choose a proper reference file is the most important. Please make sure the reference fasta and primer can get well aligned as this position is the standard to find primer binding site.   
  
The file of primers should be formatted. Primer format:"8F:AGAGTTTGATYMTGGCTCAG". Forward primer name: end with F and no R or r is allowed in name. Reverse primer name: end with R and no F or f is allowed in name.  
  
Example:  
$cp ./mock_dataset/mockinput.fasta ./ #end with .fasta  
$cp ./mock_dataset/Ref_16S2.fasta ./ #end with .fasta   
$cp ./mock_dataset/mockprimer8.txt ./  
$cp ./mock_dataset/ssu_merge_88888.silva.wang.taxonomy  ./  
$cp ./mock_dataset/ssu_merge_88888.pick.filter.fasta ./  
  
$perl ./scripts/MIPE_check_import_format.pl -i mockinput.fasta -p mockprimer8.txt -r Ref_16S2.fasta -n 88888 #Change your input file if any ERR occurs.  
$perl ./scripts/MIPE_Ref_primer_align.pl -i mockinput.fasta -p mockprimer8.txt -r Ref_16S2.fasta -n 88888  
$perl ./scripts/MIPE_get_primer_binding_sites_hmm_bigdata.pl -i mockinput.fasta -p mockprimer8.txt -r Ref_16S2.fasta -l 4 -k 01 -n 100  
$perl ./scripts/MIPE_matchtype_stat_bigdata.pl -i mockinput.fasta -p mockprimer8.txt -w ssu_merge_88888.silva.wang.taxonomy -n 88888 -r Ref_16S2.fasta  
$rm -f *88888* *tmp* 88888* *.log   

You can delete any file you want but you'd better keep ssu_merge_88888.silva.wang.taxonomy and ssu_merge_88888.pick.filter.fasta as they are the alignment and taxonomy files. The output files are in result_88888_mockinput. All *F_match_type *R_match_type *F.stat *R.stat are results of degenerate primers.  
  
2.2.1 HMMER alignment  
I suggest you to use HMMER to align as I prepared some scripts to support it. If you use HMMER to align, please use hmmreformat_nosto.pl (with stockholm file) or hmmreformat.pl (without stockholm file) to reformat it.  
$perl ./scripts/hmmreformat_nosto.pl input_align output  
$perl ./scripts/hmmreformat_nosto.pl input_sto output_align output  
  

#######################################################################################################################################
  
Citing MIPE:   
  
Zou, B., Li, J., Zhou, Q., & Quan, Z. X. (2017). MIPE: A metagenome-based community structure explorer and SSU primer evaluation tool. PloS one, 12(3), e0174609.  
  
Mao DP, Zhou Q, Chen CY, Quan ZX (2012). Coverage evaluation of universal bacterial primers using the metagenomic datasets. BMC Microbiol.3;12:66.  

