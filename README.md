# Olfactory_receptor_genes

First create a folder containg your genome file (fasta format)
This folder should also contain the following files : 

- rename_fasta.pl
- R_script_numero1.R
- R_script_numero2.R 
- Classification_outgroups.prot
- Classification_OR_multifasta.prot
- The Interproscan launcher (interproscan.sh)
- clearer_ambigous_nt.py
- All_fishes_cdhit_80.fa

Several dependencies are needed : 

- EMBOSS 6.6.0.0
- IQ-TREE 1.6.12
- samtools v1.9
- BLAST 2.6.0+
- MAFFT v7.310
- CD-HIT 4.8.1
- R with the packages  "data.table", "dplyr", "plyranges", "GenomicRanges"


How to run the pipeline ?

First step consist in finding functional OR genes from the genome. Start by launching : 

./OR_Finder_Step1.bash genome.fasta All_fishes_cdhit_80.fa path_to_unitprot_database


The results will consist in phylogenetic tree containing all the putative ORs found in the genome. You can visualise this tree
on iTOL, root using non-OR outgroup sequences and select sequences that are well clustered with OR genes. The label of good ORs (sequence name without the fasta header >) should be put in a text file, with one sequence per line. Example : 

File_good_seqs.txt :

Mysequence1

Mysequence2

Mysequence3



The second part of the pipeline concist in finding OR pseudogenes and incomplete genes :

./OR_Finder_Step2.bash genome.fasta File_good_seqs.txt


This will result in several output files : 

- RESULTS_Pseudogenes.fa containing pseudogenes
- Functionnals_Edges_Truncated_cdhit.fa containing functional and incomplete genes
- Classification_fasta.prot.aln.treefile that will help you to classify OR genes in families

One can also classify pseudogenes by performing a blastx against known fish ORs and assign pseudogenes based on their best blastx match







