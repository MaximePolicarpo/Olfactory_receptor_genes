#!/bin/bash


#This script allows to find pseudogene . It need a genome as variable 1 and the ID list of functionnals OR in a file in variable 2
#The file with ID list is extracted from iTOL. It contains leaves that are clustered with known zebrafish OR genes in a maximum likelihood tree
#The folder should contain the files : clearer_ambigous_nt.py and the interproscan launcher (interproscan.sh)

#Edges : Genes at contigs borders or next to N nucleotides
#Truncated : OR Genes with no open reading frames and without stop/fs
#Pseudogenes : OR genes with stop codons or frameshifts
#Functionnal : Intact OR genes


genome=$1
leaves_from_ML_tree=$2



#Remove sequences corresponding to zebrafish database
grep -v "Dare\|TaruOR\|DilaOR\|CyseOR\|LeocOR\|GaacOR" $leaves_from_ML_tree | sed 's/ /_/g' | grep -i "^[A-Z]\|^[0-9]" > OR_from_ML.list


#extract protein sequences of functionnals ORs
for i in `cat OR_from_ML.list` ; do samtools faidx Putative_or_multifasta.prot $i >> Fasta_OR_from_ML.prot ; done 

#extract CDS sequences of functionnals ORs
for i in `cat OR_from_ML.list` ; do samtools faidx Putative_or_multifasta.fa $i >> Functionnal_ORs_multifasta.fa ; done

#remove every sequences that are exactly identic and then remove the one that have ambiguous nucleotides

cd-hit -i Functionnal_ORs_multifasta.fa -o Functionnal_ORs_multifasta_cdhit.fa -c 1

python2 clearer_ambigous_nt.py Functionnal_ORs_multifasta_cdhit.fa 0 0 #resulting file : clear_Functionnal_ORs_multifasta_cdhit.fa

#create a protein multifasta with functionnals ORs without ambiguous nucleotides


transeq clear_Functionnal_ORs_multifasta_cdhit.fa clear_Functionnal_ORs_multifasta_cdhit.prot


IFS=$'\n'

#Put every "not clear" sequences in our pseudogene list
grep ">" clear_Functionnal_ORs_multifasta_cdhit.fa | sed 's/>//g' > clears_ID.list
grep ">" Functionnal_ORs_multifasta_cdhit.fa | sed 's/>//g' > unclear_IDs.list

for i in `cat unclear_IDs.list` 
	do 
		if grep -q "$i" clears_ID.list

			then echo "clear"
		
		else 
			echo "$i" >> Ambiguous_sequences.list

		fi
done



#tlbastn our functionnals ORs against the genome with an evalue of 1e-20 to find every functionnals/non-functionnals ORs in the genome
tblastn -evalue 1e-20 -query clear_Functionnal_ORs_multifasta_cdhit.prot -db $genome -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen' -out tblastn_functionnal_or_vs_genome.tblastn -num_threads 20


#extract the coordinates of functionnals ORs found 
grep ">" Functionnal_ORs_multifasta.fa | sed 's/>//g' | sed 's/-/	/g' > Coordinates_Functionnal_ORS.txt




#We then use a R script that will extract all best-hit sequences and remove the one that have an overlap with a already found functionnal OR (with a -100/+100 window around them)
Rscript R_script_numero2_2.R

#resulting file : Pseudo_truncated_coordinates.tsv



#extract coordinates of found pseudogenes and concatenate with our first pseudogene list (genes that had N nucleotides)


sed -i 's/-/	/g' Ambiguous_sequences.list

cat Pseudogenes_Ids.list  > OR_genes_with_N.tsv

cat Pseudo_truncated_coordinates.tsv > ALL_pseudo.tsv



#Now we are goind to make a differenciation between pseudogenes and truncated genes#
#Truncated genes are genes that are <30bp close to the beginning or end of a contig

IFS=$'\n'

samtools faidx $genome

for i in `cat ALL_pseudo.tsv` ; do
	echo "$i" > ligne_en_cours

	scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ;  gene=`cut -f4 ligne_en_cours` ; blast_start=`cut -f5 ligne_en_cours` ; blast_end=`cut -f6 ligne_en_cours`
	scaffold_length=`grep "^$scaffold	" $genome.fai | cut -f2`
	contig_start=30 ; contig_end=$((scaffold_length-30))

	samtools faidx $genome $scaffold:$start-$end > Region_Genomique.fa
	samtools faidx $genome $scaffold:$blast_start-$blast_end > Region_Genomique_2.fa
	makeblastdb -in Region_Genomique.fa -dbtype nucl
	makeblastdb -in Region_Genomique_2.fa -dbtype nucl  
	samtools faidx clear_Functionnal_ORs_multifasta_cdhit.prot $gene > Gene.fa


	number_of_frames=`tblastn -query Gene.fa -db Region_Genomique_2.fa -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe' | cut -f12 | sort | uniq | wc -l` 
	number_of_stop=`tblastn -query Gene.fa -db Region_Genomique.fa -evalue 1e-20 | grep "Sbjct" | grep -o "\*" | wc -l`

	frame=`tblastn -query Gene.fa -db Region_Genomique.fa -evalue 1e-20 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sframe' | cut -f12 | head -1`


	check_1=`if [ $end -gt $contig_end ] ; then echo "CONTIG_END" ; else echo "NA" ; fi`
	check_2=`if [ $start -lt $contig_start ] ; then echo "CONTIG_START" ; else echo "NA" ; fi`


	echo "$scaffold	$start	$end	$check_1	$check_2	$number_of_frames	$number_of_stop	$frame" >> Truncated_and_pseudo_coordinates.tsv

done 


grep -v "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 > 1 || $7 > 0)' > Pseudo_coordinates.tsv
grep "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 > 1 || $7 > 0)' >> Pseudo_coordinates.tsv
grep "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 == 1 && $7 == 0)' > Edges_coordinates.tsv
grep -v "CONTIG" Truncated_and_pseudo_coordinates.tsv | awk '($6 == 1 && $7 == 0)' > Truncated_coordinates.tsv




#Filter genes to keep only those who have a 7tm transmembrane domain recognized by interproscan (IPR017452) 


rm temp

./interproscan.sh -i clear_Functionnal_ORs_multifasta_cdhit.prot -o 7tm_checking -f TSV

grep "IPR017452" 7tm_checking | cut -f1 | sort | uniq > Genes_without_gaps_in_7tm

grep ">" clear_Functionnal_ORs_multifasta_cdhit.prot | sed 's/>//g' > list_initiale

for gene in `cat list_initiale` ; do if grep -q "$gene" Genes_without_gaps_in_7tm ; then echo "clear" ; else echo "$gene" >> Genes_with_gaps_in_7tm ; fi ; done



sed -i 's/_1$//g' Genes_without_gaps_in_7tm
sed -i 's/_1$//g' Genes_with_gaps_in_7tm



#Some genes without 7tm cna be due to contig end#

IFS=$'\n'

for line in `cat Genes_with_gaps_in_7tm` ; do
	echo "$line" > temporaryfile 

	sed -i 's/-/	/g' temporaryfile
	scaffold=`cut -f1 temporaryfile` ; start=`cut -f2 temporaryfile` ; end=`cut -f3 temporaryfile`
	scaffold_length=`grep "^$scaffold	" $genome.fai | cut -f2`
	contig_start=30 ; contig_end=$((scaffold_length-30))

	prot_to_grep=`cat temporaryfile | sed 's/	/-/g'`

	if [ $end -gt $contig_end ] || [ $start -lt $contig_start ] ; then
		samtools faidx clear_Functionnal_ORs_multifasta_cdhit.fa $prot_to_grep >> Edges_genes_toadd.fa
	else
		samtools faidx clear_Functionnal_ORs_multifasta_cdhit.fa $prot_to_grep >> Pseudo_toadd.fa
	fi
done




#Extract CDS sequences of truncated and pseudogenes and create corresponding proteins sequences file

for line in `cat Edges_coordinates.tsv` ; do echo "$line" > ligne_en_cours ; scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ; frame=`cut -f8 ligne_en_cours` ; if [ $frame -lt 0 ] ; then samtools faidx $genome $scaffold:$start-$end > temp.fa ; revseq temp.fa temp_rev.fa ; cat temp_rev.fa >> Edges_multifasta_temp.fa ; else samtools faidx $genome $scaffold:$start-$end >> Edges_multifasta_temp.fa ; fi ; done 
cat Edges_multifasta_temp.fa Edges_genes_toadd.fa > Edges_multifasta.fa
sed -i 's/>/>E_/g' Edges_multifasta.fa

for line in `cat Pseudo_coordinates.tsv` ; do echo "$line" > ligne_en_cours ; scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ; frame=`cut -f8 ligne_en_cours` ; if [ $frame -lt 0 ] ; then samtools faidx $genome $scaffold:$start-$end > temp.fa ; revseq temp.fa temp_rev.fa ; cat temp_rev.fa >> Pseudogenes_multifasta_temp.fa ; else samtools faidx $genome $scaffold:$start-$end >> Pseudogenes_multifasta_temp.fa ; fi ; done
cat Pseudogenes_multifasta_temp.fa Pseudo_toadd.fa > Pseudogenes_multifasta.fa
sed -i 's/>/>P_/g' Pseudogenes_multifasta.fa 

for line in `cat Truncated_coordinates.tsv` ; do echo "$line" > ligne_en_cours ; scaffold=`cut -f1 ligne_en_cours` ; start=`cut -f2 ligne_en_cours` ; end=`cut -f3 ligne_en_cours` ; frame=`cut -f8 ligne_en_cours` ; if [ $frame -lt 0 ] ; then samtools faidx $genome $scaffold:$start-$end > temp.fa ; revseq temp.fa temp_rev.fa ; cat temp_rev.fa >> Truncated_multifasta.fa ; else samtools faidx $genome $scaffold:$start-$end >> Truncated_multifasta.fa ; fi ; done
sed -i 's/>/>T_/g' Truncated_multifasta.fa 




for i in `cat Genes_without_gaps_in_7tm` ; do samtools faidx clear_Functionnal_ORs_multifasta_cdhit.fa $i ; done > RESULTS_Functionnal_genes.fa



mv Pseudogenes_multifasta.fa RESULTS_TEMP_Pseudogenes.fa
mv Edges_multifasta.fa RESULTS_Temp_Edges.fa



#In a final step, we are going to remove pseudogenes that are too much identical to functionnal genes which are most likely duplicate results#

mv RESULTS_TEMP_Pseudogenes.fa RESULTS_Pseudogenes.fa
mv RESULTS_Temp_Edges.fa RESULTS_Edges_genes.fa
mv Truncated_multifasta.fa RESULTS_Truncated_genes.fa

#Count file with the number of functionnal, truncated and pseudogenes

x1=`grep -c ">" RESULTS_Functionnal_genes.fa`
x2=`grep -c ">" RESULTS_Edges_genes.fa` 
x3=`grep -c ">" RESULTS_Pseudogenes.fa`
x4=`grep -c ">" RESULTS_Truncated_genes.fa` 
echo "$x1	$x2	$x3	$x4" > Count_ORs_with_redundancy.tsv



#Classify OR genes in families

transeq RESULTS_Functionnal_genes.fa RESULTS_Functionnal_genes.prot
sed -i 's/_1$//g' RESULTS_Functionnal_genes.prot
transeq RESULTS_Edges_genes.fa RESULTS_Edges_genes.prot
sed -i 's/_1$//g' RESULTS_Edges_genes.prot
transeq RESULTS_Truncated_genes.fa RESULTS_Truncated_genes.prot
sed -i 's/_1$//g' RESULTS_Truncated_genes.prot

#Remove redundancy

cat RESULTS_Functionnal_genes.fa RESULTS_Edges_genes.fa RESULTS_Truncated_genes.fa > Functionnals_Edges_Truncated.fa

sed -i 's/:/_/g' Functionnals_Edges_Truncated.fa 


python2 clearer_ambigous_nt.py Functionnals_Edges_Truncated.fa 0 0

cd-hit -i Functionnals_Edges_Truncated.fa -o Functionnals_Edges_Truncated_cdhit.fa -c 1

mv clear_Functionnals_Edges_Truncated_cdhit.fa Functionnals_Edges_Truncated_cdhit.fa

python2 clearer_ambigous_nt.py RESULTS_Pseudogenes.fa 0 0

mv clear_RESULTS_Pseudogenes.fa RESULTS_Pseudogenes.fa


x1=`grep -v ">E_\|>T_" Functionnals_Edges_Truncated_cdhit.fa | grep -c ">"`
x2=`grep -c ">E_" Functionnals_Edges_Truncated_cdhit.fa` 
x4=`grep -c ">T_" Functionnals_Edges_Truncated_cdhit.fa`
echo "$x1	$x2	$x3	$x4" > Count_ORs_without_redundancy.tsv


#reformat

sed -i 's/ .*//g' Functionnals_Edges_Truncated_cdhit.fa ; sed -i 's/_/-/g' Functionnals_Edges_Truncated_cdhit.fa

transeq Functionnals_Edges_Truncated_cdhit.fa Functionnals_Edges_Truncated_cdhit.prot

sed -i 's/_1$//g' Functionnals_Edges_Truncated_cdhit.prot

sed -i 's/ .*//g' RESULTS_Pseudogenes.fa ; sed -i 's/_/-/g' RESULTS_Pseudogenes.fa ; sed -i 's/:/-/g' RESULTS_Pseudogenes.fa


#ML tree to classify


cat Functionnals_Edges_Truncated_cdhit.prot Classification_OR_multifasta.prot Classification_outgroups.prot > Classification_fasta.prot

mafft --maxiterate 1000 --oldgenafpair --thread -40 Classification_fasta.prot > Classification_fasta.prot.aln

iqtree -s Classification_fasta.prot.aln -st AA -nt 50 -bb 1000 -m JTT+F+I+G4


#Now that OR genes are mined, we can classify them based on th maximum likelihood tree and one can classify pseudogenes based on best blastx match




