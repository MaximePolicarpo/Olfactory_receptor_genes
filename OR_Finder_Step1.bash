#!/bin/bash


#This script allows to find every functionnals olfactory receptors in a genome. It has to be run with 2 variables : a genome and a OR database
#One can apply several more filters to resulting files like a 7tm check. 


#The folder on which scripts are run should contain the file rename_fasta.pl, R_script_numero1.R, R_script_numero2.R, Classification_outgroups.prot and Classification_OR_multifasta.prot



#tblastn of all or fish against the genome of interest 

genome=$1
or_database=$2
Unitprot_database=$3

#Remove spaces in scaffold name and transform every nucleotides in uppercase
sed -i "s/ .*//g" $genome
sed -i "s/a/A/g" $genome
sed -i "s/t/T/g" $genome
sed -i "s/g/G/g" $genome
sed -i "s/c/C/g" $genome

#Prepare the genome database
makeblastdb -in $genome -dbtype nucl

#Perform a tblastn of OR database against the genome of interest
tblastn -query $or_database -db $genome -evalue 1e-10 -outfmt 6 -out OR_vs_Genome.blastn -num_threads 30

#Rscript that permit to extract Best-hits that are non-overlapping regions with tblastn match on it
Rscript R_script_numero1.R

#Prepare a genome .fai file
samtools faidx $genome

#Remove the file Putative_or_multifasta_raw.fa if a run has already be done
rm Putative_or_multifasta_raw.fa


#For each best-hits we will extract ORFs that are longer than 750 bp
IFS=$'\n'

for line in `cat Best_hits_filtered.tsv` ; do
	
	echo "$line" > test


	scaffold=`cut -f1 test`
	initial_start=`cut -f2 test` 
	initial_stop=`cut -f3 test` 

	new_3prime=$((initial_start - 1000)) #extend the best-hit 1000bp upstream an downtream to make sure we get the ATG and TAG
	new_5prime=$((initial_stop + 1000)) 

	if [ "$new_3prime" -lt '1' ] ; then  #If the extension go below 1 then we put it to 1 (scaffold beginning)
		new_3prime=1
	fi


	samtools faidx $genome $scaffold:$new_3prime-$new_5prime > sequence.fa 
	getorf -sequence sequence.fa -outseq orf_list -minsize 750 -find 3 
	sed -i "s/>/>$scaffold:/g" orf_list ; sed -i 's/ //g' orf_list
	cat orf_list >> Putative_or_multifasta_raw.fa


done


#Rename fasta headers for it to be scaffold:start-stop

sed -i 's/(REVERSESENSE)/_reverse/g' Putative_or_multifasta_raw.fa
grep ">" Putative_or_multifasta_raw.fa | sed 's/>//g' > oldfastaheaders

sed 's/:/	/g' oldfastaheaders | sed 's/-/	/g' | sed 's/_[0-9]\[/	/g' | sed 's/\]//g' | sed 's/_reverse/	reverse/g' > temporary

IFS=$'\n'

for line in `cat temporary`
	do
		echo "$line" > temp 

		scaffold=`cut -f1 temp`

		if grep -q "reverse" temp ; then
			coord_start=`awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}' temp`
			coord_end=`awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}' temp`
		else
			coord_start=`awk 'BEGIN {FS="	"} {sum+=$2+$4-1} END {print sum}' temp`
			coord_end=`awk 'BEGIN {FS="	"} {sum+=$2+$5-1} END {print sum}' temp`

		fi


		echo "$scaffold-$coord_start-$coord_end" >> newfastaheaders
done


paste -d "	" oldfastaheaders newfastaheaders > renaming_file

perl rename_fasta.pl renaming_file Putative_or_multifasta_raw.fa > Putative_or_multifasta.fa


#Due to the -1000/+1000 extension, some orf could be found twice so we remove identical sequences


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Putative_or_multifasta.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Putative_or_multifasta_uniq.fa


#transform the cds fasta in protein fasta
transeq Putative_or_multifasta_uniq.fa Putative_or_multifasta.prot  


#remove the _1 added by the programm transeq to protein names
sed -i 's/_1$//g' Putative_or_multifasta.prot 


## two step verification : blastp followed by a maximum likelihood tree #

#perform a blastp against uniprot db 
blastp -query Putative_or_multifasta.prot -db $Unitprot_database -evalue 1e-5 -outfmt "6 qseqid sseqid sscinames scomnames pident length mismatch gapopen qstart qend sstart send evalue stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads 20

#extract all sequences that match the best against on OR
grep -i "olfactory\|odorant" blastp_result | cut -f1 | sort | uniq > OR_from_blast.list

#Extract corresponding protein sequences
for i in `cat OR_from_blast.list` ; do samtools faidx Putative_or_multifasta.prot $i >> Fasta_OR_from_blast.prot ; done 

#Put resulting sequences, OR from zebrafish and one outgroup sequence (adra2da gene), align with mafft E-INS and make a maximum likelyhod tree

cat Fasta_OR_from_blast.prot Classification_OR_multifasta.prot Classification_outgroups.prot > Putative_or_plus_known_or_plus_outgroup.prot

mafft --maxiterate 1000 --oldgenafpair --thread -40 Putative_or_plus_known_or_plus_outgroup.prot > Putative_or_plus_known_or_plus_outgroup.prot.aln

iqtree -s Putative_or_plus_known_or_plus_outgroup.prot.aln -st AA -nt 50 -bb 1000 -m JTT+F+I+G4


echo "End of the initial search and time for manual verification of the maximum likelihood tree"













