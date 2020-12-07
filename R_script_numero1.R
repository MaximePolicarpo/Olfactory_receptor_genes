library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")


blast_rslt <- read.table("OR_vs_Genome.blastn", header=FALSE, sep="\t")
colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore")




blast_rslt_col_filter <- blast_rslt %>% dplyr::select(sseqid, sstart, send, evalue, length) 



blast_rslt_col_filter <- blast_rslt_col_filter %>% mutate(i_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))


blast_rslt_col_filter <- blast_rslt_col_filter %>% mutate(i_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))

blast_rslt_col_filter <- blast_rslt_col_filter %>% mutate(strand = case_when(
  sstart < send ~ "+",
  send < sstart ~ "-"
))


blast_rslt_col_filter <- blast_rslt_col_filter %>% select(sseqid, i_start, i_end, evalue, strand, length)

colnames(blast_rslt_col_filter) <- c("seqnames", "start", "end", "evalue", "strand", "length")


blast_rslt_irange <- blast_rslt_col_filter %>% as_granges()


blast_rslt_disjoin <- reduce(blast_rslt_irange,with.revmap=TRUE)


list_revmap <- as.data.frame(mcols(blast_rslt_disjoin))



filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(blast_rslt_col_filter,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$evalue)])
}


Best_hits_filtered <- slice(blast_rslt_col_filter, filtered_data)


Best_hits_filtered <- Best_hits_filtered %>% filter(length >= 100)
write.table(Best_hits_filtered, file='Best_hits_filtered.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)





