library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")

blast_rslt <- read.table("tblastn_functionnal_or_vs_genome.tblastn", header=FALSE, sep="\t")
colnames(blast_rslt) <- c("query", "sseqid", "pident", "length", "mismatch", 
                          "gapopen", "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore", "qframe", "sframe", "qlen", "slen")




#Extract regions of functionnals olfactory receptors

regions_functionnals_or <- read.table("Coordinates_Functionnal_ORS.txt", header=FALSE, sep="\t")

colnames(regions_functionnals_or) <- c("scaffold", "start", "end")


#regions_functionnals_or <- blast_rslt %>% filter(pident == 100 & qstart== 1 & qend == qlen) %>% select(sseqid, sstart, send) 

#regions_functionnals_or <- regions_functionnals_or %>% mutate(i_start = case_when(
#  sstart < send ~ sstart,
#  send < sstart ~ send
#))


#regions_functionnals_or <- regions_functionnals_or %>% mutate(i_end = case_when(
#  sstart < send ~ send,
#  send < sstart ~ sstart
#))


blast_rslt <- blast_rslt %>% mutate(i_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))


blast_rslt <- blast_rslt %>% mutate(i_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))





#to_filter <- data.frame(query=NA, sseqid=NA, pident=NA, length=NA, mismatch=NA, gapopen=NA, qstart=NA, qend=NA, sstart=NA, send=NA, evalue=NA, bitscore=NA, qframe=NA, sframe=NA, qlen=NA)[numeric(0), ]
#for(row in 1:nrow(regions_functionnals_or)){
#  to_filter <- dplyr::bind_rows(to_filter, dplyr::filter(blast_rslt, sseqid == as.character(regions_functionnals_or[row, "scaffold"]) & i_start >= (regions_functionnals_or[row, "start"] - 100) &  i_end <= (regions_functionnals_or[row, "end"] + 100)))
#}

to_filter <- data.frame(query=NA, sseqid=NA, pident=NA, length=NA, mismatch=NA, gapopen=NA, qstart=NA, qend=NA, sstart=NA, send=NA, evalue=NA, bitscore=NA, qframe=NA, sframe=NA, qlen=NA)[numeric(0), ]
for(row in 1:nrow(regions_functionnals_or)){
  to_filter <- dplyr::bind_rows(to_filter, dplyr::filter(blast_rslt, sseqid == as.character(regions_functionnals_or[row, "scaffold"]) & i_start >= (regions_functionnals_or[row, "start"] - 100) &  i_start <= (regions_functionnals_or[row, "end"] + 100)))
}


to_filter <- data.frame(query=NA, sseqid=NA, pident=NA, length=NA, mismatch=NA, gapopen=NA, qstart=NA, qend=NA, sstart=NA, send=NA, evalue=NA, bitscore=NA, qframe=NA, sframe=NA, qlen=NA)[numeric(0), ]
for(row in 1:nrow(regions_functionnals_or)){
  to_filter <- dplyr::bind_rows(to_filter, dplyr::filter(blast_rslt, sseqid == as.character(regions_functionnals_or[row, "scaffold"]) & i_end >= (regions_functionnals_or[row, "start"] - 100) &  i_end <= (regions_functionnals_or[row, "end"] + 100)))
}



#filtered_blast_rslt <- blast_rslt[!(blast_rslt$sseqid %in% to_filter$sseqid),]

all <-rbind(blast_rslt,to_filter)
filtered_blast_rslt <- all[!duplicated(all,fromLast = FALSE)&!duplicated(all,fromLast = TRUE),] 




test <- filtered_blast_rslt %>% dplyr::select(sseqid, sstart, send, evalue, length, query, slen) 



test <- test %>% mutate(i_start = case_when(
  sstart < send ~ sstart,
  send < sstart ~ send
))


test <- test %>% mutate(i_end = case_when(
  sstart < send ~ send,
  send < sstart ~ sstart
))

test <- test %>% mutate(strand = case_when(
  sstart < send ~ "+",
  send < sstart ~ "-"
))


test <- test %>% select(sseqid, i_start, i_end, evalue, strand, length, query, slen)

colnames(test) <- c("seqnames", "start", "end", "evalue", "strand", "length", "query", "slen")



test_irange <- test %>% as_granges()


test_disjoin <- reduce(test_irange,with.revmap=TRUE)


list_revmap <- as.data.frame(mcols(test_disjoin))



filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(test,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$evalue)])
}


Best_hits <- slice(test, filtered_data)


### tests

Best_hits <- Best_hits %>% mutate(i_start = case_when(
  start > 50 ~ start-50,
  start < 50 ~ as.numeric(start)
))


Best_hits <- Best_hits %>% mutate(i_end = end+50)

Best_hits <- Best_hits %>% select(seqnames, i_start, i_end, evalue, strand, length, query, slen, start, end)
colnames(Best_hits) <- c("seqnames", "start", "end", "evalue", "strand", "length", "query", "slen", "true_start", "true_end")

test_irange <- Best_hits %>% as_granges()

test_disjoin <- reduce(test_irange,with.revmap=TRUE)

list_revmap <- as.data.frame(mcols(test_disjoin))

filtered_data <- c()
for(i in 1:nrow(list_revmap)){
  filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(Best_hits,  slice(list_revmap, i) %>% unlist(use.names=FALSE))$evalue)])
}


Pseudo_truncated_coordinates <- slice(Best_hits, filtered_data)


Pseudo_truncated_coordinates <- Pseudo_truncated_coordinates %>% mutate(new_coord_start =  case_when(
  start > 501 ~ start-500,
  start < 501 ~ 1
))

Pseudo_truncated_coordinates <- Pseudo_truncated_coordinates %>% mutate(new_coord_end =  case_when(
  end < slen-501 ~ end+500,
  end > slen-501 ~ as.numeric(slen)
))

Pseudo_truncated_coordinates <- Pseudo_truncated_coordinates %>% select(seqnames, true_start, true_end, query)

write.table(Pseudo_truncated_coordinates, file='Pseudo_truncated_coordinates.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)


