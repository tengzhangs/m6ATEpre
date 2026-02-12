load("D:\\research\\m6Atranslation\\data\\single_base_level\\pos_neg_sites_DF1infor.Rdata")
pos_gene_sites <- training_m6A_sites[[1]]
neg_gene_sites <-training_m6A_sites[[2]]
library(m6ALogisticModel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(GenomicFeatures)
annotation_file <- "D:\\hg19_GTF\\genes.gtf"
get_m6A_seq <- function(target_peak_center,annotation_file){
  txdbfile <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
  exbytx_txdb <- exonsBy(txdbfile,by = "tx")
  isoform_ambiguity_method = "longest_tx"
  if(isoform_ambiguity_method == "longest_tx"){
    Longest_tx_table <- find_longest_transcript(exbytx_txdb,txdbfile)
    Kept_tx_indx <- Longest_tx_table$TXID[Longest_tx_table$longest]
    rm(Longest_tx_table)
  } else {
    Kept_tx_indx <- T
  }
  exbytx_txdb <- exbytx_txdb[Kept_tx_indx]
  exbytx_txdb <- exbytx_txdb[countOverlaps(exbytx_txdb,exbytx_txdb) == 1]
  target_peakGR <- GRanges(seqnames = as.character(target_peak_center$seqnames),
                           IRanges(start = as.numeric(as.character(target_peak_center$start)),
                                   width = 1,names=as.character(target_peak_center$peaknum)),
                           strand = as.character(target_peak_center$strand))
  target_sites_map <-   mapToTranscripts(target_peakGR, exbytx_txdb,ignore.strand=F)
  target_tx <- exbytx_txdb[target_sites_map$transcriptsHits]
  select_target_GR <- target_peakGR[target_sites_map$xHits]
  
  genome <- BSgenome.Hsapiens.UCSC.hg19
  
  
  peaks_seq <- DNAStringSet()
  for (i in 1:length(select_target_GR)) {
    if(sum(width(target_tx[[i]]))>500){
      center_start <- start(target_sites_map[i])
      center_end <- sum(width(target_tx[[i]]))-center_start+1
      one_tx <- target_tx[[i]]
      tx_seq <- getSeq(genome, one_tx)
      tx_comb_seq <- do.call("c", as.list(tx_seq))
      if((center_start>250)&(center_end>250)){
        select_seq <- tx_comb_seq[(center_start-250):(center_start+250)]
        one_selectseq <-  DNAStringSet( select_seq)
        peaks_seq <- c(peaks_seq,one_selectseq)
      }
      if((center_start>250)&(center_end==250)){
        select_seq <- tx_comb_seq[(center_start-251):(center_start+249)]
        one_selectseq <-  DNAStringSet( select_seq)
        peaks_seq <- c(peaks_seq,one_selectseq)
      }
      
      if((center_start<=250)&(center_end>=250)){
        select_seq <- tx_comb_seq[1:501]
        one_selectseq <-  DNAStringSet( select_seq)
        peaks_seq <- c(peaks_seq,one_selectseq)
      }
      if((center_start>=250)&(center_end<250)){
        select_seq <- tx_comb_seq[(length(tx_comb_seq)-250):length(tx_comb_seq)]
        if(length(select_seq)<501){
          select_seqs <- append(rep(DNAString("N"),(501-length(select_seq))),select_seq)
          
        }
        one_selectseq <-  DNAStringSet( select_seqs)
        peaks_seq <- c(peaks_seq,one_selectseq)
      }
    }
    if(sum(width(target_tx[[i]]))<=501){
      one_tx <- target_tx[[i]]
      tx_seq <- getSeq(genome, one_tx)
      tx_comb_seq <- do.call("c", as.list(tx_seq))
      if(length(tx_comb_seq)<501){
        center_start <- start(target_sites_map[i])
        center_end <- sum(width(target_tx[[i]]))-center_start+1
        if((center_start>250)&(center_end<250)){
          tx_comb_seqs <- append(tx_comb_seq,rep(DNAString("N"),(501-length(tx_comb_seq))))
          
        }
        if((center_start<=250)&(center_end>250)){
          tx_comb_seqs <- append(rep(DNAString("N"),(501-length(tx_comb_seq))),tx_comb_seq)
          
        }
        
      }
      if(length(tx_comb_seq)==501){
        tx_comb_seqs <- tx_comb_seq
      }
      one_selectseq <-  DNAStringSet( tx_comb_seqs)
      peaks_seq <- c(peaks_seq,one_selectseq)
      # print(i)
    }
    if(length(peaks_seq)<i){
      print(i)
    }
  }
  new_peakID <- names(select_target_GR)
  names(peaks_seq) <- new_peakID
  return(peaks_seq )
}
pos_peak_seq <- get_m6A_seq(target_peak_center=pos_gene_sites,annotation_file=annotation_file)
pos_peak_seqs <- as.character(pos_peak_seq)
names(pos_peak_seqs) <- NULL

neg_peak_seq <- get_m6A_seq(target_peak_center=neg_gene_sites,annotation_file=annotation_file)
neg_peak_seqs <- as.character(neg_peak_seq)
names(neg_peak_seqs) <- NULL

m6Asites_seq <- c(pos_peak_seq,neg_peak_seq)
m6Asites_seqs <- as.character(m6Asites_seq)
names(m6Asites_seqs) <- NULL

writeXStringSet(m6Asites_seq,"./allm6A_seqs.fasta")
####cdhit process
names(pos_peak_seq) <- paste0('site_',1:length(pos_peak_seq))
names(neg_peak_seq) <- paste0("site_",(length(pos_peak_seq)+1):(length(pos_peak_seq)+length(neg_peak_seq)))
all_m6A_seq <- list(pos_peak_seq=pos_peak_seq,neg_peak_seq=neg_peak_seq)
save(all_m6A_seq,file = "./cdhit_proces/proc_m6Asites_seqs.Rdata")
