library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
genome <- BSgenome.Hsapiens.UCSC.hg19
get_peak_seq <- function(peaks){
  peaks_seq <- DNAStringSet()
  for (i in 1:length(peaks)) {
    onepeak<- peaks[[i]]
    
    one_peak_seq <-getSeq(genome, onepeak)
    one_comb_seq <- do.call("c", as.list(one_peak_seq))
    one_selectseq <-  DNAStringSet( one_comb_seq)
    peaks_seq <- c(peaks_seq,one_selectseq)
  }
  ###
  motif <- c("GGACA", "GGACC", "GGACT", "AGACA", "AGACC", "AGACT", 
             "GAACA", "GAACC", "GAACT", "AAACA", "AAACC", "AAACT",
             "TGACA", "TGACC", "TGACT", "TAACA", "TAACC", "TAACT")
  cag_loc <- vmatchPattern(motif[1], peaks_seq)
  motif_start <- start(cag_loc)
  
  for(i in 2:length(motif)) {
    cag_loc0 <- vmatchPattern(motif[i], peaks_seq)
    motif_start0 <- start(cag_loc0)
    motif_start <- mapply(c, motif_start, motif_start0, SIMPLIFY=FALSE)
  }
  
  fl <- function(x){
    if(length(x)==0){
      return(0)
    }else{
      return(1)
    }
    
  }
  motif_select_peak <- which(sapply(motif_start,fl)>0)
  motif_peak_line <- as.character(names(peaks))[motif_select_peak]
  ##motif peak sequencing select
  motif_peak_seq <- peaks_seq[motif_select_peak]
  motif_m6A_peak <- peaks[motif_select_peak]
  names(motif_peak_seq) <- motif_peak_line
  # writeXStringSet(motif_peak_seq,paste0(paths,"new_motif_peak_seq.fa"))
  return(motif_peak_seq)
}
load("./peak_exon_region.Rdata")
getpeaks_seq <- get_peak_seq(peaks = obtain_peak_exon[[1]])
save(getpeaks_seq,file = "./all_peak_seq.Rdata")
