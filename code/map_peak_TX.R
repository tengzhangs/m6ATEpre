##obtain m6A peak splitted by exon
f1 <- "./m6APeak_HeLa.csv"
m6Apeak_CTL <- read.csv(f1,header = T)
fa <- "./Mod.bed"
all_m6Apeakbed <-  read.table(fa,sep="\t",header=FALSE,stringsAsFactors =FALSE)
all_m6Apeakbed$V4 <- as.character(m6Apeak_CTL$geneID)
all_m6Apeakbed$peak_ID <- as.character(m6Apeak_CTL$name)
colnames(all_m6Apeakbed) <- c("chr",         "start" ,      "end",         "name",        "score",       "strand",      "thickStart",  "thickEnd",    "itemRgb",    
                              "blockCount",  "blockSizes",  "blockStarts","peak_ID")

GTF_file <- "./hg19_gtf/genes.gtf"
library(GenomicFeatures)

map_peak_TX <- function(peak_sites_infor,annotation_file){
  ##obtain the longest transcript
  ##mapped to the longest transcript
  # a = read.table(filepath,sep="\t",header=FALSE,stringsAsFactors =FALSE)
  a1 = peak_sites_infor
  a = a1[,-ncol(a1)]
  txdbfile <- GenomicFeatures::makeTxDbFromGFF(annotation_file)
  genes_txdb <- genes(txdbfile)
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
  
  read_peak <- GRanges(seqnames = as.character(a$chr),
                       IRanges(start = as.numeric(as.character(a$start)),
                               end = as.numeric(as.character(a$end))),strand = as.character(a$strand))
  mcols(read_peak)$name <- as.character(a$name)
  mappeak_to_longtx <-  mapToTranscripts(read_peak, exbytx_txdb,ignore.strand=F)
  select_peak_label <- as.numeric(as.character(mappeak_to_longtx$xHits))
  a = a[select_peak_label,]
  a1 = a1[select_peak_label,]
  selectpeak_genomic <- mapFromTranscripts(mappeak_to_longtx,exbytx_txdb)
  select_peaks <- selectpeak_genomic[countOverlaps(selectpeak_genomic,selectpeak_genomic) == 1]
  last_select_label <- as.numeric(as.character(select_peaks$xHits))
  # # mcols_info =a[,13:length(a[1,])]
  # # a = a[last_select_label,1:12]
  a = a[last_select_label,]
  a1 = a1[last_select_label,]
  
  #################
  # get transcripts
  no_tx = length(a[,1])
  tx_id = 1:no_tx;
  # tx_name = paste("peak",1:no_tx,sep="")
  tx_name = as.character(a1$peak_ID)
  tx_chrom = a[,1]
  # tx_strand = a[,6]
  tx_strand = a[,6]
  # tx_start = a[,2]+1
  tx_start = a[,2]
  tx_end = a[,3]
  transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
  
  # get genes
  tx_name = tx_name
  gene_id = as.character(a[,4])
  gene_id[is.na(gene_id)]="NA"
  gene=data.frame(tx_name,gene_id)
  
  # no_tx
  splicing <- lapply(1:no_tx, .spliceSingleTrans, a=a, tx_start=tx_start)
  splicing <- .combineListOfSplicing(splicing)
  
  # make txdb
  peaks = suppressWarnings(
    makeTxDb(transcripts=transcripts, 
             splicings=splicing,
             genes=gene))
  
  # generate GRangesList
  tx <- exonsBy(peaks, "tx",use.names=TRUE)
  mcols(tx) <- data.frame(mcols(tx),gene_name=as.character(a$name))
  # tx_GR <- unlist(tx)
  # mappeak_to_longtx <-  mapToTranscripts(tx_GR, exbytx_txdb,ignore.strand=F)
  # select_peak_label <- as.numeric(as.character(mappeak_to_longtx$xHits))
  # tx_map_GR <- tx_GR[select_peak_label]
  # tx_peak_GR <- split(tx_map_GR,names(tx_map_GR))
  # new_peakname <- as.character(names(tx_peak_GR))
  new_peakname <- as.character(names(tx))
  new_a <- data.frame()
  for (i in 1:length(new_peakname)) {
    one_a <- a1[which(!is.na(match(a1$peak_ID,new_peakname[i]))),]
    new_a <- rbind(new_a,one_a)
  }
  # mcols(tx_peak_GR) <- data.frame(as.character(names(tx_peak_GR)),gene_name=as.character(new_a$name))
  # results <- list(tx_peak_GR,new_a)
  results <- list(tx,new_a)
  return(results)

}

.combineListOfSplicing <- function(t){
  
  a <- paste("t[[",1:length(t),"]]", sep="")
  a <- paste(a,collapse =",")
  a <- paste("rbind(",a,")",sep="")
  c <- parse(text=a)
  b <- suppressWarnings(eval(c))
  
  return(b)
}

.spliceSingleTrans <- function(i,a,tx_start) {
  tx = a[i,]
  tx_id = i
  exon_rank=1:as.integer(tx[10])
  
  # get start
  temp = as.integer(strsplit(as.character(tx[12]), ",")[[1]]) + tx_start[i]
  exon_start=temp
  
  # get end
  temp = as.integer(strsplit(as.character(tx[11]), ",")[[1]])
  temp2 = temp + exon_start - 1
  exon_end=temp2
  
  # get CDS
  cds_start = exon_start
  cds_end = exon_end
  
  # get data frame
  splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
  return(splicing_tx)
}
obtain_peak_exon <- map_peak_TX(peak_sites_infor=all_m6Apeakbed,annotation_file=GTF_file)
save(obtain_peak_exon,file = "D:\\research\\m6Atranslation\\data\\new_features\\peak_exon_region.Rdata")
