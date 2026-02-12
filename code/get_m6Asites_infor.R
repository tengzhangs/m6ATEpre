load("./peak_exon_region.Rdata")
get_m6Asites_infor <- function(high_m6Asites,peak_infor)
peak_sites <- obtain_peak_exon[[1]]
peak_gene_infor <- obtain_peak_exon[[2]]
new_peakID <- intersect(names(peak_sites),high_m6Asites$Seq_ID)
peakID_genename <- peak_gene_infor[,c(4,ncol(peak_gene_infor))]

new_peaksites <- peak_sites[which(!is.na(match(names(peak_sites),new_peakID)))]
m6A_sites_infor <- data.frame()

for (i in 1:length(new_peakID)) {
  one_match_peak <- new_peaksites[which(!is.na(match(names(new_peaksites),new_peakID[i])))]
  one_peak_width <- as.numeric(unlist(width(one_match_peak)))
  one_m6A_position <- as.numeric(high_m6Asites[high_m6Asites$Seq_ID==new_peakID[i],]$Position)
  if(length(one_peak_width)==1){
    m6A_start_sites <- as.numeric(unlist( start(one_match_peak)))+one_m6A_position
  }
  if(length(one_peak_width)>1){
    
    cumsum_width <- cumsum(one_peak_width)
    select_range_label <- vector()
    for (j in 1:length(one_m6A_position)) {
      select_range_label[j] <- which(one_m6A_position[j]<cumsum_width)[1]
    }
    
    m6A_start_sites <- as.numeric(unlist( start(one_match_peak)))[select_range_label]+one_m6A_position
  }
  seq_name <- unique(as.character(seqnames(unlist(one_match_peak))))
  seq_strand <- unique(as.character(strand(unlist(one_match_peak))))
  one_genename <- unique(peakID_genename[which(!is.na(match(peakID_genename$peak_ID,new_peakID[i]))),]$name)
  onem6A_infor <- data.frame(peaknum=as.character(names(one_match_peak)),seqnames=seq_name,
                             start=m6A_start_sites,width=1,strand=seq_strand,gene_name=one_genename)
  m6A_sites_infor <- rbind(m6A_sites_infor,onem6A_infor)
}
return(m6A_sites_infor)
}

