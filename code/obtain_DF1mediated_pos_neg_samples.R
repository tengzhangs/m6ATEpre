###load single base m6A sites with high resolution
f1 <- "./singlebase_m6Asites.txt"
Human_m6Asites <- read.delim2(f1)
HeLa_m6Asites <- Human_m6Asites[-grep("Non-m6A site",Human_m6Asites$Classification),]
high_m6Asites <- HeLa_m6Asites[grep('Very high confidence',HeLa_m6Asites$Classification),]
load("./peak_exon_region.Rdata")
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

###map YTHDF1 binding sites to m6A peak in exon 
f1 <- "./DF1_bind_gene_site.csv"
DF1bind_sites <- read.csv(f1)


DF1bind_sites$seqname<- tolower(DF1bind_sites$seqname)
Y1bind_GR <-  GRanges(seqnames = as.character(DF1bind_sites$seqname),
                      IRanges(start = as.numeric(as.character(DF1bind_sites$start)),
                              end = as.numeric(as.character(DF1bind_sites$end))),strand = as.character(DF1bind_sites$strand))

m6Asites_GR <- GRanges(seqnames = as.character(m6A_sites_peak$seqnames),
                       IRanges(start = as.numeric(as.character(m6A_sites_peak$start)),
                               end = as.numeric(as.character(m6A_sites_peak$start))),strand = as.character(m6A_sites_peak$strand))

#####
DF1_start <- as.numeric(as.character(DF1bind_sites$start))
DF1_end <- as.numeric(as.character(DF1bind_sites$end))
DF1_center <- DF1_start+round((DF1_end-DF1_start+1)/2)
DF1_center_GR <- GRanges(seqnames = as.character(DF1bind_sites$seqname),
                         IRanges(start = DF1_center,width = 1),
                         strand = as.character(DF1bind_sites$strand))
DF1_range <- resize(DF1_center_GR, 101,fix="center")
m6AY1bind_overlap <- findOverlaps(m6Asites_GR,DF1_range,type = "within")
length(unique(m6AY1bind_overlap @from))
######train sites
select_m6A_sites <- m6A_sites_peak[unique(m6AY1bind_overlap @from),]
select_DF1_sites <- DF1bind_sites[unique(m6AY1bind_overlap @to),]
##positvie gene and negative genes with m6A sites 
f2 <- "D:\\research\\m6Atranslation\\data\\Y1_TE_change.xlsx"
Y1TE_change <- readxl::read_xlsx(f2)
Y1TE_down <- Y1TE_change[which(Y1TE_change$`Translation efficiency`<(-0.5)),]
Y1KD_TEup <- Y1TE_change[which(Y1TE_change$`Translation efficiency`>(-0.5)),]

pos_gene_sites <- select_m6A_sites[which(!is.na(match(select_m6A_sites$gene_name,Y1TE_down$`Gene symbol`))),]
neg_gene_sites <- select_m6A_sites[which(!is.na(match(select_m6A_sites$gene_name,Y1KD_TEup$`Gene symbol`))),]
training_m6A_sites <- list(pos_gene_m6Asites=pos_gene_sites,neg_gene_m6Asites=neg_gene_sites,DF1_sites=select_DF1_sites)
save(training_m6A_sites,file = "./pos_neg_sites_DF1infor.Rdata")
