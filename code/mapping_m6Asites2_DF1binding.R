
m6Asites2_DF1binding_region <- function(m6A_sites,DF1binding_sites) {
DF1bind_sites$seqname<- tolower(DF1bind_sites$seqname)
m6Asites_GR <- GRanges(seqnames = as.character(m6A_sites$seqnames),
                       IRanges(start = as.numeric(as.character(m6A_sites$start)),
                               end = as.numeric(as.character(m6A_sites$start))),strand = as.character(m6A_sites$strand))

#####
DF1_start <- as.numeric(as.character(DF1binding_sites$start))
DF1_end <- as.numeric(as.character(DF1binding_sites$end))
DF1_center <- DF1_start+round((DF1_end-DF1_start+1)/2)
DF1_center_GR <- GRanges(seqnames = as.character(DF1binding_sites$seqname),
                         IRanges(start = DF1_center,width = 1),
                         strand = as.character(DF1binding_sites$strand))
DF1_range <- resize(DF1_center_GR, 101,fix="center")
m6AY1bind_overlap <- findOverlaps(m6Asites_GR,DF1_range,type = "within")
mapped_m6A_sites <- m6A_sites[unique(m6AY1bind_overlap @from),]
select_DF1_sites <- DF1bind_sites[unique(m6AY1bind_overlap @to),]
return(mapped_m6A_sites)
  
}



