###processing cdhit m6A sites
##raw m6A sites
library(Biostrings)
#seqdata <- readBStringSet("D:\\research\\m6Atranslation\\data\\single_base_level\\seq_out.fasta")
cdhit_proc_seq_sites <- function(seqdata,pos_sites_seq,neg_sites_seq,pos_sites_infor,neg_sites_infor){
  pos_peak_seq <- pos_sites_seq  
  neg_peak_seq <- neg_sites_seq 
  pos_seq <- pos_peak_seq[which(!is.na(match(names(pos_peak_seq),names(seqdata))))]
  neg_seq <- neg_peak_seq[which(!is.na(match(names(neg_peak_seq),names(seqdata))))]
  process_m6Asites_seq <- c(pos_seq,neg_seq)
  process_m6A_seq <- list(pos_m6A_seq=pos_seq,neg_m6A_seq=neg_seq)
  pos_seq_names <- names(process_m6A_seq$pos_m6A_seq)
  neg_seq_names <- names(process_m6A_seq$neg_m6A_seq)
  #load("D:\\research\\m6Atranslation\\data\\single_base_level\\pos_neg_sites_DF1infor.Rdata")
  m6A_reg_TE_sites <- pos_sites_infor
  non_m6A_reg_TE_sites <- neg_sites_infor
  pos_new_site <-  data.frame(m6A_reg_TE_sites,sites_num=paste0("sites_",1:nrow(m6A_reg_TE_sites)))
  neg_new_site <-  data.frame(non_m6A_reg_TE_sites,sites_num=paste0("sites_",(nrow(m6A_reg_TE_sites)+1):(nrow(m6A_reg_TE_sites)+nrow(non_m6A_reg_TE_sites))))
  
  last_sites <- function(seq_names,sites_infor){
    
    new_sites <- data.frame()
    for (i in 1:length(seq_names)) {
      one_sites <- sites_infor[sites_infor$sites_num==seq_names[i],]
      new_sites <- rbind(new_sites,one_sites)
    }
    #new_site <- data.frame(new_sites,sites_num=paste0("sites_",1:nrow(new_sites)))
    return(new_sites)
  }
  pos_last_sites <- last_sites(seq_names=pos_seq_names,sites_infor=pos_new_site)
  neg_last_sites <- last_sites(seq_names=neg_seq_names,sites_infor=neg_new_site)
  process_seq_siteinfor <- list(pos_site_seqs=pos_m6A_seq,pos_sites=pos_last_sites,
                            neg_site_seqs=neg_m6A_seq,neg_sites=neg_last_sites)
  return(process_seq_siteinfor)
}

