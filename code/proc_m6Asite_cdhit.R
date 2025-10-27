load("./pos_neg_sites_DF1infor.Rdata")
load("./mapLTX_pos_neg_m6Asites.Rdata")
mapped_pos_sites <- all_m6A_sites$pos_sites
names(mapped_pos_sites) <- NULL
mapped_possites <- as.data.frame(mapped_pos_sites)
mapped_neg_sites <- all_m6A_sites$neg_sites
names(mapped_neg_sites) <- NULL
mapped_negsites <- as.data.frame(mapped_neg_sites)

pos_sites <- training_m6A_sites$pos_gene_m6Asites
neg_sites <- training_m6A_sites$neg_gene_m6Asites
select_sites <- function(mapped_sites,sites_infor){
  sites_start <- as.numeric(as.character(mapped_sites$start))
  new_sites <- data.frame()
  for (i in 1:length(sites_start)) {
    one_sites <- sites_infor[sites_infor$start==sites_start[i],]
    new_sites <- rbind(new_sites,one_sites)
  }
  #new_site <- data.frame(new_sites,sites_num=paste0("sites_",1:nrow(new_sites)))
  return(new_sites)
}
pos_new_sites <- select_sites(mapped_sites=mapped_possites,sites_infor=pos_sites)
pos_new_site <-  data.frame(pos_new_sites,sites_num=paste0("sites_",1:nrow(pos_new_sites)))

neg_new_sites <- select_sites(mapped_sites=mapped_negsites,sites_infor=neg_sites)
neg_new_site <-  data.frame(neg_new_sites,sites_num=paste0("sites_",(nrow(pos_new_site)+1):(nrow(pos_new_sites)+nrow(neg_new_sites))))

####load cd-hit process pos and neg fasta file
library(Biostrings)
pos_fasta_file <- "./pos.fasta"
pos_sequences <- readDNAStringSet(pos_fasta_file)
pos_seq_names <- names(pos_sequences)
##
neg_fasta_file <- "./neg.fasta"
neg_sequences <- readDNAStringSet(neg_fasta_file)
neg_seq_names <- names(neg_sequences)

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
process_lastsites <- list(pos_sites=pos_last_sites,neg_sites=neg_last_sites)
save(process_lastsites,file = "./proc_m6Asites_infor.Rdata")
