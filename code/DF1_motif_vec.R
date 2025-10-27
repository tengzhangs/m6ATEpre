load( "./proc_m6Asites_seqs.Rdata")
pos_m6A_seq <- process_m6A_seq$pos_m6A_seq
neg_m6A_seq <- process_m6A_seq$neg_m6A_seq
get_DF1_motif_infor <- function(peak_seq){
  DF1_motif <- c("GGAC","GAAC")
  DF1_loc <- vmatchPattern(DF1_motif[1], peak_seq)
  
  DF1_num <- vcountPattern(DF1_motif[1], peak_seq)
  DFmotif_start <- start(DF1_loc)
  
  for(i in 2:length(DF1_motif)) {
    DF1_loc0 <- vmatchPattern(DF1_motif[i], peak_seq)
    DF1_num0 <- vcountPattern(DF1_motif[i], peak_seq)
    DFmotif_start0 <- start(DF1_loc0)
    DFmotif_start <- mapply(c, DFmotif_start, DFmotif_start0, SIMPLIFY=FALSE)
    DF1_num <- DF1_num+DF1_num0
  }
  DF1motif_infor <- list(DF1_num=DF1_num,DF1motif_start=DFmotif_start)
  return(DF1motif_infor)
  
}

pos_m6A_DF1 <- get_DF1_motif_infor(peak_seq=pos_m6A_seq)
neg_m6A_DF1 <- get_DF1_motif_infor(peak_seq=neg_m6A_seq)
DF1_motif_vec <- function(DF1_start){
  DF1_motif_vec <- data.frame()
  for (i in 1:length(DF1_start)) {
    oneDF1vec <- vector(length = 501-4+1)
    if(sum(unique( DF1_start[[i]]))>0){
      oneDF1_motif <- DF1_start[[i]][order(DF1_start[[i]],decreasing = F)]
      oneDF1vec[oneDF1_motif] <- 1
    }
    if(sum(unique( DF1_start[[i]]))==0){
      oneDF1vec <- rep(0,498)
    }
    DF1_motif_vec  <- rbind(DF1_motif_vec,oneDF1vec)
  }
  colnames(DF1_motif_vec) <- NULL
  return(DF1_motif_vec)
}
pos_DF1motif_vec <- DF1_motif_vec(DF1_start=pos_m6A_DF1$DF1motif_start)
neg_DF1motif_vec <- DF1_motif_vec(DF1_start=neg_m6A_DF1$DF1motif_start)

write.csv(pos_DF1motif_vec,file = "./single_base_level/cdhit_proces/pos_motif_vec.csv",row.names = F)

write.csv(neg_DF1motif_vec,file = "./single_base_level/cdhit_proces/neg_motif_vec.csv",row.names = F)
