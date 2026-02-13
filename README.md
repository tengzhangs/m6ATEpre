# m6ATEpre
## Introduction
The most prevalent mRNA modification, N6-methyladenosine (m6A) plays an important role in various RNA metabolism, including gene expression and translation. However, the selective mechanism by which m6A sites regulate mRNA translation through YTHDF1 binding remains poorly understood, due to a lack of computational methods for identifying context-specific m6A sites that regulate translation. To address this, we developed a novel computational framework named m6ATEpre, the first tool designed to predict cell-specific m6A sites with regulating translation efficiency (m6A-reg-TE). m6ATEpre integrates multi-omics data, including MeRIP-seq data, PAR-CLIP data and Ribo-seq data, introduces a novel feature representation strategy for m6A site sequences, and employs an autoencoder to effectively capture embedded feature representations.
# The step-by-step usage example
## Peak calling for MeRIP-seq data
For MeRIP-seq data, we first used the exomPeak2 software to detect the m6A peak region.
```r
library(exomePeak2)
f1 <- "./CTL_IP1.bam"
f2 <- "./CTL_IP2.bam"
f3 <- "./CTL_Input1.bam"
f4 <- "./CTL_Input2.bam"

IP_BAM <- c(f1,f2)
INPUT_BAM <- c(f3,f4)

GENE_ANNO_GTF = "/home/disk1/zhangteng/hg19_GTF/genes.gtf"
exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff_dir = GENE_ANNO_GTF,
           genome = "hg19",
           paired_end = FALSE)
```
### Mapping m6A peak to transcriptome 
```r
m6Apeak_CTL <- read.csv(f1,header = T)
fa <- "./Mod.bed"
all_m6Apeakbed <-  read.table(fa,sep="\t",header=FALSE,stringsAsFactors =FALSE)
all_m6Apeakbed$V4 <- as.character(m6Apeak_CTL$geneID)
all_m6Apeakbed$peak_ID <- as.character(m6Apeak_CTL$name)
colnames(all_m6Apeakbed) <- c("chr",         "start" ,      "end",         "name",        "score",       "strand",      "thickStart",  "thickEnd",    "itemRgb",    
                              "blockCount",  "blockSizes",  "blockStarts","peak_ID")

GTF_file <- "./hg19_gtf/genes.gtf"
library(GenomicFeatures)
obtain_peak_exon <- map_peak_TX(peak_sites_infor=all_m6Apeakbed,annotation_file=GTF_file)
```
## Obtain single-base m6A sites
### Obtain m6A peaks' sequences with DRACH motif
```r
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
genome <- BSgenome.Hsapiens.UCSC.hg19
getpeaks_seq <- get_peak_seq(peaks = obtain_peak_exon[[1]])
writeXStringSet(getpeaks_seq,"./motif_peak_seq.fa")
```
### Obtain m6A sites in single-base resulation
```
conda activate R3.6
cd ./sramp_simple/
nohup perl runsramp.pl ./motif_peak_seq.fa ./singlebase_m6Asites.txt mature &
```
```r
#Select the high confidience m6A sites in single-base resulation
f1 <- "./singlebase_m6Asites.txt"
Human_m6Asites <- read.delim2(f1)
HeLa_m6Asites <- Human_m6Asites[-grep("Non-m6A site",Human_m6Asites$Classification),]
high_m6Asites <- HeLa_m6Asites[grep('Very high confidence',HeLa_m6Asites$Classification),]
#Get the single-base m6A sites information (including the seqnames, start sites and strand)
load("./peak_exon_region.Rdata")
m6Asites_infor <- get_m6Asites_infor(high_m6Asites=high_m6Asites,peak_infor=obtain_peak_exon)
```
## Integrate Multi-omics data to select potential m6A sites with regulating translation efficiency (m6A-reg-TE)
```r
#Mapping single-base m6A sites to YTHDF1 binding region
f1 <- "./DF1_bind_gene_site.csv"
DF1binding_sites <- read.csv(f1)
m6Asites2DF1 <- m6Asites2_DF1binding_region(m6A_sites=m6Asites_infor,DF1binding_sites=DF1binding_sites)
#Select potential m6A-reg-TE sites from single-base m6A sites with DF1 binding
f2 <- "D:\\research\\m6Atranslation\\data\\Y1_TE_change.xlsx"
Y1TE_change <- readxl::read_xlsx(f2)
Y1TE_down <- Y1TE_change[which(Y1TE_change$`Translation efficiency`<(-0.5)),]
Y1KD_TEup <- Y1TE_change[which(Y1TE_change$`Translation efficiency`>(-0.5)),]
m6A_reg_TE_sites <- m6Asites2DF1[which(!is.na(match(select_m6A_sites$gene_name,Y1TE_down$`Gene symbol`))),]
non_m6A_reg_TE_sites <- m6Asites2DF1[which(!is.na(match(select_m6A_sites$gene_name,Y1KD_TEup$`Gene symbol`))),]
```
## Obtain the sequences both for m6A-reg-TE sites and non-m6A-reg-TE sites
```r
library(m6ALogisticModel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(GenomicFeatures)
annotation_file <- "D:\\hg19_GTF\\genes.gtf"
m6A_reg_TE_sites_seq <- get_m6A_seq(target_sites=m6A_reg_TE_sites,annotation_file)
non_m6A_reg_TE_sites_seq <- get_m6A_seq(target_sites=non_m6A_reg_TE_sites,annotation_file)
pos_m6Asites_seqs <- as.character(m6A_reg_TE_sites_seq)
names(pos_m6Asites_seqs) <-  paste0('site_',1:length(pos_m6Asites_seqs))
neg_m6Asites_seqs <- as.character(non_m6A_reg_TE_sites_seq)
names(neg_m6Asites_seqs) <- paste0("site_",(length(pos_m6Asites_seqs)+1):(length(pos_m6Asites_seqs)+length(neg_m6Asites_seqs)))
m6Asites_seq <- c(neg_m6Asites_seqs,neg_m6Asites_seqs)
writeXStringSet(m6Asites_seq,"./allm6A_seqs.fasta")
##CD-Hit processing for m6A sites sequences
cd cdhit-4.8.1/
./cd-hit -i ./allm6A_seqs.fasta -o ./pos_neg_samples/seq-out.fasta -c 0.9 -T 4 -n 4
```r
##select m6A sites sequence after CD-Hit processing
library(Biostrings)
seqdata <- readBStringSet("./pos_neg_samples/seq_out.fasta")
pos_seq <- m6A_reg_TE_sites_seq[which(!is.na(match(names(pos_peak_seq),names(seqdata))))]
neg_seq <- non_m6A_reg_TE_sites_seq[which(!is.na(match(names(neg_peak_seq),names(seqdata))))]
process_m6Asites_seq <- c(pos_seq,neg_seq)
process_m6A_seq <- list(pos_m6A_seq=pos_seq,neg_m6A_seq=neg_seq)
```
