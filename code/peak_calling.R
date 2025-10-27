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
