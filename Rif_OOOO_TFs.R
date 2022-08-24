library(readr)
setwd("D:/RNASeq_data_analysis/Cohen_data/cohen_data/PCIT_analysis_allgenes")
data <- read.csv("key_rif_analysis.csv")
rownames(data) <- data$ENSEMBL
setwd("D:/RNASeq_data_analysis/Cohen_data/cohen_data/Counts_final/vst_counts")
col_data <- read.csv("count_1_16.csv")
rownames(col_data) <- col_data$X
col_data <- col_data[-1]
col_data <- as.matrix.data.frame(col_data)
colnames(col_data) <- c("a1","a2","a3", "a4", "w1","w2","w3","w4")
colnames(col_data) <- paste0(colnames(col_data), c(rep("_AR", 4), rep("_WT", 4)))
tf_data <- read_csv("tf_list.csv")
sig_genes <- as.data.frame(rownames(data))
colnames(sig_genes) <- "orf"
tf_data_genes <- as.data.frame(tf_data$ORF)
colnames(tf_data_genes) <- "orf"
library(dplyr)
all_genes <- rbind(sig_genes,tf_data_genes)
library(dplyr)
all_genes <- distinct(all_genes, orf, .keep_all = TRUE)
rownames(all_genes) <- all_genes$orf
counts <- col_data[rownames(all_genes), ]
library(CeTF)
RIF_input <- as.data.frame(counts)

#identify the crucial transcription factors

library(dplyr)
TFs <- as.character(tf_data_genes$orf)
Target <- sig_genes$orf
RIF_input_1 <- RIF_input[c(Target, TFs), ]
# Performing RIF analysis
RIF_out <- RIF(input = RIF_input_1,
               nta = length(Target),
               ntf = length(TFs),
               nSamples1 = 4,
               nSamples2 = 4)

head(RIF_out)

library(org.Sc.sgd.db)
anno_1 <- AnnotationDbi::select(org.Sc.sgd.db,keys= RIF_out$TF,
                                columns= "GENENAME",
                                keytype="ENSEMBL")
RIF_out$Gene <- anno_1

#select the DE-genes and crucial transcription factors for network construction

rif1_tfs <- subset(RIF_out, RIF1 < -1.96)
rif1_tfs_1 <- subset(RIF_out, RIF1 > 1.96)
rif2_tfs <- subset(RIF_out, RIF2 < -1.96)
rif2_tfs_1 <- subset(RIF_out, RIF2 > 1.96)
rif_tfs_all <- rbind(rif1_tfs,rif1_tfs_1, rif2_tfs, rif2_tfs_1)
setwd("D:/RNASeq_data_analysis/Cohen_data/cohen_data/rif_analysis")
write.csv(rif_tfs_all, "rif_wildtype_OOOO.csv")



#PCIT Network
setwd("/data/srijith/cohen_data")
list.files()
col_data <- read.csv("count_rme1_rsf_11.csv")
rownames(col_data) <- col_data$X
col_data <- col_data[-1]
col_data <- as.matrix.data.frame(col_data)
colnames(col_data) <- c("a1","a2","a3", "a4", "w1","w2","w3","w4")
colnames(col_data) <- paste0(colnames(col_data), c(rep("_AR", 4), rep("_WT", 4)))
library(PCIT)
PCIT_input <- as.data.frame(col_data)
PCIT_input_WT <- PCIT_input[,grep("_WT", colnames(PCIT_input))]
PCIT_input_AR <- as.matrix(PCIT_input[,grep("_AR", colnames(PCIT_input))])
c <- cor(t(PCIT_input_AR))
system.time(result <- pcit(c))
signif <- idx(result) 
plotCorCoeff(c, list("PCIT Meaningful" = 
                       + signif), col=c("red")) 
nonsignif <- idxInvert(nrow(c), signif)
c[nonsignif] <- 0 
adj <- abs(c)
library(reshape2)
adj_mat <- as.matrix(adj)
coreeres <- melt(adj_mat)
corRes <- coreeres[abs(coreeres$value) > 0.95 , ]
write.csv(corRes, "edgelist_RSF1_RME1_allgenes.csv")

