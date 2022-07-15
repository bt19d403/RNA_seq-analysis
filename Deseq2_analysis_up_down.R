library(readr)
library(sva)
setwd("D:/RNASeq_data_analysis/Cohen_data/cohen_data/Counts_final")
data <- read.csv("rawcounts_grouped.csv")
row.names(data) <- data$ORF
colnames(data)
count_data <- data[,-1]
count_matrix <- as.matrix(count_data)
meta_data <- read.csv("sample_info.csv")
row.names(meta_data) = meta_data$Sample
meta_data

all(rownames(meta_data) == colnames(count_data))
#batch correction using CombatSeq
genotype <- meta_data$condition
batch <- meta_data$Replicate

# include condition (group variable)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=genotype, full_mod=TRUE)
write.csv( adjusted_counts, "bc_grouped_4replicates.csv")
raw_counts <- adjusted_counts
#write.csv( raw_counts, "bc_grouped_3replicates.csv")
library(DESeq2)
meta_data$Genotype = factor(meta_data$Genotype)
meta_data$Replicate = factor(meta_data$Replicate)

dds = DESeqDataSetFromMatrix(countData = raw_counts ,
                             colData = meta_data,
                             design = ~ Genotype)
levels(dds$Genotype)

dds$Genotype <- relevel(dds$Genotype, ref = "RME1WRSF1WIME1cWIME1ncW")


dds = DESeq(dds)
samplelist <- resultsNames(dds)
#prefiltering
#filtering the counts less than 100   
#which results in 6305 genes out of 6571 genes
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]
levels(dds$Genotype)
dds_norm <- vst(dds)
vst_dds <- assay(dds_norm)
library(dplyr)
#normalized_counts <- assay(dds_norm) %>%
 # t() # Transpose this data
#write.csv(normalized_counts, "vst_transformed_wgcna.csv")
normalized_counts <- counts(dds, normalized=TRUE)
#write.csv(normalized_counts , "normalized_bc_counts_3replicates.csv")
#normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts , "normalized_bc_counts_1.csv")
write.csv(vst_dds , "vst_bc_counts.csv")
library(dplyr)
library(tibble)
library(org.Sc.sgd.db)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")
for (j in as.character(samplelist[2:16])) {
  caseB <- j
  res <- lfcShrink(dds, coef= caseB, type="apeglm")
  results.ordered <- as.data.frame(res) %>% 
    rownames_to_column("ENSEMBL") %>% arrange(padj)
  anno_bc <- AnnotationDbi::select(org.Sc.sgd.db,keys= results.ordered$ENSEMBL,
                                   columns= "GENENAME",
                                   keytype="ENSEMBL")
  dup_ids <- anno_bc$ENSEMBL[duplicated(anno_bc$ENSEMBL)]
  filter(anno_bc, ENSEMBL %in% dup_ids) %>% arrange(ENSEMBL)
  anno_bc <- AnnotationDbi::select(org.Sc.sgd.db,keys=results.ordered$ENSEMBL,
                                   columns=c("ENSEMBL","GENENAME","ENTREZID"),
                                   keytype="ENSEMBL") %>% 
    filter(!duplicated(ENSEMBL))
  results.annotated <- left_join(results.ordered, anno_bc,by="ENSEMBL")
  assign(paste0("DEGs_",caseB), subset(results.annotated, padj < 0.1, select = c("GENENAME","padj","log2FoldChange","ENSEMBL","pvalue")))
  assign(paste0("DEGs_up",caseB), subset(results.annotated, padj < 0.1 & log2FoldChange >0.5))
  assign(paste0("DEGs_down",caseB), subset(results.annotated, padj < 0.1 & log2FoldChange < -0.5))
}
DE_list <- list( RME1_RSF1_IME1c_IME1nc = DEGs_Genotype_RME1ORSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1_RSF1_IME1c = DEGs_Genotype_RME1ORSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1_RSF1_IME1nc = DEGs_Genotype_RME1ORSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1_RSF1 = DEGs_Genotype_RME1ORSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1_IME1c_IME1nc = DEGs_Genotype_RME1ORSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1_IME1c = DEGs_Genotype_RME1ORSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1_IME1nc = DEGs_Genotype_RME1ORSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RME1 = DEGs_Genotype_RME1ORSF1WIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                 IME1nc = DEGs_Genotype_RME1WRSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                 IME1c = DEGs_Genotype_RME1WRSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                 IME1c_IME1nc = DEGs_Genotype_RME1WRSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW  ,
                 RSF1 = DEGs_Genotype_RME1WRSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW  ,
                 RSF1_IME1nc = DEGs_Genotype_RME1WRSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RSF1_IME1c = DEGs_Genotype_RME1WRSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                 RSF1_IME1c_IME1nc =  DEGs_Genotype_RME1WRSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW)


options(max.print = 1000000)
capture.output(DE_list, file = "Deseq2_DEGenes_4replicates_LFCshinkage.txt")  

Deseq2_dat <- as.data.frame(read.table("Deseq2_DEGenes_4replicates_LFCshinkage.txt", sep="\t",header=TRUE))
Deseq2_dat
options(max.print = 1000000)
write.csv(Deseq2_dat, "DEG_DEseq2_4replicates_LFCshink_apelgm.csv")



DE_length_up <- list( RME1_RSF1_IME1c_IME1nc = length(DEGs_upGenotype_RME1ORSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1_RSF1_IME1c = length(DEGs_upGenotype_RME1ORSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1_RSF1_IME1nc = length(DEGs_upGenotype_RME1ORSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1_RSF1 = length(DEGs_upGenotype_RME1ORSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1_IME1c_IME1nc = length(DEGs_upGenotype_RME1ORSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1_IME1c = length(DEGs_upGenotype_RME1ORSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1_IME1nc = length(DEGs_upGenotype_RME1ORSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RME1 = length(DEGs_upGenotype_RME1ORSF1WIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   IME1nc = length(DEGs_upGenotype_RME1WRSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   IME1c = length(DEGs_upGenotype_RME1WRSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   IME1c_IME1nc = length(DEGs_upGenotype_RME1WRSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL)  ,
                   RSF1 = length(DEGs_upGenotype_RME1WRSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL)  ,
                   RSF1_IME1nc = length(DEGs_upGenotype_RME1WRSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RSF1_IME1c = length(DEGs_upGenotype_RME1WRSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                   RSF1_IME1c_IME1nc =  length(DEGs_upGenotype_RME1WRSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL))


DE_length_down <- list( RME1_RSF1_IME1c_IME1nc = length(DEGs_downGenotype_RME1ORSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1_RSF1_IME1c = length(DEGs_downGenotype_RME1ORSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1_RSF1_IME1nc = length(DEGs_downGenotype_RME1ORSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1_RSF1 = length(DEGs_downGenotype_RME1ORSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1_IME1c_IME1nc = length(DEGs_downGenotype_RME1ORSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1_IME1c = length(DEGs_downGenotype_RME1ORSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1_IME1nc = length(DEGs_downGenotype_RME1ORSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RME1 = length(DEGs_downGenotype_RME1ORSF1WIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      IME1nc = length(DEGs_downGenotype_RME1WRSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      IME1c = length(DEGs_downGenotype_RME1WRSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      IME1c_IME1nc = length(DEGs_downGenotype_RME1WRSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL)  ,
                      RSF1 = length(DEGs_downGenotype_RME1WRSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL)  ,
                      RSF1_IME1nc = length(DEGs_downGenotype_RME1WRSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RSF1_IME1c = length(DEGs_downGenotype_RME1WRSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL) ,
                      RSF1_IME1c_IME1nc =  length(DEGs_downGenotype_RME1WRSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW$ENSEMBL))

DE_length_up_new <- as.data.frame(DE_length_up)
DE_length_up_new <- t(DE_length_up_new)

DE_length_down_new <- as.data.frame(DE_length_down)
DE_length_down_new <- t(DE_length_down_new)

De_length_all <- cbind(DE_length_up_new, DE_length_down_new)
De_length_all <- as.data.frame(De_length_all)
De_length_all$Genotype <- rownames(De_length_all)

DE_length_all_new <- De_length_all %>% mutate( total_DEGs = V1+V2 )
colnames(DE_length_all_new) <- c("upregulated_genes" , "downregulated_genes","Genotype", "total_DEGs")
library(tidyr)
library(ggplot2)

library(dplyr)
DE_length_all_new %>% 
  dplyr::select(-total_DEGs) %>% 
  gather(type, count, upregulated_genes:downregulated_genes) %>% 
  ggplot(., aes(x=Genotype, y=count, fill=forcats::fct_rev(type))) +
  geom_bar(stat="identity") +coord_flip() + theme_minimal() +theme(legend.position = "top")

#plot upset

for (j in as.character(samplelist[2:16])) {
  caseB <- j
  res <- lfcShrink(dds, coef= caseB, type="apeglm")
  res$GENEID = row.names(res)
  dup_geneid <- res$GENEID[duplicated(res$GENEID)]
  filter(data.frame(res), GENEID %in% dup_geneid)
  assign(paste0("DEGs_name",caseB),rownames(subset(data.frame(res),padj <0.1)))
  resultsNames(dds)
}


DE_length <- list( RME1_RSF1_IME1c_IME1nc = length(DEGs_nameGenotype_RME1ORSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1_RSF1_IME1c = length(DEGs_nameGenotype_RME1ORSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1_RSF1_IME1nc = length(DEGs_nameGenotype_RME1ORSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1_RSF1 = length(DEGs_nameGenotype_RME1ORSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1_IME1c_IME1nc = length(DEGs_nameGenotype_RME1ORSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1_IME1c = length(DEGs_nameGenotype_RME1ORSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1_IME1nc = length(DEGs_nameGenotype_RME1ORSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RME1 = length(DEGs_nameGenotype_RME1ORSF1WIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   IME1nc = length(DEGs_nameGenotype_RME1WRSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   IME1c = length(DEGs_nameGenotype_RME1WRSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   IME1c_IME1nc = length(DEGs_nameGenotype_RME1WRSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW)  ,
                   RSF1 = length(DEGs_nameGenotype_RME1WRSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW)  ,
                   RSF1_IME1nc = length(DEGs_nameGenotype_RME1WRSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RSF1_IME1c = length(DEGs_nameGenotype_RME1WRSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW) ,
                   RSF1_IME1c_IME1nc =  length(DEGs_nameGenotype_RME1WRSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW))

De_length <- as.data.frame(DE_length)
row.names(De_length) = "Number of DEGs"
De_length_new <- t(De_length)
De_length_new <- as.data.frame(De_length_new)
De_length_new$Genotype <- rownames(De_length_new)
library(ggpubr)
ggbarplot(De_length_new, x = "Genotype", y = "Number of DEGs", ylab = "Number of DEGs in comparison with Wildtype",
          label.pos = "out",  orientation = "horiz", fill = "steelblue", label = F,
          lab.col = "black",
          lab.size = 3,
          lab.pos = "out")


DE_list <- list( RME1_RSF1_IME1c_IME1nc = DEGs_nameGenotype_RME1ORSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW,
                   RME1_RSF1_IME1c = DEGs_nameGenotype_RME1ORSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RME1_RSF1_IME1nc = DEGs_nameGenotype_RME1ORSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RME1_RSF1 = DEGs_nameGenotype_RME1ORSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RME1_IME1c_IME1nc = DEGs_nameGenotype_RME1ORSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RME1_IME1c = DEGs_nameGenotype_RME1ORSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RME1_IME1nc = DEGs_nameGenotype_RME1ORSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RME1 = DEGs_nameGenotype_RME1ORSF1WIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                   IME1nc = DEGs_nameGenotype_RME1WRSF1WIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                   IME1c = DEGs_nameGenotype_RME1WRSF1WIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                   IME1c_IME1nc = DEGs_nameGenotype_RME1WRSF1WIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW  ,
                   RSF1 = DEGs_nameGenotype_RME1WRSF1OIME1cWIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW  ,
                   RSF1_IME1nc = DEGs_nameGenotype_RME1WRSF1OIME1cWIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RSF1_IME1c = DEGs_nameGenotype_RME1WRSF1OIME1cOIME1ncW_vs_RME1WRSF1WIME1cWIME1ncW ,
                   RSF1_IME1c_IME1nc =  DEGs_nameGenotype_RME1WRSF1OIME1cOIME1ncO_vs_RME1WRSF1WIME1cWIME1ncW)

library(UpSetR)
DE_gns <- UpSetR::fromList(DE_list)
UpSetR::upset(DE_gns)


