setwd("/data/srijith/cohen_data")
list.files()
col_data <- read.csv("count_rsf_13.csv")
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
corRes <- coreeres[abs(coreeres$value) > 0.9 , ]
write.csv(corRes, "edgelist_RSF1_allgenes.csv")
