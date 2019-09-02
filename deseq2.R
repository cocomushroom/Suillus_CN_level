library("tximport")
library("readr")
library("tximportData")
library("DESeq2")

dir <- getwd()

samples <- read.table("samples.txt", header=TRUE)

#samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$ID
#samples[,c("pop","center","run","condition")]

files <- file.path(dir,"salmon", samples$ID, "quant.sf")
names(files) <- samples$ID
tx2gene <- read_csv(file.path(dir, "tID"))

txi <- tximport(files, type="salmon",tx2gene=tx2gene)


ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

det<-DESeq(ddsTxi)

res_det<-results(det)
summary(res_det)
