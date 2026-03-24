#DEGseq2
rm(list = ls())
#library(pheatmap)
library(tidyverse)
library(DESeq2)
#STEP-1
#Processing m0 and m1 samples
a<-read.csv("clinical m0 and m1.csv")# 16 letters column name
#colnames(a)<-toupper(names(a)) # how we have to change col names to upper
#colnames(a)<-substr(colnames(a), start = 1, stop = 16)
b<-read.csv("Cancer RSEM.csv") #never use row.names = 1 for merging files
colnames(b)<-substr(colnames(b), start = 1, stop = 16) 
any(is.double(names(a)))
any(is.double(names(b)))
any(which(colnames(b)%in% colnames(a)))
which(colnames(b)%in% colnames(a))
colnames(b) %in% colnames(a) %>% table
c<-merge(a,b,all.x = TRUE,all.y = TRUE)
d<-merge(a,b,all.x=F, all.y=T)
?merge
names(a)
names(b)
e<-merge(a,b, by = intersect(names(a), names(b)),all = T, all.x = F,all.y = F)
# combining both files by using common column names or patient IDs
write.csv(c, "Cancer RSEM.2.csv")
# manuall arrangement of Cancer RSEM.2.csv and reuplode
#STEP-2
x<-read.csv("Cancer RSEM.2.CSV", row.names = 1) # column names changed to m0 & m1
y<-x[str_detect(names(x), "M0")] # selection of m0
y[1,]
colnames(y)<-y[1,]#we are changing colnames with TCGA ID along with m0
colnames(y) <- paste(colnames(y), "m0", sep = "_")
View(y)
y<-y[-1,]
z<-x[str_detect(names(x), "M1")]#selection of m1
z[1,]
colnames(z)<-z[1,] #we are changing colnames with TCGA ID along with m1
colnames(z) <- paste(colnames(z), "m1", sep = "_")
View(z)
z<-z[-1,]
m0m1<-cbind(y,z)
write.csv(m0m1, "Cancer RSEM.3.csv")

#STEP-3
countData <- read.csv("Cancer RSEM.3.csv", row.names=1)
colData <- read.csv("design matrix.csv", row.names = 1)
# converting log2(counts+1) to counts
COU <- (2^countData)-1
write.csv(COU, "THCA counts.csv")
COU<-read.csv("ACC counts.csv",header=TRUE, row.names=1)
#removing >80% which are zeros
index0=floor(dim(countData)[2]*.2); index0 # Column dimension 
a<-countData[rowSums(countData==0)<=index0,]; dim(a)

dim(colData)
dim(countData)
#any(is.na(countData))
#countData<-na.omit(countData)
class(colData)
class(countData)
#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData)) # to match names
#countData <- countData[, rownames(colData)]

#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData ,colData , design= ~ Condition)
#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
res <- results(dds)
head(res)
class(res)
write.csv(res, "DEGs DEGseq2 cesc miR.csv")
res_sig <- subset(res, padj<.05) #Sort by adjusted p-value and display
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1)
head(res_lfc)
plotMA(res)
plotCounts(dds, gene=which.min(res$padj), intgroup="Stages.pvm") #Plotting individual genes
vsd <- vst(dds) # certain plots, we need to normalize our raw count data
install.packages("pheatmap")

pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=annot_col)
