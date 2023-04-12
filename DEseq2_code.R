options(stringsAsFactors = F)
library(DESeq2)
library(limma)
library(ggpubr)
library(pheatmap)
setwd("E:\\Git_project\\exo_code")
#circRNA
circRNA <- read.table("circRNA_expression.txt",sep = "\t",header = T,check.names = F)
for (i in 1:nrow(circRNA)) {circRNA$symbol[i] <- strsplit(circRNA$symbol,"[|]")[[i]][1]
}
circRNA <- circRNA[,-1]
circRNA <- circRNA[!duplicated(circRNA$symbol),]
rownames(circRNA) <- circRNA$symbol
circRNA <- circRNA[,-1]
cluster <- as.data.frame(matrix(data = NA,ncol = 2,nrow = 15))
colnames(cluster) <- c("ID","Cluster")
cluster$ID <- colnames(circRNA)[1:15]
cluster$Cluster <- c(rep("Tumor",10),rep("Normal",5))
cluster <- as.data.frame(cluster)


group_list=as.factor(cluster$Cluster)
colData <- data.frame(row.names=cluster$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(circRNA),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Tumor","Normal"))
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "Tumor_vs_Normal_circRNA.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
##Establish group information
resOrdered$Group <-"Not-significant"
##Establish screening criteria
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(2)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-2)))] <-"Down-regulated"
table(resOrdered$Group)
#Create a column
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
#plot
p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(Padj-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-2,2), linetype = "dashed")
ggsave("Tumor_vs_Normal_circRNA.png",p)
ggsave("Tumor_vs_Normal_circRNA.pdf",p)
ggsave("Tumor_vs_Normal_circRNA.eps",p)
sig_gene <- resOrdered[resOrdered$Group == "Down-regulated"|resOrdered$Group == "Up-regulated",]
p <- pheatmap(t(scale(t(circRNA[rownames(sig_gene),]))),
              color = colorRampPalette(c("#2f5688","white","#CC0000"))(50),
              fontsize_col = 15,border_color = "white",
              cluster_cols = F,cluster_rows = T,show_rownames =F)
ggsave("circ_heatmap.png",p,width = 6,height = 8)
ggsave("circ_heatmap.pdf",p,width = 6,height = 8)
ggsave("circ_heatmap.eps",p,width = 6,height = 8)
write.csv(cicrRNA,file = "allcicRNA_marix.csv")
write.csv(sig_gene,file = "sigcicRNA_marix.csv")
#LNCRNA
lncRNA <- read.table("lncRNA_counts.txt",sep = "\t",header = T,check.names = F)
lncRNA <- lncRNA[,-1]
lncRNA <- lncRNA[!duplicated(lncRNA$Symbol),]
rownames(lncRNA) <- lncRNA$Symbol
lncRNA <- lncRNA[,-1]
cluster <- as.data.frame(matrix(data = NA,ncol = 2,nrow = 15))
colnames(cluster) <- c("ID","Cluster")
cluster$ID <- colnames(lncRNA)[1:15]
cluster$Cluster <- c(rep("Tumor",10),rep("Normal",5))
cluster <- as.data.frame(cluster)


group_list=as.factor(cluster$Cluster)
colData <- data.frame(row.names=cluster$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(lncRNA),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Tumor","Normal"))
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "Tumor_vs_Normal_lncRNA.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
#Establish group information
resOrdered$Group <-"Not-significant"
#Establish screening criteria
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(2)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-2)))] <-"Down-regulated"
table(resOrdered$Group)
#Create a column
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
#plot
p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(Padj-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-2,2), linetype = "dashed")
ggsave("Tumor_vs_Normal_lncRNA.png",p)
ggsave("Tumor_vs_Normal_lncRNA.pdf",p)
ggsave("Tumor_vs_Normal_lncRNA.eps",p)
sig_gene <- resOrdered[resOrdered$Group == "Down-regulated"|resOrdered$Group == "Up-regulated",]
p <- pheatmap(t(scale(t(lncRNA[rownames(sig_gene),]))),
              color = colorRampPalette(c("#2f5688","white","#CC0000"))(50),
              fontsize_col = 15,border_color = "white",
              cluster_cols = F,cluster_rows = T,show_rownames =F)
ggsave("lnc_heatmap.png",p,width = 6,height = 8)
ggsave("lnc_heatmap.pdf",p,width = 6,height = 8)
ggsave("lnc_heatmap.eps",p,width = 6,height = 8)
write.csv(lncRNA,file = "alllncRNA_marix.csv")
write.csv(sig_gene,file = "siglncRNA_marix.csv")

#mRNA
mRNA <- read.table("mRNA_counts.txt",sep = "\t",header = T,check.names = F)
mRNA <- mRNA[,-1]
mRNA <- mRNA[!duplicated(mRNA$Symbol),]
rownames(mRNA) <- mRNA$Symbol
mRNA <- mRNA[,-1]
cluster <- as.data.frame(matrix(data = NA,ncol = 2,nrow = 15))
colnames(cluster) <- c("ID","Cluster")
cluster$ID <- colnames(mRNA)[1:15]
cluster$Cluster <- c(rep("Tumor",10),rep("Normal",5))
cluster <- as.data.frame(cluster)


group_list=as.factor(cluster$Cluster)
colData <- data.frame(row.names=cluster$ID, group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = ceiling(mRNA),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  
resultsNames(dds2)
res <-  results(dds2, contrast=c("group_list","Tumor","Normal"))
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered = na.omit(resOrdered)
write.csv(resOrdered,file = "Tumor_vs_Normal_mRNA.csv")

head(resOrdered)
resOrdered$logP <- -log10(resOrdered$padj)
#Establish group information
resOrdered$Group <-"Not-significant"
#Establish screening criteria
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange >(2)))] <-"Up-regulated"
resOrdered$Group[which((resOrdered$padj <(0.05))&(resOrdered$log2FoldChange <(-2)))] <-"Down-regulated"
table(resOrdered$Group)
#Create a column
resOrdered$Label = ""
resOrdered <- resOrdered[order(resOrdered$padj),]
resOrdered$symbol <- rownames(resOrdered)

up_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Up-regulated")],10)
down_gene <- head(resOrdered$symbol[which(resOrdered$Group == "Down-regulated")],10)
resOrdered_top10gene <- c(as.character(up_gene),as.character(down_gene))
resOrdered$Label[match(resOrdered_top10gene,resOrdered$symbol)] <- resOrdered_top10gene
#plot
p <-ggscatter(resOrdered,x = "log2FoldChange",y = "logP",
              color = "Group",
              palette = c("#2f5688","#BBBBBB","#CC0000"),
              size = 1,
              label = resOrdered$Label,
              font.label = 8,
              repel = T,
              xlab = "log2FoldChange",
              ylab = "-log10(Padj-value)")  + geom_hline(yintercept = 1.3, linetype = "dashed")+
  geom_vline(xintercept = c(-2,2), linetype = "dashed")
ggsave("Tumor_vs_Normal_mRNA.png",p)
ggsave("Tumor_vs_Normal_mRNA.pdf",p)
ggsave("Tumor_vs_Normal_mRNA.eps",p)
sig_gene <- resOrdered[resOrdered$Group == "Down-regulated"|resOrdered$Group == "Up-regulated",]
p <- pheatmap(t(scale(t(mRNA[rownames(sig_gene),]))),
              color = colorRampPalette(c("#2f5688","white","#CC0000"))(50),
              fontsize_col = 15,border_color = "white",
              cluster_cols = F,cluster_rows = T,show_rownames =F)
ggsave("m_heatmap.png",p,width = 6,height = 8)
ggsave("m_heatmap.pdf",p,width = 6,height = 8)
ggsave("m_heatmap.eps",p,width = 6,height = 8)
write.csv(mRNA,file = "allmRNA_marix.csv")
write.csv(sig_gene,file = "sigmRNA_marix.csv")
