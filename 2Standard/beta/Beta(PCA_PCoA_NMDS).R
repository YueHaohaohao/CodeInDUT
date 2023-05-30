
library(vegan)

otu_data <- t(read.delim("AllOTU.tsv",row.names = 1,header = T, sep = "\t"))
otu_data <- otu_data[, colSums(otu_data != 0) > 0]
group_info <- read.table("Group.tsv",header = T,row.names = 1,sep = "\t")
group_info <- group_info[match(rownames(otu_data),rownames(group_info)),]

#PERMANOVA
pmnv_result <- adonis2(otu_data ~ group_info,permutations = 999, method="bray")
# write.table("\t", "PERMANOVA.tsv", append = F, quote = F, eol = "", row.names = F, col.names = F)
# write.table(PERMANOVA,"PERMANOVA.tsv",sep = "\t",row.names = F)

#多重比较
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairadonis_result <- pairwise.adonis(x=otu_data, factors=group_info)

#距离矩阵
otu_dist <- vegdist(otu_data, method="bray", binary=F)

#PCoA
otu_pcoa <- cmdscale(otu_dist, k=2, eig=T)
otu_pcoa_points <- as.data.frame(otu_pcoa$points)
colnames(otu_pcoa_points) <- c('PC1','PC2')
otu_pcoa_points$Group <- group_info
sum_eig <- sum(otu_pcoa$eig)
eig_percent <- (otu_pcoa$eig/sum_eig*100)[1:2]

#图
library(ggplot2)
library(ggforce)
library(ggrepel)
library(ggsci)

PCoA_plot <- ggplot(otu_pcoa_points, aes(x=PC1, y=PC2, color=Group, group = Group)) +
  labs(x=paste("PC1 (", round(eig_percent[1],2), "%)", sep=""),
       y=paste("PC2 (", round(eig_percent[2],2), "%)", sep="")) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none')+
  annotate(geom= "text",x=0,y=0.5, label = paste("permanova: R^2=" ,round(pmnv_result[1,'R2'],2),";","p=",round(pmnv_result[1,'Pr(>F)'],3)) )+
  geom_text_repel(aes(label=rownames(otu_pcoa_points)),size=3,max.overlaps = 20)
PCoA_plot
ggsave('PCoA.svg',PCoA_plot,width = 5,height = 5,units = 'in')

#NMDS
otu_nmds <- metaMDS(otu_dist, k = 2)
otu_nmds_points <- data.frame(otu_nmds$points)
otu_nmds_points$Group <- group_info
NMDS_stress <- round(otu_nmds$stress,2)

NMDS_plot <- ggplot(otu_nmds_points, aes(x=MDS1, y=MDS2, color=Group, group = Group)) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank())+
  geom_text_repel(aes(label=rownames(otu_nmds_points)),size=3,max.overlaps = 20)+
  annotate(geom= "text",x = -1, y = 1, label = paste("Stress=" ,NMDS_stress))
NMDS_plot
ggsave('NMDS.svg',NMDS_plot,width = 5,height = 5,units = 'in')

#PCA
#pca分析(此处使用R内置函数prcomp()函数进行分析)
otu_pca <- prcomp(otu_data,scal=TRUE)
#预览pca分析结果
otu_pca_sum <- summary(otu_pca)
#提取出PC1及PC2的坐标
otu_pca_points <- data.frame(otu_pca$x[,1:2])
otu_pca_points$Group <- group_info
#计算各主成分解释度
exv <- otu_pca_sum$importance[2,]*100

PCA_plot <- ggplot(otu_pca_points, aes(x=PC1, y=PC2, color=Group, group = Group)) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  labs(x=paste("PC1 (", round(exv[1],2), "%)", sep=""),
       y=paste("PC2 (", round(exv[2],2), "%)", sep="")) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank())+
  geom_text_repel(aes(label=rownames(otu_pca_points)),size=3,max.overlaps = 20)

PCA_plot
ggsave('PCA.svg',PCA_plot,width = 5,height = 5,units = 'in')

