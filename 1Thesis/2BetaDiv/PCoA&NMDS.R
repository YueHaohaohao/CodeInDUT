
library(tibble)
library(dplyr)
library(vegan)

df <- t(read.delim("AllOTU.tsv",row.names = 1,header = T, sep = "\t"))
grp <- read.table("Group.tsv",header = T,row.names = 1,sep = "\t")

data <-  merge(grp,df,by="row.names")%>%
          column_to_rownames("Row.names")

data$Group <- factor(data$Group,levels = c("CK","DA","AOM","AD"))
data$Day <- factor(data$Day,levels = c("Day10","Day25","Day60"))

data <- data[order(data$Day,data$Group),]

#计算距离矩阵
PERMANOVA <- data.frame()
Percent_PCoA <- data.frame()
Result_PCoA <- data.frame()
Result_NMDS <- data.frame()
Stress_NMDS <- list()


for (day in  levels(data$Day)){
  
  day_data <- subset(data,Day==day)

  #PERMANOVA
  otu.div <- adonis2(day_data[,-(1:2)] ~ Group, data = day_data, permutations = 999, method="bray")
  otu.div$Day <- day
  PERMANOVA <- rbind(PERMANOVA ,otu.div)
  
  
  #距离矩阵
  otu_dist <- vegdist(day_data[,-(1:2)], method="bray", binary=F)
  
  #PCoA
  otu_pcoa <- cmdscale(otu_dist, k=2, eig=T)
  otu_pcoa_points <- as.data.frame(otu_pcoa$points)
  sum_eig <- sum(otu_pcoa$eig)
  eig_percent <- t(as.data.frame(round(otu_pcoa$eig/sum_eig*100,1)[1:2]))
  
  colnames(eig_percent) <- paste0("PCoA", 1:2)
  rownames(eig_percent) <- day
  Percent_PCoA <- rbind(Percent_PCoA,eig_percent)
  
  colnames(otu_pcoa_points) <- paste0("PCoA", 1:2)
  otu_pcoa_result <- cbind(otu_pcoa_points,day_data[,1],day)
  Result_PCoA <- rbind(Result_PCoA,otu_pcoa_result)
  
  #NMDS
  otu_nmds <- metaMDS(otu_dist, k = 2)
  otu_nmds_result <- data.frame(otu_nmds$points)
  otu_nmds_result <- cbind(otu_nmds_result,day_data[,1],day)
  Stress_NMDS <- append(Stress_NMDS,otu_nmds$stress)
  Result_NMDS <- rbind(Result_NMDS,otu_nmds_result)

}
Percent_PCoA$Day <- row.names(Percent_PCoA)
PERMANOVA <- cbind(rownames(PERMANOVA),PERMANOVA)
colnames(Result_PCoA)[3:4] <- c("Group","Day")
colnames(Result_NMDS)[3:4] <- c("Group","Day")

# write.table(PERMANOVA,"PERMANOVA.tsv",sep = "\t",row.names = F)



#图
library(ggplot2)
library(ggalt)
library(ggrepel)
library(ggsci)

PCoA_day10 <- ggplot(Result_PCoA[which(Result_PCoA$Day=="Day10"),], aes(x=PCoA1, y=PCoA2, color=Group, group = Group)) +
  labs(x=paste("PC1 (", Percent_PCoA[1,1], "%)", sep=""),
       y=paste("PC2 (", Percent_PCoA[1,2], "%)", sep="")
  ) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F ) +
  scale_y_continuous(limits = c(-0.1,0.1),expand = c(0.02,0),breaks = seq(-0.1,0.1,by = 0.1))+
  scale_x_continuous(limits = c(-0.3,0.3),expand = c(0.02,0),breaks = seq(-0.3,0.3,by = 0.3))+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw() +
  ggtitle("Day10_PCoA")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        )+
  annotate(geom= "text",x = -0.1, y = 0.09, label = paste("PERMANOVA: R^2=" ,round(PERMANOVA[1,'R2'],2),";","p=",round(PERMANOVA[1,'Pr(>F)'],3)) )+
  geom_text_repel(aes(label=row.names(Result_PCoA[which(Result_PCoA$Day=="Day10"),])),size=3)

PCoA_day10# 5*5

#####
PCoA_day25 <- ggplot(Result_PCoA[which(Result_PCoA$Day=="Day25"),], aes(x=PCoA1, y=PCoA2, color=Group, group = Group)) +
  labs(x=paste("PC1 (", Percent_PCoA[2,1], "%)", sep=""),
       y=paste("PC2 (", Percent_PCoA[2,2], "%)", sep="")
  ) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F ) +
  scale_y_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  scale_x_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw() +
  ggtitle("Day25_PCoA")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )+
  annotate(geom= "text",x = -0.2, y = 0.45, label = paste("PERMANOVA: R^2=" ,round(PERMANOVA[4,'R2'],2),";","p=",round(PERMANOVA[4,'Pr(>F)'],3)) )+
  geom_text_repel(aes(label=row.names(Result_PCoA[which(Result_PCoA$Day=="Day25"),])),size=3)

PCoA_day25# 5*5

#####
PCoA_day60 <- ggplot(Result_PCoA[which(Result_PCoA$Day=="Day60"),], aes(x=PCoA1, y=PCoA2, color=Group, group = Group)) +
  labs(x=paste("PC1 (", Percent_PCoA[3,1], "%)", sep=""),
       y=paste("PC2 (", Percent_PCoA[3,2], "%)", sep="")
  ) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F ) +
  scale_y_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  scale_x_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw() +
  ggtitle("Day60_PCoA")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )+
  annotate(geom= "text",x = -0.2, y = 0.45, label = paste("PERMANOVA: R^2=" ,round(PERMANOVA[7,'R2'],2),";","p=",round(PERMANOVA[7,'Pr(>F)'],3)) )+
  geom_text_repel(aes(label=row.names(Result_PCoA[which(Result_PCoA$Day=="Day60"),])),size=3)

PCoA_day60# 5*5

#####
library(cowplot)
PCoA_plot <- plot_grid(PCoA_day10,PCoA_day25,PCoA_day60,rel_widths = c(1,1,1),nrow = 1,align = "h",axis = "tb")
PCoA_plot

#####
NMDS_day10 <- ggplot(Result_NMDS[which(Result_NMDS$Day=="Day10"),], aes(x=MDS1, y=MDS2, color=Group, group = Group)) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F ) +
  scale_y_continuous(limits = c(-5e-04,5e-04),expand = c(0.02,0),breaks = seq(-5e-04,5e-04,by = 5e-04))+
  scale_x_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw() +
  ggtitle("Day10_NMDS")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )+
  geom_text_repel(aes(label=row.names(Result_NMDS[which(Result_NMDS$Day=="Day10"),])),size=3)+
  annotate(geom= "text",x = -0.4, y = 4.5e-04, label = paste("Stress=" ,round(Stress_NMDS[[1]],2)))

NMDS_day10

#####
NMDS_day25 <- ggplot(Result_NMDS[which(Result_NMDS$Day=="Day25"),], aes(x=MDS1, y=MDS2, color=Group, group = Group)) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F ) +
  scale_y_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  scale_x_continuous(limits = c(-1.5,1.5),expand = c(0.02,0),breaks = seq(-1.5,1.5,by = 1.5))+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw() +
  ggtitle("Day25_NMDS")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )+
  geom_text_repel(aes(label=row.names(Result_NMDS[which(Result_NMDS$Day=="Day25"),])),size=3)+
  annotate(geom= "text",x = -1.2, y = 0.45, label = paste("Stress=" ,round(Stress_NMDS[[2]],2)))

NMDS_day25

#####
NMDS_day60 <- ggplot(Result_NMDS[which(Result_NMDS$Day=="Day60"),], aes(x=MDS1, y=MDS2, color=Group, group = Group)) +
  geom_point(size=2)+
  scale_color_d3()+
  scale_fill_d3()+
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F ) +
  scale_y_continuous(limits = c(-0.5,0.5),expand = c(0.02,0),breaks = seq(-0.5,0.5,by = 0.5))+
  scale_x_continuous(limits = c(-1,1),expand = c(0.02,0),breaks = seq(-1,1,by = 1))+
  geom_vline(aes(xintercept=0),linetype=5,col="grey")+
  geom_hline(aes(yintercept=0),linetype=5,col="grey")+
  theme_bw() +
  ggtitle("Day60_NMDS")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )+
  geom_text_repel(aes(label=row.names(Result_NMDS[which(Result_NMDS$Day=="Day60"),])),size=3)+
  annotate(geom= "text",x = -0.8, y = 0.45, label = paste("Stress=" ,round(Stress_NMDS[[3]],2)))

NMDS_day60

#####
NMDS_plot <- plot_grid(NMDS_day10,NMDS_day25,NMDS_day60,rel_widths = c(1,1,1),nrow = 1,align = "h",axis = "tb")

PCoA_NMDS <- plot_grid(PCoA_plot,NMDS_plot,nrow = 2)
PCoA_NMDS

ggsave("PCoA_NMDS_H1.svg",PCoA_NMDS,width = 12,height = 8,units = "in")


