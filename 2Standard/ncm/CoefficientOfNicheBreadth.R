
library(spaa)
library(stringr)

spp<-read.table('OtuForNet.tsv',head=T,stringsAsFactors=F,row.names=1,sep = "\t")

CK <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="CK"))
AD <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="AD"))

NW_CK <- niche.width(CK, method = 'levins')
NW_AD <- niche.width(AD, method = 'levins')

NW_All <- data.frame(t(cbind(NW_AD,NW_CK)))
colnames(NW_All) <- "niche.width"
NW_All$Group <- str_extract(rownames(NW_All),"[[:alpha:]]+")

library(ggpubr)
library(ggsci)
library(scales)
p <- ggboxplot(data = NW_All,x="Group",y="niche.width",fill = "Group")+
      stat_compare_means(method = "t.test")+
      # stat_compare_means(label.y.npc = 1,label.x.npc = 0.1)+
      scale_fill_d3(alpha = 0.5)+
      labs(x="",y="Niche Width")+
      theme_bw()+
      theme(legend.position = "none")
p  
ggsave("NicheBreadth.svg",p,width = 4,height = 4,units = "in")

