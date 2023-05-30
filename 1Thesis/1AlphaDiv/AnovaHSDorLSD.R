############################################
data <- read.table("AlphaDivIndex.tsv",header=T,sep="\t",row.names = 1)
############################################
library(stringr)
data$Group <- str_extract(rownames(data),"[[:alpha:]]+")
data$Day   <- str_extract(rownames(data),"\\d")
############################################
library(dplyr)
data <- data %>% 
  mutate(Day=case_when(Day == "1" ~ "Day10",
                       Day == "2" ~ "Day25",
                       Day == "3" ~ "Day60"),
         Group=case_when(Group == "A" ~ "AOM",
                         TRUE~as.character(Group)))
data$Group <- factor(data$Group,levels= c("CK","AD","AOM","DA"))
############################################
library(tidyr)
library(data.table)
library(agricolae)
library(ggplot2)
library(ggsci)
############################################
Index_SRE <- subset(data,select= c('Day','Group',"Shannon","Richness","Evenness")) %>% 
  pivot_longer(cols = 3:5,names_to = "Index",values_to = "value")
############################################
Sigdata <- Index_SRE %>% 
  split(~Day+Index) %>% 
  lapply(.,function(x){aov(value~Group,data=x)}) %>% 
  lapply(.,HSD.test,'Group') %>% ##TukeyHSD检验#######
##lapply(.,LSD.test,'Group') %>% ##LSD检验##
  lapply(.,function(x){merge(x[['means']],x[['groups']],by = 'row.names')}) %>% 
  rbindlist(.,use.names=F,idcol = T) %>% 
  separate(.,col=.id,sep = "[.]",into=c('Day','Index'))

Sigdata <- Sigdata %>% subset(.,select = c('Day','Index','Row.names','value.x','std','groups'))
colnames(Sigdata) <- c('Day','Index','Group','mean','std','p.sig')
Sigdata$mean <- round(Sigdata$mean,2)
Sigdata$Group <- factor(Sigdata$Group,levels = c('CK','AD','AOM','DA'))
Sigdata$Index <- factor(Sigdata$Index,levels =  c("Shannon","Richness","Evenness"))
Sigdata$Day <- factor(Sigdata$Day,levels =  c("Day10","Day25","Day60"))
Sigdata <- unite(Sigdata,col="Day_Index",Day,Index,sep='_',remove = F)
levelDaYIndex <- with(expand.grid(levels(Sigdata$Index),levels(Sigdata$Day)),paste(Var2,Var1,sep ='_' ))
Sigdata$Day_Index <- factor(Sigdata$Day_Index,levels = levelDaYIndex )
############################################
p <- ggplot(data = Sigdata,aes(x=Group,y=mean,color=Group))+
      geom_point()+
      geom_errorbar(aes(ymin=mean-std,ymax=mean+std,width=0.1))+
      scale_y_continuous(expand = c(0.4,0))+
      facet_wrap(~Day_Index,scales = "free")+
      scale_color_lancet(alpha = 0.8)+
      geom_text(aes(x=Group,y=mean+std,label=p.sig,vjust=-1))+
      theme_bw()+
      labs(x=NULL,y=NULL)+
      theme(legend.position = 'none',panel.grid = element_blank())
p
ggsave('Alpha_Errorbar.svg',p,width = 8,height = 8,units = 'in')
############################################
