
library(tidyr)
library(tibble)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

Data <- read.table("Genus_LEfSe_Input_Day10_AD.tsv",sep = '\t',header = T)
Tax <- separate(Data[-1,],1,into = c("K","P","C","O","F","G"),sep = "[|]")
TaxLDA <- read.table("lefse.LDA_Day10 .tsv",sep = '\t',header = F)
colnames(TaxLDA) <- c("Tax","Log","Group","LDA","p")

ADSigTax <- separate(subset(TaxLDA,Group=="AD"),1,into = c("K","P","C","O","F","G"),sep = "[.]")[,1:6]

Posit <- apply(ADSigTax, 1, function(x){which(x==tail(x[!is.na(x)],1))
} )
TAX <- apply(ADSigTax, 1, function(x){tail(x[!is.na(x)],1)})

TAX_Posit <- as.data.frame(t(cbind(TAX,Posit)))
colnames(TAX_Posit) <- TAX_Posit[1,]

TAX_Sig_Abun <- lapply(TAX_Posit, function(x){Tax[which(Tax[,as.numeric(x[2])]==x[1]),-(1:6)]}) %>% 
  rbindlist(.,idcol = "Feature")

##########
p__Tenericutes__c__unclassified <- subset(Tax,P=="p__Tenericutes"&C=="unclassified")[,-(1:6)] 
p__Tenericutes__c__unclassified$Feature <- "p__Tenericutes__c__unclassified"
p__Tenericutes__c__unclassified <- p__Tenericutes__c__unclassified %>%  select("Feature",everything())
p__Tenericutes__f__unclassified <- subset(Tax,P=="p__Tenericutes"&`F`=="unclassified")[,-(1:6)] 
p__Tenericutes__f__unclassified$Feature <- "p__Tenericutes__f__unclassified"
p__Tenericutes__f__unclassified <- p__Tenericutes__f__unclassified %>%  select("Feature",everything())
p__Tenericutes__o__unclassified <- subset(Tax,P=="p__Tenericutes"&O=="unclassified")[,-(1:6)] 
p__Tenericutes__o__unclassified$Feature <- "p__Tenericutes__o__unclassified"
p__Tenericutes__o__unclassified <- p__Tenericutes__o__unclassified %>%  select("Feature",everything())
##########

TAX_Sig_Abun1<- rbind(TAX_Sig_Abun,p__Tenericutes__c__unclassified,p__Tenericutes__f__unclassified,p__Tenericutes__o__unclassified)

TAX_Sig_Abun2 <- TAX_Sig_Abun1[,-1] %>% apply(.,2,as.numeric) %>% 
  cbind(TAX_Sig_Abun1[,1],.) %>% 
  group_by(Feature) %>%
  summarize_all(sum) 

plotdata <- pivot_longer(TAX_Sig_Abun2,
                         cols = 2:ncol(TAX_Sig_Abun2),
                         names_to = "SampleID",
                         values_to = "Relative_abundance"
)
plotdata$Group <- str_extract(plotdata$SampleID,"[[:alpha:]]+")
plotdata <- plotdata %>% mutate(Group=case_when(Group == "A" ~ "AOM",TRUE~as.character(Group)))
# plotdata$Group <- factor(plotdata$Group,levels = c("CK","AD","AOM"))
plotdata$Rank <- factor(str_extract(plotdata$Feature,"[[:alpha:]]"),levels = c("p","c","o","f","g"))
plotdata$Rank[which(plotdata$Feature=="p__Tenericutes__c__unclassified")] <- "c"
plotdata$Rank[which(plotdata$Feature=="p__Tenericutes__f__unclassified")] <- "f"
plotdata$Rank[which(plotdata$Feature=="p__Tenericutes__o__unclassified")] <- "o"
setorder(plotdata,Rank,-Relative_abundance)
plotdata$Feature <- factor(plotdata$Feature,levels=plotdata$Feature[!duplicated(plotdata$Feature)])
plotdata$Group <- factor(plotdata$Group,levels = c("CK","AD","AOM"))

LEfSe_abundance <- plotdata[c("Feature","Relative_abundance","Group")] %>% 
  group_by(.,Feature,Group) %>% 
  summarize_each(.,mean) %>% 
  pivot_wider(id_col = "Group",names_from = "Feature",values_from = "Relative_abundance")
write.table(LEfSe_abundance,"LEfSe_abundance_higher.tsv",row.names = F,sep = "\t")

library(ggplot2)
library(ggsci)

p <- ggplot(data = plotdata,aes(x=Group,y=Relative_abundance,fill = Group))+
  stat_summary(geom = "bar",fun = "mean",)+
  stat_summary(geom = "errorbar",width = 0.2,
               fun.max = function(x)mean(x)+sd(x),
               fun.min = function(x)mean(x)-sd(x))+
  scale_y_continuous(labels = scales::percent_format())+
  scale_fill_lancet(alpha = 0.5)+
  facet_wrap(~Feature,scales = "free_y")+
  theme_bw()+
  theme(legend.position = c(0.9,0.1))

p

ggsave("LEfSe_sig_bar_Day10.svg",p,width = 12,height = 8,units = "in")

