
library(ggplot2)
library(tidyr)
library(scales)
library(dplyr)


# df <- read.table("Genus.tsv",header = T,sep = "\t")[1:10,] %>% 
#   adorn_totals(where = "row")

df <- read.table("Genus.tsv",header = T,sep = "\t",row.names = 1)[1:10,]
df <- rbind(df,1-colSums(df))
row.names(df) <- c(row.names(df)[1:10],"Other")
df$Genus <- row.names(df)
df <- df %>% select("Genus",everything())


plotdata <- pivot_longer(df,
             cols = 2:ncol(df),
             names_to = "Sample",
             values_to = "Relative_abundance"
             )

library(stringr)

plotdata$Day <- str_extract(plotdata$Sample,"\\d")
plotdata$Group <- str_extract(plotdata$Sample,"[[:alpha:]]+")


plotdata$Genus <- factor(plotdata$Genus,levels= row.names(df))
plotdata$Sample <- factor(plotdata$Sample,levels= colnames(df))


########堆积柱状图
pltdata <- aggregate(plotdata$Relative_abundance,by=list(plotdata$Day,plotdata$Group,plotdata$Genus),FUN = mean)
colnames(pltdata) <- c("Day","Group","Genus","Relative_abundance")
pltdata <- pltdata %>% 
  mutate(Day=case_when(Day == "1" ~ "Day10",
                       Day == "2" ~ "Day25",
                       Day == "3" ~ "Day60"),
         Group=case_when(Group == "A" ~ "AOM",
                         TRUE~as.character(Group)))
pltdata$Group <- factor(pltdata$Group,levels= c("CK","AD","AOM","DA"))

library(ggsci)

p <- ggplot(data = pltdata,aes(x=Group,y=Relative_abundance,fill=Genus))+
  geom_bar(stat = "identity",position = position_fill(reverse = T))+
  scale_fill_igv(alpha = 0.9)+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_y_continuous(expand = c(0,0.01),labels = scales::percent_format())+
  theme_bw()+
  labs(x='')+
  theme(panel.grid = element_blank())+
  facet_wrap(~Day,ncol = 3,scales = "free_x")
p
ggsave("Genus.svg",p,width = 6,height = 4,units = "in")


#######显著性

library(rstatix)
library(tidyverse)
library(agricolae)
library(data.table)

Anova_Result = plotdata %>%
  split(~Day+Genus) %>%
  lapply(.,aov,formula = Relative_abundance ~ Group) %>%
  lapply(.,duncan.test,"Group") %>% 
  lapply(.,function(x){x <- merge(x[["groups"]],subset(x[["means"]],select = c("std")),by="row.names")}) %>%
  rbindlist(.,use.names = T,idcol = "Day.Genus") %>% 
  separate(.,Day.Genus,into = c("Day","Genus")) %>% 
  rename(.,Group=Row.names,label=groups) %>% 
  mutate(Day=case_when(Day == "1" ~ "Day10",
                       Day == "2" ~ "Day25",
                       Day == "3" ~ "Day60"),
         Group=case_when(Group == "A" ~ "AOM",
                         TRUE~as.character(Group)))

Anova_Result$Group <- factor(Anova_Result$Group,levels= c("CK","AD","AOM","DA"))
Anova_Result$Genus <- factor(Anova_Result$Genus,levels= row.names(df))

write.table(Anova_Result,"Genus_sig.tsv",row.names = F, sep = "\t")

library(ggsci)

p1 <- ggplot(data = Anova_Result,aes(x=Group,y=Relative_abundance,color=label))+
     geom_point()+
     geom_errorbar(aes(ymin= Relative_abundance-std,ymax= Relative_abundance+std,width = 0.2))+
     geom_text(aes(y = Relative_abundance + std, label = label), vjust = -0.2,size = 5)+
     scale_fill_d3("category20",alpha = 0.9)+
     scale_y_continuous(expand = c(0,0.01),labels = scales::percent_format())+
     theme_bw()+
     theme(panel.grid = element_blank(),
           legend.position = "none")+
     facet_grid(Day~Genus,scales = "fixed")
p1
ggsave("Genus_sig.svg",p1,width = 15,height = 5,units = "in")
