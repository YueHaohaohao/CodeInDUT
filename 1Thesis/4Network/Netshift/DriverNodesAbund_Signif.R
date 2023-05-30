
library(ggplot2)
library(tidyr)
library(scales)
library(stringr)
library(dplyr)
library(rstatix)

Drivers <- read.csv("DrivernodesData.csv",check.names = F,header = T)
AllFamilies <- read.delim("Family.tsv",check.names = F)

DriversAbund <- AllFamilies %>% subset(.,Family %in% Drivers$S_ID)
DriversAbund1 <- DriversAbund %>% .[,-ncol(.)]

DriversAbund2 <- DriversAbund1 %>% 
  pivot_longer(., cols = 2:ncol(.),
               names_to = "Sample",
               values_to = "Relative_abundance"
               ) 
DriversAbund2$Day <- str_extract(DriversAbund2$Sample,"\\d")
DriversAbund2$Group <- str_extract(DriversAbund2$Sample,"[[:alpha:]]+")
DriversAbund2 <- DriversAbund2 %>% subset(.,Group=="A"|Group=="AD")
DriversAbund2 <- DriversAbund2 %>% 
  mutate(Day=case_when(Day == "1" ~ "Day10",
                       Day == "2" ~ "Day25"),
         Group=case_when(Group == "A" ~ "AOM",
                         TRUE~as.character(Group)))
DriversAbund2$Group <- factor(DriversAbund2$Group,levels= c("AOM","AD"))

plotdata <- DriversAbund2 %>% 
  group_by(Family,Day) %>%
  dunn_test(Relative_abundance~Group)%>% 
  subset(.,select = c("Family","Day","p.adj.signif")) %>% 
  merge(DriversAbund2,.) %>% 
  subset(.,p.adj.signif=="*") %>% 
  unite(.,"Day_Family",Day,Family)

Siglabel <- subset(plotdata,select = -Sample) %>% 
  arrange(Day_Family,-Relative_abundance) %>% 
  .[!duplicated(.$Day_Family),]

library(ggpubr)
library(ggsci)

p1<- ggplot(data = plotdata,aes(x=Group,y=Relative_abundance,fill= Group))+
  stat_summary(fun = "mean", geom = "bar",position = "dodge")+
  stat_summary(geom = "errorbar",width = 0.2)+
  geom_text(data = Siglabel,aes(y=Relative_abundance,x=1.5,label=p.adj.signif))+
  geom_segment(data = Siglabel,aes(x=1,xend=2,yend=Relative_abundance))+
  scale_fill_lancet(alpha = 0.5)+
  scale_y_continuous(labels = percent_format())+
  facet_wrap(~Day_Family,scales = "free_y",nrow = 1)+
  theme_bw()+
  theme(legend.position = "none",axis.title.x = element_blank())
p1

ggsave("Driver_abundance.svg",p1,width = 14,height = 3,units = "in")

# library(ggpubr)
# library(ggsci)
# p1<- ggplot(data = plotdata,aes(x=Group,y=Family,size=Relative_abundance,color=p.adj.signif))+
#   geom_point()+
#   stat_summary(fun = "mean", geom = "point")+
#   scale_color_lancet()+
#   facet_grid(~Day)

