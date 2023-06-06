
#https://cloud.tencent.com/developer/article/1634477

library(NST)
comm=t(read.delim("OtuForNet.tsv",row.names= 1))[,1:500]
group=read.delim("GroupForNet.tsv",row.names= 1)
tnst=tNST(comm=comm, group=group, dist.method="jaccard",output.rand=TRUE,nworker=4)

# 查看结果，NST.i.bray
tnst[[2]]

# Bootstrapping test
tnstbt <- nst.boot(nst.result=tnst, group=group)
tnstbt[[2]]

# # Bootstrapping test
tnstbt <- nst.boot(nst.result=tnst, group=group)
summary <- tnstbt[['summary']] %>% filter(.,Index=='NST')
summary$Group <- factor(summary$Group,levels = c('CK','AD','AOM','DA'))
compare <- tnstbt[['compare']]

library(ggplot2)
library(ggsignif)
library(ggsci)

p <- ggplot(summary,aes(x=Group,y=mean,fill=Group))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin=mean-stdev,ymax=mean+stdev),width=0.2)+
  scale_y_continuous(limits = c(0,0.9),expand = c(0,0),labels = scales::percent_format())+
  xlab(NULL)+
  ylab('NST Index')+
  scale_fill_lancet(alpha = 0.8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none')+
  geom_signif(comparisons = list(c('CK','AD'),
                                 c('AD','AOM')),
              annotations = "*",
              y_position = c(0.75,0.75))
p
ggsave("NST1.svg",p,width = 4,height = 4,units = "in") 
