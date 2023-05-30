
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

# anova test
tnstpaov <- nst.panova(nst.result=tnst, group=group)
tnstpaov

result <- tnst[["index.grp"]][,c("group","MST.i.ruzicka")]
colnames(result) <- c("Group","Stochasticity")
result$Deterministic <- 1-result$Stochasticity
library(tidyr)
result <- result %>% pivot_longer(.,cols=c("Stochasticity","Deterministic"),names_to="Index",values_to = "Ratio")
result$Group <- factor(result$Group,levels = c('CK','AD'))
# # 导出Bootstrapping test结果
# write.table("\t", file="NST_result.tsv", append = F, quote = F, eol = "", row.names = F, col.names = F)
# write.table(result,"NST_result.tsv",sep = "\t",append = T)
library(ggplot2)
library(ggsci) 
library(scales)

p <-  ggplot(data = result,aes(x=Group,y=Ratio,fill=Index))+
      geom_col(position = position_stack())+
      scale_fill_manual(values = c("#0099b4","#FDAF91"))+
      scale_y_continuous(limits = c(0,1),expand = c(0,0),labels = percent_format())+
      geom_hline(yintercept = 0.5,lty=2)+
      xlab(NULL)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.position = 'top',
            legend.margin = margin(rep(-10,4),unit = "pt"))
p
ggsave("NST.svg",p,width = 3,height = 5,units = "in")
