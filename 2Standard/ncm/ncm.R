
library(Hmisc)
library(minpack.lm)
library(stats4)
library(stringr)


#spp: 物种或分类群的丰度表，行是分类群，列是样本
spp<-read.table('OtuForNet.tsv',head=T,stringsAsFactors=F,row.names=1,sep = "\t")

CK <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="CK"))
AD <- subset(spp,select = which(str_extract(colnames(spp),"[[:alpha:]]+")=="AD"))

NCM <- function(grp){
  
        spp<-t(grp)
        ##将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
        #用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
        N <- mean(apply(spp, 1, sum))
        p.m <- apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
        spp.bi <- 1*(spp>0)
        freq <- apply(spp.bi, 2, mean)
        freq <- freq[freq != 0]
        C <- merge(p, freq, by=0)
        C <- C[order(C[,2]),]
        C <- as.data.frame(C)
        C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
        p <- C.0[,2]
        freq <- C.0[,3]
        names(p) <- C.0[,1]
        names(freq) <- C.0[,1]
        d = 1/N
        m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
        m.fit  #获取 m 值
        m.ci <- confint(m.fit, 'm', level=0.95)
        freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
        pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
        Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
        Rsqr  #获取模型的 R2
        Rsqr <- round(Rsqr,3)
        Nm <- round(coef(m.fit)*N,3)
        m <- round(coef(m.fit),3)
        
        #输出 3 个统计结果数据表，包括各物种或分类群的平均相对丰度（p.csv）、出现频率（freq.csv）和预测的出现频率（freq.pred.csv）
        # write.csv(p, file = paste0(deparse1(substitute(grp)),"_p.csv"))
        # write.csv(freq, file = paste0(deparse1(substitute(grp)),"_freq.csv"))
        # write.csv(freq.pred, file = paste0(deparse1(substitute(grp)),"_freq.pred.csv"))
        
        #p 是平均相对丰度（mean relative abundance）
        #freq 是出现频率（occurrence frequency）的观测值
        #freq.pred 是出现频率（occurrence frequency）的预测值，即中性模型的拟合值
        
        #绘制统计图
        bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
        bacnlsALL$Group <- deparse1(substitute(grp))
        bacnlsALL$Rsqr <- Rsqr
        bacnlsALL$Nm <- Nm
        bacnlsALL$m <- m
        bacnlsALL$U_L <-  "In prediction"
        bacnlsALL$U_L[bacnlsALL$freq <= bacnlsALL$Lower] <- "Below prediction"
        bacnlsALL$U_L[bacnlsALL$freq >= bacnlsALL$Upper] <- "Above prediction"
        return(bacnlsALL)
}

NCM_CK<- NCM(CK)
NCM_AD<- NCM(AD)

Resultall <- rbind(NCM_CK,NCM_AD<- NCM(AD))
Resultall$Group <- factor(Resultall$Group,levels = c('CK','AD'))

Resultall$R_N <- paste("R^2=",Resultall$Rsqr,"     ","m=",Resultall$m)

library(ggplot2)
p  <- ggplot(data = Resultall,aes(x=log10(p),y=freq,color=U_L)) +
  geom_point(position = "identity",size = .01)+
  scale_color_manual(values = c('#A52A2A',"#29A6A6","black"))+
  labs(x='Mean Relative Abundance (log10)',y='Frequency of Occurance')+
  # annotate(geom = "text",x=-Inf,y=Inf,vjust=1.5,hjust= -0.1,label= paste("R^2=",Resultall$Rsqr,"\n","Nm=",Resultall$Nm))+
  geom_line(aes(y=freq.pred),color="blue",lwd=1,alpha=0.5)+
  geom_line(aes(y=Lower),color="blue",lwd=0.5,lty=2)+
  geom_line(aes(y=Upper),color="blue",lwd=0.5,lty=2)+
  facet_wrap(~Group+R_N,ncol = 2,scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(0,0),
        legend.position = c(0.25,0.05),
        legend.background = element_rect(linewidth = .1,linetype = 2,color = "black"),
        legend.spacing.y = unit(0,"line"),
        legend.text = element_text(size = unit(7,"pt")),
        legend.margin = margin(rep(5,4),unit = "pt"),
        legend.key.size = unit(10,"pt"))
p  
ggsave("ncm.svg",p,width = 6,height = 3,units = "in")

