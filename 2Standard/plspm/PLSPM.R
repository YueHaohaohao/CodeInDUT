#安装 plspm 包
#install.packages('devtools')
#devtools::install_github('gastonstat/plspm')
#加载 plspm 包
library(plspm)
#读取数据
dat <- read.delim('plspmdata.tsv', sep = '\t',check.names = F,row.names = 1)
dat <- scale(dat,scale = T,center = T)
#指定潜变量，在 R 中以列表（list）存储变量和潜变量的关系
#您可以直接指定列名称，或者指定列的下标都可以，我个人习惯指定列名称
dat_blocks <- list(
    DA = 'DA',
    AOM = c('TOC','TN'), 
    Env_Drivers = c('Ferrous', 'pH'), 
    Comm_Structure = c('PC1','PC2')
)
dat_blocks
#通过 0-1 矩阵描述潜变量之间的关联，其中 0 代表变量间没有关联，1 代表有关联
DA <- c(0, 0, 0, 0)
AOM <- c(0, 0, 0, 0)
Env_Drivers <- c(1, 1, 0, 0)
Comm_Structure <- c(1, 1, 1, 0)

dat_path <- rbind(DA, AOM, Env_Drivers, Comm_Structure)
colnames(dat_path) <- rownames(dat_path)
dat_path

#指定因果关系，可选 A（代表列是行的因） 或 B（代表行是列的因）
dat_modes <- rep('A', 4)
dat_modes

##一个简单的 PLS-PM，更多参数详情 ?plspm
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes, scaled = T)
dat_pls
summary(dat_pls)

#结果内容比较多，细节部分还需自行参阅 plspm 包的用户手册：
#完整版手册，235页：https://www.gastonsanchez.com/PLS_Path_Modeling_with_R.pdf
#简版手册，10页：https://rdrr.io/cran/plspm/f/inst/doc/plspm_introduction.pdf

#以下仅展示了一部分相对重要的内容

#查看路径系数的参数估计值，以及相关的统计信息
dat_pls$path_coefs
dat_pls$inner_model

#查看因果关系的路径图，详情 ?innerplot
innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

#查看作为外源潜变量和内源潜变量的状态
dat_pls$inner_summary

#查看变量间的影响状态
dat_pls$effects

#查看观测变量和潜变量关系，可通过 outerplot() 画图展示类似路径图的结构，详情 ?outerplot
dat_pls$outer_model
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')
outerplot(dat_pls, what = 'weights', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')

#goodness-of-fit 值可以帮助评估模型优度
dat_pls$gof

#查看潜变量得分，可以理解为标准化后的潜变量的值
dat_pls$scores

library(tidyr)
Effect <- dat_pls$effects
Effect <- Effect %>% separate(.,relationships,into = c("from","to"),sep = " -> ")

plotdata <- subset(Effect,to=="Comm_Structure") %>% 
  pivot_longer(.,cols=3:5,names_to="Effect_Way",values_to = "Effect_Coefficient")
plotdata$from<- factor(plotdata$from,levels = c("DA",'AOM','Env_Drivers'))

library(ggplot2)
library(ggsci)
p <- ggplot(data = plotdata,aes(x=from,y=Effect_Coefficient,fill=Effect_Way,width =0.7))+
      geom_col(position = position_dodge(width = .7),color="black")+
      # geom_bar(stat = "identity",position = position_dodge(width = .9))
      geom_hline(yintercept = 0)+
      scale_y_continuous(limits = c(-1,1),expand = c(0,0))+
      labs(x="")+
      facet_wrap(~to,scales = "free_x")+
      scale_fill_tron(alpha = 0.8)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            legend.position = c(0.5,0.9),
            legend.title = element_blank(),
            legend.direction = 'horizontal')
p
ggsave('Effect.svg',p,width = 3,height = 3,units = "in")


#输出潜变量的值
#latent <- data.frame(dat_pls$scores)
#latent <- cbind(dat$site, latent)
#write.csv(latent, 'latent.csv')
write.table(Effect,'Effect.tsv',sep = '\t',row.names = F)
# write.table(dat_pls$inner_summary,'inner_summary.tsv',sep = '\t')
# write.table(dat_pls$effects,'effects.tsv',sep = '\t',row.names = F)
