#参考https://mp.weixin.qq.com/s/125K35O6Q7HAz27XhzcDBw#
detach(ggClusterNet)
remove.packages(ggClusterNet)
devtools::install_github('taowenmicro/ggClusterNet')


#--导入所需R包#-------
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(ggsci)

#-----导入数据#-------
metadata <- read.delim('GroupForNet.tsv', row.names = 1)
otutab <- read.delim('OtuForNet.tsv', row.names = 1)
taxonomy <- read.table('Taxonomy.tsv', row.names = 1, header = T)

#-----转换数据#-------
ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)))
ps
rank_names(ps)

#-----分组网络#-------
groupNetresult = network(ps = ps,
                         N = 200,
                         layout_net = "model_Gephi.2",
                         r.threshold=0.8,
                         p.threshold=0.05,
                         method = "pearson",
                         ncol = 2,
                         label = FALSE,
                         path = '.',
                         zipi = TRUE)
groupNet = groupNetresult[[1]]
groupNet
ggsave(filename = 'groupNet.pdf',plot = groupNet,width = 20,height = 10,units = 'in')

#-----网络稳健性#------
robustresult = Robustness.Random.removal(ps = ps,
                                         Top = 200,
                                         r.threshold= 0.8,
                                         p.threshold=0.05,
                                         method = "pearson")
robustplot = robustresult[[1]]
robustplot

#--计算负相关的比例#----
negCorr = negative.correlation.ratio(ps = ps,
                                     Top = 500,
                                     degree = TRUE,
                                     zipi = FALSE,
                                     r.threshold= 0.9,
                                     p.threshold=0.01,
                                     method = "pearson")

negCorrplot = negCorr[[1]]
negCorrplot

#------网络抗毁性#------
library(tidyfst)
library(pulsar)
library(ggClusterNet)
library(phyloseq)
library(tidyverse)
naturalConresult = natural.con.microp (ps = ps,
                                       Top = 500,
                                       r.threshold= 0.9,
                                       p.threshold=0.01,
                                       method = "spearman",
                                       norm = F,
                                       end = 150,# 小于网络包含的节点数量
                                       start = 0,
                                       con.method = "pulsar")
naturalConplot = naturalConresult[[1]]
naturalConplot

# #-----模块比较#------
# mdCompresult = module.compare.net.pip(ps = ps,
#                                       Top = 500,
#                                       degree = TRUE,
#                                       zipi = FALSE,
#                                       r.threshold= 0.9,
#                                       p.threshold=0.01,
#                                       method = "spearman",
#                                       padj = F,
#                                       n = 3)

################################################################################

#------单个网络#------
psCK = ps %>%
  filter_taxa(function(x) sum(x) > 0, TRUE) %>%
  scale_micro("rela") %>%
  subset_samples.wt("Group" ,'CK') %>%
  filter_OTU_ps(500)

psAD = ps %>%
  filter_taxa(function(x) sum(x) > 0, TRUE) %>%
  scale_micro("rela") %>%
  subset_samples.wt("Group" ,'AD') %>%
  filter_OTU_ps(500)

#-----模块化展示#----
moduleResu = module_display.2(pst = psCK,
                              r.threshold= 0.9,
                              p.threshold=0.01,
                              # select.mod = c("model_1","model_2","model_3","model_4"),#选择指定模块可视化
                              Top = 500,
                              num = 20, # 模块包含OTU数量少于55个的不展示,
                              leg.col = 5)
moduleplot <- moduleResu[[3]]

moduleplot+
  scale_color_d3()+
  scale_fill_d3()+
  theme(legend.position = 'right')

#--这里我们指定三个模块，绘制这三个模块物种组成图表
select.mod = c("model_11","model_15","model_4","model_9")
#--模块信息，共有三列，第一列是OTU，第二列是吗模块的分类，第三列是模块捏的边数量
mod1 = moduleResu$mod.groups %>% filter(group %in% select.mod)
#按照属水平进行绘制
moduleCompres = module_composition(pst = psCK,mod1 = mod1,j = "Genus")
moduleCompplot = moduleCompres[[1]]
moduleCompplot
ggsave(filename = 'moduleCompplot.pdf',moduleCompplot,width = 10,height = 5,units = 'in')

################################################################################
#--sparcc方法计算相关矩阵,多组网络的比较#----
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = FALSE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.8,
  p.threshold = 0.05,
  maxnode = 2,
  method = "spearman",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 3
)

#  建议保存一下输出结果为R对象，方便之后不进行相关矩阵的运算，节约时间
saveRDS(tab.r,"network.pip.sparcc.rds")
tab.r = readRDS("./network.pip.sparcc.rds")

#-提取全部图片的存储对象
plot = tab.r[[1]]
# 提取网络图可视化结果
p0 = plot[[1]]#网络图
p0
p1 = plot[[2]]#ZiPi图
p1
p2 = plot[[3]]#随机网络幂率分布
p2

#--提取相关矩阵,这是一个list存储的相关矩阵
dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

# 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
# saveRDS(cortab,"cor.matrix.all.group.rds")
# cor = readRDS("./cor.matrix.all.group.eds")

################################################################################
#--网络显著性比较#-----
mdcmpnetData= module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 3)

mdcmpnetTab = mdcmpnetData[[1]]
head(mdcmpnetTab)
