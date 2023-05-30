#调用R包
library(dplyr)
library(linkET)
library(ggplot2)
#读取数据
speciese <- read.table("Feature.tsv",header = T,row.names = 1)#读取物种数据
env <- read.table("Envfactor.tsv",header = T,row.names = 1,check.names = F)#读取环境因子
### mantel test
mantel <- mantel_test(speciese, env,
                        spec_select = list(Genus_Abun = 1:20,
                                          Alpha_div = 24:26,
                                          Beta_div = 27:28)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
mantel
write.table(mantel,"manteltest.tsv",sep = "\t",row.names = F)

#mantel test部分设为默认配色set_default_style()
#设为其他配色
set_corrplot_style(colours = c("#0B66AA", "white", "#D3001A"))
#绘图
manteltest_plot <- qcorrplot(correlate(env, method = "spearman",engine = "Hmisc"), type = "upper", diag = FALSE) +
  #环境因子相关的填充方式
  geom_square() + 
  #显著性标记
  geom_mark( sig_thres = 0.05, only_mark = T, sig_level = 0.05, mark = "*", vjust = 0.65, size = 5, color = "grey90")+ 
  #设置mantel test与相关性热图
  geom_couple(aes(colour = pd, size = rd), data = mantel, 
              curvature = nice_curvature(0,by="to"),
)+
  # 添加对角线标签，根据出图效果调整，避免标签被遮蔽
  geom_diag_label(
    geom = "text",
    angle=0,
    nudge_x = -0.2,
    nudge_y = 0.3)+
  #只显示对角线标签，element_blank不起作用，用element_text(colour =" white")代替
  theme(
    axis.text.x.top = element_blank(),
    axis.text.y.left = element_text()
    )+
  #设置相关性热图的颜色
  #scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) + 
  scale_size_manual(values = c(0.3, 0.5, 1)) +
  #设置线的颜色
  scale_colour_manual(values = c("#ed0000","#42b540","#C6A88D")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35")),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3)),
         fill = guide_colorbar(title = "Spearman's r")
         )
  # labs(title = "Mantel test")
manteltest_plot

ggsave("mantel_test_plot1.svg",manteltest_plot,width = 8, height = 6, units = "in")

# #热图设置在左下
# #绘图
# manteltest_plot1 <- qcorrplot(correlate(env, method = "spearman",engine = "Hmisc"), type = "lower", diag = FALSE) +
#   #环境因子相关的填充方式
#   geom_square() + 
#   #显著性标记
#   geom_mark( sig_thres = 0.05, only_mark = T, sig_level = 0.05, mark = "*", vjust = 0.65, size = 5, colour = "grey90")+ 
#   #设置mantel test与相关性热图
#   geom_couple(aes(colour = pd, size = rd), data = mantel, 
#               curvature = nice_curvature(0.1,by="to"),
#               offset_x=list(Genus_Abun=-2.5,Alpha_div=-2.5,Beta_div=-3.5),
#               offset_y=list(Genus_Abun=0.5,Alpha_div=2,Beta_div=2.5))+
#   # 添加对角线标签，根据出图效果调整，避免标签被遮蔽
#   geom_diag_label(
#     geom = "text",
#     angle=0,
#     nudge_x = -0.2,
#     nudge_y = 0.3)+
#   #只显示对角线标签，element_blank不起作用，用element_text(colour =" white")代替
#   theme(
#     axis.text.x.bottom = element_blank(),
#     axis.text.y.left = element_text(),
#     legend.position = "bottom",
#     legend.justification = c(0,0)
#   )+
#   #设置相关性热图的颜色
#   #scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) + 
#   scale_size_manual(values = c(0.3, 0.5, 1)) +
#   #设置线的颜色
#   scale_colour_manual(values = c("#ed0000","#42b540","#C6A88D")) +
#   #调整图例
#   guides(size = guide_legend(title = "Mantel's r",
#                              override.aes = list(colour = "grey35"), 
#                              title.position = "top",
#                              direction = "horizontal",
#                              order = 2),
#          colour = guide_legend(title = "Mantel's p", 
#                                override.aes = list(size = 3),
#                                title.position = "top",
#                                direction = "horizontal",
#                                order = 1),
#          fill = guide_colorbar(title = "Spearman's r", order = 3,
#                                title.position = "top",
#                                direction = "horizontal",barheight = 0.2)
#   )
# # labs(title = "Mantel test")
# manteltest_plot1



