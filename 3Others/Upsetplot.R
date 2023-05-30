
library(ComplexUpset)
library(ggplot2)

data <- read.table("Upsetdata.tsv",header = T,sep = "\t",row.names = 1)

p1 <- upset(data,colnames(data),
      intersections= list(c('CK_Day10'),c('CK_Day25'),c('AD_Day10'),c('AD_Day25'),c('AOM_Day10'),c('AOM_Day25'),
                          c('CK_Day10', 'AD_Day10'),c('AOM_Day10', 'AD_Day10'),c('CK_Day25', 'AD_Day25'),c('AOM_Day25', 'AD_Day25')),
      min_degree=1,max_degree=2,mode = "inclusive_intersection", width_ratio=0.2, height_ratio = 0.4,
      base_annotations = list('intersection_size'=intersection_size(text=list(vjust=-0.2,hjust=0.5,angle=0),width = 0.6,mode = "inclusive_intersection"))
      # queries=list(upset_query(intersect=c('CK_Day60', 'AD_Day60'), color='orange',fill = "orange"),
      #              upset_query(intersect=c('CK_Day25', 'AD_Day25'), color='red',fill = "red"),
      #              upset_query(intersect=c('CK_Day10', 'AD_Day10'), color='blue',fill = "blue")
      #             )
      )+
      theme(plot.margin = margin(t = -1.5, r = -2, b = -1.5, l = -2, unit = "cm"))
p1  

ggsave('upsetplot.svg',p1,width = 6,height = 6,units = 'in')

library(ggvenn)

venndatarow <- t(data)
venndata <- list()


for(i in rownames(venndatarow)){
    venndata[[i]] <- names(venndatarow[i,])[venndatarow[i,]>0]
}

data_Day10 <- venndata[c(1:2)]
data_Day25 <- venndata[c(3:4)]
data_Day60 <- venndata[c(5:6)]


p2 <- ggvenn(data_Day10,stroke_color = "grey",set_name_size = 4,fill_color = c("#FFF5CA","#A6BEFF"))
p3 <- ggvenn(data_Day25,stroke_color = "grey",set_name_size = 4,fill_color = c("#FFF5CA","#A6BEFF"))
p4 <- ggvenn(data_Day60,stroke_color = "grey",set_name_size = 4,fill_color = c("#FFF5CA","#A6BEFF"))


# library(gridExtra)
# 
# p5 <- ggdraw() +
#         draw_plot(p2, x = 0, y = 2/3, width = 1, height = 1/3) + 
#         draw_plot(p3, x = 0, y = 1/3, width = 1, height = 1/3) + 
#         draw_plot(p4, x = 0, y = 0, width = 1, height = 1/3)
# 
# ggdraw() +
#   draw_plot(p1, x = 0, y = 0, width = 0.6, height = 1) + 
#   draw_plot(p5, x = 0.6, y = 0, width = 0.4, height = 1) +
#   draw_plot_label(label = c("a", "b"),
#                   x = c(0.05, 0.65), y = c(1, 1), size = 18)


library(cowplot)
p <- ggdraw() +
      draw_plot(p1, x = 0, y = 0, width = 0.6, height = 1) + 
      draw_plot(p2, x = 0.65, y = 2/3, width = 0.35, height = 1/3) + 
      draw_plot(p3, x = 0.65, y = 1/3, width = 0.35, height = 1/3) +
      draw_plot(p4, x = 0.65, y = 0, width = 0.35, height = 1/3) +
      draw_plot_label(label = c("a", "b","c","d"),
                    x = c(0, 0.7, 0.7, 0.7), y = c(1, 1, 2/3, 1/3), size = 20)
p      
ggsave('upsetplot1.svg',p,width = 25,height = 15,units = 'in')


