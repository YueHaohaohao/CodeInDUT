library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(multcomp)

data <- read.table('AlphaDivIndex.tsv', header = TRUE, sep = '\t', row.names = 1)
data$Group <- gsub('\\d+', '', rownames(data))

data <- data %>% 
  mutate(Group = case_when(Group == 'A' ~ 'AOM', TRUE ~ as.character(Group)),
         Group = factor(Group, levels = c('CK', 'AD', 'AOM', 'DA')))

Index_SRE <- data %>% subset(.,select= c('Group','Shannon','Richness','Evenness')) %>% 
  pivot_longer(cols = -Group, names_to = 'Index', values_to = 'value')

# 绘制箱线图
Bxp <- ggboxplot(Index_SRE, x = 'Group', y = 'value',
                 bxp.errorbar = TRUE, fill = 'Group',
                 ggtheme = theme_bw(), facet.by = 'Index',
                 scales = 'free') +
  scale_y_continuous(expand = c(0.2, 0)) +
  labs(x = NULL, y = NULL)+
  theme(legend.position = 'none',panel.grid = element_blank())+
  stat_compare_means(method = 't.test',comparisons = list(c('AD', 'CK'),
                                        c('AD', 'AOM')),
                     label = 'p.signif',tip.length = 0.01,vjust = -0.5)
Bxp
ggsave('AlphaBox1.svg', Bxp, width = 7, height = 3, units = 'in')

HSD_Reslut <- lapply(split(Index_SRE, ~ Index), function(x) {
  aov(value ~ Group, data = x) %>% 
    agricolae::HSD.test(., 'Group') %>% 
    {merge(.$means, .$groups, by = 'row.names')}
}) %>% 
  data.table::rbindlist(use.names = FALSE, idcol = TRUE) %>% 
  rename(Index = .id, Group = Row.names, mean = value.x, Label = groups) 

HSD_Reslut$mean <- round(HSD_Reslut$mean, 2)
HSD_Reslut$Group <- factor(HSD_Reslut$Group, levels = c('CK', 'AD', 'AOM', 'DA'))
HSD_Reslut$Index <- factor(HSD_Reslut$Index, levels = c('Shannon', 'Richness', 'Evenness'))

#均值+误差棒+多重比较
p <- ggplot(data = HSD_Reslut, aes(x = Group, y = mean, color = Group)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std, width = 0.2)) +
  scale_y_continuous(expand = c(0.4, 0)) +
  facet_wrap(~ Index, scales = 'free') +
  scale_color_lancet(alpha = 0.8) +
  geom_text(aes(x = Group, y = mean + std, label = Label, vjust = -1)) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'none', panel.grid = element_blank())
ggsave('AlphaBox2.svg', p, width = 7, height = 3, units = 'in')

#箱线图+多重比较
p1 <- ggplot(data = HSD_Reslut, aes(x = Group, fill = Group)) +
  geom_errorbar(aes(ymin = Min, ymax = Max, width = 0.5)) +
  geom_boxplot(aes(y = mean, ymin = Min, lower = Q25, middle = Q50, upper = Q75, ymax = Max),
               stat = 'identity') +
  scale_color_lancet(alpha = 0.8) +
  geom_text(aes(x = Group, y = Max, label = Label, vjust = -1)) +
  facet_wrap(~ Index, scales = 'free') +
  scale_y_continuous(expand = c(0.2, 0)) +
  theme_bw()+
  labs(x=NULL,y=NULL)+
  theme(legend.position = 'none',panel.grid = element_blank())
p1
ggsave('AlphaBox.svg', p1, width = 7, height = 3, units = 'in')
