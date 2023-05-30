library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggsci)
library(multcomp)
library(data.table)

# 读取数据
df <- read.table("Genus.tsv", header = TRUE, sep = "\t", row.names = 1)[1:10,]
df <- rbind(df, 1 - colSums(df))
row.names(df) <- c(row.names(df)[1:10], "Other")
df <- cbind(Genus = rownames(df), df)

# 转换数据格式
df <- pivot_longer(df, cols = -Genus, names_to = "Sample", values_to = "Relative_abundance")
df$Group <- gsub('\\d+', '', df$Sample)

# 绘制堆积柱状图
plotdata <- df %>%
  group_by(Group, Genus) %>%
  summarise(Relative_abundance = mean(Relative_abundance)) %>%
  mutate(Group = case_when(Group == 'A' ~ 'AOM', TRUE ~ as.character(Group)),
         Group = factor(Group, levels = c('CK', 'AD', 'AOM', 'DA')),
         Genus = factor(Genus, levels = unique(df$Genus)))

p <- ggplot(data = plotdata, aes(x = Group, y = Relative_abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  scale_fill_igv(alpha = 0.9) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0.01), labels = scales::percent_format()) +
  theme_bw() +
  labs(x = NULL) +
  theme(panel.grid = element_blank())
p
ggsave("Genus.svg", p, width = 4, height = 4, units = "in")

# 显著性分析
HSD_Result <- df %>%
  group_by(Genus) %>%
  do({
    aov_result <- aov(Relative_abundance ~ Group, data = .)
    hsd_result <- agricolae::HSD.test(aov_result, 'Group')
    merge(hsd_result$means, hsd_result$groups, by = 'row.names')
  }) %>%
  ungroup() %>%
  rename(Group = 'Row.names', Relative_abundance = 'Relative_abundance.x', Label = 'groups') %>%
  mutate(Group = case_when(Group == "A" ~ "AOM", TRUE ~ as.character(Group)),
         Group = factor(Group, levels = c('CK', 'AD', 'AOM', 'DA')))

p1 <- ggplot(data = HSD_Result, aes(x = Group, y = Relative_abundance, fill = Label)) +
  geom_errorbar(aes(ymin = Relative_abundance, ymax = Relative_abundance + std), width = 0.2) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(y = Relative_abundance + std, label = Label), vjust = -0.2, size = 5) +
  scale_fill_d3("category20") +
  scale_y_continuous(expand = c(0.1, 0), labels = scales::percent_format()) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  facet_wrap(~Genus, scales = "free") +
  labs(x = NULL)
p1

ggsave("Genus_sig.svg", p1, width = 8, height = 6, units = "in")
