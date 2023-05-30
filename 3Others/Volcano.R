
#https://mp.weixin.qq.com/s/6Xo71Ev7VWG5iPoTLJxFJw
#https://zhuanlan.zhihu.com/p/437712423

library(ggVolcano)
library(dplyr)

#导入数据
kodata <- read.table('allko.tsv',sep = '\t',row.names = 1,header = T)
group <- c(rep("AD", 3), rep("CK",3)) %>% factor(., levels = c("AD", "CK"), ordered = F)

#差异分析
library(limma) #加载limma包
# 导入数据
kodata <- read.table('allko.tsv', sep = '\t', row.names = 1, header = TRUE)
# 归一化处理
kodata_norm <- normalizeBetweenArrays(kodata)
# 对数转换
kodata_log <- log2(kodata_norm + 1)
# 设置分组信息
group <- c(rep("AD", 3), rep("CK", 3)) %>% factor(levels = c("AD", "CK"), ordered = FALSE)
# 创建设计矩阵
design <- model.matrix(~ group + 0)
colnames(design) <- c('AD', 'CK')
# 拟合线性模型
fit <- lmFit(kodata_log, design)
# 设置对比矩阵
komatrix <- makeContrasts(AD - CK, levels = design)
# 拟合对比
fit_contrast <- contrasts.fit(fit, komatrix)
# 使用eBayes方法进行贝叶斯分析
fit_ebayes <- eBayes(fit_contrast)
# 获取差异结果，包括logFC和调整后的p值
tempOutput <- topTable(fit_ebayes, adjust = "fdr", n = Inf)
#导出差异结果
nrDEG = na.omit(tempOutput) ## 去掉数据中有NA的行或列
diffsig <- nrDEG %>% rownames_to_column(.,'ID')

#图
data <- add_regulate(diffsig, log2FC_name = "logFC",
                     fdr_name = "adj.P.Val",log2FC = 1, fdr = 0.05)

ggvolcano(data, x = "log2FoldChange", y = "padj",
          fills = c("#269846","#b4b4d8","#e94234"),
          colors = c("#269846","#b4b4d8","#e94234"),
          label = 'ID', label_number = 20,
          legend_title = '',
          legend_position = 'UR',
          output = FALSE)
