
library(vegan)
otu = t(read.table("AllOTU.tsv", header=T, sep="\t", quote = "", row.names=1, comment.char=""))
# print(paste0("All sample rarefaction as following"))
# rowSums(otu)

##Alpha diversity
# vegan::estimateR计算obs, chao1和ACE指数
estimateR = t(estimateR(otu))[,c(1,2,4)]
colnames(estimateR) = c("Richness", "Chao1", "ACE")

# vegan::diversity计算多样性指数shannon, simpson和invsimpson
Shannon = diversity(otu, index = "shannon")
Simpson = diversity(otu, index = "simpson")
Evenness <- Shannon / log(estimateR(otu)[1, ], exp(1))

# 合并几种指数
alpha_div = cbind(estimateR, Shannon, Simpson, Evenness)
head(alpha_div, n=1)

##保存alpha多样性指数
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("SampleID\t", file="AlphaDivIndex.tsv", append = F, quote = F, eol = "", row.names = F, col.names = F)

# 保存统计结果，有waring正常
suppressWarnings(write.table(alpha_div, file="AlphaDivIndex.tsv", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
