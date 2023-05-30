#输入相对丰度表
alldata <- read.delim("Family.tsv",sep = "\t",header = T,row.names = 1)

library(dplyr)
#data选择"CK","^A(\\d+)","DA","AD"
CK <- select(alldata,matches("CK|taxonomy"))
A <- select(alldata,matches("^A(\\d+)|taxonomy"))
DA <- select(alldata,matches("DA|taxonomy"))
AD <- select(alldata,matches("AD|taxonomy"))

Major <- function(data){
  data1 <- data[,-ncol(data)]
  #计算每行0的个数
  N0 <- apply(data1 ==0,1,sum)
  #判断每行的平均值
  N1 <- rowMeans(data1)>=0.0001
  #提取0值大于80%的行号
  R0 <- which(N0 <= ncol(data1)*0.5)
  #提取平均相对丰度大于0.1%的行号
  R1 <- which(N1 == T)
  #得到行号交集
  R <- c(intersect(R0,R1))
  #得到core community
  major <- data[R,]
  write.table("Taxonomy\t", file=paste0("MajorFeature_50%_0.01%_",deparse(substitute(data)),".tsv"), append = F, quote = F, eol = "", row.names = F, col.names = F)
  suppressWarnings(write.table(major,paste0("MajorFeature_50%_0.01%_",deparse(substitute(data)),".tsv" ),quote = F, sep="\t", row.names = T, col.names = T,append = T))
  
  return(major)
}

Major(CK)
Major(A)
Major(DA)
Major(AD)
