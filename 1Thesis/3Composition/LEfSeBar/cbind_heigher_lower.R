
H <- read.table("LEfSe_abundance_higher.tsv",sep = "\t",header = T)
L <- read.table("LEfSe_abundance_lower.tsv",sep = "\t",header = T)[,-1]
HL <- cbind(H,L)
write.table(HL,"LEfSe_abundance.tsv",sep = '\t',row.names = F)
