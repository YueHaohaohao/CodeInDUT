#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
library(Hmisc)

Net <- function(Grp){data <- read.delim(paste0('MajorFeature_50%_0.01%_',deparse1(substitute(Grp)),'.tsv'), row.name = 1, check.names = FALSE)
        #以属水平丰度为例，“genus_table.txt” 
        genus <- data[,-ncol(data)]
        #计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
        genus_corr <- rcorr(t(genus), type = 'spearman')
        #阈值筛选
        #将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
        r <- genus_corr$r
        r[abs(r) < 0.4] <- 0
        #选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
        p <- genus_corr$P
        # p <- p.adjust(p, method = 'fdr')    #可选 p 值校正，这里使用 BH 法校正 p 值
        p[p>=0.01] <- -1
        p[p<0.01 & p>=0] <- 1
        p[p==-1] <- 0
        #根据上述筛选的 r 值和 p 值保留数据
        z <- r * p
        diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
        #如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
        # write.table(data.frame(z, check.names = FALSE), 'species_corr_CK.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
        ##获得网络
        library(igraph)
        #将邻接矩阵转化为 igraph 网络的邻接列表
        #构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
        g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
        g
        #自相关也可以通过该式去除
        g <- simplify(g)
        
        #孤立节点的删除（删除度为 0 的节点）
        g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
        
        #Zipi
        library(ggClusterNet)
        library(ggplot2)
        res = ZiPiPlot(igraph = g,method = "cluster_fast_greedy")
        plotzipi <- res[[1]]
        plotzipi
        ggsave(filename = paste0(deparse1(substitute(Grp)),'_zipi.svg'),plotzipi,width = 5,height = 4,units = 'in')
        netproperties <- net_properties.2(g,n.hub = T) %>%
          cbind(rownames(.),.)
        colnames(netproperties) <- c("properties",deparse1(substitute(Grp)))
        write.table(netproperties,paste0(deparse1(substitute(Grp)),'_netproperties.tsv'),sep="\t",row.names = F)

        
        
        
        #该模式下，边权重代表了相关系数
        #由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
        E(g)$correlation <- E(g)$weight
        E(g)$weight <- abs(E(g)$weight)
        E(g)$corr <- sign(E(g)$correlation)
        #为节点（微生物属）添加属性信息（界门纲目科属水平注释）
        #“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
        library(tidyr)
        tax <- separate(data["taxonomy"],1,into = c("K","P","C","O","F","G"),sep="[;]")
        tax <- tax[as.character(V(g)$name), ]
        V(g)$phylum <- tax$P
        V(g)$class <- tax$C
        V(g)$order <- tax$O
        V(g)$family <- tax$`F`
        V(g)$genus <- tax$G
        #查看网络图
        g
        plot(g)
        ##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
        #邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
        #还可以由 igraph 的邻接列表转换
        adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
        # write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix_CK.txt', col.names = NA, sep = '\t', quote = FALSE)
        #边列表
        edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表
        edge_list <- data.frame(
            source = edge[[1]],
            target = edge[[2]],
            correlation = E(g)$correlation,
            weight = E(g)$weight,
            corr = E(g)$corr
        )
        head(edge_list)
        write.table(edge_list, paste0('network.edge_list_',deparse1(substitute(Grp)),'.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
        #节点属性列表
        node_list <- data.frame(
            label = names(V(g)),
            phylum = V(g)$phylum,
            class = V(g)$class,
            order = V(g)$order,
            family = V(g)$family,
            genus = V(g)$genus
        )
        head(node_list)
        write.table(node_list, paste('network.node_list_',deparse1(substitute(Grp)),'.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
        #边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
        #此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
        #graphml 格式，可使用 gephi 软件打开并进行可视化编辑
        write.graph(g,paste0(deparse1(substitute(Grp)),'_network.graphml'), format = 'graphml')
        # write.graph(g,paste0(deparse1(substitute(Grp)),'_network.gml'), format = 'gml')
}
Net(CK)
Net(A)
Net(DA)
Net(AD)


