---
bibliography: [../ref.bib]
fignos-cleveref: True
fignos-plus-name: 图
tablenos-cleveref: True
tablenos-plus-name: 表
tablenos-caption-name: 表
---

# cMAP 药物重定位

## 1 数据与方法

### 1.1 数据描述

数据来源[@mitchell2017comparative]：
	
Mitchell KA, Zingone A, Toulabi L, Boeckelman J et al. Comparative Transcriptome Profiling Reveals Coding and Noncoding RNA Differences in NSCLC from African Americans and European Americans. Clin Cancer Res 2017 Dec 1;23(23):7412-7425. PMID: 29196495

文献URL：https://pubmed.ncbi.nlm.nih.gov/29196495/

GEO Accession：GSE101929

Tumor tissues：non-small cell lung cancer patient *T

Normal tissues：non-small cell lung cancer patient *N

这是一套非裔美国人和欧洲裔美国人非小细胞肺癌的基因表达数据。包含 66 个样本，其中 34 个非肿瘤样本，32 个肿瘤样本。用 Affymetrix 基因芯片及操作软件分析得到数据。

### 1.2 实验方法

1. 从芯片数据筛选差异表达基因，并分为上、下调两组
2. 根据 p 值排序，每组基因取前 150 个导入 Connectivity Map 查找药物

## 2 实验平台

### 2.1 硬件平台

- Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz @ 1.99GHz
- 16GB RAM + 8GB Swap

### 2.2 操作系统

- Deepin 20.2.4

### 2.3 软件工具

- RStudio 1.3.1093 with R 4.1.1 [@team2013r]
- R package: dplyr [@hadley2021dplyr]
- R package: pheatmap [@raivo2019pheatmap]
- Connectivity Map：https://clue.io/cmap [@subramanian2017next;@raivo2019pheatmap]
- DrugBank：https://go.drugbank.com/ [@wishart2018drugbank]

## 3 实验流程

### 3.1 下载数据

下载 GSE101929 的表达矩阵数据。

### 3.2 筛选差异表达基因

代码详见附录（[filterDRGs.R](###4.1-filterdrgs.r)）。p 值设置为 0.0000001，根据 logFC 筛选上调和下调基因并按 p 值升序输出。

差异表达基因热图如{@fig:heatmap}

![差异表达基因热图](https://img.limina.top/blog/差异表达基因热图.png){#fig:heatmap}

### 3.3 使用Connectivity Map查找药物

在 Connectivity Map -> Tools -> Query 如{@fig:cMAP-query}选择参数并输入各自排名前 150 的上下调基因后提交。

![GSE101929的cMAP查询](https://img.limina.top/blog/GSE101929的cMAP查询.png){#fig:cMAP-query}

查看 Detail List。

在左侧的 Perturbagen Type 下选择 Compound，因为我在定义上下调基因时将癌症组织的高表达定义为上调，cMap是对我定义的上下调基因来进行模拟和打分的，所以将表格中的 Score 升序排列，得分越负，即表示对癌症的抑制效果越好（{@fig:order_drugs}）。

![筛选排序抑癌药物](https://img.limina.top/blog/筛选排序抑癌药物.png){#fig:order_drugs}

### 3.4 调研可能的重定位药物

到 DrugBank 数据库中搜索排名前十的药物，发现排名第六的 Vorinostat 被用于治疗既往全身治疗后进行性、持续性或复发性皮肤T细胞淋巴瘤(CTCL)患者。

到谷歌学术以关键词“vorinostat lung cancer”搜索研究论文，发现两篇代表性的文章[@owonikoko2010vorinostat;@ramalingam2010carboplatin]，提到 Vorinostat 可以增加非小细胞肺癌细胞中卡铂和紫杉醇的活性。所以我认为 Vorinostat 可以重定位为非小细胞肺癌的新药。

## 4 附录

### 4.1 filterDRGs.R

```R

library(dplyr)
library(hgu133plus2.db)
library(pheatmap)

rm(list=ls())
setwd("/media/limin/Office/Study/4_Senior/Bioinformatics_Skills_Training/5_Drug-Repurposion/")

# cat -n GSE101929_series_matrix.txt|grep Sample_title > 38
type <- read.table("GSE101929_series_matrix.txt", skip = 37, nrows = 2, row.names = 1)
expmatrix <- read.delim("GSE101929_series_matrix.txt",  comment.char = '!', header = F, row.names = 1, stringsAsFactors = F)
colnames(expmatrix) <- type[1, ]

# 排序癌症和癌旁样本
data_tumor <- dplyr::select(expmatrix, ends_with("T"))
data_normal <- dplyr::select(expmatrix, ends_with("N"))
expmatrix_cluster <- cbind(data_tumor, data_normal)
colnames(expmatrix_cluster) <- expmatrix_cluster[1, ]
expmatrix_cluster <- expmatrix_cluster[-1, ]

ntumor <- ncol(data_tumor)
nnormal <- ncol(data_normal)

# 探针转GeneSymbol
probeset <- rownames(expmatrix_cluster) #提取行名为探针名
probe_symbol <- t(as.data.frame(as.list(hgu133plus2SYMBOL[probeset]), check.names = F))

symmatrix_str <- merge(probe_symbol, expmatrix_cluster, by = "row.names", all.y = T)
symmatrix_str <- na.omit(symmatrix_str)[, -1]
symmatrix <- cbind(symmatrix_str[, 1], as.data.frame(lapply(symmatrix_str[, -1], as.numeric)))
colnames(symmatrix)[1] <- "symbol"

# t检验
upgenes <- data.frame()
downgenes <- data.frame()
for (i in 1:nrow(symmatrix)) {
  if (sd(symmatrix[i, 2:(ntumor+1)]) != 0 & sd(symmatrix[i, (ntumor+2):ncol(symmatrix)]) != 0) {
    test <- t.test(symmatrix[i, 2:(ntumor+1)], symmatrix[i, (ntumor+2):ncol(symmatrix)], na.rm = TRUE)
    fold <- mean(as.matrix(symmatrix[i, 2:(ntumor+1)])) / mean(as.matrix(symmatrix[i, (ntumor+2):ncol(symmatrix)]))
    logFC <- log2(fold)
    if (test$p.value < 0.0000001 && logFC > 0) {
      upgenes <- rbind(upgenes, c(symmatrix[i, 1], test$p.value))
    }
    else if (test$p.value < 0.0000001 && logFC < 0) {
      downgenes <- rbind(downgenes, c(symmatrix[i, 1], test$p.value))
    }
  }
}

colnames(upgenes) <- c("symbol", "pvalue")
colnames(downgenes) <- c("symbol", "pvalue")

# 对p值进行BH矫正
uppadjust <- p.adjust(as.numeric(upgenes$pvalue), method = "BH", n = nrow(upgenes))
downpadjust <- p.adjust(as.numeric(downgenes$pvalue), method = "BH", n = nrow(downgenes))
upgenes$pvalue <- uppadjust
downgenes$pvalue <- downpadjust
# 排序
upgenes <- upgenes[order(as.numeric(upgenes$pvalue)), ]
downgenes <- downgenes[order(as.numeric(downgenes$pvalue)), ]
# 去重
upgenes <- upgenes[!duplicated(upgenes$symbol), ]
downgenes <- downgenes[!duplicated(downgenes$symbol), ]

# 热图数据
upgenesmatrix <- merge(upgenes[1:150, ], symmatrix, all.x = T)
downgenesmatrix <- merge(downgenes[1:150, ], symmatrix, all.x = T)
heatmapmatrix <- rbind(upgenesmatrix, downgenesmatrix)
heatmapmatrix <- heatmapmatrix[!duplicated(heatmapmatrix$symbol), ][, -2]
rownames(heatmapmatrix) <- heatmapmatrix[, 1]
heatmapmatrix <- heatmapmatrix[, -1]

pheatmap(heatmapmatrix, scale = "row", cluster_cols = TRUE, clustering_distance_rows = "correlation", margins = c(5,5), legend = TRUE,
         border = FALSE, cellheight = 9, cellwidth = 48, treeheight_row = 100, main = "DEGs_heatmap", file = 'heatmap.pdf')

write.csv(upgenes, "upgenes.csv", row.names = F)
write.csv(downgenes, "downgenes.csv", row.names = F)

```
