---
bibliography: [../ref.bib]
fignos-cleveref: True
fignos-plus-name: 图
tablenos-cleveref: True
tablenos-plus-name: 表
tablenos-caption-name: 表
---

# A549 非小细胞肺癌细胞中的 NRF2 结合位点 ChIP-Seq 生物信息学分析

## 1 数据与方法

### 1.1 数据描述

来源文献[@namani2019genome]：

Namani, Akhileshwar, et al. "Genome-wide global identification of NRF2 binding sites in A549 non-small cell lung cancer cells by ChIP-Seq reveals NRF2 regulation of genes involved in focal adhesion pathways." Aging (Albany NY) 11.24 (2019): 12600.

文献URL：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6949066/

GEO Accession：GSE141497

Experiment：NRF2_ChIPSeq SRX7268483-SRR10588628

Control：input DNA SRX7268484-SRR10588629

### 1.2 数据背景

核因子 erythroid-derived-2-like 2（NRF2）通过与其启动子区域中的抗氧化剂响应元件结合来调节其下游基因。NRF2 的过度激活导致了包括非小细胞肺癌（NSCLC）在内的各种癌症的肿瘤发生和耐药性。

文章对含 KEAP1（Kelch-like ECH-associated protein 1）的突变 A549（NSCLC的一种）细胞使用了 ChIP-Seq 技术研究了 NRF2 的全基因组结合位点，并进行了系列后续分析对 NSCLC 中 NRF2 调节的基因和途径进行了进一步研究。

### 1.3 实验方法

1. 下载作者上传的两组 SRA 数据（对照组和 NRF2 实验组）
2. 进行质量控制
3. 使用人类参考基因组 hg38 进行比对，与原文的 hg19 进行比较
4. 寻找 peaks 并注释、进行可视化
5. 分析 motif 并与原文进行比较

## 2 实验平台

### 2.1 硬件平台

- Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz @ 1.99GHz
- 16GB RAM + 8GB Swap

### 2.2 操作系统

- Deepin 20.2.3
- Deepin 20.2.4

### 2.3 软件工具

- FastQC 0.11.8[@andrews2017fastqc]
- MultiQC 1.11[@doi:10.1093/bioinformatics/btw354]
- Bowtie2 2.3.4.3[@langmead2012fast]
- MACS2 2.2.7.1[@zhang2008model]
- RStudio 1.3.1093 with R 4.1.1[@team2013r]
- ChIPseeker 1.28.3[@yu2015chipseeker]
- TxDb.Hsapiens.UCSC.hg38.knownGene 3.13.0[@team2019txdb]
- IGV Version user not_set[@thorvaldsdottir2013integrative]
- SAMtools 1.9[@li2009sequence]
- Bedtools 2.27.1[@quinlan2010bedtools]
- RSAT-peak-motifs: http://rsat.sb-roscoff.fr/peak-motifs_form.cgi [@thomas2012rsat;@thomas2012complete]

## 3 实验流程

### 3.1 软件安装

```shell
# FastQC
sudo apt-get install fastqc

# MultiQC
conda install -c bioconda multiqc

# Bowtie2
sudo apt-get install bowtie2

# MACS2
conda install macs2

# R & RStudio
略

# ChIPseeker
# 在 RStudio 中 使用 BiocManager 来安装
BiocManager::install("ChIPseeker")
# 会自动安装依赖人类参考基因组 hg19（TxDb.Hsapiens.UCSC.hg19.knownGene），如果使用 hg38，需自行安装
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# IGV
conda install igv
```

### 3.2 从 NCBI 下载 ChIP-seq reads 数据

```shell
mkdir -p GSE141497/SRA
cd GSE141497/SRA
for ((i=8; i<=9; i++));
do axel -n 100 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/02$i/SRR1058862$i/SRR1058862$i.fastq.gz;
done
```

### 3.3 质量控制及结果解读

```shell
# FastQC 进行质控
cd GSE141497/SRA
ls *gz | while read id; do fastqc -t 32 -o fastqc_out $id; done

# MultiQC 质控结果批量查看
multiqc fastqc_out
```

MultiQC 会生成一个 html 文件，从中可以获取如下信息

#### 3.3.1 General Statistics：所有样本数据基本情况统计

所有样本数据基本情况统计如{@tbl:general_statistics}

| Sample Name | % Dups | % GC | Length | % Failed | M Seqs |
| :--: | :--: | :--: | :--: | :--: | :--: |
| SRR10588628 | 47.6%  | 58%  | 51 bp  |    9%    |  30.7  |
| SRR10588629 | 27.2%  | 41%  | 51 bp  |    9%    |  29.2  |

Table:所有样本数据基本情况统计表. {#tbl:general_statistics}

- % Dups：重复 reads 比例
- M Seqs 即为 reads，单位（million）
- 平均 reads 长度51bp

#### 3.3.2 Sequence Counts：序列数量统计

![SRR10588628和SRR10588629的重复水平](https://img.limina.top/blog/SRR10588628和SRR10588629的重复水平-2021-10-04.png){#fig:sequence_counts}

由{@fig:sequence_counts}可以看出测序中的一些重复水平，这个重复水平和测序深度以及序列本身的表达情况有关。SRR10588628 样本的重复水平较高。

#### 3.3.3 Sequence Quality Histograms：每个 reads 各位置碱基的平均测序质量

![SRR10588628和SRR10588629的碱基平均质量分数](https://img.limina.top/blog/SRR10588628和SRR10588629的碱基平均质量分数-2021-10-04.png){#fig:mean_quality_scores}

从碱基平均质量分数（{@fig:mean_quality_scores}）来看，两个样本的总体质量都很好。

#### 3.3.4 Per Sequence Quality Scores：具有平均质量分数的 reads 的数量

![SRR10588628和SRR10588629具有平均质量分数的reads的数量](https://img.limina.top/blog/SRR10588628和SRR10588629具有平均质量分数的reads的数量-2021-10-04.png){#fig:per_sequence_quality_scores}

估算各颜色区域曲线下面积，可以看出低质量 reads 占整体 reads 的比例。{@fig:per_sequence_quality_scores}中大部分数值都在绿色区间，所以测序结果的质量还是比较好的。

#### 3.3.5 Per Base Sequence Content：每个 reads 各位置碱基 ATCG 的比例

![SRR10588628和SRR10588629每个reads各位置碱基的比例](https://img.limina.top/blog/SRR10588628和SRR10588629每个reads各位置碱基的比例-2021-10-04.png){#fig:per_base_sequence_content}

reads 每个位置的颜色由不同比例的4种颜色混合而成，哪一个碱基的比例大，则趋近于这个碱基所代表的颜色。正常情况下每个位置每种碱基出现的概率是相近的。如果 ATGC 在任何位置的差值大于10%，则会显示 warning（橙色），如果 ATGC 在任何位置的差值大于20%，则显示 fail（红色）。
{@fig:per_base_sequence_content}可以看出，SRR10588628 样本大约在前3bp的部分颜色非常不均匀，说明前3bp的四种碱基的比例差异比较大，这可能有过表达的序列的污染，所以该样本是 warning。

#### 3.3.6 Per Sequence GC Content：reads 的平均 GC 含量

![SRR10588628和SRR10588629的reads的平均GC含量](https://img.limina.top/blog/SRR10588628和SRR10588629的reads的平均GC含量-2021-10-04.png){#fig:per_sequence_GC_content}

正常的样本的GC含量曲线会趋近于正态分布曲线，形状接近正态但偏离理论分布的情况提示我们可能有系统偏差。偏离理论分布的 reads 超过15%时显示 warning，偏离理论分布的 reads 超过30%时显示 fail。
{@fig:per_sequence_GC_content}可以看出，SRR10588629 样本与理论分布偏离超过了15%。

#### 3.3.7 Adapter Content：接头含量

结果未发现接头含量均高于 0.1% 的样本，可能是已经去除过接头的数据。

### 3.4 从 NCBI 下载人类参考基因组

搜索 Homo sapiens

这里我下载较新的 hg38 参考基因组的 FASTA 序列（{@fig:download_hg38}）

![下载hg38参考基因组](https://img.limina.top/blog/下载hg38参考基因组-2021-10-03.png){#fig:download_hg38}

下载完成后解压（{@fig:uzip_hg38}），在 ncbi_dataset/data/GCF_000001405.39/ 文件夹下有以下文件：

- 22+X+Y+线粒体基因组的 FASTA 序列文件
- 某些染色体上的 unlocalized sequences文件：已被定为到某条染色体上，但方向或具体位置仍未确定
- 一个 unplaced sequences 文件：尚未被定位到某条染色体

![hg38参考基因组解压](https://img.limina.top/blog/hg38参考基因组解压-2021-10-04.png){#fig:uzip_hg38}

### 3.5 建立索引

不使用 unlocalized sequences 和 unplaced sequences

将22+X+Y+线粒体基因组的 FASTA 序列合并成一个文件

```shell
mkdir chrfa
for file in `ls chr?.fna chr??.fna -v`; do
        echo $file
        cat $file >> chrfa/hg38.fna
done
```

建立索引

```shell
cd chrfa
bowtie2-build hg38.fna hg38_index
```

运行约两个半小时，建立完成出来 6 个索引文件（{@fig:hg38_index}）

![hg38索引文件](https://img.limina.top/blog/hg38索引文件-2021-10-04.png){#fig:hg38_index}

### 3.6 比对到参考基因组

```shell
cd $Home/media/limin/Office/Study/4_Senior/Bioinformatics_Skills_Training/1_ChIP-Seq
# 参数说明：
# -p 32：使用32线程
# --reorder 多线程运算时，比对结果在顺序上会和文件中reads的顺序不一致，使用该选项使其一致。
# -U <fq.gz>：单末端数据比对
# -S <sam>：输出 sam 格式结果文件
# 2> <out>：屏幕回显重定向到 out 文件
for ((i=8; i<=9; i++));
do bowtie2 -p 32 --reorder 
–x Homo_sapiens_reference_genome_GRCh38.p13/ncbi_dataset/data/GCF_000001405.39/chrfa/hg38_index 
–U GSE141497/SRA/SRR1058862$i.fastq.gz 
-S bowtie2_align/SRR1058862$i.sam 
2> bowtie2_align/SRR1058862$i.out;
done
```

SRR10588628 比对了约 95mins，SRR10588629 比对了约 176mins，结果文件如{@fig:comparison_output}

![SRR10588628和SRR10588629的bowtie2比对输出文件](https://img.limina.top/blog/SRR10588628和SRR10588629的bowtie2比对输出文件-2021-10-06.png){#fig:comparison_output}

从{@fig:SRR_output_comparison}可以看到 SRR10588629.sam 比 SRR10588628.sam 文件大了约 0.7G，在两个 out 文件中也可以发现，SRR10588628 的 reads 匹配率极低，只有 6.79%；而 SRR10588629 的 reads 匹配率相对较高，达到了 52.88%。

![SRR10588628和SRR10588629的out文件对比](https://img.limina.top/blog/SRR10588628和SRR10588629的out文件对比-2021-10-05.png){#fig:SRR_output_comparison}

reads 匹配率的差异原因是 SRR10588628 是 NRF2 结合组，SRR10588629 是对照组；而两组的匹配率都比较低，可能是因为我使用了 hg38 作为参考基因组。

### 3.7 MACS2 做 peak calling

```shell
# 参数说明：
# -t <sam>：处理组文件
# -c <sam>：对照组文件
# -f <AUTO,SAM,BAM...>：文件类型
# -g <hs,mm,ce,dm>：有效基因组大小，可以选择程序自带的或自定义
# --keep-dup 1：保留一个重复，每个位置上一个 reads
# --outdir：输出文件的目录
# -n：输出文件的前缀名
# -B：以 bedGraph 格式存放 fragment pileup, control lambda, -log10pvalue 和 log10qvale
macs2 callpeak -t bowtie2_align/SRR10588628.sam -c bowtie2_align/SRR10588629.sam -f SAM -g hs --keep-dup 1 --outdir macs2_callpeak -n macs2 -B &>macs2_callpeak/MACS2.out
```

运行约 3.5mins，得到如{@fig:macs2_output}文件

![SRR10588628和SRR10588629的MACS2输出文件](https://img.limina.top/blog/SRR10588628和SRR10588629的MACS2输出文件-2021-10-06.png){#fig:macs2_output}

- .bdg：能够用 UCSC genome browser 转换成更小的 bigWig 文件；
- _model.r：可以通过$ Rscript _model.r 作图，得到基于提供数据的peak模型；
- _peaks.xls：以表格形式存放 peak 信息，和bed格式类似，但是以1为基，而bed文件是以0为基。所以 xls 的坐标都要减一才是 bed 文件的坐标；
- _summits.bed：Browser Extensible Data，记录每个 peak 的 peak summits，即记录极值点的位置。MACS 建议用该文件寻找结合位点的 motif。能够直接载入 UCSC browser，用其他软件分析时需要去掉第一行；
- MACS2.out：MACS2 运行时的信息。

查看查找出了多少 peaks：

```shell
wc -l macs2_summits.bed
```

结果：

```shell
2161 macs2_summits.bed
```

共查找出了 2161 个 peaks，与原文的 2395 个 peaks 差不多。

作图

```shell
sudo Rscript macs2_model.r 
```

### 3.8 ChIPseeker 对 peeks 注释

使用 ChIPseeker 对 peeks 进行注释

详细代码见附录（[chipseeker.R](#_41-chipseekerr)）

查看 peaks 在各个染色体上的分布（{@fig:peaks_over_chrs}）

![GSE141497的peaks在各个染色体上的分布](https://img.limina.top/blog/GSE141497的peaks在各个染色体上的分布.png){#fig:peaks_over_chrs}

查看 peaks 在所有基因的启动子附近的分布情况（{@fig:peaks_near_promoter}），这里取上下游 3kb

![GSE141497的peaks在所有基因的启动子附近的分布情况](https://img.limina.top/blog/GSE141497的peaks在所有基因的启动子附近的分布情况.png){#fig:peaks_near_promoter}

绘制信号强度曲线图（{@fig:read_count_frequency}）

![GSE141497的peaks在所有基因启动子附近分布的信号强度曲线图](https://img.limina.top/blog/GSE141497的peaks在所有基因启动子附近分布的信号强度曲线图.png){#fig:read_count_frequency}

注释结果存入 csv（{@fig:peaks_annotation}），数据展示了关联的基因以及对应的基因组区域的类别，根据此结果，可以提取关联基因进行下游的功能富集分析。

![GSE141497的peaks注释结果](https://img.limina.top/blog/GSE141497的peaks注释结果.png){#fig:peaks_annotation}

注释的结果的可视化，堆叠条形图和饼图如{@tbl:annotation_visualization}

| ![GSE141497的peaks注释结果可视化条形图](https://img.limina.top/blog/GSE141497的peaks注释结果可视化.png) | ![GSE141497的peaks注释结果的可视化饼图](https://img.limina.top/blog/GSE141497的peaks注释结果的可视化.png) |
| ------------------------------------------------------------ | ------------------------------------------------------------ |

Table:GSE141497的peaks注释结果可视化. {#tbl:annotation_visualization}

原文结果中，在所有的 peaks 中有 38.3% 位于内含子区域，12.7% 位于外显子区域，33.8% 位于基因间，0.6% 位于5'非翻译区，4.7% 位于3'非翻译区。

而我的结果显示。所有 peaks 中有 29.93% 位于内含子区域，13.55% 位于外显子区域，19.07% 位于基因间，0.69% 位于5'非翻译区，4.72% 位于3'非翻译区。

考虑到找到的 peak 数量差不多，百分比即可说明数量的差异。两结果还是比较相近的，内含子区域和基因间的 peaks 数量差异较大，这可能是由 hg38 参考基因组产生。

### 3.9 IGV 对 peaks 进行可视化

选择参考基因组为 Human(GRCh38/hg38)，导入 macs2 分析的输出结果 summits.bed 或 peaks.narrowPeak，结果展示如下，可以直观地看出 peaks 在基因组上的分布（{@fig:peaks_over_genome}）。

![GSE141497的peaks在基因组上的分布](https://img.limina.top/blog/GSE141497的peaks在基因组上的分布.png){#fig:peaks_over_genome}

### 3.10 Motif 分析

#### 3.10.1 Bedtools 转 narrowPeak 文件为 fna

```shell
# 使用 samtools 构建参考基因组的 fai 索引文件
# 此步骤不是必须的，因为 bedtools 会自动构建此文件
cd Homo_sapiens_reference_genome_GRCh38.p13/ncbi_dataset/data/GCF_000001405.39/chrfa
samtools faidx hg38.fna
# 转 narrowPeak 文件
cd ../ && mkdir motif_analysis && cd motif_analysis
bedtools getfasta -fi ../Homo_sapiens_reference_genome_GRCh38.p13/ncbi_dataset/data/GCF_000001405.39/chrfa/hg38.fna -bed ../macs2_callpeak/macs2_peaks.narrowPeak -fo macs2_peaks.fna
```

#### 3.10.2 RSAT web server 分析 motif

上传及参数设置：

- 上传 peak sequences 文件
- Reduce peak sequences 中 Cut peak sequences: +/- 设置为 0 来分析完整的数据集
- Motif Discovery parameters 中 oligomer lengths 勾选 6 和 7，勾选 Discover over-represented spaced word pairs
- Compare discovered motifs with databases 中我选择了比较新的 ENCODE-ENCODE(Human TFs)(2018-03)
- 输入邮箱后点击 GO 即可

在结果的 Sequence composition 表格中，可以看到 peaks 长度大致分布在 150 ~ 400 之间（{@fig:peaks_length_distribution}），是符合实验方法的。

![GSE141497不同长度peaks的分布](https://img.limina.top/blog/GSE141497不同长度peaks的分布.png){#fig:peaks_length_distribution}

查看 GSE141497 原始文献，文中作者使用了两种方法来根据 peaks 预测 motif。一种是使用 HOMER 软件进行已知 motif 的富集分析分析；另一种是从头预测算法。结果如{@fig:motif_in_article}。

![文献中的motif](https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=138704&s=96&r=1&c=1){#fig:motif_in_article}

{@fig:motif_in_article} A 显示的是作者使用 HOMER 对已知 motif 的富集结果，富集结果集中在 bZIP（碱性亮氨酸拉链）家族；{@fig:motif_in_article} B 显示的是作者使用从头预测算法预测 motif 的结果中排名第一的 motif，是一条包含 12bp 的一致性 NRF2 ARE 序列（ATGACTCAGCAA）的序列，在所有转录因子结合位点中占 34.47% （697/2,395）；{@fig:motif_in_article} C 显示的是作者使用 motif 比较工具 STAMP 将从头预测的 motif 与已知的 ARE motif 进行比较的结果。根据已知 motif 数据库（JASPAR）的 HOMER 查询 motif 中，排名第一的为 NRF2 TFBS，与从头预测的 NRF2 ARE 序列（TGACNNNGC）相似性最大。

使用 RSAT 分析出的前两个 motif 如{@fig:motif_analysis}

![GSE141497的motif分析](https://img.limina.top/blog/GSE141497的motif分析.png){#fig:motif_analysis}

通过对比可以发现，文献中的 NRF2 ARE 序列（TGACNNNGC）包含于 motif 1，motif 2 与 {@fig:motif_in_article} A 中的 Nrf2 motif 序列一致。

我使用人类参考基因组 hg38 虽然没有匹配到较多的 peaks，但是仍分析到了文章中已有的部分 motif，而且更加精确。

## 4 附录：R代码

### 4.1 chipseeker.R

```R

library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())
setwd("/media/limin/Office/Study/4_Senior/Bioinformatics_Skills_Training/1_ChIP-Seq/peaks_annotation/")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile("../macs2_callpeak/macs2_peaks.narrowPeak")
bed <- readPeakFile("../macs2_callpeak/macs2_summits.bed")

# ChIP Peaks over Chromosomes
png("peaks.png", width = 1300, height = 800)
covplot(peaks, weightCol="V5")
dev.off()

# Distribution of peaks near promoters of all genes

# 染色体列不是标准的chr1会报错
# newStyle <- mapSeqlevels(seqlevels(peaks), "UCSC")
# 发现UCSC和RefSeq没有映射。。需要手动改

# Fixed with:
# 排过序的才能这么干
chr <- paste0("chr", c(1:22, "X", "Y", "M"))
peakchr <- renameSeqlevels(peaks, chr)

# 定义tss上下游的距离
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peakchr, windows = promoter)

# peaks在所有基因的启动子附近的分布情况
png("Distribution of peaks near promoters of all genes.png")
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = "#aaaaff")
dev.off()

# 查看peaks在所有基因的启动子附近的分布情况，信号强度曲线图
png("The Average Profile of ChIP peaks binding to TSS region.png")
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5' -> 3')", ylab = "Read Count Frequency")
dev.off()

peakAnno <- annotatePeak(peakchr, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
peakAnno_df <- as.data.frame(peakAnno)
write.csv(peakAnno_df, "peaks_annotation.csv")

png("Pie-summarize the distribution of peaks over different type of features.png")
plotAnnoPie(peakAnno)
dev.off()

png("Bar-summarize the distribution of peaks over different type of features.png")
plotAnnoBar(peakAnno)
dev.off()

png("vennpie-summarize the distribution of peaks over different type of features.png")
vennpie(peakAnno)
dev.off()

# 查看peaks的长度分布，统计长度在1000bp以下的peaks
peaksLength = abs(peakAnno_df$end-peakAnno_df$start)
peaksLength = peaksLength[peaksLength < 500]
png("Histogram of peak length")
hist(peaksLength, breaks = 50, col = "#aaaaff", xlim = c(0, 500), xlab = "peak length", main = "Histogram of peak length")
dev.off()

```
