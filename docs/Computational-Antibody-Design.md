---
bibliography: [../ref.bib]
fignos-cleveref: True
fignos-plus-name: 图
tablenos-cleveref: True
tablenos-plus-name: 表
tablenos-caption-name: 表
---

# 计算抗体设计

## 1 数据与方法

### 1.1 数据描述

结构来源于 [@kurella2014structure]

### 1.2 实验方法

1. 使用小鼠的抗体序列建立人的同源抗体模型
2. 人源抗体模型的能量最小化优化

## 2 实验平台

### 2.1 硬件平台

- Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz @ 1.99GHz
- 16GB RAM + 8GB Swap

### 2.2 操作系统

- Deepin 20.2.4
- Windows 11

### 2.3 软件工具

- RCSB PDB: https://www.rcsb.org/ [@berman2000protein]
- SAbPred-ABodyBuilder: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/ [@leem2016abodybuilder]
- IMGT/DomainGapAlign: http://imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cg [@ehrenmann2010imgt;@ehrenmann2011imgt]
- Pymol [@PyMOL]
- Swiss-PDBViewer [@guex1997swiss]

## 3 实验流程

### 3.1 获取小鼠抗体的 FASTA 序列

我分配到了 3L7E 这一结构，在 RCSB PDB 数据库中搜索，这是一个小鼠的抗 IL-13 抗体 C836 的晶体结构（{@fig:3L7E.fasta}），下载其的 FASTA 序列。

![下载3L7E的FASTA序列](https://img.limina.top/blog/下载3L7E的FASTA序列.png){#fig:3L7E.fasta}

### 3.2 建立人的同源抗体模型

使用 SAbPred 的 ABodyBuilder 工具进行建模

将重链和轻链序列输入 ABodyBuilder，使用默认参数建模（{@fig:abodybuilder}）

![3L7E的ABodyBuilder建模](https://img.limina.top/blog/3L7E的ABodyBuilder建模.png){#fig:abodybuilder}

建模完成后下载所有的结果文件（{@fig:download_model}）

![下载ABodyBuilder的建模结果](https://img.limina.top/blog/下载ABodyBuilder的建模结果.png){#fig:download_model}

### 3.3 替换小鼠抗体蛋白框架

在 IMGT 上传序列并选择人类，获得小鼠抗体重链和轻链序列中与人最相同的序列（{@fig:IMGT}）。

![IMGT上传序列并选择人类](https://img.limina.top/blog/IMGT上传序列并选择人类.png){#fig:IMGT}

结果显示了人类抗体框架中与小鼠框架不同的残基，如{@fig:framework_difference}。

![人类与小鼠抗体框架的残基差异](https://img.limina.top/blog/人类与小鼠抗体框架的残基差异.png){#fig:framework_difference}

将残基变化整理下来，另存为 mutate.csv，如{@fig:replace_csv}。

![框架残基变化](https://img.limina.top/blog/框架残基变化.png){#fig:replace_csv}

使用 Pymol 进行残基的突变（替换）

1. 图形界面操作可用：Wizard -> Mutagenesis -> Protein
2. 由于要替换的残基过多（29+25=54），使用脚本操作。
该突变命令脚本来源于 Autolife's Blog[^1]，基于其我稍作了修改，可以让 Pymol 自动选择空间构象最好的残基位置。

[^1]:http://pymol.chenzhaoqiang.com/intro/advanceManual.html，下载地址为：http://pymol.chenzhaoqiang.com/_downloads/mutate.py

```python

#refer: https://pymolwiki.org/index.php/Rotkit

from pymol import cmd

def mutate(molecule, chain, resi, target):
    target = target.upper()
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    cmd.get_wizard().set_mode("%s" % target)
    selection = "/%s//%s/%s" % (molecule, chain, resi)
    cmd.get_wizard().do_select(selection)
    cmd.get_wizard().apply()
    # cmd.set_wizard("done")
    cmd.set_wizard()
    # cmd.refresh()

cmd.extend("mutate", mutate)

```

再写一个小脚本输出突变所有残基的命令，存入 mutate.pml。

```python
import pandas as pd

data = pd.read_csv("mutate.csv")
with open('mutate.pml', 'w') as f:
    for index, row in data.iterrows():
        # print("mutate 3L7E_fv_model, chain = row"+['chain'], ", resi = row"+['resi'], "target = row"+['target'])
        f.writelines("mutate 3L7E_fv_model, chain = " + row['chain'] + " , resi = " + str(row['resi']) + ", target = " + row['target'] + "\n")
        
```

在 Pymol 命令行运行以下命令即可。

```shell
cd /media/limin/Office/Study/4_Senior/Bioinformatics_Skills_Training/4_Computational-Antibody-Design
run mutate.py
run mutate.pml
```

{@fig:mutated}为第一次突变完成的结构

![第一次突变完成的结构](https://img.limina.top/blog/第一次突变完成的结构.png){#fig:mutated}

### 3.4 人源抗体模型的能量最小化

通过 DeepView-Swiss-PDBViewer 查看突变后结构的能量。

Tools -> Compute Energy(Force Field) 计算能量

输出的能量文件（部分）如{@fig:energy1}所示，红色标出的残基即表示能量太高不符合要求。

![残基能量（部分）](https://img.limina.top/blog/残基能量（部分）.png){#fig:energy1}

{@tbl:residues_change}整理了需要优化的残基与原结构对应残基变化。

|              | Chain H | Chain H | Chain H | Chain L | Chain L | Chain L | Chain L | Chain L |
| ------------ | ------- | ------- | ------- | ------- | ------- | ------- | ------- | ------- |
| 经突变的残基 | LYS80   | LYS84   | VAL113  | THR20   | THR88   | GLN105  | LEU107  | TYR116  |
| 原结构的残基 | LYS80   | GLY84   | VAL113  | THR20   | THR88   | GLN105  | HIS107  | TYR116  |

Table:需要优化的残基与原结构对应残基变化 {#tbl:residues_change}

{@fig:regions}展示了这八个残基所在的三个区域

![八个残基所在的三个区域](https://img.limina.top/blog/八个残基所在的三个区域.png){#fig:regions}

1. Chain H 上的 LYS84 与周围残基存在较大空间位阻（{@fig:H-LYS84} A），经尝试最后将 LYS84 突变为最优结构下空间位阻最小的 SER（{@fig:H-LYS84->SER} B）。

![H-LYS84->SER](https://img.limina.top/blog/HLYS84->SER.png){#fig:H-LYS84}

2. Chain H 上的 VAL113 和 Chain L 上的 TYR116 之间存在较大的空间位阻（{@fig:L-TYR116} A）。经尝试决定突变 Chain L 上的 TYR116 为 SER，因为 SER 与 TYR 都是极性不带电氨基酸，都含有羟基，且发现突变后无空间位阻（{@fig:L-TYR116} B）。

![L-TYR116->SER](https://img.limina.top/blog/L-TYR116->SER.png){#fig:L-TYR116}

3. Chain L 上的 LEU107 与周围残基都有较大的空间位阻（{@fig:L-LEU107} A），此残基是已经经过突变的，原残基为 HIS，是极性带正电氨基酸，经尝试三种极性带正电氨基酸的最优构象的空间位阻都较大，所以将其突变为极性不带电的 SER，几乎无空间位阻（{@fig:L-LEU107} B）。

![L-LEU107->SER->ASP](https://img.limina.top/blog/L-LEU107->SER->ASP.png){#fig:L-LEU107}

但是当再使用 DeepView-Swiss-PDBViewer 查看突变后结构的能量时发现此处的能量还是很大。所以将 LEU107 突变为了另一个最优构象下空间位阻最小的极性带负电氨基酸 ASP，然后将 GLN105 突变成无空间位阻的极性不带电的 SER（{@fig:L-LEU107} C）。

4. Chain L 上的 THR88 与临近残基有空间位阻（{@fig:L-THR88} A），将其突变为同为极性不带电且都含有羟基的 SER（{@fig:L-THR88} B）。

![L-THR88->SER](https://img.limina.top/blog/L-THR88->SER.png){#fig:L-THR88}

突变完成后再使用 DeepView-Swiss-PDBViewer 查看结构的能量，发现已经没有标红的残基了，而且整个蛋白的能量为负，符合要求（{@fig:energy2}）。

![第二次突变后的结构能量](https://img.limina.top/blog/第二次突变后的结构能量.png){#fig:energy2}

