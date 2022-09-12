---
bibliography: [../ref.bib]
fignos-cleveref: True
fignos-plus-name: 图
tablenos-cleveref: True
tablenos-plus-name: 表
tablenos-caption-name: 表
---

# 基于 Pharmmaker 的肺腺癌相关靶标蛋白的药效团模型构建

## 1 数据与方法

### 1.1 数据描述

- 原始文献[@xu2020integrative]：

  Xu, Jun-Yu, et al. "Integrative proteomic characterization of human lung adenocarcinoma." Cell 182.1 (2020): 245-261.

  文献URL：https://www.sciencedirect.com/science/article/pii/S0092867420306760

- 朱宇等[^1]同学的前期工作筛选出的744个差异蛋白网络（coreNet.txt）

[^1]:成员包括：朱宇、吴代、杨彬、刘雨、赵倩、王怡欣、陈彦洁、努尔沙依达·西依提


### 1.2 实验方法

1. 文献调研，获取肺腺癌的组学数据并筛选差异蛋白
2. 分离已知靶标蛋白与非靶标蛋白
3. 富集分析，寻找肺腺癌相关重要通路与功能的参与蛋白
4. 成药性模拟及构建药效团模型

## 2 实验平台

### 2.1 硬件平台

- Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz @ 1.99GHz
- 16GB RAM + 18GB Swap

### 2.2 操作系统

- Windows 11 Pro 21H2 22000.65 ~ 22000.194
- Deepin 20.2.3 ~ 20.3.4

### 2.3 软件工具

- Anaconda3 2020.11 with Python 3.8.5 [@van1991interactively]
- Rstudio 1.4.1103 with R 4.0.3 [@team2013r]
- Cytoscape 3.8.2 [@shannon2003cytoscape]
- Venny 2.1: https://bioinfogp.cnb.csic.es/tools/venny/ [@oliveros2018interactive]
- DrugBank: https://go.drugbank.com/ [@wishart2018drugbank]
- Therapeutic Target Database (TTD): http://db.idrblab.net/ttd/ [@wang2020therapeutic]
- KEGG: Kyoto Encyclopedia of Genes and Genomes：https://www.genome.jp/kegg/ [@kanehisa2000kegg]
- Uniprot: https://www.uniprot.org/ [@uniprot2021uniprot]
- AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/ [@jumper2021highly]
- PyMOL 2.2.0 [@PyMOL]
- CavityPlus: http://www.pkumdl.cn:8000/cavityplus/index.php [@yuan2013binding;@yuan2011ligbuilder;@xu2018cavityplus]
- DruGUI: http://prody.csb.pitt.edu/tutorials/drugui_tutorial/ [@bakan2012druggability]
- Pharmmaker: http://prody.csb.pitt.edu/pharmmaker/ [@lee2020pharmmaker]
- NAMD 2.14 [@phillips2020scalable]
- VMD 1.9.3 OpenGL [@humphrey1996vmd]
- WPS 2019 11.1.0.10702-release 正式版

## 3 实验流程

### 3.1 文献调研及差异蛋白获取

1. 阅读并了解原始文献的工作
   该文献由来自中科院上海药物所谭敏佳研究员、中国医学科学院肖汀教授、北京蛋白质组研究中心贺福初院士、汪宜研究员、上海交通大学李婧教授等合作，在 Cell 上发表。对 103 例肺腺癌的样品进行了包括蛋白质组学、磷酸化修饰组学、转录组学和全外显子测序的整合性分析。通过基于蛋白质组层面的分子特征，可以将肺腺癌分为 3 个亚型（S-I, S-II, and S-III），并具有不同的临床预后。除此之外，研究团队还发现了新的潜在药物靶点，以及血液中的 HSP 90β 可以作为肺腺癌的潜在预后标志物。该研究从机理、诊断和治疗方面，都提供了新的重要的信息。在这篇文章的支撑材料的 S4 和 S5 表中，提供了完整的蛋白质组学和磷酸组学数据，且这些数据已经经过了一定的处理，将被用于我们的实验。
2. 差异蛋白由朱宇等同学的工作获得
   差异表达蛋白相互作用网络文件：coreNet.txt，包含 744 个差异蛋白。

### 3.2 已知靶标的识别与分离

#### 3.2.1 DrugBank 数据库的下载和解析

1. 注册账号并通过审核后下载 xml 文件
   xml 文件压缩包为 140MB，目前的版本为 5.1.8，解压后约 1.38GB（{@fig:drugbank_xml}）。

![DrugBank数据库xml文件下载](https://img.limina.top/blog/DrugBank数据库xml文件下载-2021-10-12.png){#fig:drugbank_xml}

2. 处理 DrugBank 的数据
   解析代码使用 python 编写，基于 Github 上的开发者 Deshan-Zhou 贡献的代码修改，仓库地址：https://github.com/Deshan-Zhou/deal_DrugBank 。完整代码见附录 1：([deal_DrugBank.py](#deal_DrugBank.py))。使用该代码获取了DrugID，DrugName，TargetID，GeneName，Indication 以及他们之间的映射关系。
   编写 R 脚本将上述映射关系结合，获得DrugBank中所有靶标-药物-适应症的映射列表。完整代码见附录 1：([DrugBank_TDI.R](#DrugBank_TDI.R))
   获得的 DrugBank 中所有靶标-药物-适应症的映射列表如{@fig:drugbank_TDI}。

![DrugBank所有靶标-药物-适应症的映射列表](https://img.limina.top/blog/DrugBank所有靶标-药物-适应症的映射列表.png){#fig:drugbank_TDI}

#### 3.2.2 TTD 数据库的下载和解析

1. TTD 数据库可以直接下载数据
   如{@fig:drugbank_TDI}下载需要的 P1-01-TTD_target_download.txt，P1-05-Drug_disease.txt 和 P1-06-Target_disease.txt 文件[^2]。

![TTD数据库的数据下载](https://img.limina.top/blog/TTD数据库的数据下载.png){#fig:drugbank_TDI}

2. 处理 TTD 的数据
   编写 R 脚本处理 TTD 下载的数据。获得靶标与药物的映射关系 Target_Drug.csv([Target_Drug.R](#Target_Drug.R))、靶标和疾病与适应症的映射关系 Target_Indication.csv([Target_or_Drug_Indication.R]{#Target_or_Drug_Indication.R})。将上述两个 csv 文件结合，处理获得 TTD 中所有靶标-药物-适应症的映射列表。完整代码见附录 1：([TTD_TDI.R](#TTD_TDI.R))。
   获得的TTD中所有靶标-药物-适应症的映射列表如{@fig:ttd_TDI}。

![TTD所有靶标-药物-适应症的映射列表](https://img.limina.top/blog/TTD所有靶标-药物-适应症的映射列表.png){#fig:ttd_TDI}

[^2]:没有使用 Download 中提供的药物-靶标映射文件 P1-07-Drug-TargetMapping.xlsx，因为其中的映射关系都是来自于文献，信息不全。

#### 3.2.3 匹配靶标蛋白

编写 R 脚本，结合两数据库处理完成的靶标-药物-适应症的映射列表，匹配 PPI 网络中的蛋白哪些是已知的靶标蛋白[^3]。输出四个文件：PPI 网络中的靶标蛋白列表(coretarget.csv)、非靶标列表(non-coretarget.csv)、DrugBank 中匹配到的靶标-药物-适应症映射信息(DrugBank_Target.csv)和 TTD 中匹配到的靶标-药物-适应症映射信息(TTD_Target.csv)。完整代码见附录 1：([target_match.R](#target_match.R))

[^3]:此处的靶标蛋白是指数据库中记录的所有靶标蛋白，并非肺腺癌特异性的靶标蛋白，肺腺癌的特异靶标将在后面的分析中提取出。

靶标在两个数据库中的分布情况如{@fig:target_distribution}

![靶标在两个数据库中的分布情况](https://img.limina.top/blog/靶标在两个数据库中的分布情况.png){#fig:target_distribution}

### 3.3 富集分析

将所有蛋白进行 GO 功能富集分析和 KEGG 通路富集分析，完整代码见附录 1([enrichment.R](#enrichment.R))

GO 功能富集分析结果如{@fig:GO}，CC、MF、BP 各展示了排名前 10 的富集。

![744个蛋白的GO功能富集分析](https://img.limina.top/blog/744个蛋白的GO功能富集分析.png){#fig:GO}

KEGG 通路富集分析结果如{@fig:KEGG}，图中展示了 p 值前 20 的通路

![744个蛋白的KEGG通路富集分析](https://img.limina.top/blog/744个蛋白的KEGG通路富集分析.png){#fig:KEGG}

但是这些通路不是我所关心的。到 KEGG 的疾病数据库（KEGG DISEASE）中搜索 lung cancer，找到两条通路（{@fig:pathway}）：hsa05222 （小细胞肺癌）和 hsa05223 （非小细胞肺癌）

![肺癌特异性通路](https://img.limina.top/blog/肺癌特异性通路.png){#fig:pathway}

提取 744 个差异蛋白中富集到此两条通路上的蛋白，并分成已知靶标和潜在靶标两类。

在 Cytoscape 中进行可视化，把肺腺癌相关的靶标与潜在的靶标标记出来，如{@fig:network}，其中红色圆形为已知靶标，紫色三角形为潜在靶标。

![744个差异蛋白中在肺腺癌特异通路中的已知靶标与潜在靶标的分布](https://img.limina.top/blog/744个差异蛋白中在肺腺癌特异通路中的已知靶标与潜在靶标的分布.png){#fig:network}

生成子网络并重新排布一下（{@fig:new_network}）

![重排的已知靶标和潜在靶标分布](https://img.limina.top/blog/重排的已知靶标和潜在靶标分布.png){#fig:new_network}

可以看到 ITGA2B 与三个已知靶标有相互作用，所以我非常草率地选择它进行后续研究。

### 3.4 获得 ITGA2B 的结构

到 Uniprot 数据库搜索 ITGA2B，结果中的第一条即为人类的 ITGA2B，该基因编码的蛋白是整合素 α-IIb，Uniprot ID 为 P08514。

点击 Structure 查看已知的结构，该基因碱基序列共有 1039 位，已有的 pdb 结构中 6V4P（1-963）相对比较完整，所以决定优先选择其进行实验（{@fig:uniprot_6v4p}）。

![ITGA2B的pdb结构](https://img.limina.top/blog/ITGA2B的pdb结构.png){#fig:uniprot_6v4p}

到 RCSB.PDB 数据库下载 6v4p.pdb，在 PyMOL 中打开查看（{@fig:pymol_6v4p}）

![在pymol中查看6v4p](https://img.limina.top/blog/在pymol中查看6v4p.png){#fig:pymol_6v4p}

该结构有 4 条链，其中 C、D 链（图中黄色和紫色）为抗体，删除抗体和杂原子（Ca^2+^、Mg^2+^等），另存为一份纯蛋白文件 p6v4p.pdb 进行实验。

### 3.5 DruGUI 的使用

#### 3.5.1 DruGUI: VMD 插件的安装

DruGUI 当前版本为 1.1[^4]

[^4]:截止2021年10月14日，可能是作者没有更新好插件包，需要修改 drugui/pkgIndex.tcl 文件的第 11 行，将 1.0 更改为 1.1，或将 drugui/drugui.tcl 文件的第 26 行的 1.1 改为 1.0，才能在 VMD 中正常加载此插件。

```shell
conda install prody
conda install vmd

# drugui 插件下载
http://prody.csb.pitt.edu/tutorials/drugui_tutorial/drugui_plugin_files.zip

# 在 VMD 命令行输入以下命令获取 VMD 的插件目录
global env; puts $env(VMDDIR)
# 复制插件到 VMD 的插件目录
cp -r drugui /home/limin/anaconda3/lib/plugins/noarch/tcl
# 修改 loadplugins.tcl 文件
sudo gedit /home/limin/anaconda3/lib/scripts/vmd/loadplugins.tcl
# 在 186 行添加
vmd_install_extension drugui drugui_tk "Modeling/DruGUI"
```

#### 3.5.2 namd2 的安装

```shell
# 下载地址：https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD
wget https://www.ks.uiuc.edu/Research/namd/2.14/download/946183/NAMD_2.14_Linux-x86_64-multicore.tar.gz
tar -zxvf NAMD_2.14_Linux-x86_64-multicore-CUDA.tar.gz -C /home/limin/
# 在 ~/.bashrc 文件中写入
export PATH=~/namd2:$PATH
```

#### 3.5.3 建立 psf 文件

下载教程文件：

- Unix/Mac：
  - 教程PDF：http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix.pdf
  - 需要的文件：http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-files.tar.gz
- Windows：
  - 教程PDF：http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win.pdf
  - 需要的文件：http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-files.zip
- *其中需要的文件为 top_all27_prot_lipid.inp*

根据教程建立 psf 文件。

#### 3.5.4 使用 DruGUI 插件

点击 Extensions -> Modeling -> DruGUI，加载获得的 pdb 文件和 psf 文件，然后使用默认参数进行体系准备（{@fig:prepare_system}）。

![Druggability参数设置](https://img.limina.top/blog/Druggability参数设置.png){#fig:prepare_system}

使用 DruGUI 进行体系准备后生成如{@fig:prepare_system_out}文件。

![6v4p成药性模拟的输出文件](https://img.limina.top/blog/6v4p成药性模拟的输出文件.png){#fig:prepare_system_out}

系统在 VMD 中的展示如下，体系比较大，蛋白有很多缺失（{@fig:fail_system}）。

![6v4p系统在vmd中的展示](https://img.limina.top/blog/6v4p系统在vmd中的展示.png){#fig:fail_system}

运行文件夹下的`.sh`文件产生报错。查看 ubq 与 p6v4p 的 pdb 结构对比（{@fig:compare_6v4p}）。

![ubp与p6v4p的pdb结构对比](https://img.limina.top/blog/ubp与p6v4p的pdb结构对比.png){#fig:compare_6v4p}

可见 A 链几乎没有构建出来，所以删除 A 链，只使用 B 链（蓝色）进行实验。

重复步骤 3.4.2 ~ 3.4.3，{@fig:new_system}为新体系的展示。

![6v4p的B链系统展示](https://img.limina.top/blog/6v4p的B链系统展示.png){#fig:new_system}

再运行新的`.sh`文件，可正常运行。两天三夜后发现 namd 的默认步数参数似乎太大了，所以手动终止。查看 sim.log 文件，已经跑了 2680000 步，此时 sim 文件夹中有如{@fig:namd_out}文件[^5]。

![sim文件夹下的文件](https://img.limina.top/blog/sim文件夹下的文件.png){#fig:namd_out}

[^5]:输出文件的前缀名打错了，应为 p6v4pb，但是影响不大。

因为是强制终止，为了保持完整性，再跑 2000 步。修改 simrestart.conf 文件，第 16 行的 firsttimestep 改为 2682000，最后一行的 run 改为 2000，终端运行`namd2 simrestart.conf > simrestart.log`生成 sim1.dcd 文件。

#### 3.5.5 网格计算

点击 Extensions -> Modeling -> DruGUI，选择 Calculate Grids，如{@fig:calculate_grids}选择 psf、pdb 和上一步生成的 sim.dcd 文件和 sim1.dcd 文件。

![calculate_grids参数设置](https://img.limina.top/blog/calculate_grids参数设置.png){#fig:calculate_grids}

输出文件如{@fig:calculate_grids_out}所示[^6]

[^6]:Druggability options and parameters 为可选参数，如果设置了则会同时进行下一步的成药性分析（Assess Druggbility）

![grid_calculation的输出文件](https://img.limina.top/blog/grid_calculation的输出文件.png){#fig:calculate_grids_out}

#### 3.5.6 成药性分析

如果网格计算中已经设置了 Druggability options and parameters，此步骤可跳过。

点击 Extensions -> Modeling -> DruGUI，选择 Assess Druggability，选择上一步生成的四个 dx 文件[^7]，如{@fig:druggbility_analysis}。

[^7]:此步骤需要四个dx文件和上一步生成的heavyatoms.pdb都在输出文件夹中

![druggbility_analysis参数设置](https://img.limina.top/blog/druggbility_analysis参数设置.png){#fig:druggbility_analysis}

### 3.6 Pharmmaker 的使用

将 DruGUI 的结果文件夹更名为 drugui-simulation 和 drugui-analysis，新建 pharmmaker 的分析文件夹。Pharmmaker 教程[^8]可以在此处获得：http://prody.csb.pitt.edu/tutorials/pharmmaker/

[^8]:这个教程写得是真的离谱，debug详见附录。

#### 3.6.1 缩小分析残基范围

在 CavityPlus 的 Cavity 模块中 上传该结构 B 链的 pdb 文件来预测潜在的药物结合口袋（{@fig:cavityplus_B}）。

![在CavityPlus上传6v4p的B链](https://img.limina.top/blog/在CavityPlus上传6v4p的B链.png){#fig:cavityplus_B}

预测完成后对结果文件 outputcavity.txt 分析，绘制口袋残基的分布图[^9]（{@fig:pocket_resi_distribution}），代码见附录：([residues.R](#residues.R))。

[^9]:这张图上点的横向散布并没有什么意义，只是因为画在一条线上根本看不清。作此图的目的也只是为了看一下组成口袋的残基分布来缩小后续分析的残基范围。

![6v4pb预测口袋残基分布](https://img.limina.top/blog/6v4pb预测口袋残基分布){#fig:pocket_resi_distribution}

没有明显的大段非口袋残基位点，最后将残基范围定为 65 ~ 426。

#### 3.6.2 高亲和力残基分析

先挑选了一小段序列，使用探针 IPRO 测试，发现残基亲和力都较低，所以将 cutoff 定为 50。

```shell
env VMDARGS='text with blanks' vmd -dispdev text -e ../pharmmaker/CORE/highaffresid.tcl -args ../drugui-simulation/p6v4qb.pdb ../drugui-simulation/p6v4qb_sim/sim.dcd U IPRO,ACET,IPAM,ACAM 65 426 50
```
生成各探针的残基亲和力文件 out-chain-probe.dat 和残基编号文件 out-chain-probe-highaffresid.dat等。

合并各探针的结果文件 out-chain-probe.dat 并用 WPS 作{@fig:affinity}。

![6v4pb残基对四种探针的亲和力](https://img.limina.top/blog/6v4pb残基对四种探针的亲和力.png){#fig:affinity}

#### 3.6.3 高亲和力残基附近的热点

教程有误，上一步默认输出的 probe-list.dat 文件中的探针名是一行一个，回车换行的，而下述步骤需要手动修改 probe-list.dat 文件，将探针写在一行，并用空格隔开。

这里我另存为 probelist.dat，然后在终端运行以下代码。

```shell
../pharmmaker/run_hotspotsNearHighAffResids.sh chain-list.dat probelist.dat 8 ../drugui-analysis/grid_calculation/p6v4pb/p6v4pb_all_hotspots.pdb ../drugui-analysis/grid_calculation/p6v4pb_heavyatoms.pdb .
```
生成的 highAffHotspots.pdb 文件是对成药性分析的输出文件 all_hotspots.pdb 的抽取。根据 all_hotspots.pdb 文件注释可知，软件对亲和力的定义为自由能（第 11 列），并在其中升序排列。自由能越低，残基亲和力越高。

#### 3.6.4 关注感兴趣区域的热点

此部分需要对蛋白的背景知识有一定的理解，我没有，所以选择粗暴地排序亲和力来选择对几种探针都有高亲和力的残基区域作为我感兴趣的区域。

查看各探针的亲和位点文件，取一个合理的 cutoff=-1.3，筛选出文件 highAffHotspots.pdb 中自由能小于 -1.5 的探针，另存为 highAffHotspots_filter.pdb。

```python
f = open("highAffHotspots.pdb", "r")
fn = open("highAffHotspots_filter.pdb", "w")
for i in f:
    b_factor = i[61:66].strip()
    if eval(b_factor) <= -1.3:
        fn.writelines(i)
f.close()
fn.close()
```

到 PyMOL 中与 6v4p 的 B 链进行叠合。按照探针种类上色，IPRO 设置为橙色，IPAM 设置为红色，ACET 设置为蓝色。

```python
# 选择探针
select IPRO,resn IPRO
select IPAM,resn IPAM
select ACET,resn ACET
```

可以发现在红圈区域三种探针都有亲和性高的残基（{@fig:red_circle}）。

![高亲和力探针区域](https://img.limina.top/blog/高亲和力探针区域.png){#fig:red_circle}

查看高亲和力残基与探针分子的相对位置

```python
# 选择高亲和力的残基
select IPRO_resi,resi 93+125+126+129+130+137+141+173+181+182+183+187+211+212+234+261+270+272+273+274+277+278+280+282+301+312+313+316+319+320+323+338+339+342+346+360+361+362+364+387+388+406
select IPAM_resi,resi 66+71+76+91+126+129+323+336+364+365+378+406
select ACET_resi,resi 98+137+143+208+350+354+384
```

根据 highaffresid.dat 文件将探针对应的残基设置为对应颜色，即可确定此处的残基号（{@fig:probe-highaffresid}）。

![探针对应高亲和残基](https://img.limina.top/blog/探针对应高亲和残基.png){#fig:probe-highaffresid}

将此区域的高亲和热点残基和探针分别写入 interest.dat 文件和 interest.pdb 文件（或者由之前的文件删除得到）。

文件示例如下

```text
<!-- 例：highaffresid_interest.dat -->

IPRO U 312 313 316 319 320 338 339 342 346
ACET U 350 354
IPAM U 323 336

<!-- 例：highAffHotspots_interest.pdb -->
ATOM      1  C2  ACETA   9       5.000  15.500   1.000  0.99 -1.72           M  
ATOM      2  C2  ACETA   9       5.000  13.000  -1.000  1.00 -1.67           M  
ATOM      3  C2  ACETA   9       0.500  11.000   4.500  0.93 -1.44           M  
ATOM      4  C2  ACETA   9       4.000  16.500   3.000  1.00 -1.31           M  
ATOM      5  C2  ACETA   9       0.000  14.000   0.500  0.99 -1.70           M  
ATOM      6  C2  ACETA   9       4.500  15.500  -1.500  0.97 -1.45           M  
ATOM      7  C2  ACETA   9       0.500  15.500   2.500  1.00 -1.38           M  
ATOM      8  C2  IPAMA   9       1.500  15.500  -2.500  1.00 -1.50           M  
ATOM      9  C2  IPROA   9       7.000  19.000 -10.000  1.00 -1.73           M  
ATOM     10  C2  IPROA   9       1.000  21.000  -8.500  1.00 -1.52           M  
ATOM     11  C2  IPROA   9      -1.500  19.500  -6.000  0.97 -1.38           M  
ATOM     12  C2  IPROA   9       0.000  17.000  -4.000  0.66 -2.06           M  
ATOM     13  C2  IPROA   9      11.000  15.500  -5.500  1.00 -1.57           M  
ATOM     14  C2  IPROA   9       8.500  15.000  -5.000  1.00 -1.90           M  

```

#### 3.6.5 热点和带有探针的高亲和力残基处快照

接下来对 pharmmaker/run_snapshot.sh 脚本进行修改：

1. 修改 line16 探针文件为刚才的 probelist.dat
2. 修改 line132 `sed -e "s/SSTEP/$STEP/g" -e "s/CUTOFF/$CUTOFF2/g" -e "s/AAA/resname $FPROBE and chain P and not hydrogen/g"  -e "s/BBB/resname $FPROBE and chain M and resid $FF/g" $pharmmaker_dir/CORE/snapshot2.tcl > $resdir/__ligb.tcl` 中的 `chain P` （`chain M` 好像应该一般不会有问题，具体看文件）

   如何修改：

   先将 line190 注释掉，运行一次脚本，在 z.chain.resi.probe 文件夹下 找到 v-com-ok.pdb 文件，使用文本编辑器打开，找到最后的探针行，修改为探针名后的链名即可，示例如下。

   ```text
   ATOM   5849  C2  ACETX  74      -3.253  12.969  11.563  0.00  0.00      XXX  C
   ATOM   5850  O3  ACETX  74      -4.237  12.852  12.353  0.00  0.00      XXX  O
   ATOM   5851  O4  ACETX  74      -3.158  13.892  10.773  0.00  0.00      XXX  O
   ATOM   5852  C1  ACETX  74      -2.244  11.821  11.575  0.00  0.00      XXX  C
   ATOM   5853  H11 ACETX  74      -2.625  10.836  11.231  0.00  0.00      XXX  H
   ATOM   5854  H12 ACETX  74      -1.760  11.743  12.572  0.00  0.00      XXX  H
   ATOM   5855  H13 ACETX  74      -1.407  12.101  10.900  0.00  0.00      XXX  H
   ATOM      1  C2  ACETM   1       5.000  15.500   1.000  0.99 -1.72           M
   ATOM      2  C2  ACETM   2       0.000  14.000   0.500  0.99 -1.70           M
   ATOM      3  C2  ACETM   3       5.000  13.000  -1.000  1.00 -1.67           M
   ATOM      4  C2  ACETM   4       4.500  15.500  -1.500  0.97 -1.45           M
   ATOM      5  C2  ACETM   5       0.500  11.000   4.500  0.93 -1.44           M
   ATOM      6  C2  ACETM   6       0.500  15.500   2.500  1.00 -1.38           M
   ATOM      7  C2  ACETM   7       4.000  16.500   3.000  1.00 -1.31           M
   ```
3. 将 line166 和 line167 的输出文件名互换（res 和 hs），作者写反了。
4. （可选）可以先在 PyMOL 中计算一下链 X 和链 M 上的距离，据此对 sh 文件中 line20 的 CUTOFF2=1.5 作修改，我这里改成了 2.0。

修改完成后在终端运行以下代码。

```shell
../pharmmaker/run_snapshot.sh ../drugui-simulation/p6v4qb.pdb ../drugui-simulation/p6v4qb_sim/sim.dcd
```

运行完成后在 snapshot 目录下有子目录 z.chain.resi.probe，其中有如下四类文件[^10]：out-detail-hs-x.dat，out-detail-res-x.dat，outfr-x.dat 和 outfr-count。文件解读略。

[^10]:教程中提到还有 outfr.dat 文件，但是提供的代码中并没有输出此文件，而在下一步的药效团构建中，multiprobe.sh 脚本输出了此文件。

#### 3.6.6 药效团构建

对 pharmmaker/CORE/multiprobe.sh 文件进行修改：

1. 修改 line9 探针文件为之前的 probelist.dat
2. 修改 line26 `grep "$FPROBE $FCHAIN" highaffresid.dat > ____tt` 中高亲和力残基文件为之前的 highaffresid_interest.dat

修改完成后在终端运行以下代码。

```shell
../pharmmaker/CORE/multiprobe.sh
```

输出 multipleprobe1.dat 和 multipleprobe2.dat 文件示例如下

```text
<!-- multipleprobe1.dat -->
383	5
360	5
351	5
749	4
743	4
...

<!-- multipleprobe2.dat -->
snapshot/z.U.350.ACET/outfr.dat:frame 10   
snapshot/z.U.350.ACET/outfr.dat:frame 13   
snapshot/z.U.350.ACET/outfr.dat:frame 14   
snapshot/z.U.350.ACET/outfr.dat:frame 17   
snapshot/z.U.350.ACET/outfr.dat:frame 78  
...
```

multipleprobe2.dat 文件与教程中不同，先不管。根据 multipleprobe1.dat 文件可以发现排名前三的为轨迹的第 383、360、351 帧，在这三帧中，热点与高亲和残基的集合数量都是最高的 5。

“热点与高亲和残基的集合数量“我不太懂它的意思。我猜是上一步中 out-detail-hs* 和 out-detail-res* 文件中出现的探针的交集。但是这个想法似乎不正确。

查看这些帧中的残基、探针和热点信息。以 383 帧为例。

```shell
# 例：Frame 383
$ cat snapshot/z.*/out-detail-res* | grep '383 '
383 3972 GLU 312 CG  7483 IPRO 319 C1  3.47063946723938
383 3975 GLU 312 CD  7483 IPRO 319 C1  3.8982059955596924
383 3976 GLU 312 OE1  7483 IPRO 319 C1  3.4354872703552246
383 4374 ASN 339 ND2  8219 IPRO 380 C3  3.641819477081299
383 4421 GLN 342 CG  8223 IPRO 380 OH2  3.711785316467285
383 4429 GLN 342 C  8223 IPRO 380 OH2  3.948803424835205
383 4430 GLN 342 O  8223 IPRO 380 OH2  3.790900468826294
383 4489 ASP 346 CB  8223 IPRO 380 OH2  3.620969772338867
383 4493 ASP 346 OD1  8215 IPRO 380 C1  3.785810708999634
383 4545 LYS 350 CD  9869 ACET 156 O4  3.6140317916870117
383 4548 LYS 350 CE  9869 ACET 156 O4  3.312983989715576
383 4551 LYS 350 NZ  9867 ACET 156 C2  2.996321439743042
383 4551 LYS 350 NZ  9868 ACET 156 O3  2.761639356613159
383 4551 LYS 350 NZ  9869 ACET 156 O4  2.6285452842712402
383 4545 LYS 350 CD  9869 ACET 156 O4  3.6140317916870117
383 4548 LYS 350 CE  9869 ACET 156 O4  3.312983989715576
383 4551 LYS 350 NZ  9867 ACET 156 C2  2.996321439743042
383 4551 LYS 350 NZ  9868 ACET 156 O3  2.761639356613159
383 4551 LYS 350 NZ  9869 ACET 156 O4  2.6285452842712402

$ cat snapshot/z.*/out-detail-hs* | grep '383 '
383 IPRO 228 C1 IPRO 5 1.9869842529296875
383 IPRO 261 C3 IPRO 2 1.8488157987594604
383 IPRO 351 C3 IPRO 2 1.8488157987594604
383 IPRO 380 C3 IPRO 2 1.8488157987594604
383 ACET 156 C2 ACET 6 1.146291732788086
383 ACET 156 O4 ACET 6 1.8000375032424927
383 ACET 156 C1 ACET 6 1.1198018789291382
383 ACET 156 O3 ACET 7 1.989550232887268
```

三帧结构共同的探针为 IPRO228、IPRO224、IPRO261、IPRO351、IPRO380 和 ACET156。

使用 VMD，导入 psf 和 dcd 文件，取出轨迹的第 383、360、351 帧[^11]。

[^11]:注意 VMD 的 Frame 是从 0 开始编号的，所以应该取 Frame 382、359 和 350。

在 PyMOL 中叠合，可以看到只有两个探针分子是在我所感兴趣的区域内的（{@fig:3frames_align}）。

```python
remove solvent
remove chain I
select probe, chain X and resi 228+224+261+351+380+156
```

![三帧结构和探针的叠合](https://img.limina.top/blog/三帧结构和探针的叠合.png){#fig:3frames_align}

为了便于查看，只保留 383 帧。{@fig:pharmodel}所示区域可以后续用于做药效团模型的建立。

![可用于药效团建模的成药性模拟快照](https://img.limina.top/blog/可用于药效团建模的成药性模拟快照.png){#fig:pharmodel}

## 附录 1：数据库信息提取代码

### deal_DrugBank.py {#deal_DrugBank.py}

```python
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-07-27 10:33:21
 # @Description: This python script implements a class to parse the xml file 
 #  downloaded from DrugBank to obtain DrugID, DrugName, TargetID, GeneName, 
 #  Indications and their mapping relationship.

# Variable declaration:
#     dbid : DrugBank id
#     dbname : DrugBank name
#     ptid : protein id(target)
#     dise : indication
#     gname : gene name

from xml.sax.handler import ContentHandler
from xml.sax import parse
import pandas as pd

class ExtractData(ContentHandler):

    def __init__(self):
        # Mapping relationship
        self.dbid_dbname = {}
        self.dbid_ptid = {}
        self.dbid_dise = {}
        self.ptid_gname = {}
        # Current DrugID and TargetID
        self.curr_id = ""
        self.ptid = ""
        # Limitation of traversal area
        self.limit = 0

    # Get the content of the traversed tag, such as <ele>content.....</ele>
    def characters(self, content):
        if self.limit == 2:
            self.curr_id = content
            self.limit = 3

        elif self.limit == 4:
            self.dbid_dbname[self.curr_id] = content
            self.limit = 0

        elif self.limit == 5:
            self.dbid_dise[self.curr_id] = content
            self.limit = 0
            
        elif self.limit == 9:
            self.ptid_gname[self.ptid] = content
            self.limit = 0

    # Called when traversing to the beginning of the label        
    def startElement(self, name, attrs):
        if name == "drug":
            self.limit = 1
            
        if self.limit == 1 and name == "drugbank-id" and attrs:
            if attrs["primary"] == "true":
                self.limit = 2

        elif self.limit == 3 and name == "name":
            self.limit = 4

        elif name == "indication":
            self.limit = 5

        elif name == "targets":
            self.limit = 7

        elif self.limit == 7 and name == "polypeptide":
            self.dbid_ptid.setdefault(self.curr_id, set()).add(attrs["id"])
            self.ptid = attrs["id"]
            self.limit = 8
            
        elif self.limit == 8 and name == "gene-name":
            self.limit = 9
            
    # Called when the traversal to the end of the label 
    def endElement(self, name):
        if name == "drug-interactions":
            self.limit = 0
            
        elif name == "targets":
            self.limit = 0

    # Called at the end of the traversal
    def endDocument(self):
        # Mapping of DrugBank id and DrugName
        list1_key = []
        list1_val = []
        list1_columns = "DrugName",
        for key,val in self.dbid_dbname.items():
            list1_key.append(key)
            list1_val.append(val)
        file1 = pd.DataFrame(index = list1_key, columns = list1_columns, data = list1_val)
        file1.to_csv('dbid_dbname.csv')

        # Mapping of DrugBank id and Indication
        list2_key = []
        list2_val = []
        list2_columns = "Indication",
        for key,val in self.dbid_dise.items():
            list2_key.append(key)
            list2_val.append(val)
        file2 = pd.DataFrame(index = list2_key, columns = list2_columns, data = list2_val)
        file2.to_csv('dbid_dise.csv')

        # Interaction mapping between DrugBank id and TargetID
        list3_key = []
        list3_val = []
        for key,val in self.dbid_ptid.items():
            list3_key.append(key)
            list3_val.append(list(val))
        file3 = pd.DataFrame(index = list3_key, data = list3_val)
        file3.to_csv('dbid_ptid.csv')

        # Interaction mapping between TargetID and DrugBank id
        ptid_dbid = {}
        for key,val in self.dbid_ptid.items():
            for v in val:
                ptid_dbid.setdefault(v, list()).append(key)
        file4 = pd.DataFrame(index = list(ptid_dbid.keys()), data = list(ptid_dbid.values()))
        file4.to_csv('ptid_dbid.csv')

        # Mapping of TargetID and GeneName
        list5_key = []
        list5_val = []
        list5_columns = "GeneName",
        for key,val in self.ptid_gname.items():
            list5_key.append(key)
            list5_val.append(val)
        file5 = pd.DataFrame(index = list5_key, columns = list5_columns, data = list5_val)
        file5.to_csv('ptid_gname.csv')

parse('../drugbank.xml', ExtractData())

```

### DrugBank_TDI.R {#DrugBank_TDI.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-07-27 00:14:48
 # @Description: This R code file is used to obtain the final target_drug_indication form of DrugBank

library(tidyverse)
setwd("D:/Study/Project/Graduation/deal_DrugBank")

rm(list = ls())
targetid_drugid <- read.delim("ptid_dbid.csv", header = T, sep = ",", check.names = F, stringsAsFactors = F, na.strings = c("NA", ""))
drugid_drugname <- read.delim("dbid_dbname.csv", header = T, sep = ",", check.names = F, stringsAsFactors = F)
targetid_genename <- read.delim("ptid_gname.csv", header = T, sep = ",", check.names = F, stringsAsFactors = F, na.strings = c("NA", "\n"))
drugid_indication <- read.delim("dbid_dise.csv", header = T, sep = ",", check.names = F, stringsAsFactors = F, na.strings = c("NA", "\n"))

colnames(drugid_drugname)[1] <- "DrugID"
colnames(targetid_genename)[1] <- "TargetID"
colnames(drugid_indication)[1] <- "DrugID"

# Unfinished loop with higher time complexity(The second alternative method)
# targetid_drugid_single <- data.frame()
# for (i in 1:nrow(targetid_drugid)) {
#   for (j in targetid_drugid[i, 2:ncol(targetid_drugid)]){
#     if (j != ""){
#       targetid_drugid_single <- rbind(targetid_drugid_single, c(targetid_drugid[i, 1], targetid_drugid[i, j]))
#     }
#     else {
#       break
#     }
#   }
# }

drugid_targetid_single <- data.frame()
for (i in 1:nrow(targetid_drugid)) {
  ndrug <- sum(!is.na(targetid_drugid[i, ])) - 1
  druglist <- as.data.frame(t(as.data.frame(c(targetid_drugid[i, 1:ndrug+1]))))
  druglist$Target <- targetid_drugid[i, 1]
  drugid_targetid_single <- rbind(drugid_targetid_single, druglist)
}
colnames(drugid_targetid_single)[1:2] <- c("DrugID", "TargetID")

drugid_targetid_single$DrugName <- NA
drugid_targetid_single$GeneName <- NA
drugid_targetid_single$Indication <- NA

drug_cycle <- 0
for (drug in drugid_targetid_single$DrugID) {
  drug_cycle <- drug_cycle + 1
  drug_row <- which(drugid_drugname$DrugID == drug)
  if (length(drug_row) != 0) {
    drugid_targetid_single$DrugName[drug_cycle] <- drugid_drugname$DrugName[drug_row]
  }
}

target_cycle <- 0
for (target in drugid_targetid_single$TargetID) {
  target_cycle <- target_cycle + 1
  target_row <- which(targetid_genename$TargetID == target)
  if (length(target_row) != 0) {
    drugid_targetid_single$GeneName[target_cycle] <- targetid_genename$GeneName[target_row]
  }
}

indication_cycle <- 0
for (drug in drugid_targetid_single$DrugID) {
  indication_cycle <- indication_cycle + 1
  drug_row <- which(drugid_indication$DrugID == drug)
  if (length(drug_row) != 0) {
    drugid_targetid_single$Indication[indication_cycle] <- drugid_indication$Indication[drug_row]
  }
}

DrugBank_TDI <- select(drugid_targetid_single, c(2, 4, 1, 3, 5))
write.csv(DrugBank_TDI, file = "DrugBank_Target_Drug_Indication.csv", row.names = F)

```

### Target_Drug.R {#Target_Drug.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-07-27 00:18:43
 # @Description: This R code file is used to obtain the corresponding relationship 
 #  between the target and the drug from the downloaded file P1-01-TTD_target_download.txt.

library(tidyverse)
setwd("D:/Study/Project/Graduation/deal_TTD")

rm(list = ls())
# You need to comment out the previous description information line in the txt file with "#"
# This P1-01-TTD_target_download.txt file needs to be converted to csv file first because of the line break problem
target_all <- read.delim("P1-01-TTD_target_download.csv", header = F, row.names = NULL, sep = ",", 
                         check.names = F, stringsAsFactors = F, comment.char = "#", na.strings = c("NA", ""))

# For each target, not all the information is available, so for convenience, 
# I extract the TARGETID, GENENAME, and DRUGINFO lines that I need
targetid_row <- which(target_all[, 2] == "TARGETID")
genename_row <- which(target_all[, 2] == "GENENAME")
druginfo_row <- which(target_all[, 2] == "DRUGINFO")

target_all_filter <- target_all[c(targetid_row, genename_row, druginfo_row), ]
target_all_filter_order <- target_all_filter[order(as.numeric(row.names(target_all_filter))), ]

ntgid <- c()
targetid_filter_row <- which(target_all_filter_order[, 2] == "TARGETID")
genename_filter_row <- which(target_all_filter_order[, 2] == "GENENAME")
druginfo_filter_row <- which(target_all_filter_order[, 2] == "DRUGINFO")
for (i in 2:length(targetid_filter_row)) {
  ntgid <- append(ntgid, targetid_filter_row[i] - targetid_filter_row[i - 1])
}
ntgid_last <- nrow(target_all_filter_order) - targetid_filter_row[length(targetid_filter_row)] + 1
ntgid <- append(ntgid, ntgid_last)

target_locate <- 0
targetid_genename <- c()
drugid <- c()
drugname <- c()
for (n in ntgid) {
  if (n > 2) {
    target_locate = target_locate + 1
    targetid_genename <- append(targetid_genename, 
                                rep(paste(target_all_filter_order[targetid_filter_row[target_locate], 3], 
                                          target_all_filter_order[targetid_filter_row[target_locate]+1, 3], 
                                          sep = ",,"), n-2))
    
    for (j in 3:n) {
      drugid <- append(drugid, target_all_filter_order[targetid_filter_row[target_locate]+j-1, 3])
      drugname <- append(drugname, target_all_filter_order[targetid_filter_row[target_locate]+j-1, 4])
    }
  }
  else {
    target_locate = target_locate + 1
  }
}
drug <- paste(drugid, drugname, sep = ",,")

target_drug_str <- paste(targetid_genename, drug, sep = ",,")
target_drug <- as.data.frame(target_drug_str, StringsAsFactors = FALSE)

target_drug <- separate(data = target_drug, col = target_drug_str, 
                        into = c("TargetID", "GeneName", "DrugID", "DrugName"), sep = ",,")
write.csv(target_drug, file = "Target_Drug.csv", row.names = F)

```

### Target_or_Drug_Indication.R {#Target_or_Drug_Indication.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-07-27 00:18:58
 # @Description: This R code file is used to process the mapping relationship 
 #  between targets and indications or drugs and indications.

library(tidyverse)
setwd("D:/Study/Project/Graduation/deal_TTD")

rm(list = ls())
# You need to comment out the previous description information line in the txt file with "#"
target_dis <- read.delim("P1-06-Target_disease.txt", header = F, row.names = NULL, sep = "\t", 
                         check.names = F, stringsAsFactors = F, comment.char = "#", na.strings=c("NA", ""))

indication_row <- which(target_dis[, 2] == "INDICATI")
target_dis[indication_row, 3] <- target_dis[indication_row, 4]
target_dis_nogroup <- target_dis[, -ncol(target_dis)]

# This step is actually not necessary...
target_dis_nona <- na.omit(target_dis_nogroup)


# # Drug_Indication can also be treated in the same way
# target_dis_nona <- read.delim("P1-05-Drug_disease.txt", header = F, row.names = NULL, sep = "\t",
#                          check.names = F, stringsAsFactors = F, comment.char = "#", na.strings=c("NA", ""))
# ntgid <- c()
# targetid_row <- which(target_dis_nona[, 2] == "TTDDRUID")


# For Target_Indication
ntgid <- c()
targetid_row <- which(target_dis_nona[, 2] == "TARGETID")

for (i in 2:length(targetid_row)) {
  ntgid <- append(ntgid, targetid_row[i] - targetid_row[i - 1])
}
ntgid_last <- nrow(target_dis_nona) - targetid_row[length(targetid_row)] + 1
ntgid <- append(ntgid, ntgid_last)

target_locate <- 0
targetid <- c()
indication <- c()
for (n in ntgid) {
  target_locate = target_locate + 1
  targetid <- append(targetid, rep(paste(target_dis_nona[targetid_row[target_locate], 3], 
                        target_dis_nona[targetid_row[target_locate]+1, 3], 
                        sep = ",,"), n-2))
  for (j in 3:n) {
    indication <- append(indication, target_dis_nona[targetid_row[target_locate]+j-1, 3])
  }
}

target_ind_str <- paste(targetid, indication, sep = ",,")
target_ind <- as.data.frame(target_ind_str, StringsAsFactors = FALSE)
# For Target_Indication
target_ind <- separate(data = target_ind, col = target_ind_str, 
                       into = c("TargetID", "TargetName", "Indication"), sep = ",,")
write.csv(target_ind, file = "Target_Indication.csv", row.names = F)


# # For Drug_Indication
# target_ind <- separate(data = target_ind, col = target_ind_str, 
#                        into = c("DrugID", "DrugName", "Indication"), sep = ",,")
# write.csv(target_ind, file = "Drug_Indication.csv", row.names = F)

```

### TTD_TDI.R {#TTD_TDI.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-07-27 00:19:12
 # @Description: This R code file is used to obtain the final target_drug_indication form of TTD

library(tidyverse)
setwd("D:/Study/Project/Graduation/deal_TTD")

rm(list = ls())
target_drug <- read.delim("Target_Drug.csv", header = T, sep = ",", check.names = F, stringsAsFactors = F)
drug_ind <- read.delim("Drug_Indication.csv", header = T, sep = ",", check.names = F, stringsAsFactors = F)

drug_ind$TargetID <- NA
drug_ind$GeneName <- NA
drug_ind_cycle <- 0
for (drug in drug_ind$DrugID) {
  drug_ind_cycle <- drug_ind_cycle + 1
  rowid <- which(target_drug$DrugID == drug)
  if (length(rowid) != 0) {
    drug_ind$TargetID[drug_ind_cycle] <- target_drug$TargetID[rowid]
    drug_ind$GeneName[drug_ind_cycle] <- target_drug$GeneName[rowid]
  }
}

TTD_TDI <- select(drug_ind, c(4, 5, 1, 2, 3))
write.csv(TTD_TDI, file = "TTD_Target_Drug_Indication.csv", row.names = F)

```

### target_match.R {#target_match.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-07-29 15:23:30
 # @Description: This R file is used to screen the real target proteins in the 
 #  protein-protein interaction network and map them to two databases, DrugBank and TTD.

library(tidyverse)
setwd("D:/Study/Project/Graduation/target_match")

rm(list = ls())
protein_net <- read.table("../coreNet.txt", stringsAsFactors = F)
DB_info <- read.delim("../deal_DrugBank/DrugBank_Target_Drug_Indication.csv", sep = ",", header = T, stringsAsFactors = F)
TTD_info <- read.delim("../deal_TTD/TTD_Target_Drug_Indication.csv", sep = ",", header = T, stringsAsFactors = F)

protein <- distinct(as.data.frame(c(protein_net[, 1], protein_net[, 2])))
colnames(protein)[1] <- "GeneName"

DB_merge <- merge(protein, DB_info, by = "GeneName", sort = T)
TTD_merge <- merge(protein, TTD_info, by = "GeneName", sort = T)

nDB <- 0
nTTD <- 0
target_list <- c()
for (i in protein$GeneName) {
  if (i %in% DB_merge$GeneName == T) {
    nDB <- nDB + 1
    target_list <- append(target_list, i)
  }
  else if (i %in% TTD_merge$GeneName == T) {
    nTTD <- nTTD + 1
    target_list <- append(target_list, i)
  }
}

target <- as.data.frame(target_list)
colnames(target)[1] <- "GeneName"
non_target <- setdiff(protein, target)

write.csv(DB_merge, file = "DrugBank_Target.csv", row.names = F)
write.csv(TTD_merge, file = "TTD_Target.csv", row.names = F)
write.csv(target, file = "coretarget.csv", row.names = F)
write.csv(non_target, file = "non-coretarget.csv", row.names = F)

```

### enrichment.R {#enrichment.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-10-12 21:55:20
 # @Description: If using the background gene set, load the DOSE package

# library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
setwd("D:/Study/4_Senior/Bioinformatics_Skills_Training/2_Drug-Virtual-Screening/enrichment_analysis/")

rm(list = ls())
target_symbol <- read.table("../target_match/coretarget.csv", header = T, stringsAsFactors = F)
non_target_symbol <- read.table("../target_match/non-coretarget.csv", header = T, stringsAsFactors = F)
all_symbol <- rbind(target_symbol, non_target_symbol)

# Specific target and non-target of lung cancer in 744 proteins
database_lctarget <- read.table("../lungcancer_target.csv", sep = ",", header = T, stringsAsFactors = F)
lc_target_symbol <- na.omit(distinct(dplyr::intersect(dplyr::select(database_lctarget, 2), all_symbol)))
non_lc_target_symbol <- dplyr::setdiff(all_symbol, lc_target_symbol)

# data(geneList, package = "DOSE") # Background gene set for enrichment analysis
# geneid <- bitr(target_symbol$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# GO
ego <- enrichGO(gene = all_symbol$GeneName, 
                # universe = names(geneList), # Background gene set, need geneid to use this parameter
                OrgDb = org.Hs.eg.db, 
                keyType = "SYMBOL",
                ont = "all", # One of CC, BP, MF, all
                pAdjustMethod = "BH", # Correction methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                pvalueCutoff = 1, 
                qvalueCutoff = 1) # If using geneid, add a parameter: readable = TRUE

go <- as.data.frame(ego)
go_CC <- go[go$ONTOLOGY == "CC",][1:10, ]
go_MF <- go[go$ONTOLOGY == "MF",][1:10, ]
go_BP <- go[go$ONTOLOGY == "BP",][1:10, ]
go_enrich_df <- data.frame(ID = c(go_CC$ID, go_MF$ID, go_BP$ID), 
                           Description = c(go_CC$Description, go_MF$Description, go_BP$Description), 
                           GeneNumber = c(go_CC$Count, go_MF$Count, go_BP$Count), 
                           group = factor(c(rep("Cellular Component", 10), 
                                           rep("Molecular Function", 10), 
                                           rep("Biological Process", 10)), 
                                         levels = c("Cellular Component", "Molecular Function", "Biological Process")))
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

color <- c("#8DA1CB", "#FD8D62", "#66C3A5")
ggplot(data = go_enrich_df, aes(x = number, y = GeneNumber, fill = group)) + 
  geom_bar(stat = "identity", width = 0.8) + coord_flip() + scale_fill_manual(values = color) + 
  scale_x_discrete(labels = rev(go_enrich_df$Description)) + xlab("GO term") + theme_bw() + 
  theme(axis.text = element_text(face = "bold", color = "gray50", size = 12), 
        legend.title = element_text(size = 14), legend.text = element_text(size = 13)) + 
  labs(title = "The Most Enriched GO Terms of All")

# KEGG
# In KEGG enrichment analysis, symbol does not support human species (hsa) and needs to be converted to geneid
geneid <- bitr(all_symbol$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ekegg <- enrichKEGG(gene = geneid$ENTREZID, 
                    organism = "hsa", # Shorthand for human species:"hsa"
                    pvalueCutoff = 1)
barplot(ekegg, showCategory = 20, title = "The Most Enriched KEGG pathway of All")

# convert geneid to genesymbol
ekegg2 <- setReadable(pairwise_termsim(ekegg), 
                      OrgDb = "org.Hs.eg.db", 
                      keyType = "ENTREZID")

pathways <- data.frame(ekegg2@result$ID, ekegg2@result$Description)

# Look for potential target
pathresult <- ekegg2@result
lcpathway <- pathresult[which(pathways$ekegg2.result.ID == "hsa05222" | pathways$ekegg2.result.ID == "hsa05223"), ]
lcgene <- data.frame()
for (i in lcpathway$geneID) {
  lcgene_temp <- as.data.frame(strsplit(i, '/'))
  colnames(lcgene_temp)[1] <- "GeneName"
  lcgene <- rbind(lcgene, lcgene_temp)
}
lcgene <- distinct(lcgene)
known_target <- dplyr::intersect(lcgene, lc_target_symbol)
potential_target <- dplyr::intersect(lcgene, non_lc_target_symbol)

write.csv(known_target, "know_lctarget.csv", row.names = F)
write.csv(potential_target, "potential_lctarget.csv", row.names = F)

```

### residues.R {#residues.R}

```R
 # @Author: Min Li
 # @Email: mli.bio@outlook.com
 # @Last Modified by: Min Li
 # @Timestamp for Last Modification: 2021-10-21 19:29:52
 # @Description: ...

library(tidyverse)
library(stringr)
library(ggplot2)
# library(ggpubr)
library(ggrepel)

rm(list=ls())
setwd("/media/limin/Office/Study/4_Senior/Bioinformatics_Skills_Training/2_Drug-Virtual-Screening/pharmmaker/cavity_plus/")
data <- read.table("cavity/outputcavity.txt", header = F, fill = T, skip = 1, stringsAsFactors = F)
data <- separate(data, V1, into = c("cavity", "V1"), sep = "\\|")[, -1]

resnum <- numeric()
for (i in 1:nrow(data)) {
  resnum <- c(resnum, as.numeric(str_extract_all(data[i, ], "\\d+")))
}
resnum <- na.omit(resnum)

dim(resnum) <- c(281, 1)
range <- as.data.frame(resnum)
range$V2 <- "res"

png("Distribution of Residual Number", width = 960, height = 1080)
ggplot(range, aes(x = V2, y = resnum)) + 
  stat_boxplot(geom = "errorbar", width = 0.15) + 
  geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white") +  
  geom_jitter(aes(fill = V1), width = 0.3, shape = 21, size = 2.5) + 
  geom_text_repel(aes(V2, resnum, label = V1, color = "#FD8D62"), max.overlaps = 22)
dev.off()

```

## 附录 2：关于 Pharmmaker 官方教程的 debug

> 本实验临近结束快做完了我才发现，Pharmmaker 在 Github 上还有一个版本，与官网提供的脚本是不一致的。。建议后人使用 GitHub 上的版本，可能 bug 会少一点。

### DruGUI 的使用注意点

1. 需要修改 drugui/pkgIndex.tcl 文件的第 11 行，将 1.0 更改为 1.1，或将 drugui/drugui.tcl 文件的第 26 行的 1.1 改为 1.0，才能在 VMD 中正常加载此插件。根据其 GUI 界面和 log 文件，我认为正确的改法应该是后者。
2. 建立 psf 文件不能使用 VMD 自带的 Automatic PSF Builder 插件，建议根据 NAMD2 的教程进行构建。
3. 使用的拓扑文件是 top_all27_prot_lipid.inp ，有些蛋白可能会有问题，建立的 psf 文件有大段缺失，需要自主决定选择链。
4. 成药性模拟的 sim.conf 配置文件默认模拟时间太长（40ns），可以修改 sim.conf 文件参数少跑一点。

### Pharmmaker 的使用注意点

1. 高亲和力残基分析（教程中的 3.1）中可以先挑选一小段序列，使用某一探针测试来自定义合适的 cutoff。
2. 寻找高亲和力残基附近的热点时（教程中的 3.2），要将上一步默认输出的 probe-list.dat 文件中的探针名写在一行，并用空格隔开；提供的命令行有误，链文件和探针文件参数写反了。
3. 关注感兴趣区域的热点（教程中的 3.3），最好需要对蛋白的背景知识有一定的理解。
4. 热点和带有探针的高亲和力残基处快照（教程中的 3.4），对 pharmmaker/run_snapshot.sh 脚本进行修改：
   1. line16 探针文件还是上述第 2 点中的探针文件
   2. line132 `sed -e "s/SSTEP/$STEP/g" -e "s/CUTOFF/$CUTOFF2/g" -e "s/AAA/resname $FPROBE and chain P and not hydrogen/g"  -e "s/BBB/resname $FPROBE and chain M and resid $FF/g" $pharmmaker_dir/CORE/snapshot2.tcl > $resdir/__ligb.tcl` 中的 `chain P` 和 `chain M`。一般将 P 改为 X 即可，具体可试运行一次脚本，在 z.chain.resi.probe 文件夹下 找到 v-com-ok.pdb 文件，使用文本编辑器打开（也可用 PyMOL），找到最后的探针行，修改为探针名后的链名即可（注意分清 probe 和 hotspot）。
   3. 将 line166 和 line167 的输出文件名互换（res 和 hs），作者写反了。
   4. （可选）可以先在 PyMOL 中计算一下链 X 和链 M 上的距离，据此对 sh 文件中 line20 的 CUTOFF2=1.5 作修改。
5. 根据提供的脚本代码，上一步中并不会输出 outfr.dat 文件，但是教程中说有这个文件......这个文件是在下一步的药效团构建中，由 multiprobe.sh 脚本输出的。
6. 药效团构建步骤（教程中的 3.5）中：
   1. line9 探针文件还是上述第 2 点中的探针文件
   2. multipleprobe2.dat 文件与教程中展示的不同，是提供的代码的问题......实际输出的此文件没什么意义。

总之我非常不理解 Pharmmaker 的作者是如何写出此种奇怪的教程的。。
