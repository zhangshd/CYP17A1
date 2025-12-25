# **针对猪CYP17A1抑制剂的计算机辅助药物设计与虚拟筛选综合研究报告**

## **1\. 绪论：猪CYP17A1作为药物靶点的生物学意义与研究背景**

### **1.1 雄烯酮、公猪膻味与免疫去势的各种可能性**

在现代养猪业中，公猪膻味（Boar Taint）是一个长期存在的经济与动物福利问题。这种令人不悦的气味主要源于两种化合物在脂肪组织中的积累：雄烯酮（5$\\alpha$-androst-16-en-3-one）和粪臭素（3-methylindole）1。其中，雄烯酮是一种信息素，其生物合成途径主要发生在睾丸的莱迪希细胞（Leydig cells）中。该合成路径与雄激素（睾酮）的合成高度重叠，而细胞色素P450 17A1（CYP17A1，亦称P450c17）正是这一路径中的关键限速酶。

CYP17A1是一种双功能酶，催化孕烯醇酮（Pregnenolone）和孕酮（Progesterone）的17$\\alpha$-羟化反应，随后进行17,20-裂解反应，生成脱氢表雄酮（DHEA）和安雄烯二酮。这些产物不仅是睾酮的前体，也是雄烯酮合成的直接底物。因此，抑制猪睾丸中的CYP17A1活性，理论上可以阻断雄烯酮的合成，从而在不进行物理去势（手术阉割）的情况下消除公猪膻味。这种“药物去势”或“免疫去势”策略对于改善动物福利、减少手术带来的感染风险以及提高饲料转化率具有重大的现实意义 2。

### **1.2 靶点蛋白：猪源CYP17A1 (P19100)**

本研究锁定的靶点为猪（*Sus scrofa*）的CYP17A1蛋白，其UniProt登录号为**P19100** 4。该蛋白由509个氨基酸组成，定位于内质网膜。与其人类直系同源物（P05093）相比，猪CYP17A1在序列上表现出高度的保守性（序列一致性\>80%），但在底物特异性和催化动力学上可能存在细微的物种差异 5。

虽然人类CYP17A1已有多个高分辨率晶体结构被解析（如与抑制剂阿比特龙结合的3RUK），但截至目前，**猪源CYP17A1尚未有已发表的晶体结构**。这使得针对该靶点的药物设计必须依赖同源建模（Homology Modeling）技术。此外，作为典型的细胞色素P450酶，其活性中心包含一个血红素（Heme）辅基——原卟啉IX铁（Protoporphyrin IX containing Fe）。该辅基通过一个保守的半胱氨酸硫醇键轴向配位到铁原子上，这种特殊的电子环境为分子对接和分子动力学（MD）模拟带来了巨大的挑战 7。

### **1.3 本报告的研究目标与技术路线**

本报告旨在建立一套详尽的计算工作流，以解决用户提出的核心难点：如何在缺乏晶体结构且含有复杂金属辅基的情况下，针对猪CYP17A1进行有效的虚拟筛选。报告将涵盖以下关键领域：

1. **结构建模策略**：利用Swiss-Model构建高质量的猪CYP17A1三维结构，重点阐述模板选择的依据及活性口袋的定义。  
2. **血红素辅基处理**：深入探讨Heme-Fe的电子结构特性，提供在分子对接（AutoDock Vina）和动力学模拟（AMBER）中处理该辅基的具体参数与方法。  
3. **力场参数详解**：系统梳理文献中针对P450血红素铁的力场参数（特别是Shahrokh等人的参数集），并提供参数生成的实操指南。  
4. **既往研究综述**：分析现有文献中关于CYP17A1及相关P450酶的模拟工作，提炼可供借鉴的经验。

## ---

**2\. 猪CYP17A1的同源建模与活性口袋定义**

### **2.1 序列分析与关键残基比对**

在进行建模之前，必须对猪CYP17A1（P19100）与潜在模板（人源CYP17A1）进行精细的序列比对。研究表明，P450酶的活性位点具有高度的可塑性，但关键的催化残基在进化上是严格保守的。

**关键残基分析：**

* **Asn202 (人源编号) / Asn203 (猪源推测编号)**：在人类CYP17A1中，Asn202残基对于底物的定向至关重要。它与底物A环的3$\\beta$-羟基形成氢键。突变研究（N202S）显示，该位点的改变会逆转酶对孕烯醇酮与孕酮的底物偏好性 5。文献5的序列比对显示，在猪（*Sus scrofa*）序列中，该位置同样也是**天冬酰胺（Asn）**。这意味着人源抑制剂（如阿比特龙）在猪CYP17A1中的结合模式可能与在人体中高度相似，均依赖于该Asn残基提供的氢键锚定。  
* **轴向配位半胱氨酸（Cys442）**：这是P450酶的标志性残基，负责与血红素铁形成硫醇配位键。在建模过程中，必须确保该半胱氨酸的硫原子与血红素铁之间的距离保持在2.3-2.4 Å的物理合理范围内。

### **2.2 Swiss-Model建模策略与模板选择**

由于P450酶存在显著的“诱导契合”效应（Induced Fit），活性口袋的体积和形状会随结合配体的不同而发生显著变化。因此，选择一个与“抑制剂结合状态”相符的模板至关重要。

**模板筛选评估表：**

| PDB ID | 来源物种 | 结合配体 | 分辨率 | 适用性分析 |
| :---- | :---- | :---- | :---- | :---- |
| **3RUK** | *Homo sapiens* | **阿比特龙 (Abiraterone)** | 2.60 Å | **最优推荐**。该结构展示了酶与高亲和力甾体抑制剂结合时的构象，活性口袋处于“开放”且适应大分子甾体骨架的状态 9。 |
| **4NKV** | *Homo sapiens* | 阿比特龙 | 2.60 Å | 高度适用，与3RUK类似。 |
| **4NKW** | *Homo sapiens* | 孕烯醇酮 (底物) | 2.50 Å | 较差。底物结合构象可能未充分暴露适合抑制剂结合的疏水口袋。 |
| **5IRQ** | *Homo sapiens* | Orteronel (非甾体) | 2.75 Å | 适用于非甾体抑制剂筛选，但对于阿比特龙类似物，3RUK更佳。 |

**具体操作流程：**

1. **输入序列**：将P19100的FASTA序列输入Swiss-Model。  
2. **强制指定模板**：虽然系统会自动搜索，但在本研究中，建议手动指定\*\*3RUK (Chain A)\*\*作为模板，以确保模型继承抑制剂结合态的口袋特征。  
3. **辅基处理（关键步骤）**：在建模选项中，必须显式勾选包含**HEM**（血红素）配体。Swiss-Model会将模板中的血红素坐标直接迁移到新模型中。  
4. **配体占位**：建议同时勾选包含模板中的配体（阿比特龙，ABR）。虽然最终的对接需要移除它，但在建模阶段保留它可以作为“模具”，防止侧链在能量最小化过程中塌陷进活性口袋，从而保持口袋的体积 11。

### **2.3 活性口袋的定义与网格生成**

基于以3RUK为模板构建的猪CYP17A1同源模型，活性口袋的位置可以明确定义在血红素平面的上方（Distal side）。

* **口袋中心（Grid Center）**：应设定为血红素中心铁原子（Fe）的坐标。  
* **口袋尺寸（Grid Size）**：建议设置为 $22 \\times 22 \\times 22$ Å 或 $25 \\times 25 \\times 25$ Å。这一范围不仅覆盖了血红素上方的催化中心，也包含了底物进入的通道（F-G环区域）以及Asn202所在的底物锚定区 12。  
* **空间限制**：在定义口袋时，需注意I螺旋（I-helix）的位置，它是底物结合的主要疏水壁。

## ---

**3\. 核心难点攻关：血红素（Heme）辅基的处理方法**

用户查询的核心难点在于“PROTOPORPHYRIN IX CONTAINING FE”辅基在计算模拟中的处理。这是一个典型的生物无机化学计算问题：过渡金属（铁）的电子结构复杂，且存在共价配位键，常规的力场和对接软件往往无法准确描述。

### **3.1 分子对接中的血红素处理：电荷与原子类型修正**

问题描述：  
标准的AutoDock Vina打分函数主要基于范德华力、氢键和疏水作用，它将金属原子视为简单的带电球体。然而，许多CYP17A1抑制剂（如阿比特龙、酮康唑）属于“II型结合剂”，其机制是抑制剂分子上的氮原子（如吡啶氮或咪唑氮）与血红素铁形成配位键。Vina默认的范德华半径参数往往会判定这种短距离的金属-氮相互作用（约2.1-2.3 Å）为“原子重叠”（Steric Clash），从而给予极高的惩罚分，导致无法搜索到正确的结合构象 13。  
解决方案：手动修正PDBQT文件  
为了在Vina中正确模拟这种金属配位作用，必须对受体（Receptor）的PDBQT文件进行特殊处理。文献12提供了经验证的修正方案：

1. **铁原子的电荷设定**：尽管在生理状态下血红素铁通常为+3价（休止态），但在Vina的静电势模型中，直接使用+3.0的电荷可能导致过强的非特异性静电吸附。研究建议将铁原子的部分电荷（Partial Charge）手动修改为 **\+2.000** 或 **\+1.500**。这既能提供足够的静电吸引力引导氮原子靠近，又不会造成计算崩溃。  
2. **原子类型（Atom Type）**：必须确保PDBQT文件中铁原子的原子类型被标记为 Fe（区分大小写）。AutoDock Vina内部包含针对Fe的范德华参数（$R\_{ii} \\approx 1.20$ Å），这比一般的金属参数更允许配体靠近。  
3. PDBQT文件编辑示例：  
   打开生成的receptor.pdbqt文件，找到HETATM行，修改如下：  
   HETATM 4482 FE HEM A 500 5.890 24.568 \-1.058 1.00 18.80 \+2.000 Fe  
   （注意倒数第二列为电荷，最后一列为原子类型）。

进阶方案：使用GalaxyDock2-HEME  
如果AutoDock Vina的效果不佳，强烈推荐使用专门针对血红素蛋白开发的对接工具 GalaxyDock2-HEME 16。

* **优势**：该软件在评分函数中引入了专门的“金属-配体配位能”项（Metal-Ligand Coordination Energy Term），能够精确识别并奖励氮原子与铁原子之间的配位几何构型。  
* **性能**：基准测试显示，在区分铁结合型抑制剂与非结合型分子方面，GalaxyDock2-HEME显著优于AutoDock Vina和Gold 17。

### **3.2 分子动力学模拟中的血红素处理：Shahrokh力场参数**

问题描述：  
在AMBER或GROMACS等MD软件中，标准的力场（如ff14SB）不包含血红素与轴向半胱氨酸之间的键合参数。此外，铁原子的电荷分布并非简单的点电荷，而是通过卟啉环和硫原子高度离域的。直接使用通用参数会导致模拟过程中活性中心解体。  
解决方案：Shahrokh et al. (2012) 参数集  
文献8一致指向 Shahrokh 等人于2012年发表在《Journal of Computational Chemistry》上的工作，这是目前P450模拟领域的“金标准”参数。  
该参数集基于高精度的量子力学（QM/B3LYP）计算得出，涵盖了CYP催化循环中的多个状态。对于抑制剂筛选，主要关注以下两种状态：

1. **休止态（Resting State, Hexacoordinated）**：铁为六配位（配体为水），自旋态可能是低自旋（S=1/2）。  
2. **抑制剂结合态**：如果模拟的是II型抑制剂，铁与抑制剂氮原子形成六配位。

参数集关键数据提取（参考19）：  
Shahrokh参数集将血红素（HEM）和轴向半胱氨酸（CYS）视为一个整体进行电荷拟合。以下是部分关键原子的RESP电荷参考值（具体数值需下载原始补充材料确认）：

| 原子名称 | 原子类型 | 说明 | 典型电荷 (Resting State) |
| :---- | :---- | :---- | :---- |
| **FE** | Fe | 铁中心 | \+0.40 \~ \+0.60 (非+2/+3) |
| **NA/NB/NC/ND** | no | 卟啉氮 | \-0.30 \~ \-0.50 |
| **SG** | S | 轴向半胱氨酸硫 | \-0.30 \~ \-0.60 |

实操指南：如何获取并使用这些参数  
由于版权原因无法直接生成完整文件，用户需按照以下步骤操作：

1. **下载**：访问Shahrokh文章的Supporting Information或AMBER Parameter Database 20，下载 HEM.mol2（包含原子类型和电荷）和 frcmod.heme（包含力常数）。  
2. 拓扑构建（TLeap）：  
   在AMBER的tleap模块中，必须显式定义铁与半胱氨酸硫之间的键。  
   Bash  
   source leaprc.protein.ff14SB  
   source leaprc.gaff  
   loadamberparams frcmod.heme       \# 加载Shahrokh力场参数  
   loadoff HEM.lib                   \# 加载血红素库文件  
   complex \= loadpdb structure.pdb  
   bond complex.FE complex.SG        \# 关键：手动创建Fe-S键

   *注意：structure.pdb中的半胱氨酸残基名称可能需要改为CYM（去质子化态）以匹配力场。*

### **3.3 替代方案：使用MCPB.py自行生成参数**

如果无法获取Shahrokh参数，或者因为猪CYP17A1的几何构型特殊导致标准参数不适用，建议使用AmberTools中的 **MCPB.py** 工具自行生成参数 21。

**简明教程：**

1. **模型准备**：从同源模型中截取包含Heme、Cys442（侧链）和配体的小分子簇模型。  
2. **QM计算**：使用Gaussian进行几何优化和频率计算（B3LYP/6-31G\*）。  
3. **参数生成**：  
   * 运行 MCPB.py \-i input.in \-s 1 生成高斯输入文件。  
   * 运行高斯计算。  
   * 运行 MCPB.py \-i input.in \-s 2 进行Seminario力常数拟合。  
   * 运行 MCPB.py \-i input.in \-s 3 进行RESP电荷拟合。  
   * 运行 MCPB.py \-i input.in \-s 4 生成最终的 .lib 和 .frcmod 文件。

## ---

**4\. 虚拟筛选工作流与文献回顾**

### **4.1 已有的相关工作回顾**

在针对猪CYP17A1的特定研究方面，文献11报道了早在2014年就有团队基于哺乳动物模板构建了*Sus scrofa* CYP17A1的同源模型，并进行了初步分析。文献2和3则通过体外实验（微粒体孵育）研究了猪CYP17A1对孕烯醇酮和孕酮的代谢动力学，证实了其与人类酶的功能相似性。

在虚拟筛选方面，文献23描述了针对**人源**CYP17A1的大规模筛选，发现了苯并咪唑类化合物（如Compound 2）具有与阿比特龙相当的活性。这些研究为猪CYP17A1的筛选提供了重要的配体参考库。

### **4.2 推荐的虚拟筛选流程**

综合上述分析，本报告推荐以下针对猪CYP17A1的分级筛选流程：

#### **第一阶段：配体库准备**

* **来源**：ZINC20数据库（Drug-like子集）、PubChem或专门的类固醇衍生物库。  
* **预处理**：由于CYP17A1的底物具有严格的立体化学特征，必须使用RDKit或OpenBabel生成精确的3D构象，并保留手性中心。

#### **第二阶段：高通量对接（AutoDock Vina）**

* **设置**：使用前述修正电荷（Fe \+2.0）的PDBQT受体文件。  
* **参数**：Exhaustiveness \= 8。  
* **筛选标准**：按结合能排序，选取前10%-20%的化合物。同时增加几何过滤条件：要求配体中必须有一个氮原子距离血红素铁 \< 2.5 Å（针对II型抑制剂）。

#### **第三阶段：精细重打分（GalaxyDock2-HEME）**

* 对Vina筛选出的优胜分子，使用GalaxyDock2-HEME进行重对接。该步骤能有效剔除那些虽然范德华力适配但缺乏金属配位能力的假阳性分子。

#### **第四阶段：分子动力学验证（AMBER）**

* 选取前10-20个候选分子。  
* 使用**Shahrokh参数集**描述血红素，使用**GAFF2**描述小分子配体。  
* 进行100 ns的显式溶剂MD模拟。  
* **评价指标**：计算结合自由能（MM/GBSA），并监测Fe-N距离的稳定性（应稳定在2.1-2.3 Å）。同时监测Asn202与配体羟基/酮基的氢键占有率。

## ---

**5\. 结论与展望**

针对猪CYP17A1抑制剂的虚拟筛选是一项极具挑战但也充满前景的工作。通过使用人源结构（3RUK）作为模板进行高质量同源建模，并结合Shahrokh等人专门开发的血红素力场参数，可以克服该靶点无晶体结构和含金属辅基的两大技术障碍。

本报告强调，单纯依赖通用的对接软件极易导致失败，必须在原子电荷层面（针对Vina）或评分函数层面（使用GalaxyDock）对血红素铁进行特殊处理。通过这一精细化的计算策略，有望筛选出能够高效抑制猪CYP17A1、从而解决公猪膻味问题的新型先导化合物，实现“计算驱动的动物福利改善”。

### **关键资源清单**

* **同源建模**：Swiss-Model (Template: 3RUK)  
* **MD参数**：Shahrokh et al. (2012) AMBER Heme Parameters 或 MCPB.py  
* **对接软件**：AutoDock Vina (需手动修正Fe电荷) / GalaxyDock2-HEME  
* **关键位点**：Asn202 (氢键锚定), Cys442-Fe (轴向配位)

---

**参考文献引证说明**：本报告中引用的标记对应用户提供的原始研究片段。例如，8指向Shahrokh等人的参数文献，17指向GalaxyDock2-HEME的相关研究。

#### **引用的著作**

1. Pork Production with Entire Males: Directions for Control of Boar Taint \- PMC \- NIH, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC7552340/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7552340/)  
2. A Homodimer Model Can Resolve the Conundrum as to How Cytochrome P450 Oxidoreductase and Cytochrome b5 Compete for the Same Binding Site on Cytochrome P450c17 | Request PDF \- ResearchGate, 访问时间为 十二月 4, 2025， [https://www.researchgate.net/publication/311820884\_A\_Homodimer\_Model\_Can\_Resolve\_the\_Conundrum\_as\_to\_How\_Cytochrome\_P450\_Oxidoreductase\_and\_Cytochrome\_b5\_Compete\_for\_the\_Same\_Binding\_Site\_on\_Cytochrome\_P450c17](https://www.researchgate.net/publication/311820884_A_Homodimer_Model_Can_Resolve_the_Conundrum_as_to_How_Cytochrome_P450_Oxidoreductase_and_Cytochrome_b5_Compete_for_the_Same_Binding_Site_on_Cytochrome_P450c17)  
3. A Missense Mutation in the Human Cytochrome b5 Gene causes 46,XY Disorder of Sex Development due to True Isolated 17,20 Lyase Deficiency, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC3388247/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3388247/)  
4. CYP17A1 \- Steroid 17-alpha-hydroxylase/17,20 lyase \- Sus scrofa (Pig) | UniProtKB, 访问时间为 十二月 4, 2025， [https://www.uniprot.org/uniprotkb/P19100/entry](https://www.uniprot.org/uniprotkb/P19100/entry)  
5. Human P450 CYP17A1: Control of Substrate Preference by Asparagine 202 \- PMC \- NIH, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC5801141/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5801141/)  
6. Identification and Analysis of CYP7A1, CYP17A1, CYP20A1, 访问时间为 十二月 4, 2025， [https://www.jstage.jst.go.jp/article/jvms/advpub/0/advpub\_14-0313/\_article/-char/en](https://www.jstage.jst.go.jp/article/jvms/advpub/0/advpub_14-0313/_article/-char/en)  
7. New AMBER force field parameters of heme iron for cytochrome P450s determined by quantum chemical calculations of simplified models \- PubMed, 访问时间为 十二月 4, 2025， [https://pubmed.ncbi.nlm.nih.gov/15812779/](https://pubmed.ncbi.nlm.nih.gov/15812779/)  
8. Quantum mechanically derived AMBER-compatible heme parameters for various states of the cytochrome P450 catalytic cycle \- PubMed, 访问时间为 十二月 4, 2025， [https://pubmed.ncbi.nlm.nih.gov/21997754/](https://pubmed.ncbi.nlm.nih.gov/21997754/)  
9. 3RUK: Human Cytochrome P450 CYP17A1 in complex with Abiraterone \- RCSB PDB, 访问时间为 十二月 4, 2025， [https://www.rcsb.org/structure/3RUK](https://www.rcsb.org/structure/3RUK)  
10. 6WR1: Human steroidogenic cytochrome P450 17A1 mutant N52Y with inhibitor abiraterone \- RCSB PDB, 访问时间为 十二月 4, 2025， [https://www.rcsb.org/structure/6WR1](https://www.rcsb.org/structure/6WR1)  
11. Comparative Modeling and Molecular Dynamics Simulation of Substrate Binding in Human Fatty Acid Synthase \- PubMed Central, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC4394236/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4394236/)  
12. Tutorial: Molecular Docking intro, 访问时间为 十二月 4, 2025， [https://www.kfc.upol.cz/wp-content/uploads/2023/12/KFC\_DD\_Tutorial\_docking\_vz08.pdf](https://www.kfc.upol.cz/wp-content/uploads/2023/12/KFC_DD_Tutorial_docking_vz08.pdf)  
13. Please can you tell me how can I add charge \+3 to Fe ion using AutoDock tool or any relevant program? | ResearchGate, 访问时间为 十二月 4, 2025， [https://www.researchgate.net/post/Please\_can\_you\_tell\_me\_how\_can\_I\_add\_charge\_3\_to\_Fe\_ion\_using\_AutoDock\_tool\_or\_any\_relevant\_program](https://www.researchgate.net/post/Please_can_you_tell_me_how_can_I_add_charge_3_to_Fe_ion_using_AutoDock_tool_or_any_relevant_program)  
14. How do I dock antifungal drugs to a HEME-containing protein? \- ResearchGate, 访问时间为 十二月 4, 2025， [https://www.researchgate.net/post/How-do-I-dock-antifungal-drugs-to-a-HEME-containing-protein](https://www.researchgate.net/post/How-do-I-dock-antifungal-drugs-to-a-HEME-containing-protein)  
15. Combined Machine Learning and Molecular Modelling Workflow for the Recognition of Potentially Novel Fungicides \- MDPI, 访问时间为 十二月 4, 2025， [https://www.mdpi.com/1420-3049/25/9/2198](https://www.mdpi.com/1420-3049/25/9/2198)  
16. GalaxyDock2-HEME: Protein-ligand docking for heme proteins \- PubMed, 访问时间为 十二月 4, 2025， [https://pubmed.ncbi.nlm.nih.gov/36809651/](https://pubmed.ncbi.nlm.nih.gov/36809651/)  
17. GalaxyDock2‐HEME : P rotein–ligand docking for heme proteins \- ResearchGate, 访问时间为 十二月 4, 2025， [https://www.researchgate.net/publication/368698736\_GalaxyDock2-HEME\_P\_rotein-ligand\_docking\_for\_heme\_proteins](https://www.researchgate.net/publication/368698736_GalaxyDock2-HEME_P_rotein-ligand_docking_for_heme_proteins)  
18. (PDF) Quantum Mechanically Derived AMBER-Compatible Heme Parameters for Various States of the Cytochrome P450 Catalytic Cycle \- ResearchGate, 访问时间为 十二月 4, 2025， [https://www.researchgate.net/publication/51715715\_Quantum\_Mechanically\_Derived\_AMBER-Compatible\_Heme\_Parameters\_for\_Various\_States\_of\_the\_Cytochrome\_P450\_Catalytic\_Cycle](https://www.researchgate.net/publication/51715715_Quantum_Mechanically_Derived_AMBER-Compatible_Heme_Parameters_for_Various_States_of_the_Cytochrome_P450_Catalytic_Cycle)  
19. Quantum mechanically derived AMBER-compatible heme parameters for various states of the cytochrome P450 catalytic cycle \- PMC \- NIH, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC3242737/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3242737/)  
20. AMBER parameter database (Bryce Group: Computational Biophysics and Drug Design \- University of Manchester), 访问时间为 十二月 4, 2025， [http://amber.manchester.ac.uk/](http://amber.manchester.ac.uk/)  
21. MCPB.py: A Python Based Metal Center Parameter Builder | Request PDF \- ResearchGate, 访问时间为 十二月 4, 2025， [https://www.researchgate.net/publication/295910556\_MCPBpy\_A\_Python\_Based\_Metal\_Center\_Parameter\_Builder](https://www.researchgate.net/publication/295910556_MCPBpy_A_Python_Based_Metal_Center_Parameter_Builder)  
22. Molecular Dynamics Simulations of a Cytochrome P450 from Tepidiphilus thermophilus (P450-TT) Reveal How Its Substrate-Binding Channel Opens \- PubMed Central, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC8231624/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8231624/)  
23. Non-steroidal CYP17A1 Inhibitors: Discovery and Assessment \- PMC \- PubMed Central, 访问时间为 十二月 4, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC10226049/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10226049/)  
24. Synthesis and Structure–Activity Relationships of Novel Non-Steroidal CYP17A1 Inhibitors as Potential Prostate Cancer Agents \- MDPI, 访问时间为 十二月 4, 2025， [https://www.mdpi.com/2218-273X/12/2/165?type=check\_update\&version=1](https://www.mdpi.com/2218-273X/12/2/165?type=check_update&version=1)