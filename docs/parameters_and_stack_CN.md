# 平台的技术栈与每个参数的来龙去脉（只给 Daqian 看，2026-11-05 晚）

两部分：先讲整个东西是怎么搭起来的（前端、后端、分析引擎、数据怎么流），再逐页讲每个参数：它是什么、用在哪一步、为什么默认这样、改了会怎样、什么时候该改。所有默认值都能在界面上改，并且写进每次运行记录的 `parameters.json`。

---

## 第一部分：技术栈

### 一句话
一个 **R 包**（`labproteomics`），里面有一组普通的 R 函数（分析引擎）和一个 **Shiny** 应用（界面）。运行时只有一个 R 进程在本机跑，浏览器只是显示器。没有 Python，没有数据库，没有服务器，不联网。

### 前端
- **Shiny**（R 包）：R 代码描述界面，Shiny 在浏览器里渲染成网页，并把按钮、下拉框的值传回 R。你在页面上点一下，R 里的一个函数就被调用。
- **bslib**：Bootstrap 5 的布局与主题（六页导航、侧栏、卡片、瓷砖）。配色、字体是一段 CSS（`R/app.R` 里的 `app_css()`），没有下载网络字体，用的是系统字体栈。
- **DT**：表格用 DataTables（JavaScript 库），所以能排序、搜索、翻页、选行。
- 图不是前端画的：ggplot2 在 R 里画成 PNG（ragg 渲染，Arial），浏览器只显示图片。火山图的点击是 Shiny 的 `plotOutput(click=)` 把坐标传回 R，R 用 `nearPoints` 找最近的点。
- 唯一的一小段 JavaScript：Rescue 页隐藏后给导航重新编号。

### 后端
- 就是那个 R 进程（`shiny::runApp()`）。每个页面是一个 Shiny module（`R/mod_*.R`），模块之间通过一个共享的 `reactiveValues`（`state`）传递数据集、QC 结果、对比结果。
- 计算是**同步**的：点「Run QC」后 R 把整个 QC 算完再把结果交给页面，期间显示进度条（`withProgress`）。单用户、每步半分钟左右，所以暂时不需要异步。
- 所有统计都在 `R/` 的引擎函数里（`import_quant`、`attach_annotation`、`design_check`、`run_scorecard`、`marker_scores`、`qc_global`、`run_contrast`、`rescue_classify`、`ranked_genes`、`gsea_collection`、`ora_collection`、`export_adapters`、`write_run_record`）。Shiny 代码里不写统计，这是 CLAUDE.md 里的硬规则，目的是引擎可以脱离界面在脚本里复用，也方便测试。

### 分析引擎依赖的包（版本记录在每次运行记录的 `package_versions.json`）
| 环节 | 包 | 用到的函数 / 作用 |
|---|---|---|
| 读 DIA-NN 报告 | **msdap** 1.2.2（Koopmans 等，GitHub）、**arrow**（读 parquet） | `import_dataset_diann()`：读 `report.parquet`，按 q 值标记检出，每个修饰序列只留一个电荷态，歧义肽段删除；基因名由本包从报告的 `Genes` 列补上（msdap 不读它） |
| FASTA（可选） | msdap | `import_fasta()`：蛋白描述与基因名 |
| 过滤、归一化、蛋白定量 | msdap | `filter_dataset()`：肽段过滤 → **vwmb** 归一化 → **MaxLFQ** 蛋白 roll-up |
| 差异分析 | msdap + **DEqMS** + **limma** | `dea(dea_algorithm = "deqms")`：limma `lmFit` + `eBayes`，再用 DEqMS `spectraCounteBayes` 按肽段数校正方差；p 用 `sca.p`，q 用 BH |
| limma 对照与 treat | limma | `eBayes(trend = TRUE)` 做一致性检查；`treat(lfc = log2(1.2))` 直接检验倍数阈值 |
| 留一法 | 本包 `deqms_fit()`（镜像 msdap 的 DEqMS 调用） | 在固定的蛋白矩阵上逐一去掉一个 subject / run 重拟合，不重做归一化 |
| 检出差异 | R 自带 `fisher.test`、`p.adjust` | 检出次数的 Fisher 精确检验 |
| 设计检查与成绩单 | R 自带 `kruskal.test`、`mad`、`median` | 批次 × 组列联表、批次内 MAD 标记、批次间 / 组间 Kruskal-Wallis |
| 队列 QC | R 自带 `cor`、`prcomp`、`hclust` | 相关矩阵、PCA、聚类排序 |
| 富集 | **fgsea**（`fgseaMultilevel`、`fora`）、**msigdbr** + msigdbdf（离线基因集） | GSEA 按 DEqMS t 统计量排序；ORA 用本对比定量到的蛋白作背景 |
| 直系同源映射 | 本包 `mapped_symbols()` | 豚鼠 accession → 人类 symbol，一对一行，重复 symbol 取 |t| 最大 |
| 图 | **ggplot2**、ggrepel、patchwork；**ragg**、**systemfonts**（Arial）；cairo（PDF 嵌字体） | 主题在 `R/theme.R` |
| 记录与报告 | jsonlite、digest（MD5）、**rmarkdown + knitr + pandoc**（HTML 报告）、**chromote**（Chrome 打印 PDF）、**rmzqc**（HUPO-PSI mzQC） | `write_run_record()`、`render_qc_report()`、`write_mzqc()` |
| 测试 | testthat | 单元测试 + golden 测试（论文数字必须复现） |

### 数据怎么流
1. `import_quant()` → msdap 数据集 `ds`：`ds$peptides`（每行一个肽段 × run：log2 强度、RT、电荷、置信度、是否检出）、`ds$proteins`、`ds$samples`。
2. `attach_annotation()` → 注释表并进 `ds$samples`（group、subject_id、batch、run_order、额外列），标记 exclude 的 run 直接移除。
3. QC：`qc_global()` 对**所有** run 做一次宽松过滤 + 归一化 + roll-up，得到 QC 用的蛋白矩阵；成绩单读 `report.stats.tsv`；marker 面板用**未归一化**的肽段强度。
4. 每个对比：`run_contrast()` 只取这两组的 run，按参数过滤 + 归一化 + roll-up，再 DEqMS；然后加 treat 列、检出差异、留一法。
5. Rescue、富集在对比结果表上做。
6. `write_run_record()` 把一切写进 `runs/<时间戳>_<label>/`。

### 为什么是 R 而不是 Python
PhD 的锁定流程是 R（MS-DAP、DEqMS、limma、fgsea 都是 R / Bioconductor 生态），直接搬最保真；Shiny 让界面和引擎在同一个进程、同一种语言里，不用维护两套环境；同事装一个 R 加几个包就能跑。以前 PhD 项目里那个 Python 的元数据检查，已经用 R 重写在 `validate_annotation()` 里。

### 打包与安装
R 包的标准结构：`DESCRIPTION`（依赖清单）、`NAMESPACE`（导出函数）、`R/`、`inst/report/`（报告模板）、`tests/`。`R CMD INSTALL .` 或 `remotes::install_local()` 装进用户的 R 库；`launch_mac.command` 只是调用 `Rscript scripts/launch_app.R`。给同事的安装：`scripts/install_dependencies.R` 一次装齐依赖（CRAN + Bioconductor + GitHub 上的 msdap），再装本包。

---

## 第二部分：参数逐个讲

格式：**参数（默认值）**。是什么 / 用在哪一步 / 为什么这样设 / 改了会怎样 / 什么时候改。

### 页面 1 Data

**Format（auto）**。输入文件的来源软件。auto 按扩展名和表头判断。只影响用哪个 msdap 导入函数；Spectronaut、FragPipe 没用真实数据测过。

**FASTA（可选）**。给蛋白加描述和基因名。有就用；没有时 DIA-NN 报告自带的 `Genes`、`Protein.Names` 列补上，marker 面板和富集都靠基因名，所以二者至少要有一个。豚鼠数据必须给 FASTA，因为报告里的基因名不全（RHO、GNAT1 没有）。

**Precursor confidence, q-value（0.01）**。msdap 的 `confidence_threshold`：precursor 的 q 值 ≤ 0.01 才算「检出」（detect）。用在导入那一步，决定 `ds$peptides$detect`，后面「每组至少检出 1 次」的过滤靠它。DIA-NN 输出本身已按 1% FDR 过滤，所以 0.01 等于「信任 DIA-NN 的 FDR」；调低（如 0.001）会让检出更严、肽段更少、缺失更多；调高没有意义。论文用 0.01。

**注释表的列**。`sample_id` 必须和报告里的 run 名完全一致；`group` 决定能比什么；`subject_id` 决定能不能配对（同一动物的两只眼、同一供体的两种处理）；`batch` 决定 QC 里的分层（成绩单按批次内标记、PCA 用形状区分）和设计检查里的混杂判定；`run_order` 用在漂移图和「顺序是否分块」检验；`exclude = TRUE` 在导入后立刻移除该 run，写进记录。

**Design check 的 QC 识别规则（`pool|qc|hela|standard|blank`）**。按组名或样品名识别 pooled QC / 标准品，识别到的 run 不进对比。改它只在你的命名不同时。

### 页面 2 Quality control

**Minimum peptides per protein（QC 矩阵，1）**。QC 这一步的蛋白矩阵只要求 ≥ 1 个肽段。为什么比差异分析（2）松：QC 想看到尽可能多的蛋白来判断 run 的状态，缺失率、相关性在更大的矩阵上更稳定；差异分析才需要用 2 个肽段保证定量可靠。改成 2 会让 QC 矩阵和 DEA 矩阵一致，缺失率会变低一点，结论一般不变。

**Normalisation（vwmb）**。msdap 的「variation within, mode between」：先在组内让样本间的变异最小，再按组间的众数对齐。这是论文锁定的方法，也是 msdap 对多组设计的推荐。备选：vsn（方差稳定）、loess、rlr（稳健线性回归）、median、msempire。QC 页改它只影响 QC 矩阵；DEA 页有自己的同名参数。什么时候改：组间蛋白总量确实不同（例如组织不一样）时 vwmb 的「组间众数对齐」假设成立；若两组大部分蛋白都变了，所有基于「大多数不变」的归一化都有问题，这时要先讨论。

**Outlier flag threshold, MADs（3）**。三处都用 3 个 MAD：队列层三个指标（缺失率高、相关性低、强度偏移）在组内标记（当每组 ≥ 6 个 run，否则全体一起）；成绩单在批次内（当每个批次 ≥ 4 个 run）；marker 面板全体一起。MAD 乘 1.4826 使之与正态分布的 SD 可比，3 个 MAD 大致对应 3 SD，是 QC 里惯用的「明显离群」尺度。调成 2 会标出更多 run（更敏感、更多假警报），调成 4 更保守。队列层「两项以上才算候选」是为了避免单一指标误判。

**成绩单的两个固定常数**。鉴定数低于中位数的 70%（`id_fraction = 0.7`）：鉴定数是最敏感的综合指标，用相对比例而不是 MAD，因为鉴定数的分布常常偏斜；70% 是经验值。漏切率上限 0.25（`missed_cleavage_max`）：只有当这一层的中位数本身低于 0.25 时才用它标单个 run；整个队列都高于 0.25 时改报「全队列消化不完全」。HUVEC 公开数据中位数 0.33 就是这种情况，你的数据在 0.25 以下。

**Marker 面板：Tissue（按数据选）与 min_share（0.1%）**。组织决定每个集合的角色：expected（低了标）、neighbour / contamination（高了标）、check（只看）。`min_share` 是一个地板：邻近 / 污染集合的份额至少达到总信号的 0.1% 才可能被标，防止 0.001% 这种噪声级别的波动被 MAD 放大成警报。自定义 CSV（set, gene, role）会替换内置名单。

**Exclude runs**。手动选择；在对比里生效，QC 图永远显示所有 run；写入记录的 `excluded_after_qc`。规矩是先看证据再排除、排除后不回头。

### 页面 3 Differential analysis

**Case / Control / Pair runs by**。`case` 是分子，log2FC > 0 表示 case 高。配对列（通常 subject_id）进入 msdap 的模型作为区组因子（design = `~ subject + condition`），只用两组都有该 subject 的运行；没有配对就是普通两组比较。配对的作用是把动物 / 供体间的差异从误差里扣掉，HUVEC 数据配对比非配对多出 175 个 hit 就是这部分。

**Fold-change threshold（1.2，即 |log2FC| ≥ 0.263）**。论文锁定值。它是「生物学上值得看」的最低倍数，不是统计量；1.2 对蛋白组是常见的宽松阈值（很多论文用 1.5 或 2）。用在三处：hit 的判据；limma treat 的原假设（|log2FC| ≤ 0.263）；rescue 的 Tier 2 残差带。调高会减少 hit，treat 会更严，Tier 2 的带更宽。

**Significance threshold（0.05）**与 **Gate（p / q / treat）**。p 门：未校正 p < 0.05 且过倍数阈值，论文的主判据；q 门：BH 校正后的 q < 0.05，控制 FDR，在 n = 10 配对时常常一个都不剩；treat 门：limma treat 的 q，检验的是「倍数是否显著大于阈值」，比 p 门严得多。三种 hit 列同时都在表里，门槛只决定哪一列叫 hit，并决定留一法和 rescue 用哪一列。

**Minimum peptides per protein（2）**。msdap `filter_min_peptide_per_prot`：一个蛋白至少 2 个肽段才定量，单肽段蛋白的 MaxLFQ 值不可靠。论文锁定。改成 1 会多几百个蛋白但方差更大。

**Minimum detection fraction per group（0.5）**。msdap `filter_fraction_quant`：一个肽段要在每组至少 50% 的 run 里有定量值才保留（另外每组至少 1 次检出、1 次定量）。这是 MS-DAP 的默认思路：缺失太多的肽段会把 MaxLFQ 的估计拖坏。调低（0.3）保留更多低丰度肽段、缺失更多；调高（0.8）更干净、蛋白更少。注意它是按组算的，所以组内样本少时颗粒很粗（n = 6 时 0.5 就是 3 个）。

**Normalisation / Model**。同上；Model 可选 DEqMS 或 limma eBayes（不按肽段数校正方差）。DEqMS 是论文锁定方法，它用「肽段数多的蛋白方差小」这一经验关系调整每个蛋白的先验方差，多肽段蛋白的 p 更可信。

**Also run limma-trend（TRUE）**。在同一个蛋白矩阵上再跑一次 limma `eBayes(trend = TRUE)`，给出 log2FC 的相关系数（concordance）和 treat 列。它不改 hit，只提供第二引擎的印证；treat 门槛需要它，所以选 treat 时会强制打开。关掉能省几秒。

**Leave-one-out stability（TRUE，上限 60 个单位）**。配对时按 subject、非配对时按 run 逐一去掉一个，在同一个蛋白矩阵上重拟合 DEqMS（或 treat）；全样本重拟合与 msdap 的 p 一致到 1e-13，这是自检。hit 在任一次重拟合里掉出判据 = fragile。它不重做归一化和 roll-up（成本和意义都不够），所以叫「模型稳定性」而不是「全流程稳定性」。单位超过 60 会拒绝，防止跑太久；可以关掉。

**Detection-only（present ≥ 0.8，absent ≤ 0.2）**。「一组 ≥ 80% 的 run 检出、另一组 ≤ 20%」的蛋白没有可靠的倍数（一边几乎没值），所以不进火山图，用检出次数的 Fisher 精确检验单独列。两个比例是惯用的「几乎全有 / 几乎全无」；n = 6 时 0.8 就是 5 个、0.2 就是 1 个。Fisher 检验忽略配对、功效低，只当线索。

### 页面 4 Rescue

参数全部继承页面 3（倍数、显著性、门槛）；定义：Tier 1 = 疾病对比显著 且 治疗对比显著且方向相反；Tier 2 = 疾病对比显著 且 治疗对比不显著 且 |疾病 log2FC + 治疗 log2FC| < log2(1.2)；Disease-only = 疾病对比显著但两者都不满足；Not perturbed = 疾病对比不显著。drug-only 标记 = 药物单独作用显著 且 方向与治疗作用相同，只标注不过滤。`disease_fragile` / `treatment_fragile` 来自各自对比的留一法。这些定义就是论文 rescue_v2 的 2-case 框架，一个字没改。

### 页面 5 Enrichment

**Method（GSEA / ORA）**。GSEA 用全部定量蛋白按 DEqMS 的 t 统计量排序（`fgseaMultilevel`，`eps = 0`，随机种子 1 保证可复现）；ORA 只取 hit，用 `fora` 做超几何检验，背景 = 本对比定量到且能映射到 symbol 的蛋白。背景用整个基因组是富集分析最常见的错误（Wijesooriya 等 2022），因为蛋白组只测到一部分基因。

**Collection**。Hallmark（50 个，概览）、GO BP/CC/MF（细但冗余）、Reactome、KEGG MEDICUS、C8 细胞类型。都是 msigdbr 自带、离线。

**Hits to test（both / up / down）**，只对 ORA：上调和下调一起测会稀释方向相反的两个通路，分开测更清楚。

**Ortholog table / Species**。豚鼠必须映射到人类 symbol（一对一行，重复取 |t| 最大）；人、小鼠数据直接用自己的 symbol，物种选对应的 msigdbr 集合。

**Minimum / maximum set size（15 / 500；ORA 下限 5）**。太小的集合几个基因就能显著，太大的集合（如「代谢」）没有解释力；15 和 500 是 fgsea 的常规默认。ORA 的下限放到 5 是因为交集后的集合本来就小。

### 页面 6 Run record

**Adapters（TRUE）**：IPA、iPathwayGuide、STRING、Cytoscape 上传表。**HTML report（TRUE）**：需要 pandoc（RStudio 自带）。**QC report（TRUE）**：独立的 QC 报告 HTML + PDF（PDF 需要 Chrome）和 mzQC 文件（需要 rmzqc 包）。

### QC 报告与 mzQC 里有什么、没什么
PDF 报告 = 设计检查 + 三层 QC 的图表 + 每个 run 的判定表 + 排除决定 + 参数与版本。mzQC（HUPO-PSI 的 QC 交换格式，JSON）只能装受控词表里有的指标，所以目前写入每个 run 的：鉴定 precursor 数（MS:4000103）、鉴定 precursor 的平均电荷（MS:4000174）、污染蛋白份额（MS:4000177，取角色为 contamination 的 marker 集合）、批次标签（MS:4000088，元数据）、输入文件与 DIA-NN 版本；其余指标在同一文件夹的 CSV 里。MsQuality 这个包要读谱图（mzML），平台的输入是搜库报告，所以没有接；以后若把 mzML 也交给平台，可以把它的指标并进成绩单和 mzQC。
