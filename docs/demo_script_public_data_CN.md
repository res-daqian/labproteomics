# 公开数据版 demo 脚本（备用；PXD079891，第三版界面）

## 演示前检查（周四晚 / 周五早）
- [ ] 终端里 `cd` 到本目录，运行 `Rscript scripts/launch_app.R`（或双击 `launch_mac.command`），确认浏览器能打开应用。
- [ ] 确认 `demo_data/report.parquet` 在（164024858 字节）。不在的话运行 `Rscript scripts/get_demo_data.R`（从 PRIDE 下载约 164 MB）。
- [ ] 用 demo 数据从头到尾点一遍（下面的流程），确认每页都出图。
- [ ] 备用：现成的运行记录 `runs/<时间戳>_demo_backup_v3/report.html`（卡住时直接打开）。
- [ ] 万一第三版出问题：`git checkout demo-2026-11-06-safe` 然后 `R CMD INSTALL --no-docs --no-multiarch .`，回到周三的第一版（demo 后 `git checkout main`）。注意第一版的 demo 按钮指向的是论文数据，现在已移到 `golden_data/`，回退时要把 `golden_data/` 的文件复制回 `demo_data/`。
- [ ] 关掉占内存的程序；浏览器缩放 100–110%；关闭通知。
- [ ] 准备一句话回答「数据会不会上传」：不会，全部在本机 R 进程里运行，浏览器只是界面。

## 数据背景（开场时一句话带过）
公开数据集 PRIDE PXD079891（Zhang 等，J Physiol 2026）：人脐静脉内皮细胞（HUVEC），8 位供体（4 女 4 男），每位供体的细胞分别转染对照 ASO 和针对 lncRNA MANCR 的 ASO；16 个 run，Orbitrap Exploris 480，DIA-NN 1.9。对比 `MANCR_KD vs Control` 按供体配对。作者在差异分析前剔除了供体 SZF3（PCA 离群）；平台的 QC 会自己找到这一对。

## 开场（1 分钟）
- 导师的需求：把 PhD 的分析链打包，让不写代码的同事也能用同一套方法。
- 它是什么：本地运行的 Shiny 应用；输入是 DIA-NN（或 Spectronaut、FragPipe）的输出；统计引擎是一组普通 R 函数，界面只是调用它们。
- 和现成工具的区别：QC 分三层（设计、单个 run、队列）而且先于统计；peptide 级 MS-DAP 路径（vwmb + MaxLFQ + DEqMS）；配对 / 区组设计；每个 hit 带留一法稳定性；每次分析自动生成可复现的运行记录和外部工具上传表。
- demo 用公开数据，不是实验室自己的数据；同事拿到的软件里也只有这套公开数据。

## 页面 1 Data（2 分钟）
1. 点 **Fill in the demo dataset** → **Import**（约 4 秒）。
2. 看顶部瓷砖：16 runs、109301 个肽段、8573 个蛋白；右下角「Search engine settings」自动读出 DIA-NN 1.9 的版本和设置。没有 FASTA 时基因名和蛋白名直接取自 DIA-NN 报告。
3. 点 **Template with the run names**：模板已预填 run 名，同事只需填 group / subject_id / batch / run_order。
4. 点 **Validate and attach** → 左边是校验，右边是 **Design check**：
   - 两组各 8 个 run；batch（这里是供体性别，女 1–8 号、男 9–16 号两段采集）与组别交叉，可以分离；
   - 8 位供体都有两组的 run，可以配对；
   - 没有 pooled QC（橙色提醒；以后实验要带）。
   - 页面只有五页：Rescue 需要「疾病 vs 对照」和「治疗 vs 疾病」两个对比，这套数据只有两组，所以那一页自动隐藏；换成四组数据它会出现。

## 页面 2 Quality control（3 分钟）
1. 组织已自动选成 Cultured endothelial cells (HUVEC)，点 **Run QC**（约 35 秒）。等待时讲三层：run 成绩单（仪器）、marker 面板（解剖 / 培养物组成）、队列（归一化之后）。
2. **Run scorecard**：14 个指标，每个 run 一个点，批次内 3 MAD 标记。SZF3NC 和 SZF3KD 在鉴定数、信号量、归一化稳定性、肽段长度上同时出界（两个以上标记）；整个队列漏切率中位数 0.33，平台把它报成「全队列消化不完全」而不是逐个 run 标红；批次间系统差异：鉴定数、MS1 质量误差。
3. **Marker panel**：内皮 marker 占总信号 0.5–0.9%（expected）；角蛋白在 SZF4KD 冲到 1.3%（handling 污染）；Y 染色体蛋白只在男性供体的 run 里出现（性别核对，可抓样品调换）。说明不同组织（视网膜、RPE/脉络膜、巩膜、泪液）有各自的角色设置，名单可以一起改。
4. 队列层：缺失率（SZF3 两个 run 约 26%，其他 7–12%）、居中、PCA（SZF3NC 在 PC1 上孤悬，SZF3KD 次之；其余按组分开）、相关矩阵、离群筛查（两项以上标记 = SZF3KD、SZF3NC）、漂移、批次、同一供体一致性。
5. 演示排除：在 **Exclude runs** 选 SZF3NC 和 SZF3KD，说明这和作者的决定一致，而且是在看差异分析之前、凭 QC 证据决定的；排除会写进运行记录。

## 页面 3 Differential analysis（3 分钟）
1. 加两个对比：
   - MANCR_KD vs Control，pair = subject_id（同一供体配对，主分析）
   - MANCR_KD vs Control，pair = none，label 填 `MANCR_KD_vs_Control_unpaired`（对照：不配对会损失多少）
2. 打开 **Parameters**：默认值就是论文锁定的设置；门槛可以选 p、q 或 limma treat；留一法可关。
3. 点 **Run all contrasts**（每个约 35 秒）。
4. Summary（排除 SZF3 后 7 对）：配对 5978 tested，1324 hits（q < 0.05 的 987；treat 31；fragile 441；detection-only 38）；非配对 1149 hits（q 下 429）。配对比非配对多出的就是供体间差异被吸收掉的部分。
5. 火山图：空心点 = 去掉一位供体就不显著的 hit。点最上面的 MLKL 或 MASTL → 右侧 per-run 图（配对连线、定量次数、留一法结果）。
6. 下面的 **Detection-only differences**：一组有、一组没有的蛋白（38 个），不混进火山图，用 Fisher 检验单独列。

## 页面 4 Enrichment（1.5 分钟）
1. 选 MANCR_KD_vs_Control；GSEA + Hallmark，物种已是 Homo sapiens、不用直系同源表，点 **Run enrichment**（约 1 秒）→ 9 个集合 padj < 0.05：干扰素 α / γ 反应、MYC targets、EMT、TNFA/NF-κB、炎症反应、UPR 全部在敲低后下调；排序用模型 t 统计量。
2. 切到 ORA + Reactome，再点一次：背景集 = 本对比定量到的蛋白，不是整个基因组（这是富集分析最常见的错误）。GO、KEGG、细胞类型集合都在下拉里。豚鼠数据则勾选 ortholog 表。

## 页面 5 Run record（1 分钟）
1. 填 label（如 `labmeeting_demo`）和 notes，点 **Write run record**（约 10 秒）。
2. 点 **Open HTML report**：输入 MD5（和 PRIDE 上的文件一致）、搜库设置、设计检查、成绩单、marker、QC、每个对比（含 treat / 留一法 / 检出差异）、富集、上传表清单都在一个文件夹里；`adapters/` 里是 IPA / iPathwayGuide / STRING / Cytoscape 的上传表。

## 收尾（1 分钟）
- 边界：DIA-NN 导入已用两套真实数据验证（公开数据；论文数据的 golden 测试与论文数字一致）；Spectronaut / FragPipe 导入已接入但没有数据测试，成绩单目前只认 DIA-NN 的 stats 文件；计算是同步的；留一法只重拟合模型，不重做归一化。
- 下一步：请同事用自己的数据试；每种组织的 marker 名单一起定；表型关联；Spectronaut 的公开数据；Windows 安装；名字征集。
- 问大家：你们的数据类型、常做的比较、有没有 pooled QC？

## 常见问题的回答
- 「为什么 p 不是 q？」：默认沿用论文（配对 n = 7–10 时 q 过严；q 列和 treat 列都在，可一键切换门槛）。
- 「配对是怎么处理的？」：subject_id 进入模型作为区组因子；只用两组都有的供体。
- 「留一法说 441 个 fragile，是不是假阳性？」：不是假阳性的判定，是证据等级：n = 7 时一个 hit 依赖单个供体很常见，报告时要写明。
- 「Rescue 页去哪了？」：需要三组以上（疾病、对照、治疗）才有意义，两组数据自动隐藏；引擎和测试都还在。
- 「这个 demo 数据会随软件给大家吗？」：会，它是公开的（PRIDE CC0），`scripts/get_demo_data.R` 一条命令下载。实验室的数据不在软件里。
- 「数据安全？」：本机运行，不联网（下载 demo 数据除外），不上传。
