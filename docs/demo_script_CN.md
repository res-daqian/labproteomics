# Lab meeting demo 脚本（2026-11-06，约 12 分钟，第三版界面，用论文的 Long-term 数据）

Daqian 2026-11-05 晚决定：周五用自己的数据演示。数据在 `golden_data/`（只在本机，不进 git），用专门的启动入口打开。给同事的出厂 demo 仍是公开数据 PXD079891，那套的脚本在 `docs/demo_script_public_data_CN.md`。

## 演示前检查（周四晚 / 周五早）
- [ ] 双击 `launch_mac_thesis_data.command`（第一次：右键 > 打开），确认浏览器打开应用；左上角空状态提示写着 "fills in the Long-term retina data of the PhD project"。终端等价命令：`Rscript scripts/launch_app.R golden_data`。
- [ ] 用下面的流程从头到尾点一遍，确认每页都出图。
- [ ] 备用：现成的运行记录 `runs/<时间戳>_demo_backup_v3_thesis/report.html`（卡住时直接打开）。
- [ ] 万一第三版出问题：回到周三的第一版。那一版的 demo 按钮读 `demo_data/`，所以要先把 `golden_data/` 里的文件复制到 `demo_data/`（公开数据的文件可以先挪开），再执行下面两行；demo 后 `git checkout main` 并恢复文件。

```
git checkout demo-2026-11-06-safe
R CMD INSTALL --no-docs --no-multiarch .
```

- [ ] 关掉占内存的程序（整个流程峰值约 6 GB）；浏览器缩放 100–110%；关闭通知。
- [ ] 准备一句话回答「数据会不会上传」：不会，全部在本机 R 进程里运行，浏览器只是界面。

## 开场（1 分钟）
- 导师的需求：把 PhD 的分析链打包，让不写代码的同事也能用同一套方法。
- 它是什么：本地运行的 Shiny 应用；输入是 DIA-NN（或 Spectronaut、FragPipe）的输出；统计引擎是一组普通 R 函数，界面只是调用它们。
- 和现成工具的区别：QC 分三层（设计、单个 run、队列）而且先于统计；peptide 级 MS-DAP 路径（vwmb + MaxLFQ + DEqMS）；配对 / 区组设计；每个 hit 带留一法稳定性；每次分析自动生成可复现的运行记录和外部工具上传表。
- 今天用我自己的 Long-term 数据演示；给大家的软件里带的是一套公开数据（人内皮细胞，8 位供体配对），一条命令下载。

## 页面 1 Data（2 分钟）
1. 点 **Fill in the demo dataset** → **Import**（约 7 秒）。
2. 看顶部瓷砖：36 runs、肽段数、蛋白数；右下角「Search engine settings」自动读出 DIA-NN 2.3.1 的版本和设置。
3. 点 **Template with the run names**：模板已预填 run 名，同事只需填 group / subject_id / batch / run_order。
4. 点 **Validate and attach** → 左边是校验，右边是 **Design check**：
   - 红色一条：四个跨批次对比 confounded（CNT/CNTA 在 Date_0721，LIM/LIMA 在 Date_0828）→ 这就是论文里 LIM vs CNT 只能当 supportive 的原因，平台在算任何统计之前就说出来；
   - 两个配对对比在同一批次内；18 只动物可配对；没有 pooled QC（橙色提醒；以后实验要带）。
   - 附上注释后导航多出 **4 Rescue**：这页只在设计有三组以上时出现，两组数据会自动隐藏。

## 页面 2 Quality control（3 分钟）
1. 组织已自动选成 Retina，点 **Run QC**（约 35 秒）。等待时讲三层：run 成绩单（仪器）、marker 面板（解剖）、队列（归一化之后）。
2. **Run scorecard**：14 个指标，每个 run 一个点，批次内 3 MAD 标记。标记的 CNTA_L_4L、LIMA_L_89L 与论文 QC 一致；展开「Flag rules and systematic differences」：两个批次在肽段长度、电荷、保留时间偏差上系统不同，是混杂的直接证据。
3. **Marker panel**：晶状体 0.002%、RPE 0.008%、血红蛋白 0.07% → 解剖干净；只有 CNT_L_4R 骨骼肌偏高。说明不同组织（RPE/脉络膜、巩膜、泪液、培养细胞）有各自的角色设置，名单可以一起改。
4. 队列层：缺失率、居中、PCA（PC1 同时是组别和批次）、相关矩阵（冷色阶）、离群筛查（CNTA_L_4L、LIM_L_1R、LIM_L_3L）、漂移、批次、双眼一致性。
5. 演示排除：在 **Exclude runs** 选 L_CTRL_4 与 L_CTRL_BLK 两只动物的四个 run，说明排除写进运行记录，而且是在看差异分析之前决定的。

## 页面 3 Differential analysis（3 分钟）
1. 加三个对比：
   - LIM_L vs CNT_L，pair = none（非配对，批次混杂 → supportive）
   - LIMA_L vs LIM_L，pair = subject_id（同一动物双眼配对，confirmatory）
   - CNTA_L vs CNT_L，pair = subject_id（药物单独作用）
2. 打开 **Parameters**：默认值就是论文锁定的设置；门槛可以选 p、q 或 limma treat；留一法可关。
3. 点 **Run all contrasts**（约 1 分钟）。
4. Summary：tested / hits 与论文一致（1985；47；303）；多出三列：treat hit、fragile、detection-only。
5. 火山图（青 = 下调，红 = 上调）：空心点 = 去掉一只动物就不显著的 hit（LIMA_L vs LIM_L 47 个里 28 个）。点一个点或选一行 → 右侧 per-run 图（配对连线、定量次数、留一法结果）。
6. 下面的 **Detection-only differences**：一组有、一组没有的蛋白（LIM_L vs CNT_L 83 个），不混进火山图，用 Fisher 检验单独列。

## 页面 4 Rescue（1.5 分钟）
1. disease = LIM_L vs CNT_L；treatment = LIMA_L vs LIM_L；drug-only = CNTA_L vs CNT_L。
2. 点 **Classify** → Tier 1 = 15，Tier 2 = 180（与论文一致）；瓷砖上标出 15 个 Tier 1 里 11 个在某个对比里 fragile；定义印在旁边；drug-only 只标注不过滤。

## 页面 5 Enrichment（1.5 分钟）
1. 选 LIM_L vs CNT_L；GSEA + Hallmark；直系同源表已自动填好并勾选（豚鼠 → 人），物种 Homo sapiens；点 **Run enrichment**（约 2 秒）→ NES 柱图；排序用模型 t 统计量。
2. 切到 ORA + Reactome，再点一次：背景集 = 本对比定量到的 5965 个蛋白，不是整个基因组（这是富集分析最常见的错误）。GO、KEGG、细胞类型集合都在下拉里。

## 页面 6 Run record（1 分钟）
1. 填 label（如 `labmeeting_demo`）和 notes，点 **Write run record**（约 10 秒）。
2. 点 **Open HTML report**：输入 MD5、搜库设置、设计检查、成绩单、marker、QC、每个对比（含 treat / 留一法 / 检出差异）、rescue、富集、上传表清单都在一个文件夹里；`adapters/` 里是 IPA / iPathwayGuide / STRING / Cytoscape 的上传表。

## 收尾（1 分钟）
- 边界：DIA-NN 导入已用两套真实数据验证（本数据的数字与论文一致，公开数据全流程跑通）；Spectronaut / FragPipe 导入已接入但没有数据测试，成绩单目前只认 DIA-NN 的 stats 文件；计算是同步的；留一法只重拟合模型，不重做归一化。
- 下一步：请同事用自己的数据试；每种组织的 marker 名单一起定；表型关联；Spectronaut 的公开数据；Windows 安装；名字征集。
- 问大家：你们的数据类型、常做的比较、有没有 pooled QC？

## 常见问题的回答
- 「为什么 p 不是 q？」：默认沿用论文（n = 10 配对，q 过严；q 列和 treat 列都在，可一键切换门槛）。
- 「配对是怎么处理的？」：subject_id 进入模型作为区组因子；只用两组都有的动物。
- 「留一法说 28 个 fragile，是不是假阳性？」：不是假阳性的判定，是证据等级：n = 10 时一个 hit 依赖单只动物很常见，报告时要写明。
- 「结果和论文一样吗？」：一样。三个对比和 rescue 的数字逐项核对过（golden test）；新加的列只描述 hit，不改变 hit。
- 「这个 demo 数据会随软件给大家吗？」：不会。软件里带的是公开数据集（PRIDE CC0）。
- 「数据安全？」：本机运行，不联网（下载公开 demo 数据除外），不上传。
