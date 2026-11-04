# PROGRESS — 当前状态（先读这页）

Last updated: 2026-11-05 凌晨（第二版：QC 三层、检出差异、treat、留一法、ORA、导出适配、新界面）

## 1. 目标与范围
- 导师要求：把 PhD 的 MS 蛋白组分析流程打包成 **本地运行的 Shiny 应用**，给不写代码的同事用。
- 第一次 demo：**本周五 2026-11-06 lab meeting**（Daqian：做不完也没关系）。
- 平台：Mac + Windows。上游软件：DIA-NN、Spectronaut 为主；ProteinPilot（SWATH）、FragPipe 也有。
- demo 数据：PhD Part II **Long-term 阶段的 36 个样本**（CNT_L 8、CNTA_L 8、LIM_L 10、LIMA_L 10），由 `scripts/make_demo_subset.R` 从 PhD 原始 `report.parquet` 生成（只读、gitignored）；另加同一 36 个 run 的 `report.stats.tsv` 和完整的 `report.log.txt`（同样 gitignored）。
- 长期计划（Daqian 2026-11-04 晚）：给同事的出厂 demo 改用公开数据集（DIA-NN / Spectronaut 输出），demo 后执行（候选见 memory）。
- Daqian 2026-11-05 回答：pooled QC 以后尽量带，旧数据基本没有；设计各种都有，每组 n 5–10；组织以视网膜、RPE/脉络膜、巩膜、泪液为主；小鸡直系同源表现在没有；同事基本留着 full report。

## 2. 现在做到哪
- **安全点**：git tag `demo-2026-11-06-safe` = 周三晚可演示的第一版（commit 474627d）。周五早上如果第二版有任何问题，可以回到它：`git checkout demo-2026-11-06-safe` 后重新 `R CMD INSTALL .`（看完 demo 再 `git checkout main`）。
- **第二版（2026-11-05 凌晨）已 headless 全流程跑通**，console 零错误，golden 测试全部通过（17 项断言），单元测试 73 项通过。
- 新功能（都已进界面和运行记录）：
  1. 页面 1 **设计检查**：batch × group 列联表；每对组别判定 confounded / crossed / 同一批次；run order 的 Kruskal-Wallis；配对可行性；pooled QC 识别。demo：四个跨批次对比全部 confounded（CNT/CNTA 在 Date_0721，LIM/LIMA 在 Date_0828），两个配对对比在批次内。另读取 DIA-NN 的 `report.log.txt`（版本、命令行、设置）。
  2. 页面 2 **Run 成绩单**：解析 `report.stats.tsv` 的 14 个指标，批次内 3 MAD 标记，另给批次间 / 组间的系统差异（Kruskal-Wallis）。demo：批次内标记 CNTA_L_4L（总量、MS2 信号）和 LIMA_L_89L（MS1 质量误差、漏切）——与论文 76 样本 QC 的候选离群一致；两个批次在肽段长度、电荷、RT 偏差、峰宽上系统不同（p ≈ 1e-6），是混杂的直接证据。
  3. 页面 2 **组织纯度 / 污染 marker 面板**：12 个集合（光感受器、RPE、晶状体、血细胞、血浆、巩膜/脉络膜 ECM、平滑肌、骨骼肌、皮肤角蛋白、眼表上皮、泪液、Y 染色体），角色随组织（视网膜 / RPE-脉络膜 / 巩膜 / 泪液 / 其他）变化，可上传自定义 CSV。demo：解剖很干净（晶状体 0.002%、RPE 0.008%、血红蛋白 0.07%），只有 CNT_L_4R 骨骼肌比例偏高。
  4. 页面 3 **检出差异（detection-only）**：一组 ≥ 80% 检出、另一组 ≤ 20% 的蛋白单独列表，Fisher 精确检验（不考虑配对，标注为探索性）。demo：LIM_L vs CNT_L 83 个，CNTA_L vs CNT_L 24 个，LIMA_L vs LIM_L 3 个。
  5. 页面 3 **limma treat** 列（treat_p / treat_q）与第三种门槛「treat q」；默认门槛不变。demo：LIM_L vs CNT_L treat q < 0.05 有 687 个；两个配对对比 0 个。
  6. 页面 3 **留一法稳定性**：在固定蛋白矩阵上逐一去掉一只动物（配对）或一个 run（非配对）重拟合 DEqMS（全样本重拟合与 msdap 的 p 值一致到 1e-13）。hit 在任一次重拟合中丢失 = fragile，火山图上空心点。demo：LIMA_L vs LIM_L 47 个 hit 中 28 个 fragile；CNTA_L vs CNT_L 303 中 217；LIM_L vs CNT_L 1985 中 384；15 个 Tier 1 里 11 个在两个对比之一 fragile。
  7. 页面 5 **GSEA 扩展到 Hallmark / GO BP、CC、MF / Reactome / KEGG MEDICUS / C8 细胞类型**（msigdbr 自带数据，离线可用；会话内缓存）；**ORA** 用本对比定量到的蛋白作背景（fgsea::fora），可选只测上调或下调。
  8. 页面 6 **导出适配**：每个对比写 IPA、iPathwayGuide、STRING（hit 列表 + 带值表）、Cytoscape 节点表（含 rescue tier）和 manifest。
  9. **界面重做**：白底、细边框、小写字距标题、统计瓷砖、步骤编号导航、状态胶囊、每页一句说明、空状态提示；黑-灰-红不变。
- 实测耗时（M4 Pro，36 样本）：导入 4 s；QC（成绩单 + marker + 归一化）约 35 s；单个对比 12–27 s（含留一法；可在参数里关闭）；GSEA Hallmark 1 s、Reactome 3 s、GO BP 5 s；ORA < 1 s；运行记录 + 报告约 10 s。
- 备用运行记录：`runs/<时间戳>_demo_backup_v2/report.html`（第二版，headless 生成）；第一版的在 `runs/2026-11-04_201047_demo_backup/`。
- 启动方式不变：双击 `launch_mac.command`，或 `Rscript scripts/launch_app.R`。
- 演示脚本：`docs/demo_script_CN.md`（已按第二版更新）。QC 讲义：`docs/qc_primer_CN.md`。
- 未做：表型关联 / WGCNA（需要表型表和配对感知模型，demo 后讨论）；msqrob2 第二引擎（需装 Bioconductor 包，demo 前不动）；小鸡直系同源表；Windows 实测；renv；异步计算；Spectronaut / FragPipe 实测数据（成绩单目前只认 DIA-NN 的 stats 文件）；公开数据集出厂 demo。

## 3. Golden 数字（必须复现 PhD 锁定结果）
| 对比 | tested | hits p | hits q | 来源 |
|---|---|---|---|---|
| LIM_L vs CNT_L（非配对，排除 L_CTRL_4 / L_CTRL_BLK） | 6044 | 1985 | 1915 | PhD `dea_enhanced_summary.csv` |
| LIMA_L vs LIM_L（配对，同排除） | 6193 | 47 | 0 | 同上 |
| CNTA_L vs CNT_L（配对，6 对） | 6223 | 303 | 0 | 同上 |
| Long-term rescue | Tier 1 = 15，Tier 2 = 180，drug flag = 2 | | | PhD `rescue_v2_summary.md` |
| LIMA_E vs LIM_E（配对；只在完整 76 样本数据上检查） | 5963 | 36 | 0 | PhD `dea_enhanced_summary.csv` |

第二版新增的回归值（不是论文数字，只防止无意改动）：留一法全样本重拟合与 msdap p 值差 < 1e-8；LIMA_L vs LIM_L fragile + robust = 47；LIM_L vs CNT_L 检出差异 83、treat hit 687；LIMA_L vs LIM_L treat hit 0。

## 4. 风险 / 已知问题
- msdap 的函数假定其 Depends 已 attach（`is_character` 无命名空间）→ 包内用 `ensure_msdap()` 先 attach；DESCRIPTION 把 msdap 放在 Depends。
- **msdap 导入有两处「全数据集依赖」**（2026-11-04，细节见 `PROVENANCE.md`「Demo data」）：① 每个修饰序列只保留一个电荷态、歧义序列整条删除，都按整个数据集决定；② 归一化 / rollup 对输入行顺序数值敏感。`make_demo_subset.R` 两点都处理了。对论文的含义：结果可复现，但 p 值贴近 0.05 的边缘 hit 本来就不稳——第二版的留一法把这一点量化了（见上）。
- 留一法只重拟合统计模型，不重做归一化和 rollup；界面和报告都写明。
- 检出差异的 Fisher 检验忽略配对；n = 6–10 时功效很低，只能当线索。
- marker 面板的基因名单是起点；每种组织的名单要和 Daqian 一起定（RPE/脉络膜、巩膜、泪液尤其）。豚鼠 FASTA 里 RHO、GNAT1 没有基因名，光感受器集合靠其他 33 个蛋白。
- 成绩单只认 DIA-NN 的 `report.stats.tsv`；Spectronaut 的 pipeline report 还没接。
- renv 暂未启用；版本由 `DESCRIPTION` + 运行记录的 `package_versions.json` 记录。
- 计算同步 + `withProgress`；异步是后续升级。
- GUI 的队列 QC 离群标记分层用「组」（每组 ≥ 6 个样本时），PhD 用 stage × pair_family；36 样本下两项以上标记的 run = CNTA_L_4L、LIM_L_1R、LIM_L_3L（CNTA_L_4L 与论文一致），这是有意的泛化。

## 5. 当前 vs 已作废
| 项 | 状态 |
|---|---|
| 第二版界面 + QC 三层 + DEA 三项 + ORA + 导出（2026-11-05 凌晨） | 当前 |
| 第一版六页界面（2026-11-04） | 被第二版取代；tag `demo-2026-11-06-safe` 保留作为回退点 |
| demo 数据 = Long-term 36 样本（`make_demo_subset.R`） | 当前；「全部 76 样本」「Early 40 样本子集」作废 |
| 出厂 demo 改公开数据集 | 计划中，demo 后执行 |

## 6. 最近变更
- 2026-11-04 下午：项目建立；demo 数据克隆；engine + 界面初版；全流程 headless 检查通过；golden 测试通过；首次 commit。
- 2026-11-04 晚：demo 数据改为 Long-term 36 样本；msdap 的电荷态 / 歧义序列选择与行顺序依赖处理；文档同步；备用运行记录重建。
- 2026-11-05 凌晨：QC 讨论后的第二版（见第 2 节 1–9）；新单元测试 4 个文件；golden 测试扩展；README / PROVENANCE / CLAUDE.md / 演示脚本 / QC 讲义更新；两次 commit（engine；界面 + 文档）。
