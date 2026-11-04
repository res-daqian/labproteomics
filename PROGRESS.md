# PROGRESS — 当前状态（先读这页）

Last updated: 2026-11-05 凌晨（第三版：公开数据 demo、新配色与 Arial 图、新界面字体与配色、Rescue 页按设计自动隐藏）

## 1. 目标与范围
- 导师要求：把 PhD 的 MS 蛋白组分析流程打包成 **本地运行的 Shiny 应用**，给不写代码的同事用。
- 第一次 demo：**本周五 2026-11-06 lab meeting**（Daqian：做不完也没关系）。
- 平台：Mac + Windows。上游软件：DIA-NN、Spectronaut 为主；ProteinPilot（SWATH）、FragPipe 也有。
- **demo 数据（2026-11-05 起）：公开数据集 PRIDE PXD079891**（Zhang 等，J Physiol 2026；CC0）：HUVEC，8 位供体（4 女 4 男）配对，对照 ASO vs MANCR ASO，16 个 run，DIA-NN 1.9。`scripts/get_demo_data.R` 从 PRIDE 下载 `report.parquet`（164024858 字节）、stats、log、SDRF，并生成注释表（`demo_data/README.md`）。Daqian 2026-11-05：「放弃用我自己的数据」。
- **论文数据只留给 golden 测试**：Long-term 36 个样本的子集移到 `golden_data/`（整个目录 gitignored，只读），由 `scripts/make_golden_subset.R` + `make_golden_annotation.R` 重建；不随软件分发。
- Daqian 2026-11-05 回答：pooled QC 以后尽量带，旧数据基本没有；设计各种都有，每组 n 5–10；组织以视网膜、RPE/脉络膜、巩膜、泪液为主；小鸡直系同源表现在没有；同事基本留着 full report；界面要更美观、字体要换；demo 数据不支持 rescue 时可暂时去掉该模块。

## 2. 现在做到哪
- **安全点**：git tag `demo-2026-11-06-safe` = 周三晚可演示的第一版（commit 474627d）。回退：`git checkout demo-2026-11-06-safe` 后 `R CMD INSTALL --no-docs --no-multiarch .`（看完 demo 再 `git checkout main`）。注意第一版的 demo 按钮指向论文数据，回退时要把 `golden_data/` 里的文件复制回 `demo_data/`。
- **第三版（2026-11-05 凌晨）已在公开数据上 headless 全流程跑通**：console 零错误；单元测试 73 项通过；golden 测试 17 项通过（论文数字不变）。
- 第三版改动：
  1. **demo 数据换成 PXD079891**（见第 1 节）。注释表：`group` = Control / MANCR_KD；`subject_id` = 供体；`batch` = 供体性别（女 1–8 号、男 9–16 号两段采集，是唯一能当批次的分组）；`sex` 同值；无预先排除。
  2. **导入器**：没有 FASTA 时，基因名和蛋白名直接取自 DIA-NN 报告的 `Genes` / `Protein.Names` 列（msdap 不读这两列）。人类 Swiss-Prot FASTA 从 UniProt / EBI 三条路都下载失败（SSL 错误），为可选项；有网时 `get_demo_data.R` 会再试。
  3. **图形主题**（`R/theme.R`）：Daqian 给的四套配色（qualitative、diverging、sequential 冷 / 暖）；Arial（systemfonts 检测，ragg 画 PNG，cairo 嵌字体写 PDF）；只留左、下两条细轴，无网格、无边框；柱宽 0.5–0.55；样品点用空心圆；图例在上方。语义：红 = 病例 / 上调 / 标记，青 = 对照 / 下调 / 主系列。
  4. **界面主题**（`R/app.R`）：正文 Avenir Next / Segoe UI（系统字体栈，不下载网络字体），标题与瓷砖数字用 Charter / Palatino / Georgia 衬线；工作色青 #015069，红只用于警示和 hit；瓷砖顶部青色细线；消息点按级别着色（绿 / 蓝 / 橙 / 红）；导航不换行，状态合并成一枚胶囊。
  5. **Rescue 页按设计自动隐藏**：附上注释后，组别 < 3 时 `nav_hide`，≥ 3 时 `nav_show`；步骤编号随可见页面重排（小段 JS）。引擎 `rescue.R` 和 golden 测试不变。
  6. **成绩单规则**：漏切率的绝对限 0.25 只在该层中位数低于 0.25 时标单个 run；整个队列都高时改报「全队列消化不完全」（本数据中位数 0.33）。
  7. **marker 面板**：新增内皮细胞集合；新组织选项「Cultured endothelial cells (HUVEC)」；demo 按钮通过 `demo_settings.json` 自动选组织和物种。
  8. 报告里 p / q 值改为三位有效数字（原来四位小数把 1e-7 显示成 0）。
- **公开数据上的结果（排除 SZF3 供体后，7 对）**：成绩单两项以上标记 = SZF3NC、SZF3KD（作者也剔除了这位供体）；批次间系统差异：鉴定数、MS1 质量误差；marker：内皮 0.5–0.9%，SZF4KD 角蛋白 1.3%，Y 染色体蛋白只在男性 run；队列 QC 候选离群 = SZF3KD、SZF3NC；配对对比 5978 tested、1324 hits（q 987、treat 31、fragile 441、detection-only 38），非配对 1149 hits（q 429）；GSEA Hallmark 9 个集合显著（干扰素 α/γ、MYC targets、EMT、TNFA/NF-κB、炎症、UPR 下调）；ORA Reactome 3 个。
- 实测耗时（M4 Pro，16 样本）：导入 4 s；QC 约 35 s；单个对比约 35 s（含留一法）；GSEA 1–2 s；运行记录 + 报告约 10 s。
- 备用运行记录：`runs/<时间戳>_demo_backup_v3/report.html`（第三版，headless 生成，含排除 SZF3）。
- 启动方式不变：双击 `launch_mac.command`，或 `Rscript scripts/launch_app.R`。演示脚本：`docs/demo_script_CN.md`（已按第三版重写）。QC 讲义：`docs/qc_primer_CN.md`（例子仍来自论文数据，开头已注明）。
- 第二版（2026-11-05 凌晨早些时候）的内容仍在：设计检查、run 成绩单、marker 面板、检出差异、limma treat 门槛、留一法、GSEA/ORA 多集合、导出适配、运行记录与报告。
- 未做：表型关联 / WGCNA；msqrob2 第二引擎；小鸡直系同源表；Windows 实测；renv；异步计算；Spectronaut / FragPipe 实测数据（成绩单只认 DIA-NN 的 stats 文件；候选 Spectronaut 公开数据 PXD033653，775 MB）；人类 FASTA 的下载；若要 Windows 上字体完全一致，需打包一款 OFL 网络字体（需 Daqian 同意下载）。

## 3. Golden 数字（必须复现 PhD 锁定结果；数据在 `golden_data/`）
| 对比 | tested | hits p | hits q | 来源 |
|---|---|---|---|---|
| LIM_L vs CNT_L（非配对，排除 L_CTRL_4 / L_CTRL_BLK） | 6044 | 1985 | 1915 | PhD `dea_enhanced_summary.csv` |
| LIMA_L vs LIM_L（配对，同排除） | 6193 | 47 | 0 | 同上 |
| CNTA_L vs CNT_L（配对，6 对） | 6223 | 303 | 0 | 同上 |
| Long-term rescue | Tier 1 = 15，Tier 2 = 180，drug flag = 2 | | | PhD `rescue_v2_summary.md` |
| LIMA_E vs LIM_E（配对；只在完整 76 样本数据上检查） | 5963 | 36 | 0 | PhD `dea_enhanced_summary.csv` |

第二版新增的回归值（不是论文数字）：留一法全样本重拟合与 msdap p 值差 < 1e-8；LIMA_L vs LIM_L fragile + robust = 47；LIM_L vs CNT_L 检出差异 83、treat hit 687；LIMA_L vs LIM_L treat hit 0。运行：`NOT_CRAN=true Rscript -e 'testthat::test_file("tests/testthat/test-golden.R")'`（先 `pkgload::load_all()`）。

## 4. 风险 / 已知问题
- msdap 的函数假定其 Depends 已 attach（`is_character` 无命名空间）→ 包内用 `ensure_msdap()` 先 attach；DESCRIPTION 把 msdap 放在 Depends。
- **msdap 导入有两处「全数据集依赖」**（细节见 `PROVENANCE.md`「Golden data」）：① 每个修饰序列只保留一个电荷态、歧义序列整条删除，都按整个数据集决定；② 归一化 / rollup 对输入行顺序数值敏感。golden 子集两点都处理了；公开数据的 parquet 原样使用（运行记录里的 MD5 就是 PRIDE 文件的 MD5）。
- 留一法只重拟合统计模型，不重做归一化和 rollup；界面和报告都写明。
- 检出差异的 Fisher 检验忽略配对；n = 6–10 时功效很低，只能当线索。
- marker 面板的基因名单是起点；每种组织的名单要和 Daqian 一起定（RPE/脉络膜、巩膜、泪液尤其）。
- 成绩单只认 DIA-NN 的 `report.stats.tsv`；Spectronaut 的 pipeline report 还没接。
- demo 注释里把供体性别放在 `batch` 列（README 里说明了理由）；严格说那不是处理批次。
- 字体：Mac 上 Avenir Next / Charter，Windows 上退到 Segoe UI / Palatino Linotype；图中 Arial 两边都有。
- renv 暂未启用；版本由 `DESCRIPTION` + 运行记录的 `package_versions.json` 记录。计算同步 + `withProgress`；异步是后续升级。

## 5. 当前 vs 已作废
| 项 | 状态 |
|---|---|
| 第三版：公开数据 demo + 新配色 / 字体 + Rescue 自动隐藏（2026-11-05 凌晨） | 当前 |
| 第二版：QC 三层 + DEA 三项 + ORA + 导出 + 白底界面（2026-11-05 凌晨） | 内容保留，界面与 demo 数据被第三版取代 |
| 第一版六页界面（2026-11-04） | 被取代；tag `demo-2026-11-06-safe` 保留作为回退点 |
| demo 数据 = 论文 Long-term 36 样本 | 作废（移到 `golden_data/`，只作 golden 测试） |
| 黑-灰-红配色 | 作废，换成 Daqian 的配色 |

## 6. 最近变更
- 2026-11-04 下午：项目建立；demo 数据克隆；engine + 界面初版；全流程 headless 检查通过；golden 测试通过；首次 commit。
- 2026-11-04 晚：demo 数据改为 Long-term 36 样本；msdap 的电荷态 / 歧义序列选择与行顺序依赖处理；文档同步；备用运行记录重建。
- 2026-11-05 凌晨（早）：QC 讨论后的第二版；新单元测试 4 个文件；golden 测试扩展；两次 commit（engine；界面 + 文档）。
- 2026-11-05 凌晨（晚）：第三版：公开数据 PXD079891 作 demo（论文数据移到 golden_data）；导入器从 DIA-NN 报告取基因名；新配色 + Arial 图形主题；新界面字体与配色；Rescue 页自动隐藏；漏切率规则改为队列级；内皮 marker 集合；报告 p 值格式；README / PROVENANCE / CLAUDE.md / 演示脚本 / QC 讲义更新；三次 commit（数据与引擎；图形与界面；文档）。
