# PROGRESS — 当前状态（先读这页）

Last updated: 2026-11-04 晚（demo 数据改为 Long-term 36 样本，数字与论文一致）

## 1. 目标与范围
- 导师要求：把 PhD 的 MS 蛋白组分析流程打包成 **本地运行的 Shiny 应用**，给不写代码的同事用。
- 第一次 demo：**本周五 2026-11-06 lab meeting**（Daqian：做不完也没关系）。
- 平台：Mac + Windows。上游软件：DIA-NN、Spectronaut 为主；ProteinPilot（SWATH）、FragPipe 也有。
- demo 数据：PhD Part II **Long-term 阶段的 36 个样本**（CNT_L 8、CNTA_L 8、LIM_L 10、LIMA_L 10），由 `scripts/make_demo_subset.R` 从 PhD 原始 `report.parquet` 生成（只读、gitignored）。Early 阶段不再放在 demo 里。
- 长期计划（Daqian 2026-11-04 晚）：给同事的出厂 demo 改用公开数据集（DIA-NN / Spectronaut 输出），不再带 Daqian 自己的数据；候选正在调研，demo 后执行。

## 2. 现在做到哪
- **可演示版本已就绪（2026-11-04）**：包 `labproteomics` 已安装；六页界面用 headless Chrome 按 demo 流程逐页点过、截图检查过（Data → QC → 排除两只动物 → 三个 Long-term 对比 → Rescue → Hallmark GSEA → 运行记录 + HTML 报告），浏览器 console 零错误。
- 冒烟测试 + golden testthat 全部通过（8 项；Early 块只在完整 76 样本数据存在时检查）：三个 Long-term 对比与论文数字完全一致；rescue Tier 1 = 15、Tier 2 = 180；drug flag 在两层内 = 2。
- 实测耗时（M4 Pro，36 样本）：导入 4 s；全量 QC 约 30 s；单个对比 8–21 s；GSEA 1.4 s；运行记录 + 报告约 4 s；整个流程峰值内存约 6 GB。
- 备用运行记录：`runs/2026-11-04_201047_demo_backup/report.html`。
- 启动方式：双击 `launch_mac.command`，或 `Rscript scripts/launch_app.R`，或 R 里 `labproteomics::run_app(demo_dir = "demo_data", runs_dir = "runs")`。
- 演示脚本：`docs/demo_script_CN.md`。
- 未做：Windows 机器实测（launcher 已写，未验证）；renv；异步计算；Spectronaut / FragPipe 实测数据；GO / Reactome；公开数据集出厂 demo；名字。

## 3. Golden 数字（必须复现 PhD 锁定结果）
| 对比 | tested | hits p | hits q | 来源 |
|---|---|---|---|---|
| LIM_L vs CNT_L（非配对，排除 L_CTRL_4 / L_CTRL_BLK） | 6044 | 1985 | 1915 | PhD `dea_enhanced_summary.csv` |
| LIMA_L vs LIM_L（配对，同排除） | 6193 | 47 | 0 | 同上 |
| CNTA_L vs CNT_L（配对，6 对） | 6223 | 303 | 0 | 同上 |
| Long-term rescue | Tier 1 = 15，Tier 2 = 180，drug flag = 2 | | | PhD `rescue_v2_summary.md` |
| LIMA_E vs LIM_E（配对；只在完整 76 样本数据上检查） | 5963 | 36 | 0 | PhD `dea_enhanced_summary.csv` |

## 4. 风险 / 已知问题
- msdap 的函数假定其 Depends 已 attach（`is_character` 无命名空间）→ 包内用 `ensure_msdap()` 先 attach；DESCRIPTION 把 msdap 放在 Depends。
- **msdap 导入有两处「全数据集依赖」**（2026-11-04 发现，细节见 `PROVENANCE.md`「Demo data」）：① 每个修饰序列只保留一个电荷态、歧义序列整条删除，都按整个数据集决定，所以子集不能直接按行过滤；② 归一化 / rollup 对输入行顺序数值敏感（同样的肽段、打乱行顺序 → CNTA_L vs CNT_L 有 8 个 p 在 0.02–0.056 之间的边缘蛋白改变状态）。`make_demo_subset.R` 两点都处理了。对论文的含义：结果可复现（输入锁定、顺序固定），但 p 值贴近 0.05 的边缘 hit 本来就不稳，不要单独解读。
- renv 暂未启用（两天内不冒安装风险），版本由 `DESCRIPTION` + 运行记录的 `package_versions.json` 记录；demo 后再做。
- 计算同步 + `withProgress`（QC 约 30 秒会阻塞界面）；异步是后续升级。
- Spectronaut / FragPipe 导入只是对 msdap 导入器的薄包装，**没有数据测试过**，界面上要如实标注。
- GUI 的 QC 离群标记分层用「组」（每组 ≥ 6 个样本时），PhD 用 stage × pair_family；36 样本下两项以上标记的 run = CNTA_L_4L、LIM_L_1R、LIM_L_3L（CNTA_L_4L 与论文一致），这是有意的泛化。

## 5. 当前 vs 已作废
| 项 | 状态 |
|---|---|
| demo 数据 = Long-term 36 样本（`make_demo_subset.R`，保留肽段选择与行顺序） | 当前（2026-11-04 晚）；「全部 76 样本」（2026-11-04 下午）作废；「Early 40 样本子集」作废 |
| 出厂 demo 改公开数据集 | 计划中（候选调研），demo 后执行 |

## 6. 最近变更
- 2026-11-04 下午：项目建立；demo 数据克隆；engine + 界面初版；冒烟测试修复 msdap attach 问题与 group_colours 正则；界面改为非填充布局（图高度正确）、QC 摘要改计数；全流程 headless 检查通过；golden 测试通过；首次 commit。
- 2026-11-04 晚：demo 数据改为 Long-term 36 样本；发现并处理 msdap 的电荷态 / 歧义序列选择与行顺序依赖；golden 测试按数据自动跳过 Early；界面文字、README、演示脚本、CLAUDE.md、PROVENANCE 同步；备用运行记录重建。
