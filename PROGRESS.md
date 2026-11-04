# PROGRESS — 当前状态（先读这页）

Last updated: 2026-11-04（第一天）

## 1. 目标与范围
- 导师要求：把 PhD 的 MS 蛋白组分析流程打包成 **本地运行的 Shiny 应用**，给不写代码的同事用。
- 第一次 demo：**本周五 2026-11-06 lab meeting**（Daqian：做不完也没关系）。
- 平台：Mac + Windows。上游软件：DIA-NN、Spectronaut 为主；ProteinPilot（SWATH）、FragPipe 也有。
- demo 数据：PhD Part II 全部 76 个样本（只读克隆到 `demo_data/`），演示走 Long-term 对比。

## 2. 现在做到哪
- 包骨架 `labproteomics`（临时名）+ engine 九个文件 + 六页 Shiny 界面已写出。
- 冒烟测试：导入 7 s；注释校验通过；golden 对比核对中（见 §3）。
- 未做：安装成包并用浏览器逐页检查；README / 安装说明；launcher；commit。

## 3. Golden 数字（必须复现 PhD 锁定结果）
| 对比 | tested | hits p | hits q | 来源 |
|---|---|---|---|---|
| LIMA_E vs LIM_E（配对） | 5963 | 36 | 0 | PhD `dea_enhanced_summary.csv` |
| LIM_L vs CNT_L（非配对，排除 L_CTRL_4 / L_CTRL_BLK） | 6044 | 1985 | 1915 | 同上 |
| LIMA_L vs LIM_L（配对，同排除） | 6193 | 47 | 0 | 同上 |
| CNTA_L vs CNT_L（配对，6 对） | 6223 | 303 | 0 | 同上 |
| Long-term rescue | Tier 1 = 15，Tier 2 = 180，drug flag = 2 | | | PhD `rescue_v2_summary.md` |

## 4. 风险 / 已知问题
- msdap 的函数假定其 Depends 已 attach（`is_character` 无命名空间）→ 包内用 `ensure_msdap()` 先 attach；DESCRIPTION 把 msdap 放在 Depends。
- renv 暂未启用（两天内不冒安装风险），版本由 `DESCRIPTION` + 运行记录的 `package_versions.json` 记录；demo 后再做。
- 计算同步 + `withProgress`（全量 QC 约 1 分钟会阻塞界面）；异步是后续升级。
- Spectronaut / FragPipe 导入只是对 msdap 导入器的薄包装，**没有数据测试过**，界面上要如实标注。
- GUI 的 QC 离群标记分层用「组」（每组 ≥ 6 个样本时），PhD 用 stage × pair_family；结果可能略有差异，这是有意的泛化。

## 5. 当前 vs 已作废
| 项 | 状态 |
|---|---|
| demo 数据 = 全部 76 样本，演示 Long-term | 当前（2026-11-04）；「Early 40 样本子集」作废 |

## 6. 最近变更
- 2026-11-04：项目建立；demo 数据克隆；engine + 界面初版；冒烟测试修复 msdap attach 问题。
