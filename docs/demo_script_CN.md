# Lab meeting demo 脚本（2026-11-06，约 10 分钟）

## 演示前检查（周四晚 / 周五早）
- [ ] 终端里 `cd` 到本目录，运行 `Rscript scripts/launch_app.R`，确认浏览器能打开应用。
- [ ] 用 demo 数据从头到尾点一遍（下面的流程），确认每页都出图。
- [ ] 备用：现成的运行记录 `runs/2026-11-04_201047_demo_backup/report.html`（卡住时直接打开）；周四也可以自己再写一份。
- [ ] 关掉占内存的程序（整个流程峰值约 6 GB 内存）；浏览器缩放 110%；关闭通知。
- [ ] 准备一句话回答「数据会不会上传」：不会，全部在本机 R 进程里运行，浏览器只是界面。

## 开场（1 分钟）
- 导师的需求：把 PhD 的分析链打包，让不写代码的同事也能用同一套方法。
- 它是什么：本地运行的 Shiny 应用；输入是 DIA-NN（或 Spectronaut、FragPipe）的输出；统计引擎是一组普通 R 函数，界面只是调用它们。
- 和现成工具（FragPipe-Analyst 等）的区别：peptide 级 MS-DAP 路径（vwmb + MaxLFQ + DEqMS）、配对 / 区组设计、每次分析自动生成可复现的运行记录、物种直系同源映射、rescue 分层。

## 页面 1 Data（1.5 分钟）
1. 点 **Fill in the demo dataset** → **Import**（约 7 秒）。
2. 看顶部数字：36 runs（Long-term 阶段四组）、肽段数、蛋白数；格式 DIA-NN。
3. 点 **Template with the run names**：说明模板已预填 run 名，同事只需填 group / subject_id。
4. 点 **Validate and attach** → 右侧出现每组 run 数与动物数；提到校验会拦住常见错误（run 名对不上、组为空、重复）。

## 页面 2 Quality control（2 分钟）
1. 点 **Run QC**（约 30 秒）。等待时讲图：缺失率、强度居中、PCA、相关矩阵、上样顺序漂移、批次、同一动物双眼一致性、离群筛查。
2. 结果：离群筛查列出候选 run（两项以上标记：CNTA_L_4L、LIM_L_1R、LIM_L_3L）；CNTA_L_4L 与论文一致；候选只是提示，不会自动排除。
3. 演示排除：在 **Exclude runs** 选 L_CTRL_4 与 L_CTRL_BLK 两只动物的四个 run（CNT_L_4L、CNTA_L_4R 等），说明排除会写进运行记录。

## 页面 3 Differential analysis（3 分钟）
1. 加三个对比：
   - LIM_L vs CNT_L，pair = none（非配对，消化批次混杂 → 只能当 supportive）
   - LIMA_L vs LIM_L，pair = subject_id（同一动物双眼配对，confirmatory）
   - CNTA_L vs CNT_L，pair = subject_id（药物单独作用）
2. 打开 **Parameters**：默认值就是论文锁定的设置（2 肽段、50% 检出、vwmb、MaxLFQ、DEqMS、FC 1.2、p < 0.05），全部可改且会记录。
3. 点 **Run all contrasts**（约 45 秒）。
4. 看 summary：tested / hits 与论文一致（LIM_L vs CNT_L 1985；LIMA_L vs LIM_L 47；CNTA_L vs CNT_L 303）；limma 一致性 r 约 1。
5. 看 volcano，点一个红点 → 右侧 per-run 图（配对对比有连线）；表格可只看 hits，可下载 CSV。

## 页面 4 Rescue（1.5 分钟）
1. disease = LIM_L vs CNT_L；treatment = LIMA_L vs LIM_L；drug-only = CNTA_L vs CNT_L。
2. 点 **Classify** → Tier 1 = 15，Tier 2 = 180（与论文一致）；定义就印在旁边；drug-only 只标注不过滤。
3. 图：横轴疾病效应、纵轴治疗效应，虚线是 FC 阈值，点线是完全逆转。

## 页面 5 Enrichment（1 分钟）
1. 选 LIM_L vs CNT_L；直系同源表已预填（豚鼠 → 人）。
2. 点 **Run Hallmark GSEA**（约 10 秒）→ NES 柱图；提一句排序用的是模型 t 统计量。

## 页面 6 Run record（1 分钟）
1. 填 label（如 `labmeeting_demo`）和 notes，点 **Write run record**（约 30 秒）。
2. 点 **Open HTML report**：输入文件 MD5、参数、版本、所有表和图都在一个文件夹里，任何人可以复查。

## 收尾（1 分钟）
- 第一版的边界：DIA-NN 导入已用真实数据验证（demo 是论文 Long-term 的 36 个 run，数字与论文一致）；Spectronaut / FragPipe 导入已接入但没有数据测试；计算是同步的（QC 约 30 秒会等待）；富集目前只有 Hallmark。
- 下一步：请同事用自己的数据试；Windows 安装说明；异步计算；GO / Reactome；Spectronaut 实测数据；名字征集。
- 问大家：你们各自的数据类型和最常做的比较是什么？（决定第二版的导入器和功能优先级）

## 常见问题的回答
- 「为什么 p 不是 q？」：默认沿用论文（n = 10 配对，q 过严；q 列保留，可一键切换门槛）。
- 「配对是怎么处理的？」：subject_id 进入模型作为区组因子；只用两组都有的动物。
- 「结果和论文一样吗？」：一样。demo 数据是论文 Long-term 的 36 个 run，三个对比和 rescue 的数字与论文逐项核对过（golden test）；完整 76 样本数据在开发时也核对过。
- 「这个 demo 数据会随软件给大家吗？」：不会。demo 数据不在代码仓库里；给同事的版本计划换成公开数据集。
- 「数据安全？」：本机运行，不联网，不上传。
