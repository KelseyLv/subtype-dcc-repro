# Subtype-DCC 复现会话纪要（压缩版）

> 便于下次打开仓库后直接接续工作。项目根：`subtype-dcc-repro`。

---

## 〇、项目目录结构（速览）与近期进度摘要

### 〇.1 顶层目录在干什么（「乱」多来自 `vendor/` + `repro_out/`）

```
subtype-dcc-repro/
├── .venv/                 # Python 虚拟环境（始终在仓库根使用）
├── repro/                 # 本仓库「复现与作图」脚本（评测、SNF 对接、KIRC 论文图等）
├── repro_out/             # 运行产物：指标 JSON、SNF TSV、KIRC 图与缓存（易膨胀，可定期扫）
├── vendor/Subtype-DCC/    # 上游 Subtype-DCC 训练代码与 results/*.dcc、*.fea（体积极大）
├── data/clinical/         # 临床随访等（如 cBioportal 风格 TSV）
├── benchmarks/            # R baseline（SNF 为主）：模板与最小安装脚本
├── docs/                  # 论文 PDF、本会话纪要 SESSION_HANDOFF
├── scripts/               # 数据下载 / Windows 环境辅助
├── external/              # 可选外部仓库快照（如 Subtype-GAN）
├── requirements.txt
├── README.md
└── subtype-dcc-repro-plan.md
```


| 目录                        | 建议心智模型                                                          |
| ------------------------- | --------------------------------------------------------------- |
| `**repro/**`              | 唯一需要经常打开的「自己的代码」；Python 入口几乎都在此。                                |
| `**repro_out/**`          | **一次性输出**：可整夹备份/清理；**不要**当源码改。子目录含义见下表。                         |
| `**vendor/Subtype-DCC/`** | **上游子项目**：训练、checkpoint、`results/*.dcc`；batch 训练时 `cwd` 常设在此树下。 |
| `**data/`**               | 输入数据；与 `repro/data_io.py` 的加载路径一致。                              |
| `**docs/`**               | 人读文档 + 论文 PDF；**SESSION_HANDOFF** 为接续工作的总目录。                    |


### 〇.2 `repro_out/` 子结构（最容易堆杂）


| 路径                                 | 内容                                                                                                                 |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `*.metrics.json`                   | 各癌种 Subtype-DCC（或自指定标签）的 `evaluate.py` 结果。                                                                         |
| `batch_metrics_summary.json`       | `batch_evaluate.py` 汇总。                                                                                            |
| `*.snf.tsv` / `*.snf.metrics.json` | R SNF 输出的簇标签与评测；`snf_table1_compare.json` 为与论文 Table 1 SNF 行对照。                                                    |
| `kirc_manuscript_cache/`           | **KIRC 论文 Fig 4–7 的中间结果**（`Z.npy`、`go_fig5.csv`、`kegg_fig6.csv` 等）；**先 `kirc_manuscript_compute_cache.py` 再单图重绘**。 |
| `kirc_manuscript_figs/`            | **当前主用的 KIRC 论文风格图**（`KIRC_fig4A_tsne.png` … `fig7`）。                                                              |
| `kirc_paper_figs/`                 | 早期试验输出（CSV、Enrichr 日志等）；**根目录重复 PNG 已删**，若不再需要可整夹删除。                                                               |
| `kirc_manuscript_cache_test/`      | 若存在，多为调试残留，**确认无用可删**。                                                                                             |
| `KIRC_train.log` 等 `*.log`         | 训练或批处理日志。                                                                                                          |


### 〇.3 近期工作进度摘要（跨会话，便于对齐「做到哪了」）

- **评测与临床**：修复 `data_io` 合并临床表时的索引对齐问题后，KIRC 等可正常出 **log-rank** 与 **6 项临床显著数**；与论文 Table 1 的数值对齐策略见 **§三、§七**。
- **SNF baseline**：以 CRAN `SNFtool` 为主路径（`benchmarks/snf_template.R` + `install_snf_only.R`）；九癌种批量跑通，汇总 `**repro_out/snf_table1_compare.json`**（§六–七）。
- **Subtype-DCC 训练**：批量入口 `repro/run_subtype_dcc_batch.py`；产物在 `**vendor/Subtype-DCC/Subtype-DCC/results/`** 与 `save/`（§九）。若部分癌种未完成，可只传未完成的 `--cancers` 续跑后再 `batch_evaluate.py`。
- **KIRC 论文风格 Fig 4–7（本地 Subtype-DCC 结果）**：实现 **t-SNE、GO（BP/CC/MF 三色条形）、KEGG 点图、miRNA–通路热图（相关法代理）**；**Fig 4B** 为与 4A 同坐标的 2×2 标志基因散点；代码拆为 `**kirc_manuscript_lib.py` + `compute_cache` + 各 `plot_kirc_fig*.py`**（§十一），避免每次改图都重跑 Enrichr。
- **仓库整洁**：已删除 `repro_out` 下与 `**kirc_manuscript_figs`** 重复的根目录 `**figure4_KIRC.png`** 及 `**kirc_paper_figs/*.png`**；保留 `**kirc_manuscript_figs`** 为当前主图输出。
- **队友换数据集**：见 `**docs/CUSTOM_DATASET_SUBTYPE_DCC.md`**（数据源建议、工作流、`fea`/临床格式、仓库建议补充项）。仍用 TCGA、癌种未定：见 `**docs/TCGA_SUBTYPE_DCC_GENERIC.md`**。从 **UCSC Xena `.gz` 转 `*.fea`** 的清洗与甲基化特征选择：`**docs/XENA_TO_FEA_PREPROCESSING.md`**。
- **论文 Table 1 / Fig.2 vs 本仓库箱线图（`repro/plot_figure1_table1_boxplots.py`）**：原文性能对比在 **Fig.2**（**非 Fig.1**），为 **Subtype-DCC + 十种** 基线的 **(A) −log₁₀ log-rank p、(B) 临床显著项个数**；**Table 1 加粗**表示 **该癌种该格该方法最优**。Results/Discussion 据此写多数据集上的优越性，**不是**「任意五法箱线图的中位数必最高」。本复现箱线仅为 **五法 × 九癌种** + **GerkeLab 随访、单次训练**（论文对波动方法写 **5 次平均**），故 **中位数看不出 Subtype-DCC 优势并不否定论文表述**；对齐叙事可做 **逐癌种夺冠** 或 **与 Fig.2 同口径方法集**。更短摘录：`docs/handoff_sessions.md`。
- **LIHC 扩展 Subtype-DCC（数据 / Xena 预处理 / 训练 cwd / loss·NaN 排障 / 超参策略）**：见 `**docs/LIHC_EXTENSION_SUBTYPE_DCC.md`**（一站式；预处理细节仍链到 `XENA_TO_FEA_PREPROCESSING.md`）。

---

## 一、背景与目标

- **课题**：健康管理课程相关，复现论文 *Subtype-DCC: decoupled contrastive clustering…*（Briefings in Bioinformatics, 2023），论文 PDF 位于上级目录课程资料中。
- **数据**：多组学特征来自 Subtype-GAN / Subtype-DCC 管线下的 `*.fea`（通过 junction 指向 `vendor/Subtype-DCC/subtype_file/fea`）；临床随访使用 `**data/clinical/cBioportal_data.tsv`**（GerkeLab / cBioPortal 风格）。
- **环境**：仓库根目录 `.venv`，依赖见 `requirements.txt`。
- **训练性质**：`train.py` 为 **无监督**，无 train/test 划分；对某癌种全部样本训练后导出 `**results/<CANCER>.dcc`**（簇标签）与 `**<CANCER>.fea`**（嵌入）。

---

## 二、要复现 / 常用的代码模块


| 模块                     | 路径                                                          | 作用                                                                                                                   |
| ---------------------- | ----------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| 训练                     | `vendor/Subtype-DCC/Subtype-DCC/train.py`                   | 单癌种训练；checkpoint：`save/<CANCER>/checkpoint_last.tar`；输出：`results/<CANCER>.dcc`、`.fea`                                |
| 评测（论文 Table 1 指标）      | `repro/evaluate.py`                                         | −log₁₀(log-rank *p*)、6 项临床检验显著个数                                                                                     |
| 批量评测（9 癌种）             | `repro/batch_evaluate.py`                                   | 循环跑 `evaluate.py` 逻辑；汇总 `repro_out/batch_metrics_summary.json`                                                       |
| **Table 1 全方法评测（5×9）** | `repro/batch_evaluate_table1_methods.py`                    | Subtype-DCC `.dcc` + `repro_out/*.snf.tsv` + kmeans/spectral/nmf TSV；写出 `repro_out/table1_reproduction_metrics.json` |
| **与论文 Table 1 对照**     | `repro/table1_compare_to_paper.py`                          | 读 `data/paper_table1_reference.json`（SNF 行已录入；其余可从 PDF 补全）；生成 `repro_out/table1_reproduction_vs_paper.{json,md}`     |
| 临床与 omics IO           | `repro/data_io.py`                                          | 加载临床表（legacy PANCAN / GerkeLab merged）；`resolve_fea_root`                                                            |
| Kaplan–Meier           | `repro/plot_survival.py`                                    | 按簇的生存曲线                                                                                                              |
| Figure 4 风格图（旧版）       | `repro/plot_figure4.py`                                     | **A**：t-SNE；**B**：离散散点（红/蓝）；与 **KIRC manuscript** 的 2×2 t-SNE 标志物图（`plot_kirc_fig4b.py`）是两套产物，勿混用                    |
| KIRC manuscript 图      | `repro/kirc_manuscript_*.py`、`plot_kirc_fig*.py`（§〇、§十一）    | 论文 Fig 4–7 风格；**先缓存再单图重绘**                                                                                           |
| SNF vs Table 1         | `repro/compare_snf_table1.py`                               | 读取 `repro_out/<CANCER>.snf.metrics.json`，写出 `snf_table1_compare.json`                                                |
| R：SNF baseline         | `benchmarks/snf_template.R`、`benchmarks/install_snf_only.R` | CRAN `SNFtool`；输出与 `evaluate.py` 兼容的 TSV                                                                             |
| Subtype-DCC 批量训练       | `repro/run_subtype_dcc_batch.py`                            | 按列表 **串行** 调用 `train.py`；`cwd` 固定为 `vendor/Subtype-DCC/Subtype-DCC`                                                  |


### 常用命令（在仓库根目录）

```powershell
cd "D:\LateInUniversity\Term10\健康管理中的智能数据分析\subtype-dcc-repro"

# 训练（**必须先 cd 到 train.py 所在目录**，否则相对路径 `config/config.yaml` 会报 FileNotFoundError）
Set-Location "D:\LateInUniversity\Term10\健康管理中的智能数据分析\subtype-dcc-repro\vendor\Subtype-DCC\Subtype-DCC"
..\..\..\.venv\Scripts\python.exe train.py -c KIRC --workers 4

# 新癌种示例 LIHC（不在 cancer_dict 内，必须传 --cluster_number）：
# ..\..\..\.venv\Scripts\python.exe train.py -c LIHC --cluster_number 3 --workers 4

# 等价：在仓库根用批量脚本（内部已 cwd 到 Subtype-DCC\Subtype-DCC）：
# Set-Location "D:\LateInUniversity\Term10\健康管理中的智能数据分析\subtype-dcc-repro"
# .\.venv\Scripts\python.exe repro\run_subtype_dcc_batch.py --cancers LIHC --extra-args "--cluster_number 3 --workers 4"

# 单癌种评测
.\.venv\Scripts\python.exe repro\evaluate.py `
  --labels vendor\Subtype-DCC\Subtype-DCC\results\KIRC.dcc `
  --clinical data\clinical\cBioportal_data.tsv -c KIRC `
  --out-json repro_out\KIRC.metrics.json

# 九癌种 Subtype-DCC 批量评测（仅 .dcc）
.\.venv\Scripts\python.exe repro\batch_evaluate.py

# 九癌种 × 五种标签源（Subtype-DCC + SNF + sklearn 三法）→ 与论文 Table 1 对照表
.\.venv\Scripts\python.exe repro\batch_evaluate_table1_methods.py
.\.venv\Scripts\python.exe repro\table1_compare_to_paper.py

# Figure 4 风格图（KIRC，旧版 plot_figure4；输出路径可自定）
.\.venv\Scripts\python.exe repro\plot_figure4.py -c KIRC --out repro_out\figure4_KIRC.png

# KIRC 论文 Fig 4–7（推荐）：见 §十一（先 compute_cache 再 plot_kirc_fig*）
```

### `train.py` 增强（若已合并进 vendor）

- `**--log_every**`：每隔 N epoch 打印 loss（默认 10）。
- **早停**：plateau / patience 相关参数；推理默认加载 `**checkpoint_last.tar`**。
- Loss 为 **epoch 内 batch loss 之和**，绝对数值随样本数变化，主要看趋势。

---

## 三、问题排查记录

### 1. PowerShell 找不到 `.venv\Scripts\python.exe`；或 `train.py` 报 `config/config.yaml` 找不到

- **venv 路径**：在 `vendor\Subtype-DCC\Subtype-DCC` 下需 **向上三级** 才到仓库根，故为 `..\..\..\.venv\Scripts\python.exe`（误用 `..\..\` 只到 `vendor`，找不到 `.venv`）。
- `**config/config.yaml`**：`train.py` 用 相对路径 读配置与 `dataloader` 的 `../subtype_file/fea`，工作目录必须是 `vendor\Subtype-DCC\Subtype-DCC`。若在仓库根执行 `python vendor\...\train.py`，会触发 FileNotFoundError: config/config.yaml。正确做法见 §二「常用命令」（先 `Set-Location` 到该子目录再跑），或用 `**repro\run_subtype_dcc_batch.py`**（内部已 `cwd` 到该目录）。

### 2. 评测结果 `n_significant_clinical: 0` 且六项 `p` 全为 `1.0`

- **原因**：`repro/data_io.py` 中 `_load_clinical_merged` 构建 `out` 时，把带 **原始整数索引** 的 `Series` 放进 `pd.DataFrame(..., index=patient_barcode)`，pandas **按索引对齐**导致性别/AJCC 等列 **整列 NaN**；`days`/`status` 因使用 `.values` 仍正常，故 log-rank 可算、临床项不可算。
- **修复**：增加 `**_col_vals_aligned`**，用 `**Series.to_numpy()`** 按 **行顺序** 写入，与 `patient_barcode` 索引一一对应。
- **修复后 KIRC 示例**：`neg_log10_logrank_p ≈ 6.24`，`**n_significant_clinical = 6`**（与论文「6 项」计数一致）；论文 Table 1 KIRC 的 **−log₁₀ p = 8.79** 仍可能因 **随访定义/数据源/PyTorch 随机性** 不完全一致。

### 3. 与论文 Table 1 对比时要注意

- **临床项数**：对齐 Gerke 表并修复 IO 后可复现「几项显著」。
- **log-rank 数值**：需与论文同款随访表与预处理才可能接近。

### 4. Figure 4B 初版用热图显得「连续型」

- **期望**：论文 Figure 4B 为 **离散红/蓝** 展示「本簇高、他簇低」。
- **修复**：`plot_figure4.py` 改为 **两张散点图**（x = 簇标签 + 抖动，y = RNA 表达；红 = 宿主亚型，蓝 = 其他），布局为左 t-SNE + 右上下各一子图。

### 5. 其他癌种是否还会遇到「临床全 NaN」？

- **训练**：不涉及临床表，无此问题。
- **评测**：使用 Gerke merged 路径时，**当前对齐修复对所有癌种通用**；若某表缺列，仍可能单项不可用，但不会是「索引对齐」类整列 NaN。

---

## 四、关键产物路径速查


| 产物               | 路径                                                                                                 |
| ---------------- | -------------------------------------------------------------------------------------------------- |
| 模型 checkpoint    | `vendor/Subtype-DCC/Subtype-DCC/save/<CANCER>/checkpoint_last.tar`                                 |
| 簇标签 / 嵌入         | `vendor/Subtype-DCC/Subtype-DCC/results/<CANCER>.dcc`、`.fea`                                       |
| 评测 JSON          | `repro_out/<CANCER>.metrics.json`、`batch_metrics_summary.json`                                     |
| Figure 4（旧散点脚本）  | `repro/plot_figure4.py` 的 `--out` 路径自定（曾用 `repro_out/figure4_<CANCER>.png`）                        |
| KIRC 论文图 Fig 4–7 | `repro_out/kirc_manuscript_figs/KIRC_fig*.png`；中间结果见 `repro_out/kirc_manuscript_cache/`（§〇、§十一）    |
| 多组学输入            | `vendor/Subtype-DCC/subtype_file/fea/<CANCER>/*.fea`                                               |
| SNF 标签 / 指标      | `repro_out/<CANCER>.snf.tsv`、`<CANCER>.snf.metrics.json`；九癌种汇总 `repro_out/snf_table1_compare.json` |


---

## 五、后续可做事项（可选）

- 补齐 9 癌种 `train.py` 输出后运行 `repro/batch_evaluate.py`，与论文 Table 1 逐项对照。
- 如需更接近论文 log-rank，换用论文补充材料中的 **PANCAN 随访表** 并走 `data_io` 的 **legacy** 分支（列名需齐全）。
- 性能优化：GPU、batch、`num_workers` 等（训练脚本层面，非本会话修改重点）。

---

## 六、R baseline 方法选择（SNF）与安装踩坑

- **为何选 SNF**：在 NEMO（GitHub + `remotes`）、PINS（`PINSPlus` 依赖多）、SNF（**纯 CRAN `SNFtool`**）三者中，SNF **最易稳定安装**，且与论文 Table 1 中 **similarity-based** 对照一致。
- **权威脚本**：`benchmarks/snf_template.R`（读 `vendor/Subtype-DCC/subtype_file/fea/<CANCER>/*.fea` → SNF → 写 `repro_out/<CANCER>.snf.tsv`）；**仅安装依赖**：`benchmarks/install_snf_only.R`。
- `**benchmarks/r/install_deps.R`**：曾尝试一键安装多包 + GitHub，易因 **CRAN 下载超时** 失败（后台任务 exit 1）；已改为 **废弃占位**，日常请只用 `install_snf_only.R`。
- `**nemo_template.R` / `pins_template.R`**：保留为说明入口，需要时再按 `snf_template.R` 的 **TSV 输出格式** 自行接 API。

### 批量 SNF + 评测（九癌种，仓库根 PowerShell 示例）

```powershell
cd "D:\LateInUniversity\Term10\健康管理中的智能数据分析\subtype-dcc-repro"
$R = "C:\Program Files\R\R-4.3.2\bin\Rscript.exe"
$py = ".\.venv\Scripts\python.exe"
$cancers = "BRCA","BLCA","KIRC","LUAD","PAAD","SKCM","STAD","UCEC","UVM"
foreach ($c in $cancers) {
  & $R "benchmarks\snf_template.R" $c "vendor\Subtype-DCC\subtype_file\fea\$c" "repro_out\$c.snf.tsv"
  & $py repro\evaluate.py --labels "repro_out\$c.snf.tsv" `
    --clinical data\clinical\cBioportal_data.tsv -c $c --out-json "repro_out\$c.snf.metrics.json"
}
& $py repro\compare_snf_table1.py   # 生成 repro_out/snf_table1_compare.json
```

---

## 七、论文 Table 1「SNF」行 vs 本仓库复现（`cBioportal_data.tsv` + 当前 SNF 默认超参）

汇总文件：`**repro_out/snf_table1_compare.json**`（由 `repro/compare_snf_table1.py` 生成）。论文数值来自 Subtype-DCC 原文 Table 1 中 **SNF** 行（格式：−log₁₀ *p* / 显著临床项个数）。


| 癌种   | 论文 SNF   | 本复现 SNF      | 备注                                                           |
| ---- | -------- | ------------ | ------------------------------------------------------------ |
| BRCA | 0.93 / 5 | 2.25 / 6     | 生存 *p* 与临床项数均有偏差；随访表与簇划分与原文不完全一致时可预期。                        |
| BLCA | 1.31 / 6 | 0.64 / 5     | 同上。                                                          |
| KIRC | 8.19 / 6 | **9.70 / 6** | 临床项数 **一致**；−log₁₀ *p* 同量级略高，可接受。                            |
| LUAD | 2.23 / 4 | 2.30 / 3     | −log₁₀ *p* 接近；临床少 1 项显著。                                     |
| PAAD | 3.24 / 3 | 2.78 / 0     | 小样本 + 表型列易不显著；与论文差距较大时优先核对临床列缺失与簇标签。                         |
| SKCM | 5.27 / 4 | **5.63 / 4** | 临床项数 **一致**；−log₁₀ *p* 接近。评测 **n=420**（标签 446 与临床内交后减少，属正常）。 |
| STAD | 0.72 / 2 | **0.66 / 2** | **两项指标均接近**，可作为「流程正确」的强对照之一。                                 |
| UCEC | 5 / 1    | 7.18 / 1     | 显著临床 **项数一致**；log-rank 更强，与生存列定义/簇划分有关。                      |
| UVM  | 2.77 / 0 | 4.30 / 0     | 显著临床 **项数一致**（均为 0）；小样本 *p* 波动大。                             |


**结论（是否「复现无误」）**：

- **工程层面**：九癌种均已生成 `**repro_out/<CANCER>.snf.tsv`** 与 `***.snf.metrics.json`**，`evaluate.py` 与 `**data_io` 临床合并逻辑** 工作正常；**STAD、KIRC、SKCM、UCEC、UVM** 等在 **至少一项指标上与论文接近或一致**，说明 **SNF 管线 + 评测脚本与论文设定对齐**，无静默失败。
- **数值逐格与 Table 1 完全一致**：**不应作为硬性标准**。论文写明对结果波动的方法 **重复 5 次取均值**；本复现为 **单次 SNF**；且随访使用 **GerkeLab/cBioPortal 合并 OS**，与 Subtype-DCC 原文常用的 **legacy PANCAN** 路径未必相同——**log-rank 与 χ² 的绝对数值会系统性偏移**，此前 Subtype-DCC KIRC 亦有类似讨论（见第三节）。

若课程或审稿要求 **更贴近 Table 1 数字**：在固定 **SNF 超参（K、alpha、T）与随机种子** 的前提下，换 **与论文一致的随访表** 并重跑 `evaluate.py`。

---

## 八、训练与资源策略（Subtype-DCC + 与 baseline 并行）

- **Subtype-DCC**：`train.py` 默认 **早停**（plateau）、`epochs` 上限 600；本机若 **KIRC 全 epoch 约 <2h**，则九癌种 Subtype-DCC 可在 **约 1 天量级** 内规划完成（大癌种如 BRCA 预留略多时间）。
- **与 sklearn / R 并行**：若 Subtype-DCC 占用 **GPU**，**SNF（CPU）** 与 **sklearn baseline** 可并行跑；若全部为 **CPU**，则与深度训练 **竞争算力**，建议 **错峰**（训练空档跑 SNF/sklearn）。
- **记录**：建议记录每癌种 **墙钟时间、epoch 数、评测样本数**（`evaluate.py` 的 `n_samples_eval`），便于报告附录。

---

## 九、Subtype-DCC 批量训练（当前终端任务 & 明日接续）

### 9.1 你正在跑的命令（8 癌种，不含 KIRC）

在 **仓库根** `subtype-dcc-repro` 下执行（`KIRC` 已在此前训过，本批 **刻意跳过**；**BRCA 在本批中重训**）：

```powershell
cd "D:\LateInUniversity\Term10\健康管理中的智能数据分析\subtype-dcc-repro"

.\.venv\Scripts\python.exe repro\run_subtype_dcc_batch.py `
  --cancers BRCA BLCA LUAD PAAD SKCM STAD UCEC UVM `
  --extra-args "--workers 4"
```

- `**run_subtype_dcc_batch.py` 行为**：对 `--cancers` 列表 **按顺序逐个** `subprocess.run(..., cwd=vendor/Subtype-DCC/Subtype-DCC)`；等价于在该子目录依次执行  
`python train.py -c <CANCER> --workers 4`。  
- `**--extra-args`**：按空格拆成 argv 片段追加到 `train.py` 后；若要改 epoch 上限可写 `--extra-args "--workers 4 --epochs 600"`（默认配置里 `epochs` 已是 600，一般不必重复）。  
- `**--workers 4`**：传给 `train.py` / DataLoader；可按机器改为 `2` 或 `8`。

### 9.1b 续跑未完成癌种（示例：仅剩 STAD、UCEC、UVM）

与 §9.1 相同：**先在仓库根 `cd`，再用 `.\.venv\Scripts\python.exe`**（勿用系统 `python`，勿在 `Subtype-DCC\Subtype-DCC` 子目录里用错的相对路径去找 `.venv`，见 §二、§9.4）。

```powershell
cd "D:\LateInUniversity\Term10\健康管理中的智能数据分析\subtype-dcc-repro"

.\.venv\Scripts\python.exe repro\run_subtype_dcc_batch.py `
  --cancers STAD UCEC UVM `
  --extra-args "--workers 4"
```

### 9.2 训练产物存哪儿（每训完一个癌种就应出现）


| 类型                    | 路径（相对仓库根）                                                               |
| --------------------- | ----------------------------------------------------------------------- |
| 簇标签（评测入口）             | `vendor/Subtype-DCC/Subtype-DCC/results/<CANCER>.dcc`                   |
| 嵌入（t-SNE 等）           | `vendor/Subtype-DCC/Subtype-DCC/results/<CANCER>.fea`                   |
| Loss 曲线图              | `vendor/Subtype-DCC/Subtype-DCC/results/<CANCER>Train_loss.png`         |
| Checkpoint（推理用）       | `vendor/Subtype-DCC/Subtype-DCC/save/<CANCER>/checkpoint_last.tar`      |
| 每 epoch 存档（体积极大，注意磁盘） | `vendor/Subtype-DCC/Subtype-DCC/save/<CANCER>/checkpoint_*.tar`（上游默认行为） |
| 全局配置 / 随机种子           | `vendor/Subtype-DCC/Subtype-DCC/config/config.yaml`（如 `seed: 21`）       |


**预设簇数 M**（`train.py` 中 `cancer_dict`，勿与 SNF 邻居参数混淆）：BRCA=5，BLCA=5，KIRC=4，LUAD=3，PAAD=2，SKCM=4，STAD=3，UCEC=4，UVM=4。

### 9.3 明日开工检查清单

1. **看终端是否报错退出**：任一中途 `check=True` 失败会终止整批后续癌种；需从失败癌种重跑或改 `--cancers` 只列未完成的。
2. **逐癌种确认文件存在**：至少应有 `**results/<CANCER>.dcc`** 与 `**save/<CANCER>/checkpoint_last.tar`**。
3. **（推荐）全量 Subtype-DCC 评测**：本批 8 个 + 已有 **KIRC** 一并汇总：

```powershell
.\.venv\Scripts\python.exe repro\batch_evaluate.py
```

默认读取 `**vendor/Subtype-DCC/Subtype-DCC/results/*.dcc**` 与 `**data/clinical/cBioportal_data.tsv**`，写出 `**repro_out/batch_metrics_summary.json**`；单癌种也可用：

```powershell
.\.venv\Scripts\python.exe repro\evaluate.py `
  --labels vendor\Subtype-DCC\Subtype-DCC\results\LUAD.dcc `
  --clinical data\clinical\cBioportal_data.tsv -c LUAD `
  --out-json repro_out\LUAD.metrics.json
```

1. **（可选）记录训练用时**：在仓库根对单癌种包一层 `Measure-Command { ... }`，或手写表格：**癌种、起止时间、总秒数、最终 epoch**（见终端 log / `checkpoint_*` 编号）。

### 9.4 其他必要信息（防踩坑）

- **Python 路径**：务必在 **仓库根** 用 `.\.venv\Scripts\python.exe`；勿在 `vendor\Subtype-DCC\Subtype-DCC` 下用错误的 `..\..\` 去找 venv（见第三节）。  
- **数据**：`train.py` 读 `../subtype_file/fea/<CANCER>/` 四张 `*.fea`；本仓库通过 `**vendor/Subtype-DCC/subtype_file/fea`**（多为 junction 到 Subtype-GAN）提供。缺文件会直接报错。  
- **与论文 Table 1 数值**：即使训练完美，`**cBioportal_data.tsv`** 与论文 **PANCAN/legacy** 随访不一致时，**−log₁₀ log-rank *p*** 仍可能与原文表不完全一致；**临床显著项个数**通常更稳。  
- **磁盘**：若 `save/<CANCER>/` 下按 epoch 存了大量 `checkpoint_N.tar`，且空间紧张，可在确认 `**checkpoint_last.tar` + results 可用** 后 **自行清理中间 epoch 文件**（仅操作建议，非脚本默认）。  
- **接续跑**：若本批只完成部分癌种，下次可只传未完成的，例如：  
`--cancers STAD UCEC UVM --extra-args "--workers 4"`。

---

## 十、Subtype-DCC 超参与当前训练策略（如何选）

### 10.1 超参一览（按「官方默认 / 论文设定」归类）


| 类别                    | 名称                                                             | 典型取值 / 来源                                                         | 说明                                                           |
| --------------------- | -------------------------------------------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------ |
| **簇数**                | 每癌种簇数 **M**                                                    | `train.py` 内 `cancer_dict`（如 KIRC=4，BRCA=5）                       | 与论文 Table 1 各数据集预设子型数一致；可用 `--cluster_number` 覆盖。            |
| **优化**                | `learning_rate`                                                | `3e-4`（`config/config.yaml`）                                      | 论文 *Experimental settings* 与仓库默认一致。                          |
|                       | `weight_decay`                                                 | `0`                                                               | 默认无 L2。                                                      |
|                       | `batch_size`                                                   | `64`（`train.py` 默认参数）                                             | 论文与官方 README 一致；影响 DCL 与 batch 内对比。                          |
| **对比学习温度**            | 实例级（DCL）                                                       | 代码中 `**DCL(temperature=0.5)` 写死**（`train.py` `train_one_epoch`）   | 与配置里 `instance_temperature: 0.5` **数值一致**；改超参时需同时改代码或接受二者分离。 |
|                       | 簇级                                                             | `cluster_temperature`（默认 **1.0**，yaml）                            | 传入 `ClusterLoss`。                                            |
| **表示维度**              | `feature_dim`                                                  | 默认 **128**（yaml）                                                  | 决定嵌入维度及后续 ICH/CCH 头规模（见 `modules/network.py`）。               |
| **训练轮数**              | `epochs`                                                       | 默认 **600**（yaml）                                                  | 上限；配合早停通常不必跑满。                                               |
|                       | `start_epoch`                                                  | 0                                                                 | 一般不改。                                                        |
| **随机性**               | `seed`                                                         | 默认 **21**（yaml）                                                   | `torch` / `numpy` 种子；换 seed 会改变簇划分与指标。                       |
| **数据加载**              | `workers`                                                      | yaml 默认 **8**；你当前批次用 `**--workers 4`**                            | 仅影响 DataLoader 并行，**不改变**损失定义；按 CPU 核数与内存调。                  |
| **增广**                | 训练视图                                                           | 对输入加 **高斯噪声** `N(0,1)` 得两视图 `x_i, x_j`                            | 论文写明 Gaussian noise 增广；与 Subtype-DCC 仓库实现一致。                 |
| **网络宽度**              | MLP 瓶颈                                                         | `ae.py` 默认 `**[5000,2000,1000,256]`** 逐层（输入维由 `get_input_dim` 推断） | 论文 *Experimental settings* 与官方结构一致；**消融**时由补丁按输入维自动收窄。       |
| **早停（本仓库 vendor 补丁）** | `min_epochs`                                                   | 默认 **100**                                                        | 前 100 epoch 不触发 plateau 停。                                   |
|                       | `plateau_window` / `plateau_rel_tol` / `plateau_stable_epochs` | 30 / 0.005 / 40                                                   | 损失窗口相对起伏足够小时计为「平」，连续若干 epoch 则停。                             |
|                       | `improvement_patience`                                         | 默认 **0**（关闭）                                                      | 若设为 >0，则按 `improvement_min_delta` 判无提升停训。                    |
|                       | `--no_early_stop`                                              | 默认 **不**加                                                         | 加上则强制跑满 `--epochs`。                                          |
| **日志**                | `--log_every`                                                  | 默认 **10**                                                         | 控制 loss 打印间隔。                                                |
| **消融**                | `--exclude_omics`                                              | 如 `rna` 或 `CN,meth,...`                                           | 非论文主表默认；用于 Figure S1 类实验。                                    |


配置文件路径：`**vendor/Subtype-DCC/Subtype-DCC/config/config.yaml*`*。其中键均通过 `yaml_config_hook` 暴露为 `**train.py` 的命令行 `--键 值`**，例如 `--learning_rate 1e-4`（一般不必改）。

### 10.2 当前训练策略里「超参是怎么选的」

- **总原则**：复现以 **官方 Subtype-DCC 默认 + 论文 Experimental settings** 为准，**不做网格搜索**；与原文 Table 1 对比时，优先保证 **流程、簇数 M、增广类型、学习率/批大小/温度与论文叙述一致**，而非在本地重新调参刷表。  
- **你正在跑的 8 癌种 batch**：在默认 yaml 与 `cancer_dict` 不变的前提下，仅通过 `--extra-args "--workers 4"` **降低 DataLoader 线程数**，属于 **工程上的 I/O/CPU 调优**，**不是**改变算法超参。  
- **早停**：保持 **开启**（默认），在 **不少于 `min_epochs`** 的前提下，损失平台期后提前结束，节省墙钟时间；若需与「固定 600 epoch」的第三方结果逐 epoch 对齐，再加 `--no_early_stop`。  
- **若将来要系统性调参**：论文写明六类主超参（实例/簇温度、batch、表示维、学习率、epoch）曾按顺序调优（见原文 **Table S1** 叙述）；本地可固定其它五项、每次改一项，并 **固定 seed** 记录 `results/*.dcc` 与 `repro_out/*.metrics.json`。

---

## 十一、KIRC 论文风格 Fig 4–7（先算缓存，再单独出图）

避免每次改图都重跑 **t-SNE + Enrichr**：先写磁盘缓存，再按图调用绘图脚本（**不依赖 pyarrow**，GO/KEGG 表为 `*.csv`）。


| 角色                                              | 路径                                                                                                          |
| ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| 共享逻辑（加载、富集、t-SNE、各 figure 绘制函数）                 | `repro/kirc_manuscript_lib.py`                                                                              |
| **一次性计算**（t-SNE、标志基因、GO BP/CC/MF、KEGG、Fig7 矩阵等） | `repro/kirc_manuscript_compute_cache.py`                                                                    |
| 缓存读写                                            | `repro/kirc_manuscript_cache_io.py`                                                                         |
| 默认缓存目录                                          | `repro_out/kirc_manuscript_cache/`（`meta.json`、`Z.npy`、`go_fig5.csv`、`kegg_fig6.csv` 等）                     |
| 单图脚本                                            | `repro/plot_kirc_fig4a.py`、`plot_kirc_fig4b.py`、`plot_kirc_fig5.py`、`plot_kirc_fig6.py`、`plot_kirc_fig7.py` |
| 编排入口（可选：只画其中几张）                                 | `repro/kirc_manuscript_figures.py`                                                                          |


**推荐命令（仓库根）：**

```powershell
# 1) 只算一次（需联网 Enrichr；约 1–2 分钟量级，视网络而定）
.\.venv\Scripts\python.exe repro\kirc_manuscript_compute_cache.py --cache-dir repro_out\kirc_manuscript_cache

# 2) 只重画 Figure 5（读缓存，不再请求 Enrichr）
.\.venv\Scripts\python.exe repro\plot_kirc_fig5.py --cache-dir repro_out\kirc_manuscript_cache --out-dir repro_out\kirc_manuscript_figs

# 等价：编排脚本只选 fig5
.\.venv\Scripts\python.exe repro\kirc_manuscript_figures.py --from-cache --figures fig5 --cache-dir repro_out\kirc_manuscript_cache --out-dir repro_out\kirc_manuscript_figs

# 无 --from-cache：会先跑完整 compute 再画 --figures 所列全部图（慢）
.\.venv\Scripts\python.exe repro\kirc_manuscript_figures.py --figures fig4a,fig4b,fig5,fig6,fig7
```

---

*文档由会话压缩生成；若路径与本机不一致，请以当前仓库根为准替换命令中的盘符。*