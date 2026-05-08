# 给队友：在新癌种上扩展 Subtype-DCC

LIHC 已由仓库维护者完成 **Xena→`*.fea`→训练→`evaluate`→基础可视化**；你要做**其他癌种**时，把 LIHC 当作「一条完整样例流水线」对照即可。

## 建议阅读顺序（由具体到泛化）

1. **[`LIHC_EXTENSION_SUBTYPE_DCC.md`](LIHC_EXTENSION_SUBTYPE_DCC.md)**（**首选**）  
   结论、数据与 `*.fea`、训练 cwd 与超参、`evaluate` 口径、组图与汇报话术、附录里的命令与排障。扩展新癌种时：**换癌种缩写、路径与 `train.py --cluster_number`**，结构沿用本文。

2. **[`XENA_TO_FEA_PREPROCESSING.md`](XENA_TO_FEA_PREPROCESSING.md)**  
   `repro/xena_gz_to_subtype_dcc_fea.py` 的细节（缺失填补、方差截断、甲基化 top‑K 等）。

3. **[`CUSTOM_DATASET_SUBTYPE_DCC.md`](CUSTOM_DATASET_SUBTYPE_DCC.md)**、**[`TCGA_SUBTYPE_DCC_GENERIC.md`](TCGA_SUBTYPE_DCC_GENERIC.md)**  
   非论文九癌种或自定义队列时的注意点。

4. **[`SESSION_HANDOFF.md`](SESSION_HANDOFF.md)**  
   Windows / venv 路径、`train.py` 工作目录、常见报错与评测入口速查。

5. **[`README.md`](../README.md)**（仓库根）  
   九癌种 Table 1 复现、`batch_evaluate_table1_methods.py`、sklearn 基线与作图命令。

## 代码与命令示例在哪里

| 目的 | 优先看哪里 |
|------|------------|
| **一条癌种从数据到图（LIHC 模板）** | **[`LIHC_EXTENSION_SUBTYPE_DCC.md`](LIHC_EXTENSION_SUBTYPE_DCC.md)** 正文 §2–§5 与 **附录 A**（训练 cwd、完整/评测/出图命令汇总）。 |
| **Xena `.gz` → `*.fea`** | 脚本 **[`repro/xena_gz_to_subtype_dcc_fea.py`](../repro/xena_gz_to_subtype_dcc_fea.py)**（`python repro/xena_gz_to_subtype_dcc_fea.py -h`）；说明见 [`XENA_TO_FEA_PREPROCESSING.md`](XENA_TO_FEA_PREPROCESSING.md)。 |
| **Subtype-DCC 训练** | **[`vendor/Subtype-DCC/Subtype-DCC/train.py`](../vendor/Subtype-DCC/Subtype-DCC/train.py)**；LIHC 最小命令见 LIHC 文档 §4；九癌种批量 **[`repro/run_subtype_dcc_batch.py`](../repro/run_subtype_dcc_batch.py)**。 |
| **单样本 / 批量评测（log-rank、临床六项）** | **[`repro/evaluate.py`](../repro/evaluate.py)**；九癌种×五方法 **[`repro/batch_evaluate_table1_methods.py`](../repro/batch_evaluate_table1_methods.py)**；与论文表对照 **[`repro/table1_compare_to_paper.py`](../repro/table1_compare_to_paper.py)**；入口示例见根目录 **README §4**。 |
| **K-means / Spectral / NMF 基线** | **[`repro/baselines_sklearn.py`](../repro/baselines_sklearn.py)**；九癌种批跑 **[`scripts/run_sklearn_baselines_batch.ps1`](../scripts/run_sklearn_baselines_batch.ps1)**（README §5）。 |
| **生存 / t-SNE / bench 条形图** | **[`repro/plot_survival.py`](../repro/plot_survival.py)**、**[`repro/plot_tsne.py`](../repro/plot_tsne.py)**、**[`repro/run_benchmark_eval.py`](../repro/run_benchmark_eval.py)** + **`plot_benchmark_bars.py`**；示例见 **README §6**。 |
| **Manuscript 风格组图（KIRC 脚本名、LIHC 数据）** | **[`repro/kirc_manuscript_compute_cache.py`](../repro/kirc_manuscript_compute_cache.py)** + **`repro/plot_kirc_fig4a.py` … `plot_kirc_fig7.py`**；命令见 LIHC 文档 **附录 A.4**。 |
| **Table 1 风格箱线图** | **[`repro/plot_figure1_table1_boxplots.py`](../repro/plot_figure1_table1_boxplots.py)**（README §4 小节）。 |
| **数据：Subtype-GAN + 临床 + junction** | **[`scripts/run_data_setup.ps1`](../scripts/run_data_setup.ps1)**、**[`scripts/link_fea_junction.ps1`](../scripts/link_fea_junction.ps1)**、**[`scripts/download_clinical.py`](../scripts/download_clinical.py)**；说明见 **README §2**。 |

Python 脚本均可 **`python repro/<脚本>.py -h`** 查看参数；训练类命令务必在文档标明的 **工作目录**（多为 `vendor/Subtype-DCC/Subtype-DCC`）下执行。

## 仓库里与 LIHC 相关的快照（无需先重跑即可浏览）

- **[`docs/lihc_artifacts/`](../lihc_artifacts/)**：`LIHC.metrics.json`、loss 与 manuscript 风格 PNG。  
- **`repro_out/`**（若已随仓库提交）：本机全量输出目录，含 LIHC 与其他实验的中间结果；体积较大时克隆会久一些。

## 你扩展新癌种时至少要自备 / 自跑的部分

- **`vendor/Subtype-DCC/subtype_file/fea/<CANCER>/` 下四个 `*.fea`**（或按文档从 Xena `.gz` 生成；**`.fea` / `.gz` 通常仍不入 Git**）。  
- **`vendor/Subtype-DCC/Subtype-DCC/results/<CANCER>.dcc`**（训练产出；`results/` 仍被 ignore 时需本地训练或向维护者索取）。  
- **临床表**：本仓库已提供 **`data/clinical/cBioportal_data.tsv`** 时，可直接用于 `repro/evaluate.py`；若与论文 PANCAN 表不一致，报告中需声明。

## 一句话任务描述（可复制给队友）

> Clone 后先看 **`docs/EXTENSION_FOR_NEW_CANCER.md`**（含「看哪些文档」+「**代码与命令示例在哪里**」表格）。扩展新癌种以 **`docs/LIHC_EXTENSION_SUBTYPE_DCC.md`** 为主模板（正文 §2–§5 + **附录 A** 命令）；预处理看 **`docs/XENA_TO_FEA_PREPROCESSING.md`**；排障看 **`docs/SESSION_HANDOFF.md`**；九癌种 Table 1 / 通用命令见根目录 **`README.md`**。临床在 **`data/clinical/`**，输出参考 **`repro_out/`**；环境按 **README §1** 本地建 **`.venv`**。
