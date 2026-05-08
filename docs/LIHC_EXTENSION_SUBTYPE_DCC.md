# LIHC：Subtype-DCC 扩展

主文顺序：**§0 结论** → **§1–3**（背景、数据、预处理）→ **§4 超参** → **§5 训练与评测结果**（图、指标、汇报与替代方案）。其余命令、排障、代码修复等一律见 **附录**。

泛化「新癌种」仍见 `docs/CUSTOM_DATASET_SUBTYPE_DCC.md`、`docs/TCGA_SUBTYPE_DCC_GENERIC.md`；预处理细节全文见 `docs/XENA_TO_FEA_PREPROCESSING.md`。**给队友的推荐阅读顺序**见 [`docs/EXTENSION_FOR_NEW_CANCER.md`](EXTENSION_FOR_NEW_CANCER.md)。

---

## 0. 核心结论

**在本仓库当前设定下，Subtype-DCC 应用于 TCGA LIHC 的预后验证为阴性（OS 多组 log-rank 不显著）。**


| 结论项            | 内容                                                                                                                                                                                         |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **定量**         | `repro/evaluate.py` 归档：**仓库快照** [`docs/lihc_artifacts/LIHC.metrics.json`](../lihc_artifacts/LIHC.metrics.json)（与 **`repro_out/LIHC.metrics.json`** 本机重跑输出一致）— **n_samples_eval = 357**，**neg_log10_logrank_p ≈ 0.127**（*p* ≈ 0.75，不显著），**n_significant_clinical = 2**（年龄、性别与簇相关；分期/T/N 不显著，M 边缘）。 |
| **可汇报**        | 完成 **Xena→fea→训练→与论文 Table 1 同口径评测**；**不在本队列上主张 OS 优越性**；并列 **临床协变量（年龄、性别）** 与 **manuscript 级组图**（§10）说明 **组学/表型维度的可解释性**。                                                                 |
| **若需「阳性预后」叙事** | 更稳妥的是 **换癌种/换队列**（如论文其他癌种中 OS 信号较强的设定），在 **同一评测脚本** 上重跑；或 **外部队列验证**。**不推荐** 用同批数据生存去筛基因再聚类（泄漏）。                                                                                           |


---

## 1. 背景与目标

- **目标**：在 **LIHC** 上训练 Subtype-DCC，得到 `vendor/Subtype-DCC/Subtype-DCC/results/LIHC.dcc`、`LIHC.fea`，并与本仓库临床表做 `evaluate.py` 评测。
- **与论文差异**：`train.py` 的 `cancer_dict` **不含 LIHC**，必须显式 `**--cluster_number`**；临床为 `data/clinical/cBioportal_data.tsv`（与原文 PANCAN 表可能不一致，报告中需声明）。

---

## 2. 数据来源

- **组学**：UCSC Xena 四个 `.gz`（RNA、miRNA、450K 甲基化、GISTIC2 CN），置于 `vendor/Subtype-DCC/subtype_file/fea/LIHC/`（或 junction）。
- **临床**：`data/clinical/cBioportal_data.tsv`；样本 ID 经 `patient_barcode`（TCGA **前 12 位**）与标签合并。可预检：`repro/check_omics_clinical_overlap.py -c LIHC`。

---

## 3. 预处理与 `*.fea`

- **脚本（仓库根）**：`repro/xena_gz_to_subtype_dcc_fea.py` → 产出 `rna.fea`、`miRNA.fea`、`meth.fea`、`CN.fea`（逗号 CSV，**行=特征、列=样本**）。
- **样本**：四表 **列名交集**，列顺序为 **字典序排序** 后的统一列表。
- **缺失**：行内 **中位数**；`±inf`→NaN 后同样处理。
- **甲基化**：流式读入，按行 **nanvar（ddof=0）** 保留 top **K** 探针（默认 **3139**，`--meth-max-probes`）。
- **RNA / miRNA / CN**：填补后按行方差截断至默认 **3500 / 450 / 3500**（`--rna-max-genes`、`--mirna-max`、`--cn-max-genes`）。
- **本仓库 LIHC 样本列数**：**357**（以 `fea/LIHC` 及 `xena_to_fea_meta.txt` 为准）。

**训练时特征块顺序与归一化**（与 Subtype-DCC 一致）：拼接顺序 **CN → meth → miRNA → RNA**；对拼接矩阵 `**MinMaxScaler` 按列** 一次缩放（`dataloader.load_scaled_features`）。详见附录 **A.3**。

---

## 4. 超参与训练策略

**工作目录**：必须在 `vendor/Subtype-DCC/Subtype-DCC` 再执行 `train.py`，否则找不到 `config/config.yaml`（详见附录 **A.1**）。

**推荐一组（与论文常用 batch 量级一致、便于监控）**：


| 项        | 建议                            | 参数                                                              |
| -------- | ----------------------------- | --------------------------------------------------------------- |
| 簇数 M     | 3 或 4（报告写清）                   | `--cluster_number`                                              |
| 学习率      | 1e-4；不稳再试 5e-5                | `--learning_rate`                                               |
| batch    | **32** 或 64（避免与 16+极小 lr 双保守） | `--batch_size`                                                  |
| seed     | 42 等                          | `--seed`                                                        |
| workers  | 视机器                           | `--workers`                                                     |
| 簇 / 生存日志 | 每 10 epoch                    | `--cluster_log_every 10`、`--survival_logrank_every 10`（后者 0 可关） |


**注意**：`epoch` 打印的 loss 是 **当 epoch 内各 batch loss 之和**，**勿与不同 batch_size 的数值直接比绝对大小**（附录 **A.2**）。

**最小训练命令**（路径按你本机改）：

```powershell
cd <repo>\vendor\Subtype-DCC\Subtype-DCC
Remove-Item -Recurse -Force .\save\LIHC\ -ErrorAction SilentlyContinue
..\..\..\.venv\Scripts\python.exe train.py -c LIHC --cluster_number 4 `
  --learning_rate 0.0001 --seed 42 --batch_size 32 --workers 4 `
  --cluster_log_every 10 --survival_logrank_every 10
```

重训、批量脚本、早停等：**附录 A.1**。

---

## 5. 训练与评测结果（图 + 指标 + 汇报 + 替代）

### 5.1 Loss 曲线

训练结束后由 `train.py` 写出（文件名随 `cancer_type`）：

`vendor/Subtype-DCC/Subtype-DCC/results/LIHCTrain_loss.png`

在 Markdown 预览或汇报 PPT 中 **直接插入该 PNG** 即可。从本仓库文档预览时，可使用已提交的副本（与 `results/LIHCTrain_loss.png` 一致）：

![LIHC train loss](../lihc_artifacts/LIHCTrain_loss.png)

*若尚未生成，请确认已完成一次完整训练且 `results/` 目录可写；生成后可复制一份到 `docs/lihc_artifacts/` 以便随仓库展示。*

---

### 5.2 论文口径评测（`LIHC.metrics.json`）

**仓库内归档（便于在 GitHub 上直接打开）**：[`docs/lihc_artifacts/LIHC.metrics.json`](../lihc_artifacts/LIHC.metrics.json)。**本机重跑** 时 `--out-json` 默认写到 `repro_out/LIHC.metrics.json`（`repro_out/` 被 `.gitignore` 忽略，不入版控；数值应与 `docs/lihc_artifacts/` 快照一致，除非你更换了标签或临床表）。

**命令（仓库根）**：

```powershell
.\.venv\Scripts\python.exe repro\evaluate.py `
  --labels vendor\Subtype-DCC\Subtype-DCC\results\LIHC.dcc `
  --clinical data\clinical\cBioportal_data.tsv -c LIHC `
  --out-json repro_out\LIHC.metrics.json
```

**归档数值（与 JSON 一致；若你重训后请更新下表）**：


| 字段                                | 值                                |
| --------------------------------- | -------------------------------- |
| n_samples_eval                    | **357**                          |
| neg_log10_logrank_p               | **≈ 0.127**（OS **不显著**，*p*≈0.75） |
| n_significant_clinical（α=0.05，六项） | **2**                            |


**六项临床 *p*（摘要）**：**年龄**、**性别** 与簇显著；**分期/T/N** 不显著；**M** 约 0.051（边缘，未计入「显著个数」）。

---

### 5.3 Manuscript 风格组图（`docs/lihc_artifacts/` 与 `repro_out/lihc_manuscript_figs/`）

由 `repro/kirc_manuscript_compute_cache.py --cancer LIHC` + `repro/plot_kirc_fig4a.py` … `plot_kirc_fig7.py` 生成。**仓库内 PNG 副本**（便于 GitHub 预览）：[`docs/lihc_artifacts/`](../lihc_artifacts/) 下表所列 `LIHC_*.png`。**本机重跑** 时默认在 `repro_out/lihc_manuscript_figs/`（可用 `--out-dir` 改路径）；汇报建议使用 **`LIHC_*.png`** 命名（与 `KIRC_fig*.png` 同内容时可复制重命名）。


| 文件                                    | 内容                                                |
| ------------------------------------- | ------------------------------------------------- |
| Fig4A                                 | 嵌入 t-SNE（`docs/lihc_artifacts/LIHC_fig4A_tsne.png`；本机同名为 `repro_out/lihc_manuscript_figs/LIHC_fig4A_tsne.png`）                   |
| Fig4B（标题的癌种亚型是按照KIRC的癌种亚型命名的，有误，后续可改） | 标志基因网格（`LIHC_fig4B_marker_heatmap.png`，同上两处目录）           |
| Fig5                                  | GO 条形（`LIHC_fig5_GO_bars.png`）                    |
| Fig6                                  | KEGG 点图（`LIHC_fig6_KEGG_dot.png`）                 |
| Fig7                                  | miRNA–通路热图（`LIHC_fig7_miRNA_pathway_heatmap.png`） |


**重跑缓存/出图**（需 `gseapy` 与网络访问 Enrichr）：附录 **A.4**。

---

### 5.4 结果解读

- **LIHC 与方差筛的局限（预后阴性的原因之一）**：**LIHC 肿瘤异质性强**，组学上主导的主轴常来自 **免疫/基质、病因学、纯度或批次相关变异** 等，**未必沿「总生存」方向排序**。本管线对 RNA/miRNA/CN/甲基化探针均采用 **跨样本方差（无监督）** 保留高变特征，**会优先留下「跨样本波动大」的维度，而不是「与 OS 最相关」的维度**；因此聚类更易反映 **组学主导轴**，与 **随访 OS 弱对齐甚至脱钩** 是合理现象，可作为 **预后阴性的一种解释**（与实现错误无关，全链路见附录 **A.5**）。  
- **OS**：无监督目标 **未优化 log-rank**；上节所述 **方差筛 ≠ 预后特征**；加之 *n* 与事件量限制检验力 → **阴性可理解**。  
- **非 OS**：年龄、性别与簇相关 → 汇报中可做 **「表型维补充」**；**KM + 本组 t-SNE/通路图** 形成 **「生存一节 + 机制/表型一节」**，避免只有一张 KM。  
- **实现**：训练内 `neg_log10_p` 与 `evaluate.py` **同一函数**。

---

### 5.5 汇报方案（可直接改编 PPT）

1. **主句**：LIHC 扩展完成；**Table 1 同口径下 OS log-rank 不显著**（给出 −log₁₀ *p*）；**不主张本队列预后优越性**。
2. **承接**：**簇与年龄、性别相关** + **嵌入/通路图（§5.3）** → **组学—表型可讨论**。
3. **局限**：单队列、随访表与论文可能不同、无监督未做多重校正等。
4. **后续**：换 **M/seed** 敏感性；**换癌种** 或外部队列；**sklearn/SNF 基线**（附录 A.4）。

---

### 5.6 替代方案（不引入生存泄漏）


| 做法                  | 入口                                                                    |
| ------------------- | --------------------------------------------------------------------- |
| KM + log-rank 标题    | `repro/plot_survival.py`（附录 A.4）                                      |
| 旧版 Fig4 风格          | `repro/plot_figure4.py -c LIHC`                                       |
| 训练中快照 checkpoint    | `repro/monitor_subtype_dcc_checkpoint.py`                             |
| 基线聚类 + 同口径 evaluate | `repro/baselines_sklearn.py`、`repro/batch_evaluate_table1_methods.py` |


---

## 附录

### A.1 训练 cwd、完整命令与批量

- **必须**：`cd vendor/Subtype-DCC/Subtype-DCC` 再运行 `train.py`，否则 `config/config.yaml` 报错。  
- **venv**：从该目录向上 **三级** 到仓库根再调用 `.venv\Scripts\python.exe`（见 `docs/SESSION_HANDOFF.md` §二、§三）。  
- **批量**：`repro/run_subtype_dcc_batch.py --cancers LIHC --extra-args "..."`（内部已设 cwd）。  
- **早停、plateau、improvement_patience**：均为 `train.py` 已有参数，见 argparse 帮助或历史长版文档备份。

### A.2 Loss 与 batch_size

- `train_one_epoch` 返回 **本 epoch 所有 batch 的 loss 之和**，非均值。  
- `batch_size` 变小 → 每 epoch batch 数变多 → **打印的 epoch loss 总和通常变大**，不宜与另一 batch 设置直接比绝对值。

### A.3 特征拼接顺序与默认维数


| 顺序  | 块         |
| --- | --------- |
| 1   | CN.fea    |
| 2   | meth.fea  |
| 3   | miRNA.fea |
| 4   | rna.fea   |


默认行数上限：**CN 3500、meth 3139 探针、miRNA 450、RNA 3500**（均由 `xena_gz_to_subtype_dcc_fea.py` 参数控制）。`dataloader.load_omics_blocks` 对后续块 `**reindex` 与首张 CN 列序一致**。

### A.4 评测与可视化命令汇总

```powershell
# evaluate（仓库根）
.\.venv\Scripts\python.exe repro\evaluate.py --labels vendor\Subtype-DCC\Subtype-DCC\results\LIHC.dcc `
  --clinical data\clinical\cBioportal_data.tsv -c LIHC --out-json repro_out\LIHC.metrics.json

# KM
.\.venv\Scripts\python.exe repro\plot_survival.py --labels vendor\Subtype-DCC\Subtype-DCC\results\LIHC.dcc `
  --clinical data\clinical\cBioportal_data.tsv -c LIHC --out repro_out\LIHC_survival_km.png

# manuscript 缓存 + 出图（需网络）
.\.venv\Scripts\python.exe repro\kirc_manuscript_compute_cache.py --cancer LIHC `
  --labels vendor\Subtype-DCC\Subtype-DCC\results\LIHC.dcc `
  --embedding vendor\Subtype-DCC\Subtype-DCC\results\LIHC.fea `
  --cache-dir repro_out\lihc_manuscript_cache
$cd="repro_out\lihc_manuscript_cache"; $od="repro_out\lihc_manuscript_figs"
foreach ($s in "plot_kirc_fig4a","plot_kirc_fig4b","plot_kirc_fig5","plot_kirc_fig6","plot_kirc_fig7") {
  .\.venv\Scripts\python.exe "repro\$s.py" --cache-dir $cd --out-dir $od
}
# 可选：复制为 LIHC_ 前缀便于汇报
# Copy-Item repro_out\lihc_manuscript_figs\KIRC_*.png -Destination {对应 LIHC_ 名}
```

### A.5 排障与代码修复（历史记录摘要）

- **列顺序错配**：历史 LIHC 四表 `**list(columns)` 曾不一致**；已在 `vendor/Subtype-DCC/Subtype-DCC/dataloader.py` `**reindex`**；`repro/data_io.load_concat_matrix` 已对齐同一逻辑。  
- **簇损失 NaN**：`contrastive_loss.py` 中 `**torch.xlogy`** + `**clamp_min`**。  
- `**.dcc` 行序**：`train.py` 使用 `load_omics_blocks` 返回的 `**sample_names`**（与训练矩阵行一致）。  
- **训练内 log-rank**：调用 `repro/evaluate.survival_neg_log10_p`，与 `**evaluate.py`** 一致。  
- **各模态数值范围**：RNA 多为 log 型表达（不必强行 0–10）；甲基化 β∈[0,1]；CNV 为 GISTIC2 离散档 {-2,…,2}。**全链路复查结论**：拼接顺序与 `evaluate` 口径正确；阴性更宜归为生物学/统计设定，详见 **§10.4** 与 `docs/XENA_TO_FEA_PREPROCESSING.md`。

### A.6 文档索引


| 主题          | 文档                                   |
| ----------- | ------------------------------------ |
| Xena→fea 细节 | `docs/XENA_TO_FEA_PREPROCESSING.md`  |
| 新癌种通用       | `docs/CUSTOM_DATASET_SUBTYPE_DCC.md` |
| TCGA 未定癌种   | `docs/TCGA_SUBTYPE_DCC_GENERIC.md`   |
| 终端与 cwd     | `docs/SESSION_HANDOFF.md`            |


---

