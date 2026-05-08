# TCGA 通用说明（尚未选定具体癌种时）

本文约定：**数据源仍为 TCGA**，但**不必事先锁定某一个癌种**。队友确定癌种后，按下列「同一套流程」换 `TUMOR_TYPE` 即可。

---

## 1. 为什么仍推荐 TCGA

- Subtype-DCC 论文 Table 1 的设定与 **TCGA 多组学 + 临床** 一致；本仓库 `fea/` 与 `train.py` 的 **`cancer_dict` 簇数 M** 已按论文常见癌种写好若干 **TCGA 缩写**（如 BRCA、KIRC…）。
- 矩阵级数据（表达 / miRNA / 甲基化 / CNV）便于在国内通过 **NGDC 镜像 / 学术网络** 或 **UCSC Xena** 等渠道获取，**无需自测测序**即可跑通方法。

---

## 2. 未定癌种时，建议先决定的 3 件事（与「癌种名」一一对应）

| 决策项 | 说明 |
|--------|------|
| **TCGA 项目缩写 `TUMOR_TYPE`** | 如 `BRCA`、`LIHC`、`PRAD`… 即 TCGA 的 study 缩写；确定后所有路径、`-c` 参数都用它。 |
| **簇数 M** | 若 `TUMOR_TYPE` 已在 `vendor/Subtype-DCC/Subtype-DCC/train.py` 的 `cancer_dict` 中，可直接用默认 M；**若不在字典中**，必须传 **`--cluster_number`**（文献亚型数、或临床先验、或探索性网格——需在报告里写清）。 |
| **样本集合** | 常用 **Primary Tumor** 的 barcode 子集；是否排除某种质控失败样本，要在 `DATA_PROVENANCE` 类文档里固定下来。 |

未定癌种时，可先在 [TCGA GDC 项目列表](https://portal.gdc.cancer.gov/projects) 或 Xena 的 cohort 列表里 **按样本量、多组学完整度** 筛 2～3 个候选，再定 `TUMOR_TYPE`。

**已定癌种示例（LIHC）**：从 Xena 下载到 `*.fea`、训练 cwd、`--cluster_number`、loss/NaN 与超参调整等，见 **`docs/LIHC_EXTENSION_SUBTYPE_DCC.md`**。

---

## 3. TCGA 样本 ID 与本仓库的对齐（全癌种通用）

- **TCGA barcode**：列名常为 `TCGA-XX-XXXX-01A` 等；本仓库评测侧常用 **前 12 位** 作为患者级 ID（见 `repro/data_io.py` 的 `patient_barcode()`）。
- **四个 `*.fea` 文件**：列名需为 **同一次交集后的 aliquot / sample 列名**（与你们下载矩阵一致即可），且 **四张表列顺序一致**。
- **临床表**：`evaluate.py` 使用的临床 TSV 中，患者列应能与 `*.dcc` 里的 `sample_name` **匹配**（同一套 ID 规则）；若临床只有 patient 级、标签是 aliquot 级，需事先在预处理里约定映射（并写进文档）。

---

## 4. 多组学矩阵从哪里来（思路级，不定癌种）

以下均为科研中常用的 **「矩阵已算好」** 来源（具体下载界面随时间可能微调，以官网为准）：

| 模态 | 常见获取方式（TCGA 通用） |
|------|---------------------------|
| **RNA** | Xena TOIL TPM / 或 GDC gene expression quantification 汇总为 **gene × sample** 矩阵。 |
| **miRNA** | Xena miRNA Hub / GDC isoform or mature matrix。 |
| **甲基化** | 450K/EPIC probe 或 **gene-level β** 矩阵（需固定一套 lift 或注释版本）。 |
| **拷贝数 CN** | GISTIC2 thresholded **gene × sample** 或 log2 segment-level 再映射到 gene（与论文管线一致即可）。 |

下载后统一：**log/过滤规则固定一次**，再写成四个 `*.fea`（见 `docs/CUSTOM_DATASET_SUBTYPE_DCC.md` 与 `repro/data_io.BLOCK_ORDER`）。

若已从 Xena 下载 **四个 `.gz`（制表符矩阵）**，可在仓库根运行 **`repro/xena_gz_to_subtype_dcc_fea.py`** 自动生成 `rna.fea` / `miRNA.fea` / `meth.fea` / `CN.fea`（甲基化会 **按方差抽 top 探针**，避免 45 万行直接训练爆内存；参数可调）。**预处理与甲基化特征选择的完整表述**见 **`docs/XENA_TO_FEA_PREPROCESSING.md`**。

---

## 5. 定好 `TUMOR_TYPE` 之后的固定命令形态（与癌种无关）

在 **`vendor/Subtype-DCC/subtype_file/fea/<TUMOR_TYPE>/`** 已放好四个 `*.fea` 的前提下，于仓库根：

```powershell
# 训练（若 TUMOR_TYPE 不在 cancer_dict，必须加 --cluster_number M）
.\.venv\Scripts\python.exe vendor\Subtype-DCC\Subtype-DCC\train.py -c <TUMOR_TYPE> --workers 4

# 评测（临床表路径按你们实际文件）
.\.venv\Scripts\python.exe repro\evaluate.py `
  --labels vendor\Subtype-DCC\Subtype-DCC\results\<TUMOR_TYPE>.dcc `
  --clinical data\clinical\cBioportal_data.tsv -c <TUMOR_TYPE> `
  --out-json repro_out\<TUMOR_TYPE>.metrics.json
```

批量训练多个未定好的候选癌种时，可用：

`repro/run_subtype_dcc_batch.py --cancers BRCA LIHC PRAD ... --extra-args "--workers 4"`

---

## 6. 临床表仍用当前仓库的 `cBioportal_data.tsv` 时

- 该表为 **多癌种合并**；`-c <TUMOR_TYPE>` 会在评测逻辑里 **按癌种筛行**，只要 **barcode 与列名** 与 TCGA 习惯一致即可。
- 若队友 **只跑论文以外的 TCGA 癌种**，需确认该 TSV **是否含该癌种患者**；若不含，需 **扩展临床表** 或换用 **仅含该癌种的临床文件**（列格式仍建议兼容 `data_io`，或小幅改 `data_io` 的列名候选）。

---

## 7. 相关文档

- **通用换数据集 / 格式契约**：`docs/CUSTOM_DATASET_SUBTYPE_DCC.md`
- **仓库目录与进度**：`docs/SESSION_HANDOFF.md`（§〇、§九）
