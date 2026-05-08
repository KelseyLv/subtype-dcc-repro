# 新癌种 / 自建数据集：Subtype-DCC 使用说明（给队友）

本文说明：**在中国大陆网络环境下**，如何选取医学多组学数据、整理成本仓库可训练格式，以及仓库里还应具备哪些信息才能尽量「一键重跑」。

若 **数据源固定为 TCGA、但尚未选定具体癌种**：请直接读 `**docs/TCGA_SUBTYPE_DCC_GENERIC.md`**（样本 ID、四组学矩阵来源思路、定好缩写后的命令模板）。

**TCGA LIHC 扩展（Xena→fea、训练排障、稳定超参）**：`**docs/LIHC_EXTENSION_SUBTYPE_DCC.md**`。

---

## 1. 方法对数据的硬性要求（先对照再选数据）

Subtype-DCC 训练脚本（`vendor/Subtype-DCC/Subtype-DCC/train.py`）通过 `dataloader` 读取 `**subtype_file/fea/<癌种>/` 下四个矩阵文件**（与论文管线一致）：


| 文件          | 含义                              | 矩阵形状约定                    |
| ----------- | ------------------------------- | ------------------------- |
| `rna.fea`   | 基因表达（行=特征/基因，列=样本）              | 与其它三个 **列名（样本 ID）顺序完全一致** |
| `miRNA.fea` | miRNA 表达                        | 同上                        |
| `meth.fea`  | DNA 甲基化（如 probe 或 gene-level β） | 同上                        |
| `CN.fea`    | 拷贝数（segment/gene-level）         | 同上                        |


- 格式：**CSV**，`read_csv(..., index_col=0, sep=",")` 即 **第一列为特征名索引**，首行为样本 ID。
- **样本对齐**：四个文件中 **列名集合与顺序必须一致**；缺失某模态时不要用空文件糊弄，应使用官方支持的 `**--exclude_omics`**（见下文）。
- **簇数 M**：默认写在 `train.py` 的 `cancer_dict` 里；**新癌种**若未写入字典，必须在命令行传 `**--cluster_number`**（需文献/先验或内部约定）。

评测脚本（`repro/evaluate.py`）还需要 **临床表**（生存 + 可选临床协变量），列名需能被 `repro/data_io.py` 识别（与当前 `cBioportal_data.tsv` 风格兼容 easiest）。

---

## 2. 中国大陆较易获取、且医学相关的数据源（推荐思路）

下列来源 **不保证** 每条链路始终可直连，但多为国内科研常用、或提供国内镜像/合规下载途径；**优先选「已有样本对齐好的多组学矩阵」**，可极大减少预处理工作量。

### 2.1 国际公共数据 + 国内镜像 / 学术网络（最贴近论文设定）


| 数据资源                  | 医学相关性                | 与 Subtype-DCC 的匹配度                          | 大陆访问提示                                                                  |
| --------------------- | -------------------- | ------------------------------------------- | ----------------------------------------------------------------------- |
| **TCGA**（多癌种肿瘤）       | 肿瘤基因组学金标准            | **最高**：论文 Table 1 即基于 TCGA 式多组学             | 可通过 **国家生物信息中心 NGDC** 等提供的 **GDC / TCGA 数据镜像或索引** 检索下载；教育网/机构 VPN 常更稳定。 |
| **ICGC** / PCAWG 子集   | 肿瘤                   | 高（需自行对齐到统一 ID）                              | 部分项目数据可在 **NGDC / CNGBdb** 查到登记与下载说明。                                   |
| **GEO**（ArrayExpress） | 疾病与肿瘤大量表达芯片 / 少量-seq | **中**：多为 **单组学** 或双组学；需自己拼多组学或做 **模态缺失** 实验 | 多数时间可访问；大文件建议断点续传。                                                      |


### 2.2 国内门户（检索入口，具体队列需点进项目说明）

- **国家生物信息中心（CNCB）/ NGDC**：[https://ngdc.cncb.ac.cn/](https://ngdc.cncb.ac.cn/) — 检索 TCGA 备份、注册研究、组学归档（**注意数据使用协议与伦理**）。
- **国家基因组科学数据中心 / GSA 等**：关注是否提供 **矩阵级 release**（不仅是原始 fastq），否则预处理成本高。

### 2.3 自建临床队列（医院/多中心合作）

- **适用**：若有 **同一批患者的 RNA +（miRNA / 甲基化 / CNV）中至少两种** 且样本量足够（建议至少与论文同量级讨论，几十例以上更稳）。
- **要点**：脱敏 **患者主键**、统一 **ID 命名**、在文档中写清 **各组学实验平台与版本**。

### 2.4 若实在无法凑齐四种组学

可使用训练参数 `**--exclude_omics`**（逗号分隔：`CN,meth,miRNA,rna` 之一或多个），与仓库内 `repro/run_ablation.py` 思路一致；**簇数 M 与论文表不可直接对比**，需在报告里单独说明。

---

## 3. 推荐工作流（队友 A 换癌种 / 换数据集）

### 阶段 A：定题与下载

1. 选定 **癌种 + 队列**（例如某 TCGA 癌种子集或国内公开注册研究）。
2. 下载或生成四张矩阵：**RNA、miRNA、甲基化、CNV**，颗粒度可为 **gene-level**（与当前 `fea` 一致即可，不必与 TCGA 官方 probe 完全一致，但要 **全流程固定**）。
3. 准备 **临床随访表**（至少：患者 ID、生存时间、删失指示；可选：分期、性别、年龄等用于 `evaluate.py`）。

### 阶段 B：清洗与对齐（最关键）

1. **样本 ID 统一**：全部转为同一套 barcode（如 TCGA 用前 12 位 `patient_barcode` 规则可在 `data_io.patient_barcode` 参考；非 TCGA 则全队统一字符串）。
2. **取交集**：只对 **四种矩阵同时存在** 的样本建列；剔除列名重复、空列。
3. **特征侧处理**（常见做法，可按队列调整）：
   - RNA：log2(count+1) 或 log2(TPM+1)；剔除方差为 0 或缺失过高的基因（设阈值写进 README）。
   - miRNA：log 变换；过滤低表达。
   - 甲基化：β 或 M-value；剔除低方差 probe。
   - CNV：log2 ratio 或离散化后的连续值；与论文一致即可，**关键是跨样本可比**。
4. **不要在矩阵里混入临床表**；临床单独 TSV。

**若从 UCSC Xena 四张 `.gz` 经本仓库脚本转 `*.fea`**：清洗、缺失填补与甲基化 **两阶段高方差探针选择** 的正式说明见 **`docs/XENA_TO_FEA_PREPROCESSING.md`**（可直接用于复现报告「预处理」小节）。

### 阶段 C：写入本仓库约定路径

1. 在 `**vendor/Subtype-DCC/subtype_file/fea/<YOUR_COHORT>/`** 下放置：`rna.fea`、`miRNA.fea`、`meth.fea`、`CN.fea`（或删减模态 + `--exclude_omics`）。
2. 确认 `**scripts/link_fea_junction.ps1**`（若你们用 junction）仍指向该 `fea` 根目录；或保证 `**repro/data_io.resolve_fea_root()**` 能解析到该目录（当前逻辑优先 `vendor/Subtype-DCC/subtype_file/fea`）。

### 阶段 D：配置簇数并训练

在仓库根：

```powershell
cd <repo_root>
.\.venv\Scripts\python.exe vendor\Subtype-DCC\Subtype-DCC\train.py -c YOUR_COHORT --cluster_number 4 --workers 4
```

- `YOUR_COHORT` 为文件夹名（与 `fea/YOUR_COHORT/` 一致）。
- `**--cluster_number**`：新癌种 **必给**，除非你们已把该键加入 `train.py` 的 `cancer_dict`。

产物：

- `vendor/Subtype-DCC/Subtype-DCC/results/YOUR_COHORT.dcc`
- `vendor/Subtype-DCC/Subtype-DCC/results/YOUR_COHORT.fea`

### 阶段 E：生存与临床评测

将新队列临床行并入或单独指定为 `**--clinical`** 可读的 TSV（列名需与 `repro/data_io.py` 中 `_pick_col` 候选名兼容，或小幅扩展 `data_io`）。

```powershell
.\.venv\Scripts\python.exe repro\evaluate.py `
  --labels vendor\Subtype-DCC\Subtype-DCC\results\YOUR_COHORT.dcc `
  --clinical data\clinical\<your_clinical>.tsv -c YOUR_COHORT `
  --out-json repro_out\YOUR_COHORT.metrics.json
```

---

## 4. 仓库里建议补充的内容（方便队友「少踩坑、尽量一键」）


| 补充项                                                           | 目的                                                                             |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `**docs/DATA_PROVENANCE.md**`                                 | 写明：数据来自哪个数据库/文章、下载日期、版本号、过滤规则（可审稿）。                                            |
| `**docs/FEA_FORMAT.md**`                                      | 四个 `*.fea` 的 **示例路径、最小行列数示例、列名规则**（可复制粘贴当契约）。                                  |
| `**data/clinical/README.md`**                                 | 临床表 **必填列**、与 `sample_name` 对齐方式、OS 列名示例。                                      |
| **（可选）`repro/build_fea_<COHORT>.py` 或 notebook**              | 从 Xena / 你们原始表 **一键生成四张 .fea**；这是真正的「一键」核心，强烈建议至少给 **一个癌种的完整可运行示例**。           |
| `**vendor/Subtype-DCC/Subtype-DCC/train.py` 中 `cancer_dict`** | 若长期维护固定队列，可把 `YOUR_COHORT: M` 写进字典，避免忘传 `--cluster_number`。                    |
| `**.env` 或 `config/local_paths.yaml`（勿提交密钥）**                 | 若数据在共享盘，用配置指向绝对路径，避免每人改代码。                                                     |
| `**requirements.txt` 已锁定版本**                                  | 队友 `pip install -r requirements.txt` 后与你们一致；若有 R/SNF，见 `benchmarks/README.md`。 |
| `**README.md` 顶部「队友 Quickstart」**                             | 三行命令：`setup` → `train -c XXX` → `evaluate.py`；链接到本文。                           |


---

## 5. 伦理与合规（医学数据必提）

- 仅使用 **已去识别、已公开发表或已注册** 的数据；院内数据需 **伦理批件** 与脱敏流程，**勿上传 Git**。
- Git 仓库中 **不要提交** 原始可识别患者信息；大矩阵建议 **Git-LFS** 或 **网盘/对象存储 + README 链接**。

---

## 6. 与当前仓库的衔接

- 训练当前工作目录应为 `**vendor/Subtype-DCC/Subtype-DCC`**（`run_subtype_dcc_batch.py` 已处理）；队友若在子目录手跑，注意 `python` 与 `fea` 相对路径。
- 更详细的会话级说明见 `**docs/SESSION_HANDOFF.md**`（§〇 目录结构、§九 批量训练、§十 超参）。

若你们确定一个 **具体新癌种 ID + 临床表样例**，可在 Issues 或附录中贴 **5×5 迷你矩阵** 作单元测试，便于 CI 验证格式。