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

## 仓库里与 LIHC 相关的快照（无需先重跑即可浏览）

- **[`docs/lihc_artifacts/`](../lihc_artifacts/)**：`LIHC.metrics.json`、loss 与 manuscript 风格 PNG。  
- **`repro_out/`**（若已随仓库提交）：本机全量输出目录，含 LIHC 与其他实验的中间结果；体积较大时克隆会久一些。

## 你扩展新癌种时至少要自备 / 自跑的部分

- **`vendor/Subtype-DCC/subtype_file/fea/<CANCER>/` 下四个 `*.fea`**（或按文档从 Xena `.gz` 生成；**`.fea` / `.gz` 通常仍不入 Git**）。  
- **`vendor/Subtype-DCC/Subtype-DCC/results/<CANCER>.dcc`**（训练产出；`results/` 仍被 ignore 时需本地训练或向维护者索取）。  
- **临床表**：本仓库已提供 **`data/clinical/cBioportal_data.tsv`** 时，可直接用于 `repro/evaluate.py`；若与论文 PANCAN 表不一致，报告中需声明。

## 一句话任务描述（可复制给队友）

> Clone 本仓库后先看 **`docs/LIHC_EXTENSION_SUBTYPE_DCC.md`**，按其中流程把 `LIHC` 换成你的癌种；预处理与特征矩阵看 **`docs/XENA_TO_FEA_PREPROCESSING.md`**；环境与命令细节查 **`docs/SESSION_HANDOFF.md`** 和根目录 **`README.md`**。临床表在 **`data/clinical/`**；全量运行输出在 **`repro_out/`**（若已提交）。
