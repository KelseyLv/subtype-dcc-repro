# Subtype-DCC paper reproduction kit

This folder vendors a slightly patched [Subtype-DCC](https://github.com/zhaojingo/Subtype-DCC) training stack and adds **evaluation**, **sklearn baselines**, **plotting**, and **ablation** runners aligned with *Briefings in Bioinformatics* (2023), without modifying your attached plan file.

## Layout

| Path | Purpose |
|------|---------|
| [vendor/Subtype-DCC/Subtype-DCC/](vendor/Subtype-DCC/Subtype-DCC/) | Training code (`train.py`, `dataloader.py`, `modules/`, `config/`) |
| [vendor/Subtype-DCC/subtype_file/fea/](vendor/Subtype-DCC/subtype_file/README.txt) | **You populate**: four `*.fea` matrices per cancer (see Subtype-GAN) |
| [repro/](repro/) | Metrics, plots, sklearn baselines, batch helpers |
| [benchmarks/](benchmarks/) | Notes + R templates for SNF / NEMO / PINS / iCluster / MCCA / LRACluster |

## 1. Environment

```bash
cd subtype-dcc-repro
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
# Install PyTorch with CUDA from https://pytorch.org if training on GPU
```

## 2. Data

This repository may already include **`data/clinical/cBioportal_data.tsv`** (team-shared GerkeLab extract, ~9 MB). If that file is present after `git clone`, you can **skip** the clinical download step in Option A and only run the Subtype-GAN junction / `fea` setup you still need.

### Option A — one PowerShell script (recommended)

From the repo root, run:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1
```

This will: **clone Subtype-GAN** into `external/Subtype-GAN`, create a **junction**  
`vendor/Subtype-DCC/subtype_file/fea` → `external/Subtype-GAN/fea`, and **download**  
GerkeLab `data/clinical/cBioportal_data.tsv` (~9 MB) for use with `repro/evaluate.py`  
(alternative to legacy `clinical_PANCAN_patient_with_followup.tsv`).  
A transcript is written to `data_setup_transcript.log`.

Clone URLs: **defaults try GitHub first**, then `gitclone.com` mirror (mirrors often return **502**; if both fail, use VPN or download the repo ZIP manually into `external/Subtype-GAN`). Optional flags:

```powershell
# Prefer mirror first
powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1 -MirrorFirst

# Single URL only
powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1 -CloneUrl https://github.com/haiyang1986/Subtype-GAN.git
```

If you **manually downloaded** GitHub ZIP (`Subtype-GAN-master`), create the junction only:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\link_fea_junction.ps1
```

### Option B — manual

1. Clone `https://github.com/haiyang1986/Subtype-GAN` (large `fea/` matrices).
2. Link `vendor/Subtype-DCC/subtype_file/fea` → `Subtype-GAN/fea`:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\setup_windows.ps1 -SubtypeGanFeaPath "D:\path\to\Subtype-GAN\fea"
```

3. Clinical table: either **`clinical_PANCAN_patient_with_followup.tsv`** (legacy PANCAN; Subtype-GAN default) or run  
   `python scripts/download_clinical.py --source gerke --out data/clinical/cBioportal_data.tsv`.  
   Pass the path with `--clinical` below.

## 3. Train Subtype-DCC (nine cancers)

```bash
python repro/run_subtype_dcc_batch.py --extra-args "--epochs 600"
# Smoke test:
# python repro/run_subtype_dcc_batch.py --extra-args "--epochs 2" --cancers BRCA
```

Checkpoints and outputs land in `vendor/Subtype-DCC/Subtype-DCC/save/` and `vendor/Subtype-DCC/Subtype-DCC/results/` (`*.dcc`, `*.fea`).

## 4. Paper metrics (log-rank −log10 p, six clinical tests)

```bash
python repro/evaluate.py --labels vendor/Subtype-DCC/Subtype-DCC/results/BRCA.dcc --clinical PATH\clinical_PANCAN_patient_with_followup.tsv --cancer BRCA --out-json repro_out\BRCA.metrics.json
```

**Batch all five label sources × nine cancers** (Subtype-DCC `.dcc`, `repro_out/*.snf.tsv`, kmeans/spectral/nmf TSV), then compare to bundled paper Table 1 (SNF row filled; add Subtype-DCC / sklearn rows from the PDF into `data/paper_table1_reference.json` if you want full Δ columns):

```bash
.venv/Scripts/python.exe repro/batch_evaluate_table1_methods.py
.venv/Scripts/python.exe repro/table1_compare_to_paper.py
```

Outputs: `repro_out/table1_reproduction_metrics.json`, `repro_out/table1_reproduction_vs_paper.{json,md}`, per-method slims under `repro_out/table1_eval/`.

### Git-friendly snapshot (`docs/table1_artifacts/`)

In this repository, **`repro_out/` may be tracked in full** for team sharing. You can still maintain **`docs/table1_artifacts/`** as a compact mirror (same Table 1 JSON/TSV/PNGs) for stable documentation links. After you finish the five methods × nine cancers pipeline (batch evaluate, `table1_compare_to_paper.py`, and optional Figure 1 plots), refresh that mirror with:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\sync_table1_artifacts_to_docs.ps1
```

That script fills **`docs/table1_artifacts/`** with: `table1_reproduction_metrics.json`, `table1_reproduction_vs_paper.{json,md}`, `table1_eval/*.json`, `figure1_*.png`, **`label_tsv/`** (nine cohorts × SNF + kmeans + spectral + nmf), and **`batch_metrics_summary.json`** (Subtype-DCC `.dcc` batch from `batch_evaluate.py`). **Training checkpoints (`save/`), `vendor/.../results/*.dcc`, `*.fea`, Xena `.gz`, full clinical TSV, and large manuscript cache trees are not copied** — others still train or place `.dcc` locally, then re-run `batch_evaluate_table1_methods.py` if they need to refresh numbers.

Regenerate Figure 1 from the committed summary JSON:

```bash
.venv/Scripts/python.exe repro/plot_figure1_table1_boxplots.py --metrics-json docs/table1_artifacts/table1_reproduction_metrics.json --out-dir docs/table1_artifacts
```

## 5. Sklearn baselines (K-means, spectral, NMF)

Uses the **same MinMax-scaled concatenated matrix** as Subtype-DCC (`repro/data_io.py`).

**Nine cohorts × three methods (PowerShell, from repo root):** `scripts/run_sklearn_baselines_batch.ps1` (uses `.\.venv\Scripts\python.exe` and Table 1 `K` per cancer). Optional `-Evaluate` appends `repro/evaluate.py` per run; `-DryRun` prints commands only.

```bash
python repro/baselines_sklearn.py -c BRCA --method kmeans --k 5 --out repro_out/BRCA.kmeans.tsv
python repro/baselines_sklearn.py -c BRCA --method spectral --k 5 --out repro_out/BRCA.spectral.tsv
python repro/baselines_sklearn.py -c BRCA --method nmf --k 5 --out repro_out/BRCA.nmf.tsv
```

Evaluate each TSV with `repro/evaluate.py` (`cluster` column).

## 6. Figures

```bash
python repro/plot_survival.py --labels vendor/Subtype-DCC/Subtype-DCC/results/KIRC.dcc --clinical PATH\clinical.tsv --cancer KIRC --out figures/KIRC_survival.png
python repro/plot_tsne.py --fea vendor/Subtype-DCC/Subtype-DCC/results/KIRC.fea --labels vendor/Subtype-DCC/Subtype-DCC/results/KIRC.dcc --out figures/KIRC_tsne.png
python repro/run_benchmark_eval.py --clinical PATH\clinical.tsv --cancer BRCA --labels SubtypeDCC=vendor/Subtype-DCC/Subtype-DCC/results/BRCA.dcc KMeans=repro_out/BRCA.kmeans.tsv --out repro_out/BRCA.summary.json
python repro/plot_benchmark_bars.py --summary-json repro_out/BRCA.summary.json --out-prefix figures/BRCA_bench
```

## 7. Ablation (drop one omics)

```bash
python repro/run_ablation.py -c KIRC --drop rna
# Or run all four single-omic removals:
python repro/run_ablation.py -c KIRC
```

## 8. Other baselines (Subtype-GAN, SNF, NEMO, …)

See [benchmarks/README.md](benchmarks/README.md). Subtype-GAN uses a legacy TensorFlow 1 environment—keep it separate from this PyTorch env.

## Notes on patches vs upstream

- `matplotlib` uses **Agg** backend; loss curves save without blocking on `plt.show()`.
- `Encoder` **input dimension is inferred** from data so ablations change width correctly.
- `dataloader.py` supports `--exclude_omics` for ablation.
- Training saves `checkpoint_{epochs}.tar` each epoch (upstream behaviour).
