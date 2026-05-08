# Subtype-DCC paper reproduction kit

Vendors a lightly patched [Subtype-DCC](https://github.com/zhaojingo/Subtype-DCC) training stack and adds **evaluation**, **sklearn baselines**, **SNF / Table 1 batching**, **plotting**, and **ablation** helpers aligned with *Briefings in Bioinformatics* (2023). This fork also carries **team-facing snapshots**: bundled clinical table, versioned **`repro_out/`**, and curated **`docs/*_artifacts/`** for GitHub-friendly links.

## Repository layout

| Path | Purpose |
|------|---------|
| [vendor/Subtype-DCC/Subtype-DCC/](vendor/Subtype-DCC/Subtype-DCC/) | Training code (`train.py`, `dataloader.py`, `modules/`, `config/`) |
| [vendor/Subtype-DCC/subtype_file/fea/](vendor/Subtype-DCC/subtype_file/README.txt) | **You link or generate**: four `*.fea` matrices per cancer (Subtype-GAN junction or `repro/xena_gz_to_subtype_dcc_fea.py`) |
| [repro/](repro/) | Metrics, plots, sklearn baselines, batch Table 1, LIHC / manuscript helpers |
| [repro_out/](repro_out/) | **Tracked in this fork**: metrics JSON, label TSVs, Table 1 outputs, manuscript caches/PNGs, logs |
| [data/clinical/](data/clinical/) | **Tracked**: `cBioportal_data.tsv` (GerkeLab extract) for `repro/evaluate.py` |
| [docs/lihc_artifacts/](docs/lihc_artifacts/) | Small LIHC snapshot (metrics + figures) for documentation |
| [docs/table1_artifacts/](docs/table1_artifacts/) | Compact mirror of nine-cohort × five-method Table 1 artifacts |
| [docs/EXTENSION_FOR_NEW_CANCER.md](docs/EXTENSION_FOR_NEW_CANCER.md) | **Onboarding** for teammates extending to new cancers |
| [benchmarks/](benchmarks/) | R templates for SNF / NEMO / PINS / iCluster / MCCA / LRACluster |

**Still not in Git (too large or machine-local):** `vendor/.../Subtype-DCC/save/` (checkpoints), `vendor/.../results/*.dcc` / `*.fea`, Xena `*.gz` under `fea/`. Train or place `.dcc` locally, then evaluate.

## Quick start after `git clone`

1. **Python env:** §1 below.  
2. **Clinical:** use **`data/clinical/cBioportal_data.tsv`** in all `repro/evaluate.py` examples (already in the repo).  
3. **Features (`*.fea`):** obtain via §2 (Subtype-GAN junction or your own preprocessing).  
4. **New cancer extension:** read **[docs/EXTENSION_FOR_NEW_CANCER.md](docs/EXTENSION_FOR_NEW_CANCER.md)** then **[docs/LIHC_EXTENSION_SUBTYPE_DCC.md](docs/LIHC_EXTENSION_SUBTYPE_DCC.md)** as the end-to-end template.

## 1. Environment

```bash
cd subtype-dcc-repro
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
# PyTorch with CUDA: https://pytorch.org if training on GPU
```

## 2. Data (`fea` + clinical)

### Option A — PowerShell one-shot (Subtype-GAN + junction + clinical)

From the repo root:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1
```

This **clones Subtype-GAN** into `external/Subtype-GAN`, creates a **junction**  
`vendor/Subtype-DCC/subtype_file/fea` → `external/Subtype-GAN/fea`, and runs **`scripts/download_clinical.py`** so `data/clinical/cBioportal_data.tsv` exists (or is refreshed from the same GerkeLab source as the committed file). A transcript is written to `data_setup_transcript.log`.

If **`fea\` is already present** under `external/Subtype-GAN`, the script skips cloning. If you **only** need the junction (e.g. you already cloned Subtype-GAN manually):

```powershell
powershell -ExecutionPolicy Bypass -File scripts\link_fea_junction.ps1
```

Clone URLs default to **GitHub**, then `gitclone.com` (mirrors often **502**). Optional flags:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1 -MirrorFirst
powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1 -CloneUrl https://github.com/haiyang1986/Subtype-GAN.git
```

### Option B — manual

1. Clone `https://github.com/haiyang1986/Subtype-GAN` (large `fea/` matrices).  
2. Link `vendor/Subtype-DCC/subtype_file/fea` → `Subtype-GAN/fea`:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\setup_windows.ps1 -SubtypeGanFeaPath "D:\path\to\Subtype-GAN\fea"
```

3. Clinical: either rely on **`data/clinical/cBioportal_data.tsv`** from this repo or regenerate with  
   `python scripts/download_clinical.py --source gerke --out data/clinical/cBioportal_data.tsv`.  
   Pass **`--clinical data/clinical/cBioportal_data.tsv`** (or your PANCAN path) to `repro/evaluate.py`.

## 3. Train Subtype-DCC (nine cancers)

```bash
python repro/run_subtype_dcc_batch.py --extra-args "--epochs 600"
# Smoke test:
# python repro/run_subtype_dcc_batch.py --extra-args "--epochs 2" --cancers BRCA
```

Checkpoints and training outputs go to `vendor/Subtype-DCC/Subtype-DCC/save/` and `vendor/Subtype-DCC/Subtype-DCC/results/` (`*.dcc`, training-time `*.fea` where applicable). Those paths remain **gitignored**.

## 4. Paper metrics (log-rank −log10 p, six clinical tests)

```bash
python repro/evaluate.py --labels vendor/Subtype-DCC/Subtype-DCC/results/BRCA.dcc --clinical data/clinical/cBioportal_data.tsv --cancer BRCA --out-json repro_out/BRCA.metrics.json
```

**Batch all five label sources × nine cancers** (Subtype-DCC `.dcc`, `repro_out/*.snf.tsv`, kmeans/spectral/nmf TSV), then compare to bundled paper Table 1 (SNF row filled; extend `data/paper_table1_reference.json` from the PDF if you want full Δ columns):

```bash
.venv/Scripts/python.exe repro/batch_evaluate_table1_methods.py
.venv/Scripts/python.exe repro/table1_compare_to_paper.py
```

Outputs: `repro_out/table1_reproduction_metrics.json`, `repro_out/table1_reproduction_vs_paper.{json,md}`, slims under `repro_out/table1_eval/`.

### `docs/table1_artifacts/` mirror (optional)

`scripts/sync_table1_artifacts_to_docs.ps1` copies a **curated** Table 1 subset into **`docs/table1_artifacts/`** (metrics, `table1_eval` JSON, Figure 1 PNGs, `label_tsv/`, `batch_metrics_summary.json`) so documentation can link to stable paths under `docs/`.

**Not copied into that mirror** (by design): training checkpoints, `.dcc`, `.fea`, Xena `.gz`, and large `repro_out/*_manuscript_cache` trees—those stay under `repro_out/` in this fork when present, or only on your machine in a slimmer checkout.

Refresh the mirror after re-running the batch + compare + Figure 1 plots:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\sync_table1_artifacts_to_docs.ps1
```

Regenerate Figure 1 from a metrics JSON:

```bash
.venv/Scripts/python.exe repro/plot_figure1_table1_boxplots.py --metrics-json docs/table1_artifacts/table1_reproduction_metrics.json --out-dir docs/table1_artifacts
```

## 5. Sklearn baselines (K-means, spectral, NMF)

Same MinMax-scaled concatenated matrix as Subtype-DCC (`repro/data_io.py`).

**Nine cohorts × three methods (PowerShell, repo root):** `scripts/run_sklearn_baselines_batch.ps1` (Table 1 `K` per cancer). `-Evaluate` runs `repro/evaluate.py` per TSV; `-DryRun` prints commands only.

```bash
python repro/baselines_sklearn.py -c BRCA --method kmeans --k 5 --out repro_out/BRCA.kmeans.tsv
python repro/baselines_sklearn.py -c BRCA --method spectral --k 5 --out repro_out/BRCA.spectral.tsv
python repro/baselines_sklearn.py -c BRCA --method nmf --k 5 --out repro_out/BRCA.nmf.tsv
```

Evaluate each TSV with `repro/evaluate.py` (`cluster` column).

## 6. Figures

`figures/` is gitignored; write PNGs there or under `repro_out/`.

```bash
python repro/plot_survival.py --labels vendor/Subtype-DCC/Subtype-DCC/results/KIRC.dcc --clinical data/clinical/cBioportal_data.tsv --cancer KIRC --out figures/KIRC_survival.png
python repro/plot_tsne.py --fea vendor/Subtype-DCC/Subtype-DCC/results/KIRC.fea --labels vendor/Subtype-DCC/Subtype-DCC/results/KIRC.dcc --out figures/KIRC_tsne.png
python repro/run_benchmark_eval.py --clinical data/clinical/cBioportal_data.tsv --cancer BRCA --labels SubtypeDCC=vendor/Subtype-DCC/Subtype-DCC/results/BRCA.dcc KMeans=repro_out/BRCA.kmeans.tsv --out repro_out/BRCA.summary.json
python repro/plot_benchmark_bars.py --summary-json repro_out/BRCA.summary.json --out-prefix figures/BRCA_bench
```

## 7. Ablation (drop one omics)

```bash
python repro/run_ablation.py -c KIRC --drop rna
python repro/run_ablation.py -c KIRC
```

## 8. Other baselines (Subtype-GAN, SNF, NEMO, …)

See [benchmarks/README.md](benchmarks/README.md). Subtype-GAN expects a legacy TensorFlow 1 stack—keep it separate from this PyTorch env.

## Notes on patches vs upstream

- `matplotlib` uses **Agg**; loss curves save without blocking on `plt.show()`.
- `Encoder` **input dimension is inferred** from data so ablations resize correctly.
- `dataloader.py` supports `--exclude_omics` for ablation.
- Training saves `checkpoint_{epochs}.tar` each epoch (upstream behaviour).
