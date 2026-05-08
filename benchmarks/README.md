# Baselines not shipped as Python in Subtype-DCC

These methods match the paper Table 1 listing; run them in **R** with the same sample-aligned feature matrices as in `subtype_file/fea/{CANCER}/`.

| Method    | Suggested implementation |
|-----------|---------------------------|
| Subtype-GAN | [haiyang1986/Subtype-GAN](https://github.com/haiyang1986/Subtype-GAN) TensorFlow 1.x env; `-m SubtypeGAN` |
| SNF       | R package `SNFtool` — `SNF()` then spectral clustering with fixed `K` |
| NEMO      | R package `NEMO` |
| PINS      | `PINSPlus` / legacy PINS; see `pins_template.R` |
| MCCA      | `PMA` or `RGCCA`; cluster latent scores with fixed `K` |
| iCluster  | Bioconductor `iClusterPlus` |
| LRACluster| Author [LRAcluster](https://www.bioconductor.org/packages/release/bioc/html/LRAcluster.html) or CRAN equivalent |

**Fair comparison:** use the same **K** as `train.py` `cancer_dict`, identical patient ordering as in `rna.fea` columns, and repeat stochastic methods **5×** (paper) then average metrics from `repro/evaluate.py`.

Templates in this folder sketch data paths only; install Bioconductor dependencies before sourcing.
