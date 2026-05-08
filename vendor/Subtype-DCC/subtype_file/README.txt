Place multi-omics feature matrices here:

  fea/BRCA/CN.fea
  fea/BRCA/meth.fea
  fea/BRCA/miRNA.fea
  fea/BRCA/rna.fea
  ... (nine cancers)

Easiest source: clone https://github.com/haiyang1986/Subtype-GAN and junction-link this `fea`
folder to Subtype-GAN's `fea` directory (see scripts/setup_windows.ps1).

train.py resolves paths as: ../subtype_file/fea/ from vendor/Subtype-DCC/Subtype-DCC/
