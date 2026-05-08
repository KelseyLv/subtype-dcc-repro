<#
.SYNOPSIS
  Copy nine-cohort × five-method Table 1 reproduction outputs from repro_out/ into docs/table1_artifacts/ for git.

.DESCRIPTION
  Reads from repro_out/ (this fork may version repro_out/ in full; the script still copies a curated Table 1 subset into docs/ for stable links). Mirrors only small, reviewable artifacts:
  - table1_reproduction_metrics.json, table1_reproduction_vs_paper.{json,md}
  - repro_out/table1_eval/*.json
  - figure1_* boxplot PNGs
  - Per-cohort label TSVs for SNF + sklearn (snf, kmeans, spectral, nmf) for the paper nine cancers
  - Optional: batch_metrics_summary.json (Subtype-DCC .dcc batch only)

  Does NOT copy: *.fea, *.dcc, manuscript caches, LIHC-only trees, or other large repro_out content.

  Usage (repo root):
    powershell -ExecutionPolicy Bypass -File scripts\sync_table1_artifacts_to_docs.ps1
#>
[CmdletBinding()]
param(
    [switch] $DryRun,
    [switch] $SkipBatchDccSummary
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Src = Join-Path $RepoRoot "repro_out"
$Dst = Join-Path $RepoRoot "docs\table1_artifacts"

$Cancers = @(
    "BRCA", "BLCA", "KIRC", "LUAD", "PAAD", "SKCM", "STAD", "UCEC", "UVM"
)
$SklearnLike = @("snf", "kmeans", "spectral", "nmf")

function Copy-One([string]$Rel, [string]$DestRel = $null) {
    $from = Join-Path $RepoRoot $Rel
    if (-not (Test-Path -LiteralPath $from)) {
        Write-Warning "Missing (skip): $Rel"
        return
    }
    $targetRel = if ($DestRel) { $DestRel } else { $Rel -replace "^repro_out\\", "docs\table1_artifacts\" }
    $to = Join-Path $RepoRoot $targetRel
    $dir = Split-Path -Parent $to
    if ($DryRun) {
        Write-Host "COPY $Rel -> $targetRel"
        return
    }
    if (-not (Test-Path -LiteralPath $dir)) {
        New-Item -ItemType Directory -Path $dir -Force | Out-Null
    }
    Copy-Item -LiteralPath $from -Destination $to -Force
    Write-Host "OK $targetRel"
}

if (-not (Test-Path -LiteralPath $Src)) {
    Write-Error "repro_out not found at $Src — run batch_evaluate_table1_methods.py first."
}

if (-not $DryRun) {
    New-Item -ItemType Directory -Path $Dst -Force | Out-Null
}

Copy-One "repro_out\table1_reproduction_metrics.json" "docs\table1_artifacts\table1_reproduction_metrics.json"
Copy-One "repro_out\table1_reproduction_vs_paper.json" "docs\table1_artifacts\table1_reproduction_vs_paper.json"
Copy-One "repro_out\table1_reproduction_vs_paper.md" "docs\table1_artifacts\table1_reproduction_vs_paper.md"

$evalSrc = Join-Path $Src "table1_eval"
$evalDst = Join-Path $Dst "table1_eval"
if (Test-Path -LiteralPath $evalSrc) {
    if ($DryRun) {
        Write-Host "COPY repro_out\table1_eval\* -> docs\table1_artifacts\table1_eval\"
    }
    else {
        New-Item -ItemType Directory -Path $evalDst -Force | Out-Null
        Copy-Item (Join-Path $evalSrc "*") -Destination $evalDst -Force -Recurse:$false
        Write-Host "OK docs\table1_artifacts\table1_eval\ ($((Get-ChildItem $evalDst -File).Count) files)"
    }
}
else {
    Write-Warning "Missing repro_out\table1_eval — run batch_evaluate_table1_methods.py"
}

foreach ($pat in @("figure1_neg_log10_logrank_p_boxplot.png", "figure1_n_significant_clinical_boxplot.png", "figure1_table1_boxplots_combined.png")) {
    Copy-One "repro_out\$pat" "docs\table1_artifacts\$pat"
}

foreach ($c in $Cancers) {
    foreach ($m in $SklearnLike) {
        $fn = "$c.$m.tsv"
        Copy-One "repro_out\$fn" "docs\table1_artifacts\label_tsv\$fn"
    }
}

if (-not $SkipBatchDccSummary) {
    Copy-One "repro_out\batch_metrics_summary.json" "docs\table1_artifacts\batch_metrics_summary.json"
}

Write-Host "`nDone. Next: git add docs/table1_artifacts && git commit"
