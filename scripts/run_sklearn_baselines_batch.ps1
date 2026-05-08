<#
.SYNOPSIS
  Batch K-means, Spectral clustering, and NMF+argmax on nine TCGA cohorts (paper Table 1 set).

.DESCRIPTION
  Runs repro/baselines_sklearn.py with the same K as train.py / SNF (vendor Subtype-DCC cancer_dict).
  Requires: repo root .venv, fea under vendor/Subtype-DCC/subtype_file/fea/<CANCER>/.

  Usage (from anywhere):
    powershell -ExecutionPolicy Bypass -File scripts\run_sklearn_baselines_batch.ps1

  From repo root:
    .\scripts\run_sklearn_baselines_batch.ps1

  Optional:
    -DryRun    Print commands only
    -Evaluate  After each TSV, run repro/evaluate.py -> repro_out/<CANCER>.<method>.metrics.json

  See docs/SESSION_HANDOFF.md §二、§九 (venv path); compare bars: repro/run_benchmark_eval.py
#>
[CmdletBinding()]
param(
    [switch] $DryRun,
    [switch] $Evaluate
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $RepoRoot

$Py = Join-Path $RepoRoot ".venv\Scripts\python.exe"
if (-not (Test-Path -LiteralPath $Py)) {
    Write-Error "Missing venv Python: $Py  (cd to repo root; create .venv per README)"
}

# Order: paper nine cohorts; K matches train.py / benchmarks/snf_template.R CANCER_K
$Cancers = @(
    @{ Name = "BRCA"; K = 5 },
    @{ Name = "BLCA"; K = 5 },
    @{ Name = "KIRC"; K = 4 },
    @{ Name = "LUAD"; K = 3 },
    @{ Name = "PAAD"; K = 2 },
    @{ Name = "SKCM"; K = 4 },
    @{ Name = "STAD"; K = 3 },
    @{ Name = "UCEC"; K = 4 },
    @{ Name = "UVM";  K = 4 }
)
$Methods = @("kmeans", "spectral", "nmf")
$Clinical = "data\clinical\cBioportal_data.tsv"

$n = 0
$total = $Cancers.Count * $Methods.Count
foreach ($row in $Cancers) {
    $c = $row.Name
    $k = $row.K
    foreach ($m in $Methods) {
        $n++
        $outTsv = Join-Path $RepoRoot "repro_out\$c.$m.tsv"
        $args = @(
            "repro\baselines_sklearn.py",
            "-c", $c,
            "--method", $m,
            "--k", "$k",
            "--out", $outTsv
        )
        Write-Host "[$n/$total] $c $m (K=$k) -> repro_out\$c.$m.tsv"
        if ($DryRun) {
            Write-Host "  DRY: & `"$Py`" $($args -join ' ')"
            continue
        }
        & $Py @args
        if ($LASTEXITCODE -ne 0) {
            throw "baselines_sklearn failed for $c $m (exit $LASTEXITCODE)"
        }
        if ($Evaluate) {
            $outJson = Join-Path $RepoRoot "repro_out\$c.$m.metrics.json"
            $ev = @(
                "repro\evaluate.py",
                "--labels", $outTsv,
                "--clinical", (Join-Path $RepoRoot $Clinical),
                "-c", $c,
                "--out-json", $outJson
            )
            & $Py @ev
            if ($LASTEXITCODE -ne 0) {
                throw "evaluate failed for $c $m (exit $LASTEXITCODE)"
            }
        }
    }
}

Write-Host "Done. Wrote $($Cancers.Count * $Methods.Count) label TSVs under repro_out\"
if (-not $Evaluate) {
    Write-Host "Optional: re-run with -Evaluate for metrics JSON, or: .\.venv\Scripts\python.exe repro\batch_evaluate.py (Subtype-DCC .dcc only)"
}
