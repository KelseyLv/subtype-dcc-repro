# Create directory junction: vendor\Subtype-DCC\subtype_file\fea -> your Subtype-GAN fea folder
# Run from repo root:
#   powershell -ExecutionPolicy Bypass -File scripts\link_fea_junction.ps1
# Or specify source (ZIP 解压常见为 Subtype-GAN-master)：
#   powershell -ExecutionPolicy Bypass -File scripts\link_fea_junction.ps1 -FeaSource "external\Subtype-GAN-master\fea"
param(
    [string]$FeaSource = ""
)
$ErrorActionPreference = 'Stop'
$base = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)

if (-not $FeaSource) {
    $candidates = @(
        "external\Subtype-GAN\fea",
        "external\Subtype-GAN-master\fea"
    )
    foreach ($rel in $candidates) {
        $p = Join-Path $base $rel
        if (Test-Path $p) {
            $FeaSource = $p
            break
        }
    }
} else {
    if (-not [System.IO.Path]::IsPathRooted($FeaSource)) {
        $FeaSource = Join-Path $base $FeaSource
    }
}

if (-not $FeaSource -or -not (Test-Path $FeaSource)) {
    throw "Could not find fea folder. Pass -FeaSource path to folder containing BRCA\CN.fea etc."
}

$subtypeRoot = Join-Path $base 'vendor\Subtype-DCC\subtype_file'
New-Item -ItemType Directory -Path $subtypeRoot -Force | Out-Null
$feaLink = Join-Path $subtypeRoot 'fea'

if (Test-Path $feaLink) {
    Write-Host "Removing existing fea link/dir: $feaLink"
    cmd /c "rmdir `"$feaLink`""
}

Write-Host "Junction:"
Write-Host "  $feaLink"
Write-Host "  -> $FeaSource"
cmd /c "mklink /J `"$feaLink`" `"$FeaSource`""

if (-not (Test-Path (Join-Path $feaLink 'BRCA\CN.fea'))) {
    Write-Warning "Probe failed: $($feaLink)\BRCA\CN.fea not visible — check paths."
} else {
    Write-Host "OK: BRCA\CN.fea visible through junction."
}
