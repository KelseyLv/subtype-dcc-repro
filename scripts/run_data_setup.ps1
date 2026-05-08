# One-shot: clone Subtype-GAN, junction-link fea for Subtype-DCC, download GerkeLab clinical TSV.
#
# 用法：
#   powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1
#   （默认：先试 GitHub，再试 gitclone 镜像；任一成功即可）
#
#   只走镜像：  -MirrorFirst
#   指定单个 URL（不重试）：  -CloneUrl "https://github.com/haiyang1986/Subtype-GAN.git"
#
param(
    [string]$CloneUrl = "",
    [switch]$MirrorFirst
)

$ErrorActionPreference = 'Stop'

$defaultGithub = "https://github.com/haiyang1986/Subtype-GAN.git"
$defaultMirror = "https://gitclone.com/github.com/haiyang1986/Subtype-GAN.git"

if ($CloneUrl) {
    $candidates = @($CloneUrl)
} elseif ($MirrorFirst) {
    $candidates = @($defaultMirror, $defaultGithub)
} else {
    # gitclone 常见 502，默认优先官方 GitHub
    $candidates = @($defaultGithub, $defaultMirror)
}

$base = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$log = Join-Path $base 'data_setup_transcript.log'
Start-Transcript -Path $log -Force

try {
    Set-Location $base
    Write-Host "Repo root: $base"

    $sg = Join-Path $base 'external\Subtype-GAN'
    $feaTarget = Join-Path $sg 'fea'

    if (-not (Test-Path $feaTarget)) {
        $cloneOk = $false
        foreach ($url in $candidates) {
            if (Test-Path $sg) {
                Write-Host "Removing incomplete folder: $sg"
                Remove-Item -LiteralPath $sg -Recurse -Force -ErrorAction Stop
            }
            Write-Host ""
            Write-Host ">>> git clone --depth 1"
            Write-Host ">>> $url"
            & git clone --depth 1 $url $sg 2>&1 | ForEach-Object { Write-Host $_ }

            if (($LASTEXITCODE -eq 0) -and (Test-Path $feaTarget)) {
                $cloneOk = $true
                Write-Host "Clone OK (fea folder present)."
                break
            }
            Write-Warning "Failed or incomplete for this URL (exit $LASTEXITCODE). Trying next candidate if any..."
            if (Test-Path $sg) {
                Remove-Item -LiteralPath $sg -Recurse -Force -ErrorAction SilentlyContinue
            }
        }

        if (-not $cloneOk) {
            throw @"
Subtype-GAN clone failed — no fea folder after trying:
  $($candidates -join "`n  ")

Common fixes:
  - Retry later (mirror 502 / network).
  - Manual clone into external\Subtype-GAN using VPN / browser ZIP from GitHub.
  - Force GitHub only:
      powershell -ExecutionPolicy Bypass -File scripts\run_data_setup.ps1 -CloneUrl https://github.com/haiyang1986/Subtype-GAN.git
"@
        }
    } else {
        Write-Host "Subtype-GAN with fea\ already present: $sg"
    }

    $subtypeRoot = Join-Path $base 'vendor\Subtype-DCC\subtype_file'
    New-Item -ItemType Directory -Path $subtypeRoot -Force | Out-Null
    $feaLink = Join-Path $subtypeRoot 'fea'
    if (Test-Path $feaLink) {
        cmd /c "rmdir `"$feaLink`""
    }
    cmd /c "mklink /J `"$feaLink`" `"$feaTarget`""

    $clinicalOut = Join-Path $base 'data\clinical\cBioportal_data.tsv'
    python (Join-Path $base 'scripts\download_clinical.py') --source gerke --out $clinicalOut

    Write-Host "Clinical header (first line):"
    Get-Content $clinicalOut -TotalCount 1
    Write-Host "Transcript saved to $log"
}
finally {
    Stop-Transcript
}
