# Run from repo root (subtype-dcc-repro). Creates a directory junction so
# vendor/Subtype-DCC/subtype_file/fea points to Subtype-GAN fea folder.
param(
  [Parameter(Mandatory=$true)]
  [string]$SubtypeGanFeaPath
)
$target = Resolve-Path $SubtypeGanFeaPath
$link = Join-Path $PSScriptRoot "..\vendor\Subtype-DCC\subtype_file\fea"
New-Item -ItemType Directory -Force -Path (Split-Path $link) | Out-Null
if (Test-Path $link) { cmd /c rmdir $link }
cmd /c mklink /J "$link" "$target"
Write-Host "Linked $link -> $target"
