param(
    [ValidateSet("release", "debug")]
    [string]$Profile = "release",
    [switch]$SkipCargoBuild
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $RepoRoot

try {
    # Locate WiX v4 CLI
    $wixCmd = Get-Command wix -ErrorAction SilentlyContinue
    if ($wixCmd) {
        $wixExe = $wixCmd.Source
    } elseif (Test-Path "$env:USERPROFILE\.dotnet\tools\wix.exe") {
        $wixExe = "$env:USERPROFILE\.dotnet\tools\wix.exe"
    } else {
        throw "WiX Toolset v4 ble ikke funnet. Installer med: dotnet tool install --global wix && wix extension add -g WixToolset.UI.wixext"
    }

    if (-not $SkipCargoBuild) {
        Write-Host "==> Bygger cluster_overview.exe ($Profile)..." -ForegroundColor Cyan
        if ($Profile -eq "release") {
            cargo build --release -j 2
        } else {
            cargo build -j 2
        }
        if ($LASTEXITCODE -ne 0) {
            throw "Cargo build feilet med kode $LASTEXITCODE"
        }
    }

    $binPath = Join-Path $RepoRoot "target\$Profile\cluster_overview.exe"
    if (-not (Test-Path $binPath)) {
        throw "Fant ikke kjørbar fil: $binPath"
    }

    $outDir = Join-Path $RepoRoot "wix\out"
    if (-not (Test-Path $outDir)) {
        New-Item -ItemType Directory -Path $outDir | Out-Null
    }
    $outMsi = Join-Path $outDir "ClusterOverview-1.0.0-x64.msi"

    Write-Host "==> Bygger WiX v4 MSI-installasjonsprogram..." -ForegroundColor Cyan
    & $wixExe build "wix\Package.wxs" `
        -ext WixToolset.UI.wixext `
        -arch x64 `
        -d "BinPath=$binPath" `
        -o $outMsi

    if ($LASTEXITCODE -ne 0) {
        throw "WiX build feilet med kode $LASTEXITCODE"
    }

    $msiInfo = Get-Item $outMsi
    $sizeMb = [Math]::Round($msiInfo.Length / 1MB, 2)
    Write-Host ""
    Write-Host "FERDIG! MSI-installasjonsfil bygget:" -ForegroundColor Green
    Write-Host "  Fil:       $($msiInfo.FullName)"
    Write-Host "  Størrelse: $sizeMb MB"
    Write-Host "  Opphavsrett: Copyright (c) 2026 Robert Larsen"
} finally {
    Pop-Location
}
