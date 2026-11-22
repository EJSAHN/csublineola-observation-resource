param(
    [Parameter(Mandatory=$true)][string]$Data,
    [Parameter(Mandatory=$true)][string]$Out,
    [string]$Python = 'python',
    [string]$SD13 = '',
    [string]$MaskData = ''
)
$ErrorActionPreference = 'Stop'
$runArgs = @((Join-Path $PSScriptRoot 'run.py'), '--data', $Data, '--out', $Out)
if ($SD13) { $runArgs += @('--sd13', $SD13) }
if ($MaskData) { $runArgs += @('--maskdata', $MaskData) }
& $Python @runArgs
if ($LASTEXITCODE -ne 0) { throw "Replay failed with exit code $LASTEXITCODE" }
