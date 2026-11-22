param(
    [Parameter(Mandatory=$true)][string]$Data,
    [Parameter(Mandatory=$true)][string]$Out,
    [string]$Python = 'python',
    [Alias('SD13')][string]$Followup = '',
    [Alias('MaskData')][string]$Masks = ''
)
$ErrorActionPreference = 'Stop'
$runArgs = @((Join-Path $PSScriptRoot 'run.py'), '--data', $Data, '--out', $Out)
if ($Followup) { $runArgs += @('--followup', $Followup) }
if ($Masks) { $runArgs += @('--masks', $Masks) }
& $Python @runArgs
if ($LASTEXITCODE -ne 0) { throw "Verification failed with exit code $LASTEXITCODE" }
