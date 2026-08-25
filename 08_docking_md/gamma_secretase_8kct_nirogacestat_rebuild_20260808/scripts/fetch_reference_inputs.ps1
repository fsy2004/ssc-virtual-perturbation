[CmdletBinding()]
param(
    [string]$Destination
)

$ErrorActionPreference = 'Stop'
$scriptRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$packageRoot = Split-Path -Parent $scriptRoot
if (-not $Destination) {
    $Destination = Join-Path $packageRoot 'inputs\reference'
}

$resolvedPackage = [System.IO.Path]::GetFullPath($packageRoot)
$resolvedDestination = [System.IO.Path]::GetFullPath($Destination)
if (-not $resolvedDestination.StartsWith($resolvedPackage, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Destination must remain inside the rebuild package: $resolvedPackage"
}

New-Item -ItemType Directory -Path $resolvedDestination -Force | Out-Null
$downloads = @(
    @{ Name = '8KCT.cif'; Url = 'https://files.rcsb.org/download/8KCT.cif' },
    @{ Name = '8KCT.pdb'; Url = 'https://files.rcsb.org/download/8KCT.pdb' },
    @{ Name = 'O6U.cif'; Url = 'https://files.rcsb.org/ligands/download/O6U.cif' },
    @{ Name = 'O6U_ideal.sdf'; Url = 'https://files.rcsb.org/ligands/download/O6U_ideal.sdf' },
    @{ Name = '8kct_full_validation.pdf.gz'; Url = 'https://files.rcsb.org/pub/pdb/validation_reports/kc/8kct/8kct_full_validation.pdf.gz' }
)

foreach ($item in $downloads) {
    $target = Join-Path $resolvedDestination $item.Name
    if (Test-Path -LiteralPath $target) {
        throw "Refusing to overwrite an existing frozen input: $target"
    }
    Invoke-WebRequest -Uri $item.Url -OutFile $target -UseBasicParsing
    if ((Get-Item -LiteralPath $target).Length -eq 0) {
        throw "Downloaded file is empty: $target"
    }
}

$compressedValidation = Join-Path $resolvedDestination '8kct_full_validation.pdf.gz'
$validationPdf = Join-Path $resolvedDestination '8kct_full_validation.pdf'
$inputStream = [System.IO.File]::OpenRead($compressedValidation)
try {
    $gzipStream = [System.IO.Compression.GZipStream]::new(
        $inputStream,
        [System.IO.Compression.CompressionMode]::Decompress
    )
    try {
        $outputStream = [System.IO.File]::Create($validationPdf)
        try { $gzipStream.CopyTo($outputStream) } finally { $outputStream.Dispose() }
    }
    finally { $gzipStream.Dispose() }
}
finally { $inputStream.Dispose() }

$checksumPath = Join-Path $resolvedDestination 'SHA256SUMS.txt'
$checksumTargets = @($downloads | ForEach-Object { Join-Path $resolvedDestination $_.Name }) + @($validationPdf)
$lines = foreach ($target in $checksumTargets) {
    $hash = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant()
    "$hash  $([System.IO.Path]::GetFileName($target))"
}
Set-Content -LiteralPath $checksumPath -Value $lines -Encoding utf8
Write-Output "Downloaded and checksummed reference inputs in $resolvedDestination"

