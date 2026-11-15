# apply-backup.ps1 - unpack a backup made on admin.html ("Download full backup").
#
#   site/  -> data/ of this repository (site edits baked into the txt/json files,
#             hidden poems and texts left out). Review with git diff, then commit.
#   full/  -> _private/backup/latest and _private/backup/<archive name>  (gitignored)
#
# Usage: double-click tools/apply-backup.cmd (takes the newest
# thorntales-backup-*.zip from Downloads), or
#   powershell -ExecutionPolicy Bypass -File tools/apply-backup.ps1 [path\to\backup.zip] [-BackupOnly]
param(
  [string]$Zip,
  [switch]$BackupOnly
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

if (-not $Zip) {
  $downloads = Join-Path $env:USERPROFILE 'Downloads'
  $found = Get-ChildItem -Path $downloads -Filter 'thorntales-backup-*.zip' -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending | Select-Object -First 1
  if (-not $found) { throw "No thorntales-backup-*.zip in $downloads. Download one on /admin.html first." }
  $Zip = $found.FullName
}
$Zip = (Resolve-Path $Zip).Path
$name = [IO.Path]::GetFileNameWithoutExtension($Zip)
Write-Host "Backup: $Zip"

Add-Type -AssemblyName System.IO.Compression.FileSystem
$tmp = Join-Path ([IO.Path]::GetTempPath()) ("thorntales-" + [guid]::NewGuid())
[IO.Compression.ZipFile]::ExtractToDirectory($Zip, $tmp)

try {
  foreach ($part in 'site', 'full') {
    if (-not (Test-Path (Join-Path $tmp $part))) { throw "The archive has no '$part' folder: is it a Thorntales backup?" }
  }

  # 1. Private full copy: a dated folder and "latest"
  $private = Join-Path $repo '_private\backup'
  New-Item -ItemType Directory -Force $private | Out-Null
  $dated = Join-Path $private $name
  $latest = Join-Path $private 'latest'
  foreach ($dest in $dated, $latest) {
    if (Test-Path $dest) { Remove-Item -Recurse -Force $dest }
    Copy-Item -Recurse (Join-Path $tmp 'full') $dest
  }
  Copy-Item $Zip (Join-Path $private ($name + '.zip')) -Force
  Write-Host "Full backup saved to _private\backup\$name and _private\backup\latest"

  # 2. Public files over data/ (only files that are in the archive are replaced)
  if (-not $BackupOnly) {
    $src = Join-Path $tmp 'site\data'
    $count = 0
    Get-ChildItem -Recurse -File $src | ForEach-Object {
      $rel = $_.FullName.Substring($src.Length + 1)
      $dest = Join-Path (Join-Path $repo 'data') $rel
      New-Item -ItemType Directory -Force (Split-Path -Parent $dest) | Out-Null
      Copy-Item $_.FullName $dest -Force
      $count++
    }
    Write-Host "$count files copied into data\"

    # A poem that became hidden must leave the public files: a translation file that has
    # no public poems any more is not in site/, so it is removed here
    Get-ChildItem (Join-Path $repo 'data\poetry') -Filter '*.txt' | ForEach-Object {
      if (-not (Test-Path (Join-Path $src ('poetry\' + $_.Name)))) {
        if (Test-Path (Join-Path $tmp ('full\data\poetry\' + $_.Name))) {
          Remove-Item $_.FullName
          Write-Host "removed data\poetry\$($_.Name) (only hidden poems in it)"
        }
      }
    }

    Push-Location $repo
    try { git status --short -- data } finally { Pop-Location }
    Write-Host ''
    Write-Host 'Review the changes (git diff), then commit and push to publish.'
  }
} finally {
  Remove-Item -Recurse -Force $tmp
}
