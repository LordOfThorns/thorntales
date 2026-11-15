@echo off
rem Unpack the newest thorntales-backup-*.zip from Downloads: see apply-backup.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0apply-backup.ps1" %*
pause
