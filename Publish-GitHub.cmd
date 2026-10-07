@echo off
setlocal
set "ATLAS_PUBLIC_ROOT=%~dp0"
powershell.exe -NoLogo -NoProfile -Command "& ([scriptblock]::Create([System.IO.File]::ReadAllText((Join-Path $env:ATLAS_PUBLIC_ROOT 'tools\publish_github.ps1'))))"
if errorlevel 1 (
  echo GitHub publication did not complete.
  pause
  exit /b 1
)
pause
