@echo off
rem Double-click to start the Lab Proteomics Workbench (Windows).
cd /d "%~dp0"
set "RS="
where Rscript >nul 2>&1 && set "RS=Rscript"
if "%RS%"=="" (
  for /d %%d in ("%ProgramFiles%\R\R-*") do set "RS=%%d\bin\Rscript.exe"
)
if "%RS%"=="" (
  echo Rscript was not found. Install R from https://cran.r-project.org and try again.
  pause
  exit /b 1
)
"%RS%" scripts\launch_app.R
pause
