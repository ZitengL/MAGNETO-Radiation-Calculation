@echo off
rem Example run at D = 0.5 m: combined, dipole only, undulator only. The D scans are run from
rem matlab/interference_und_dipole/run_scan.m with the same command lines (D=... on the command line).

if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)
set D=0.5

cd /d "%~dp0"

"%GPTBIN%\gpt.exe" -o result_combined_D0p5.gdf combined.in D=%D% GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error
"%GPTBIN%\gpt.exe" -o result_dipole_only_D0p5.gdf dipole_only.in D=%D% GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error
"%GPTBIN%\gpt.exe" -o result_undulator_only.gdf undulator_only.in D=%D% GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error

echo Done.
goto :eof

:error
echo GPT failed with errorlevel %errorlevel%.
exit /b 1