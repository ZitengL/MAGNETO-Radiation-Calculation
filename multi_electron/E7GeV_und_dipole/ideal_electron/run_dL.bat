@echo off
rem Ideal electron, undulator turned over, dipole moved downstream by 0.3 m (check of the translation in dipole_distance_scan.m).
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)
cd /d "%~dp0"
"%GPTBIN%\gpt.exe" -o result_sm_dL0p3_dz60.gdf single_electron_dL.in dz=60e-6 und_s=-1 dL=0.3 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 (echo GPT failed & exit /b 1)
echo Done.
