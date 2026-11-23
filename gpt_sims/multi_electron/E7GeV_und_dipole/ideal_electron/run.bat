@echo off
rem Ideal electron through undulator + dipole: undulator orientation und_s = +1 / -1, dz = 60 um and 30 um (tout check).
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)
cd /d "%~dp0"
"%GPTBIN%\gpt.exe" -o result_sp_dz60.gdf single_electron.in dz=60e-6 und_s=1 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 (echo GPT failed & exit /b 1)
"%GPTBIN%\gpt.exe" -o result_sm_dz60.gdf single_electron.in dz=60e-6 und_s=-1 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 (echo GPT failed & exit /b 1)
"%GPTBIN%\gpt.exe" -o result_sm_dz30.gdf single_electron.in dz=30e-6 und_s=-1 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 (echo GPT failed & exit /b 1)
echo Done.
