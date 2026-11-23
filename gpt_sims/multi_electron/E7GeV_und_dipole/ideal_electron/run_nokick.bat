@echo off
rem Ideal electron, undulator turned over (und_s = -1), launched with x' = -4.547e-7 to cancel the undulator's residual kick.
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)
cd /d "%~dp0"
"%GPTBIN%\gpt.exe" -o result_sm_nokick_dz60.gdf single_electron_nokick.in dz=60e-6 und_s=-1 xp0=-4.547e-7 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 (echo GPT failed & exit /b 1)
echo Done.
