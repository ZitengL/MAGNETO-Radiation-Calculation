@echo off
rem Convergence copy of ..\combined.in (output steps 30 um / 1.875 um instead of 60 um / 3.75 um), D = 0.5 m.
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)
cd /d "%~dp0"
"%GPTBIN%\gpt.exe" -o result_combined_fine_D0p5.gdf combined_fine.in D=0.5 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 exit /b 1
echo Done.
