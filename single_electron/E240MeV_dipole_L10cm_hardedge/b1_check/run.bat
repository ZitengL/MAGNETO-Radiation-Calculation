@echo off
rem Hard-edge check: copy of ..\single_electron_dipole_hardedge.in with b1 = 1e6 1/m instead of 1e5 1/m.
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)

cd /d "%~dp0"

"%GPTBIN%\gpt.exe" -o result_b1e6.gdf single_electron_dipole_hardedge_b1e6.in GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error
echo Done: result_b1e6.gdf written.
goto :eof

:error
echo GPT failed with errorlevel %errorlevel%.
exit /b 1
