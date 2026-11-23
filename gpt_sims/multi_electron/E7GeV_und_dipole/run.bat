@echo off
rem One chunk of 2000 macro-particles (chunk_beam.gdf) -> chunk_result.gdf. Argument: output step in um (default 120).
rem Called by run_combined.m; the result is deleted after the MAGNETO step. Works from any directory.
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)
cd /d "%~dp0"
set DZ=%1
if "%DZ%"=="" set DZ=120
"%GPTBIN%\gpt.exe" -o chunk_result.gdf multi_electron.in dz=%DZ%e-6 GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error
goto :eof
:error
echo GPT failed with errorlevel %errorlevel%.
exit /b 1
