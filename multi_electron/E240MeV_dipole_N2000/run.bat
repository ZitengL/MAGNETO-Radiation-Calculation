@echo off
rem Run the 2000-electron GPT simulation (240 MeV, short dipole) and build stats.gdf. Works from any directory.

if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)

cd /d "%~dp0"

"%GPTBIN%\gpt.exe" -o result.gdf multi_electron_dipole.in GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error

"%GPTBIN%\gdfa.exe" -o stats.gdf result.gdf time numpar avgx avgy avgz stdx stdy stdz avgBx avgBy avgBz avgG stdG avgfBx avgfBy avgfBz
if errorlevel 1 goto :error

echo Done: result.gdf and stats.gdf written.
goto :eof

:error
echo GPT/gdfa failed with errorlevel %errorlevel%.
exit /b 1