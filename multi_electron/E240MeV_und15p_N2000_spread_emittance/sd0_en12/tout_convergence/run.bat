@echo off
rem tout convergence copy (30 um) of ..\multi_electron.in. Run this case (2000 electrons, 15-period undulator, energy spread / emittance) and build stats.gdf.
rem initial_distribution.gdf is made by ..\make_distributions.m. Works from any directory.

if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)

cd /d "%~dp0"

"%GPTBIN%\gpt.exe" -o result_dz30um.gdf multi_electron_dz30um.in GPTLICENSE=%GPTLICENSE%
if errorlevel 1 goto :error

"%GPTBIN%\gdfa.exe" -o stats_dz30um.gdf result_dz30um.gdf time numpar avgx avgy avgz stdx stdy stdz avgBx avgBy avgBz avgG stdG avgfBx avgfBy avgfBz
if errorlevel 1 goto :error

echo Done: result_dz30um.gdf and stats_dz30um.gdf written.
goto :eof

:error
echo GPT/gdfa failed with errorlevel %errorlevel%.
exit /b 1