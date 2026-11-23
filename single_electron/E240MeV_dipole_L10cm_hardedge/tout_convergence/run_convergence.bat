@echo off
rem tout-step convergence runs: copies of ..\single_electron_dipole_hardedge.in with only dz changed (7.5, 1.875 um).
if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)

cd /d "%~dp0"

for %%D in (7p5 1p875) do (
  "%GPTBIN%\gpt.exe" -o result_dz%%Dum.gdf single_electron_dipole_hardedge_dz%%Dum.in GPTLICENSE=%GPTLICENSE%
  if errorlevel 1 goto :error
)

echo Done: result_dz7p5um.gdf and result_dz1p875um.gdf written.
goto :eof

:error
echo GPT failed with errorlevel %errorlevel%.
exit /b 1
