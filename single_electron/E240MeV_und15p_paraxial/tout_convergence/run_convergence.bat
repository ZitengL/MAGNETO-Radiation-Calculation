@echo off
rem tout-step convergence runs: copies of ..\single_electron_und15p.in with only dz changed (30, 15 um)
rem and the field-map path one level deeper.

if not defined GPTBIN (echo Set the environment variable GPTBIN to the GPT bin directory. & exit /b 1)
if not defined GPTLICENSE (echo Set the environment variable GPTLICENSE to your GPT license number. & exit /b 1)

cd /d "%~dp0"

for %%D in (30 15) do (
  "%GPTBIN%\gpt.exe" -o result_dz%%Dum.gdf single_electron_und15p_dz%%Dum.in GPTLICENSE=%GPTLICENSE%
  if errorlevel 1 goto :error
)

echo Done: result_dz30um.gdf and result_dz15um.gdf written.
goto :eof

:error
echo GPT failed with errorlevel %errorlevel%.
exit /b 1