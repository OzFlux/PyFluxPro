@echo off
setlocal EnableDelayedExpansion

:: Name of the environment to activate
set "ENV_NAME=pfp_env"

:: Directory where this batch script lives
set "SCRIPT_DIR=%~dp0"
set "APP_ENTRY=%SCRIPT_DIR%PyFluxPro.py"

:: -------------------------------------------------------------
:: 1. Check for micromamba in PATH or standard user directories
:: -------------------------------------------------------------
where micromamba >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo Using micromamba ^(environment: %ENV_NAME%^)...
    for /f "delims=" %%i in ('micromamba shell hook -s cmd.exe 2^>nul') do %%i
    call micromamba activate "%ENV_NAME%"
    python "%APP_ENTRY%" %*
    goto :end
)

if exist "%USERPROFILE%\.local\bin\micromamba.exe" (
    echo Using micromamba from %USERPROFILE%\.local\bin...
    set "PATH=%USERPROFILE%\.local\bin;%PATH%"
    for /f "delims=" %%i in ('"%USERPROFILE%\.local\bin\micromamba.exe" shell hook -s cmd.exe 2^>nul') do %%i
    call micromamba activate "%ENV_NAME%"
    python "%APP_ENTRY%" %*
    goto :end
)

if exist "%LOCALAPPDATA%\micromamba\micromamba.exe" (
    echo Using micromamba from %LOCALAPPDATA%\micromamba...
    set "PATH=%LOCALAPPDATA%\micromamba;%PATH%"
    for /f "delims=" %%i in ('"%LOCALAPPDATA%\micromamba\micromamba.exe" shell hook -s cmd.exe 2^>nul') do %%i
    call micromamba activate "%ENV_NAME%"
    python "%APP_ENTRY%" %*
    goto :end
)

:: -------------------------------------------------------------
:: 2. Check for mamba in PATH
:: -------------------------------------------------------------
where mamba >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo Using mamba ^(environment: %ENV_NAME%^)...
    for /f "delims=" %%i in ('mamba shell hook 2^>nul') do %%i
    call mamba activate "%ENV_NAME%"
    python "%APP_ENTRY%" %*
    goto :end
)

:: -------------------------------------------------------------
:: 3. Check for conda in PATH or standard install paths
:: -------------------------------------------------------------
where conda >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo Using conda ^(environment: %ENV_NAME%^)...
    for /f "delims=" %%i in ('conda shell.cmd.exe hook 2^>nul') do %%i
    call conda activate "%ENV_NAME%"
    python "%APP_ENTRY%" %*
    goto :end
)

:: Fallback check common Conda / Miniforge install directories on Windows
for %%P in (
    "%USERPROFILE%\miniforge3"
    "%USERPROFILE%\miniconda3"
    "%USERPROFILE%\anaconda3"
    "%LOCALAPPDATA%\miniforge3"
    "%LOCALAPPDATA%\miniconda3"
    "%ProgramData%\miniforge3"
    "%ProgramData%\miniconda3"
    "%ProgramData%\anaconda3"
) do (
    if exist "%%~P\condabin\conda.bat" (
        echo Using conda from %%~P...
        call "%%~P\condabin\conda.bat" activate "%ENV_NAME%"
        python "%APP_ENTRY%" %*
        goto :end
    )
)

:: -------------------------------------------------------------
:: 4. Fallback: Neither found
:: -------------------------------------------------------------
echo Error: Neither micromamba nor conda was found on this system. >&2
echo Please ensure micromamba or conda is installed and available in PATH. >&2
exit /b 1

:end
endlocal