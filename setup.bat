@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul

echo ============================================================
echo   Integrative miRNA Annotation Pipeline - Setup (Windows)
echo ============================================================
echo.

:: ── 1. Check Python ──────────────────────────────────────────
echo [1/4] Checking Python installation...

python --version >nul 2>&1
if %ERRORLEVEL% == 0 (
    for /f "tokens=*" %%i in ('python --version 2^>^&1') do set PY_VER=%%i
    echo       Found: !PY_VER!
    goto :check_pip
)

py --version >nul 2>&1
if %ERRORLEVEL% == 0 (
    for /f "tokens=*" %%i in ('py --version 2^>^&1') do set PY_VER=%%i
    echo       Found (py launcher): !PY_VER!
    set PYTHON_CMD=py
    goto :check_pip
)

echo       Python not found. Downloading installer...
echo.

:: Download Python 3.11.9 installer (stable, widely compatible)
set PY_INSTALLER=%TEMP%\python_installer.exe
set PY_URL=https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe

:: Try with PowerShell (available on all modern Windows)
powershell -Command "& {[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri '%PY_URL%' -OutFile '%PY_INSTALLER%'}" 2>nul
if not exist "%PY_INSTALLER%" (
    echo [ERROR] Could not download Python automatically.
    echo         Please install Python 3.11+ manually from:
    echo         https://www.python.org/downloads/
    echo         Make sure to check "Add Python to PATH" during install.
    pause
    exit /b 1
)

echo       Installing Python 3.11.9 (this may take a minute)...
"%PY_INSTALLER%" /quiet InstallAllUsers=0 PrependPath=1 Include_test=0
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Python installation failed.
    echo         Please install manually from https://www.python.org/downloads/
    pause
    exit /b 1
)
del "%PY_INSTALLER%" >nul 2>&1

:: Refresh PATH so we can find python immediately
call RefreshEnv.cmd >nul 2>&1
set "PATH=%LOCALAPPDATA%\Programs\Python\Python311;%LOCALAPPDATA%\Programs\Python\Python311\Scripts;%PATH%"

python --version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo.
    echo [IMPORTANT] Python was installed but the PATH was not updated yet.
    echo             Please CLOSE this window, open a new CMD, and run setup.bat again.
    pause
    exit /b 0
)
echo       Python installed successfully.

:check_pip
:: ── 2. Check / upgrade pip ───────────────────────────────────
echo.
echo [2/4] Checking pip...
if not defined PYTHON_CMD set PYTHON_CMD=python

%PYTHON_CMD% -m pip --version >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo       pip not found, installing...
    %PYTHON_CMD% -m ensurepip --upgrade
)
%PYTHON_CMD% -m pip install --upgrade pip --quiet
echo       pip OK.

:: ── 3. Create virtual environment ────────────────────────────
echo.
echo [3/4] Creating virtual environment (miRNA_env)...

if exist miRNA_env (
    echo       miRNA_env already exists, skipping creation.
) else (
    %PYTHON_CMD% -m venv miRNA_env
    if %ERRORLEVEL% neq 0 (
        echo [ERROR] Could not create virtual environment.
        pause
        exit /b 1
    )
    echo       Virtual environment created.
)

:: ── 4. Install dependencies ───────────────────────────────────
echo.
echo [4/4] Installing Python dependencies...

call miRNA_env\Scripts\activate.bat

if not exist requirements.txt (
    echo [ERROR] requirements.txt not found in current directory.
    echo         Make sure you are running this script from the project root folder.
    pause
    exit /b 1
)

pip install -r requirements.txt --quiet
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Dependency installation failed.
    echo         Check your internet connection and try again.
    pause
    exit /b 1
)

:: ── Done ─────────────────────────────────────────────────────
echo.
echo ============================================================
echo   Setup complete!
echo ============================================================
echo.
echo   To activate the environment in the future:
echo     miRNA_env\Scripts\activate.bat
echo.
echo   To run the pipeline:
echo     python scripts\mirna_super_pipeline.py
echo.
pause
