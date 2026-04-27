#Requires -Version 5.0
<#
.SYNOPSIS
    Setup script for the Integrative miRNA Annotation Pipeline (Windows PowerShell)
.DESCRIPTION
    Checks for Python, installs it if missing, creates a virtual environment,
    and installs all required dependencies.
#>

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "  Integrative miRNA Annotation Pipeline - Setup (Windows)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""

# ── Helper ────────────────────────────────────────────────────
function Find-Python {
    foreach ($cmd in @("python", "python3", "py")) {
        try {
            $ver = & $cmd --version 2>&1
            if ($LASTEXITCODE -eq 0) { return $cmd }
        } catch {}
    }
    return $null
}

# ── 1. Check Python ───────────────────────────────────────────
Write-Host "[1/4] Checking Python installation..." -ForegroundColor Yellow

$pythonCmd = Find-Python

if ($null -eq $pythonCmd) {
    Write-Host "      Python not found. Downloading Python 3.11.9..." -ForegroundColor Yellow

    $installerPath = "$env:TEMP\python_installer.exe"
    $pythonUrl     = "https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe"

    try {
        [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
        Invoke-WebRequest -Uri $pythonUrl -OutFile $installerPath -UseBasicParsing
    } catch {
        Write-Host ""
        Write-Host "[ERROR] Could not download Python automatically." -ForegroundColor Red
        Write-Host "        Please install Python 3.11+ manually from:"
        Write-Host "        https://www.python.org/downloads/"
        Write-Host "        Make sure to check 'Add Python to PATH' during install."
        Read-Host "Press Enter to exit"
        exit 1
    }

    Write-Host "      Installing Python 3.11.9 (this may take a minute)..." -ForegroundColor Yellow
    $proc = Start-Process -FilePath $installerPath `
        -ArgumentList "/quiet", "InstallAllUsers=0", "PrependPath=1", "Include_test=0" `
        -Wait -PassThru
    Remove-Item $installerPath -Force -ErrorAction SilentlyContinue

    if ($proc.ExitCode -ne 0) {
        Write-Host "[ERROR] Python installation failed (exit code $($proc.ExitCode))." -ForegroundColor Red
        Write-Host "        Please install manually from https://www.python.org/downloads/"
        Read-Host "Press Enter to exit"
        exit 1
    }

    # Refresh PATH for current session
    $env:PATH = [System.Environment]::GetEnvironmentVariable("PATH","Machine") + ";" +
                [System.Environment]::GetEnvironmentVariable("PATH","User")

    $pythonCmd = Find-Python
    if ($null -eq $pythonCmd) {
        Write-Host ""
        Write-Host "[IMPORTANT] Python was installed but PATH is not updated yet." -ForegroundColor Yellow
        Write-Host "            Please CLOSE this window, open a new PowerShell, and run setup.ps1 again."
        Read-Host "Press Enter to exit"
        exit 0
    }

    Write-Host "      Python installed successfully." -ForegroundColor Green
} else {
    $ver = & $pythonCmd --version 2>&1
    Write-Host "      Found: $ver" -ForegroundColor Green
}

# ── 2. Check / upgrade pip ────────────────────────────────────
Write-Host ""
Write-Host "[2/4] Checking pip..." -ForegroundColor Yellow
& $pythonCmd -m pip --version | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "      pip not found, installing..." -ForegroundColor Yellow
    & $pythonCmd -m ensurepip --upgrade
}
& $pythonCmd -m pip install --upgrade pip --quiet
Write-Host "      pip OK." -ForegroundColor Green

# ── 3. Create virtual environment ─────────────────────────────
Write-Host ""
Write-Host "[3/4] Creating virtual environment (miRNA_env)..." -ForegroundColor Yellow

if (Test-Path "miRNA_env") {
    Write-Host "      miRNA_env already exists, skipping creation." -ForegroundColor Green
} else {
    & $pythonCmd -m venv miRNA_env
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERROR] Could not create virtual environment." -ForegroundColor Red
        Read-Host "Press Enter to exit"
        exit 1
    }
    Write-Host "      Virtual environment created." -ForegroundColor Green
}

# ── 4. Install dependencies ───────────────────────────────────
Write-Host ""
Write-Host "[4/4] Installing Python dependencies..." -ForegroundColor Yellow

& .\miRNA_env\Scripts\Activate.ps1

if (-not (Test-Path "requirements.txt")) {
    Write-Host "[ERROR] requirements.txt not found." -ForegroundColor Red
    Write-Host "        Make sure you are running this script from the project root folder."
    Read-Host "Press Enter to exit"
    exit 1
}

pip install -r requirements.txt --quiet
if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Dependency installation failed." -ForegroundColor Red
    Write-Host "        Check your internet connection and try again."
    Read-Host "Press Enter to exit"
    exit 1
}

# ── Done ──────────────────────────────────────────────────────
Write-Host ""
Write-Host "============================================================" -ForegroundColor Green
Write-Host "  Setup complete!" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Green
Write-Host ""
Write-Host "  To activate the environment in the future:"
Write-Host "    .\miRNA_env\Scripts\Activate.ps1" -ForegroundColor Cyan
Write-Host ""
Write-Host "  To run the pipeline:"
Write-Host "    python scripts\mirna_super_pipeline.py" -ForegroundColor Cyan
Write-Host ""
Read-Host "Press Enter to exit"

