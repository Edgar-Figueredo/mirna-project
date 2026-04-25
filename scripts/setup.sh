#!/bin/bash

echo "[SETUP] Creating virtual environment..."

python3 -m venv miRNA_env

echo "[SETUP] Activating environment..."
source miRNA_env/bin/activate

echo "[SETUP] Upgrading pip..."
pip install --upgrade pip

echo "[SETUP] Installing dependencies..."
pip install -r env/requirements.txt

echo "[DONE] Environment ready ✔"
echo " Run: source miRNA_env/bin/activate"
echo " Then: python scripts/mirna_super_pipeline.py"
