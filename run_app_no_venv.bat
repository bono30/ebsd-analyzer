@echo off
REM ── EBSD Analyzer — abrir o app usando o Python do sistema (sem .venv) ────
cd /d "%~dp0"
python -m pip install -r requirements.txt
streamlit run app.py
pause
