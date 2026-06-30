@echo off
REM ── EBSD Analyzer — abrir o app usando o ambiente virtual (.venv) ─────────
cd /d "%~dp0"
if not exist ".venv\Scripts\activate.bat" (
    echo Ambiente virtual nao encontrado. Execute install.bat primeiro.
    pause
    exit /b 1
)
call .venv\Scripts\activate.bat
streamlit run app.py
pause
