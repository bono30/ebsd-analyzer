@echo off
REM ── EBSD Analyzer — instalador para Windows ──────────────────────────────
REM Cria um ambiente virtual (.venv) e instala as dependencias.
cd /d "%~dp0"
echo Criando ambiente virtual...
python -m venv .venv
if errorlevel 1 (
    echo ERRO: Python nao encontrado. Instale em https://www.python.org/downloads/
    echo Marque a opcao "Add Python to PATH" durante a instalacao.
    pause
    exit /b 1
)
echo Instalando dependencias...
call .venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r requirements.txt
echo.
echo Pronto! Agora execute run_app.bat para abrir o aplicativo.
pause
