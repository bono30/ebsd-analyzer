#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_public.sh — executa o EBSD Analyzer atrás de um proxy em /port/5000.
#
# Usado pelo caminho de publicação (Perplexity publish_website), em que a página
# estática public/index.html aponta para /port/5000. Os flags de CORS/XSRF são
# desativados porque o app é servido atrás de um proxy reverso.
#
# NÃO afeta o Streamlit Community Cloud: lá o deploy usa o próprio comando
# (streamlit run app.py) e este script simplesmente não é chamado.
#
# Uso:  bash run_public.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "$0")"

exec streamlit run app.py \
  --server.address=0.0.0.0 \
  --server.port=5000 \
  --server.headless=true \
  --server.enableCORS=false \
  --server.enableXsrfProtection=false
