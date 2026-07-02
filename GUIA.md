# Guia do EBSD Analyzer (em Português)

Guia rápido para qualquer pessoa abrir, usar e publicar o **EBSD Analyzer** sem
precisar de suporte técnico.

---

## 1. O que é

O EBSD Analyzer é um aplicativo web (feito em [Streamlit](https://streamlit.io))
para analisar dados de EBSD (Electron Backscatter Diffraction). Você envia um
arquivo exportado do microscópio e o app gera automaticamente:

- Distribuição de **tamanho de grão** (histograma log-normal, CDF, teste de normalidade)
- Histograma de **desorientação** (LAGB/HAGB, curva de Mackenzie)
- Resumo de **textura** (ângulos de Euler, frações de orientações ideais)
- Detecção de **outliers** (IQR / Z-score)
- **KAM** e estimativa de densidade de discordâncias (GND)
- Figuras prontas para publicação (PNG/SVG/PDF, 300 dpi)
- Visualizações **IPF** (mapa, triângulo, cubo 3D) e **figuras de polo (PF)**

---

## 2. Link público

> O link público é gerado ao publicar no **Streamlit Community Cloud** (gratuito).
> Veja a Seção 6. Depois de publicar, anote o endereço aqui, por exemplo:
>
> **App online:** `https://<seu-usuario>-ebsd-analyzer.streamlit.app`
>
> Quem receber esse link só precisa abri-lo no navegador — não instala nada.

---

## 3. Arquivos necessários (o que precisa estar na pasta)

Para o app rodar (local ou na nuvem), estes arquivos são **obrigatórios** e já
estão no repositório:

| Arquivo | Função |
|---|---|
| `app.py` | Aplicativo principal (ponto de entrada) |
| `file_readers.py` | Leitura de arquivos CTF / BCF / CSV |
| `ctf_processing.py` | Segmentação de grãos e cálculos (KAM, estatísticas) |
| `ipf_plots.py` | Gráficos IPF (densidade, triângulo, cubo 3D) |
| `pf_plots.py` | Figuras de polo e mapa IPF 2D |
| `excel_reference.py` | Leitura opcional de planilha de referência (.xlsx/.xlsm) |
| `requirements.txt` | Lista de dependências Python |

Pasta de exemplos (opcional, mas recomendada para testar):

| Arquivo | Função |
|---|---|
| `sample_data/exemplo_graos.csv` | Exemplo por grão (CSV) — uso imediato |
| `sample_data/exemplo_pixels.ctf` | Exemplo por pixel (CTF) — testa a segmentação |
| `generate_sample.py` | Script que regenera os exemplos acima |

> **Ponto de entrada (entrypoint):** `app.py`. É o arquivo que o Streamlit executa.

---

## 4. Exemplo mínimo para teste

> ⚠️ **Os arquivos de exemplo são SINTÉTICOS** (gerados por computador pelo
> `generate_sample.py`). Servem apenas para demonstrar o app. **Não são dados
> experimentais reais** e não devem ser usados como resultado de pesquisa.

Passo a passo:

1. Abra o app (online pelo link, ou localmente — veja a Seção 5).
2. Clique em **Browse files** (Procurar arquivos) na barra lateral.
3. Selecione `sample_data/exemplo_graos.csv`.
4. As abas de análise (tamanho de grão, desorientação, textura, outliers) são
   preenchidas automaticamente.
5. Para testar a **segmentação automática de grãos** a partir de pixels, envie
   `sample_data/exemplo_pixels.ctf` no lugar.

Para regenerar os exemplos (opcional):

```bash
python generate_sample.py
```

---

## 5. Como rodar localmente

### Windows (mais simples)

1. Instale o Python em <https://www.python.org/downloads/> e marque
   **"Add Python to PATH"**.
2. Dê duplo-clique em **`install.bat`** (cria o ambiente e instala tudo).
3. Dê duplo-clique em **`run_app.bat`** (abre o navegador em `http://localhost:8501`).

### macOS / Linux (terminal)

```bash
cd ebsd-analyzer
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

O navegador abre em `http://localhost:8501`.

---

## 6. Como publicar online (link público gratuito)

A forma recomendada de hospedar é o **Streamlit Community Cloud**, que conecta
direto ao repositório do GitHub:

1. Garanta que o código está no GitHub: <https://github.com/bono30/ebsd-analyzer>.
2. Acesse <https://share.streamlit.io> e faça login com a conta do GitHub.
3. Clique em **New app** e escolha:
   - **Repository:** `bono30/ebsd-analyzer`
   - **Branch:** `main`
   - **Main file path:** `app.py`
4. Clique em **Deploy**. Em poucos minutos o app fica online em um link
   `https://...streamlit.app` que pode ser compartilhado com qualquer pessoa.

> Não é preciso configurar servidor: o Streamlit Cloud lê o `requirements.txt`
> e instala as dependências sozinho.

> ⚠️ **Privacidade (instância pública):** em uma demonstração pública, os
> arquivos enviados são processados no servidor (em memória) e não são
> armazenados permanentemente, mas **não envie dados confidenciais ou ainda não
> publicados** para uma instância pública. Para dados sigilosos, rode o app
> localmente (Seção 5).

---

## 7. Guia curto: atualizar a pasta depois (manual)

Sempre que você mudar arquivos (código, exemplos, README) e quiser que o app
online reflita a mudança, atualize o GitHub. O Streamlit Cloud **redeploy
automaticamente** ao detectar um novo commit no branch `main`.

```bash
# dentro da pasta do projeto
git add .
git commit -m "Descreva a mudanca aqui"
git push origin main
```

Se for a primeira vez configurando o repositório remoto:

```bash
git remote add origin https://github.com/bono30/ebsd-analyzer.git
git branch -M main
git push -u origin main
```

Checklist antes do push:

- [ ] O app abre localmente com `streamlit run app.py`?
- [ ] O exemplo `sample_data/exemplo_graos.csv` carrega sem erro?
- [ ] `requirements.txt` continua com as dependências corretas?

---

## 8. Formatos de arquivo aceitos

| Formato | Extensão | Origem |
|---|---|---|
| CTF | `.ctf` | Oxford Instruments / HKL Channel 5 (dados por pixel) |
| BCF | `.bcf` | Bruker Esprit (se a leitura falhar, exporte como CSV) |
| CSV / TXT | `.csv`, `.txt` | OIM, AztecCrystal, MTEX, etc. (dados por grão) |
| Excel (referência) | `.xlsx`, `.xlsm` | Exportação AztecCrystal / ESPRIT — **opcional**, só para calibração/conferência (ver Seção 9) |

Se os nomes das colunas forem diferentes, use o painel **Column mapping**
dentro do app para indicar manualmente quais colunas usar.

---

## 9. Planilha de referência Excel (.xlsx/.xlsm) — opcional

O app aceita, **de forma opcional**, uma planilha de *exportação/pós-processamento*
de EBSD gerada por softwares comerciais (por exemplo, **Oxford AztecCrystal** ou
**Bruker ESPRIT**). Não é o mapa EBSD bruto — é uma planilha com resumos já
calculados: aba `Overview`, `Grain List`, `Boundary Statistics`, histogramas de
grão, `PolePlotData...` e `MackenziePlot...`.

### Para que serve

Essa planilha **não substitui** a análise do app. Ela é usada para
**calibração e conferência**:

- **Metadados de aquisição** (`Overview`): *step size* (µm), *pixel count*,
  *raster*, *hit rate*, *zero solution count* e frações de fase.
- **Calibração do step size:** se o arquivo EBSD enviado **não** tiver o passo
  (step size), o app usa o valor da planilha para a segmentação de grãos. Se
  ambos existirem, o app mostra os dois lado a lado para comparação.
- **Resumo da `Grain List`:** nº de grãos, contagem por fase e estatísticas de
  Área, ECD, Max Feret, *Mean/Maximum Orientation Spread*.
- **`Boundary Statistics`:** comprimento e fração de LAGB (2–10°) e HAGB (>10°).
- **Textura (`PolePlotData`):** pico e média de m.u.d. e o ângulo do pico.
- **`MackenziePlot`:** desorientação média medida e fração de LAGB (<15°).

### Como usar (teste com uma planilha real)

1. Abra o app (`streamlit run app.py`).
2. Na barra lateral, em **Reference workbook (optional)**, clique em
   **Browse files** e selecione sua planilha `.xlsx` ou `.xlsm`.
3. Um painel **“Reference workbook summary (Excel)”** aparece no topo com as
   tabelas de metadados, grãos, contornos, textura e Mackenzie.
4. Envie também o arquivo EBSD (CTF/CSV) correspondente na parte de cima da
   barra lateral. Na aba **KAM / Band Contrast**, o painel de GND passa a mostrar
   o comparativo de step size e o diagnóstico de discrepância.

### Densidade de discordâncias (GND) — honestidade científica

O app calcula a densidade de GND **a partir do KAM** do mapa EBSD, com a fórmula

```
ρ_GND = 2·θ_KAM / (α · u · b)
```

onde θ_KAM está em **radianos**, `u` = step size em **metros** (µm × 1e-6),
`b` = vetor de Burgers em **metros** (nm × 1e-9) e `α` = fator de método.

- **α = 1** é a forma padrão (Kubin & Mortensen 2003; Calcagnotto et al. 2010).
- Versões anteriores deste app usavam **α ≈ 1,86**, o que **reduzia ρ nesse fator**
  — provável causa de valores de GND discrepantes (mais baixos) em relação a
  outras ferramentas. Agora `α` é um campo ajustável e o app mostra os dois
  valores no diagnóstico.
- O app **não** substitui seu ρ_GND por um valor da planilha, porque essas
  planilhas de exportação normalmente **não trazem** um valor de GND. O *Mean
  Orientation Spread* da `Grain List` é exibido apenas como **referência de
  espalhamento intragranular** — ele **não** é igual a KAM nem a GND.

### Limitações

- É uma estimativa de **limite inferior** de GND (só as GNDs resolvidas pelo
  KAM naquele step; ignora SSDs).
- Nomes de abas truncados pelo Excel são reconhecidos por prefixo; abas muito
  fora do padrão dos exportadores podem não ser lidas (o app apenas ignora as
  seções ausentes, sem quebrar).
- A leitura é **somente de valores** (não executa macros VBA do `.xlsm`).
- Caminhos absolutos privados presentes na planilha (ex.: `C:\Users\...`) **não**
  são exibidos — o app mostra apenas o nome do arquivo de origem.

---

## 10. Problemas comuns

| Problema | Solução |
|---|---|
| `Python not found` | Reinstale o Python marcando "Add to PATH" |
| Navegador não abre | Acesse `http://localhost:8501` manualmente |
| Colunas não detectadas | Use o painel **Column mapping** dentro do app |
| CSV europeu (vírgula decimal) | Ajuste "Decimal separator" para `,` na barra lateral |
| App online não atualiza | Confirme que o `git push` foi para o branch `main` |
