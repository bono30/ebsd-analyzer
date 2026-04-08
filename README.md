# EBSD Analyzer

A local Streamlit app for EBSD (Electron Backscatter Diffraction) data analysis. Drop your CSV exports and get grain-size distributions, misorientation histograms, texture summaries, outlier detection, KAM-based GND density estimates, and publication-ready plots — all in English.

---

## Features

| Module | What it does |
|---|---|
| **Grain Size** | Log-normal histogram + CDF, Shapiro-Wilk normality test, grain size by phase |
| **Misorientation** | LAGB/HAGB frequency bar chart, Mackenzie random-texture reference curve, fraction 15°–65° |
| **Texture** | Euler angle distributions, ideal orientation fractions (Cube, Goss, Brass, Gamma fiber), Φ vs φ₂ ODF section |
| **Outliers** | IQR / Z-score / Modified Z-score detection with boxplots and outlier row table |
| **KAM / IQ** | KAM distribution, Image Quality histogram, GND density estimate (Kubin–Mortensen method) |
| **Export** | All figures as ZIP (PNG/SVG/PDF at 300 dpi) + individual downloads + stats CSV |

---

## Requirements

- **Windows 10/11** (also works on macOS/Linux)
- **Python 3.9 or newer**

---

## Quick Start (Windows)

### Step 1 — Install Python (if not installed)

Download from [python.org/downloads](https://www.python.org/downloads/).
During installation, check **"Add Python to PATH"**.

### Step 2 — Install the app

Double-click `install.bat`.
This creates a virtual environment and installs all dependencies automatically.

### Step 3 — Run the app

Double-click `run_app.bat`.
Your browser opens at `http://localhost:8501`.

---

## Supported File Formats

| Format | Extension | Source software | Notes |
|---|---|---|---|
| **CTF** | `.ctf` | Oxford Instruments / HKL Channel 5 | Direct upload — full header parsed automatically |
| **BCF** | `.bcf` | Bruker Esprit EBSD | Binary container auto-extracted |
| **CSV / TXT** | `.csv`, `.txt` | OIM, AztecCrystal, MTEX, any export | Configurable separator and decimal |

### CTF column mapping (auto-detected)

| CTF column | App column | Description |
|---|---|---|
| `Phase` | Phase | Phase index → name (from header) |
| `X`, `Y` | X, Y | Position in µm |
| `Euler1`, `Euler2`, `Euler3` | φ₁, Φ, φ₂ | Euler angles in degrees |
| `MAD` | MAD | Mean Angular Deviation (°) |
| `BC` | Band Contrast | 0–255 |
| `BS` | Band Slope | 0–255 |
| `Bands` | Bands | Number of detected bands |
| `Error` | Error / CI | Indexing quality |

### BCF — if extraction fails

Bruker BCF EBSD uses a proprietary binary container. If your BCF file fails to load:
1. Open the project in **Bruker ESPRIT**
2. Go to **File → Export → Text/CSV**
3. Select columns: X, Y, Euler1, Euler2, Euler3, BC, MAD, Phase, KAM
4. Upload the resulting `.csv` file

### CSV / other exports

| Column | Auto-detected names |
|---|---|
| Grain diameter (ECD) | `Grain Diameter`, `ECD`, `Diameter` |
| Misorientation | `Misorientation`, `Misorientation Angle`, `MAD` |
| Grain area | `Area`, `Grain Area` |
| KAM | `KAM`, `Kernel Average Misorientation` |
| Euler angles | `Euler1`, `phi1`, `Phi`, `phi2` |
| Phase | `Phase`, `Phase Name` |
| CI | `CI`, `Confidence Index`, `Error` |
| IQ | `IQ`, `Image Quality`, `Fit` |

If column names differ, use the **Column mapping** panel inside the app.

---

## Sample Datasets

Two sample datasets are included in `sample_data/`:

| File | Format | Description |
|---|---|---|
| `AISI444_EBSD_sample.csv` | CSV | 3,000 grains, full data including KAM |
| `AISI444_EBSD_sample.ctf` | CTF | 2,000 pixel points, Oxford CTF format |

Both simulate **AISI 444 ferritic stainless steel** deep drawing experiments with γ-fiber (111) texture.

| Property | Value |
|---|---|
| Total points | 3,000 |
| Mean grain diameter | ~25 µm |
| LAGB fraction | ~30% |
| HAGB fraction | ~70% |
| Dominant texture | γ-fiber (111) |
| Phase | 95% Ferrite (BCC) + 5% Austenite |
| Outliers injected | 25 artificial outliers |

To use it: open the app, click **Browse files**, and select `sample_data/AISI444_EBSD_sample.csv`.

---

## Sidebar Controls

| Control | Effect |
|---|---|
| Column separator | `,` `;` `\t` or space |
| Decimal separator | `.` or `,` (for European locales) |
| Outlier method | IQR (default), Z-score, Modified Z-score |
| Remove outliers | Toggle to exclude detected outliers from plots |
| Histogram bins | Control bin count for grain size and misorientation plots |
| CI threshold | Slider to filter out low-confidence EBSD points |
| Plot format | PNG (default), SVG (vector), or PDF |

---

## GND Density Formula

The **Kernel Average Misorientation (KAM)** method uses the Kubin–Mortensen (2003) relation:

$$\rho_{GND} = \frac{2\,\theta_{KAM}}{u \cdot b}$$

Where:
- θ_KAM = local misorientation in radians
- u = 1.86 × step size (m)
- b = Burgers vector (m): Fe BCC ≈ 0.248 nm, FCC ≈ 0.254 nm

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `Python not found` | Reinstall Python and tick "Add to PATH" |
| Browser does not open | Manually go to `http://localhost:8501` |
| Columns not detected | Use the Column mapping dropdowns inside the app |
| European CSV (`,` as decimal) | Set "Decimal separator" to `,` in the sidebar |
| Streamlit not found | Run `install.bat` first to create the virtual environment |

---

## File Structure

```
ebsd_analyzer/
├── app.py                    ← Main Streamlit application
├── requirements.txt          ← Python package list
├── install.bat               ← Windows installer (creates venv)
├── run_app.bat               ← Windows launcher (uses venv)
├── run_app_no_venv.bat       ← Launcher without venv (system Python)
├── generate_sample.py        ← Script to regenerate the sample dataset
├── sample_data/
│   └── AISI444_EBSD_sample.csv   ← Ready-to-use sample dataset
└── README.md                 ← This file
```

---

## Citation

If this tool assists your research, you may acknowledge it as:

> EBSD Analyzer (2025). Local Streamlit application for grain microstructure analysis. Available at: [your institution or repository].

---

## License

MIT License — free to use, modify, and distribute.
