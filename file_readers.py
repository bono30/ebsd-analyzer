"""
EBSD file readers for CTF (Oxford/HKL) and BCF (Bruker) formats.
No external EBSD libraries required — pure Python/pandas/numpy.
"""

import io
import re
import struct
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


def _coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert columns to numeric where possible, leaving non-numeric columns
    (e.g. phase names) untouched. Replaces the removed pandas
    `to_numeric(errors="ignore")` behaviour (gone in pandas 3.0).
    """
    for c in df.columns:
        conv = pd.to_numeric(df[c], errors="coerce")
        # keep the conversion only if it did not turn valid values into NaN
        if conv.notna().sum() == df[c].notna().sum():
            df[c] = conv
    return df


# ══════════════════════════════════════════════════════════════════════════════
#  CTF READER  (Oxford Instruments / HKL Channel 5)
# ══════════════════════════════════════════════════════════════════════════════
# Standard CTF column order (tab-separated after the header block):
#   Phase  X  Y  Bands  Error  Euler1  Euler2  Euler3  MAD  BC  BS
# Some exports add extra columns; we handle those too.

CTF_STANDARD_COLS = [
    "Phase", "X", "Y", "Bands", "Error",
    "Euler1 (phi1)", "Euler2 (Phi)", "Euler3 (phi2)",
    "MAD", "BC", "BS",
]

CTF_EXTRA_COLS = ["Euler1 (phi1)", "Euler2 (Phi)", "Euler3 (phi2)", "MAD", "BC", "BS"]


def read_ctf(file_bytes: bytes) -> tuple[pd.DataFrame, dict]:
    """
    Parse a CTF file (Oxford/HKL) and return (DataFrame, metadata_dict).

    The CTF format is a plain-text file with:
      - A header block (lines starting with keywords like 'Channel', 'Phases', etc.)
      - A 'Phases' section with crystallographic data
      - A column-header line: 'Phase\\tX\\tY\\tBands\\t...'
      - Data rows (tab-separated)
    """
    try:
        text = file_bytes.decode("utf-8", errors="replace")
    except Exception:
        text = file_bytes.decode("latin-1", errors="replace")

    lines = text.splitlines()
    metadata = {}
    phase_names = {}
    header_end = 0
    col_header_line = -1

    # ── Parse header ──────────────────────────────────────────────────────────
    in_phases = False
    phase_count = 0
    current_phase_idx = 0

    for i, line in enumerate(lines):
        stripped = line.strip()

        # Detect end of header / start of data
        if stripped.lower().startswith("phase\t") or stripped.lower().startswith("phase "):
            # This is the column header row
            col_header_line = i
            header_end = i + 1
            break

        # Key-value header lines
        if "\t" in stripped:
            parts = stripped.split("\t", 1)
            key = parts[0].strip()
            val = parts[1].strip() if len(parts) > 1 else ""
            metadata[key] = val

            if key == "Phases":
                try:
                    phase_count = int(val)
                except ValueError:
                    pass
        elif stripped.startswith("Phases"):
            try:
                phase_count = int(stripped.split()[-1])
            except ValueError:
                pass

        # Phase name lines follow semicolons in the header
        # Format: <a>;<b>;<c>;<alpha>;<beta>;<gamma>;SG;<name>;...
        if phase_count > 0 and ";" in stripped and not stripped.startswith("Channel"):
            # Heuristic: lines with semicolons after 'Phases' header are phase defs
            parts_sc = stripped.split(";")
            if len(parts_sc) >= 8:
                try:
                    float(parts_sc[0])  # first field should be a float (lattice param)
                    current_phase_idx += 1
                    pname = parts_sc[7].strip() if len(parts_sc) > 7 else f"Phase{current_phase_idx}"
                    phase_names[current_phase_idx] = pname
                except ValueError:
                    pass

    # ── Parse column headers ──────────────────────────────────────────────────
    if col_header_line == -1:
        # Fallback: find any line that starts with 'Phase'
        for i, line in enumerate(lines):
            if re.match(r"Phase[\t ]", line.strip(), re.IGNORECASE):
                col_header_line = i
                header_end = i + 1
                break

    if col_header_line == -1:
        raise ValueError(
            "Could not locate the data header row in the CTF file. "
            "Expected a line starting with 'Phase\\tX\\tY...'."
        )

    raw_cols = re.split(r"[\t ]+", lines[col_header_line].strip())

    # ── Rename columns to friendly names ─────────────────────────────────────
    col_rename = {
        "phase": "Phase",
        "x": "X",
        "y": "Y",
        "bands": "Bands",
        "error": "Error",
        "euler1": "Euler1 (phi1)",
        "euler2": "Euler2 (Phi)",
        "euler3": "Euler3 (phi2)",
        "mad": "MAD",
        "bc": "BC",
        "bs": "BS",
        "phi1": "Euler1 (phi1)",
        "phi":  "Euler2 (Phi)",
        "phi2": "Euler3 (phi2)",
    }
    renamed_cols = [col_rename.get(c.lower(), c) for c in raw_cols]

    # ── Read data block ───────────────────────────────────────────────────────
    data_lines = lines[header_end:]
    # Remove comment/blank lines
    data_lines = [l for l in data_lines if l.strip() and not l.strip().startswith(";")]

    data_text = "\n".join(data_lines)
    try:
        df = pd.read_csv(
            io.StringIO(data_text),
            sep=r"\s+",
            header=None,
            names=renamed_cols,
            engine="python",
            on_bad_lines="skip",
        )
    except Exception as e:
        raise ValueError(f"Failed to parse CTF data block: {e}")

    # Convert numerics
    for col in df.columns:
        if col != "Phase":
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Map phase indices to names if available
    if phase_names and "Phase" in df.columns:
        df["Phase"] = df["Phase"].map(
            lambda x: phase_names.get(int(x), f"Phase {int(x)}")
            if pd.notna(x) else x
        )

    # Compute grain diameter from step size if X/Y available
    if "X" in df.columns and "Y" in df.columns:
        xs = df["X"].dropna().unique()
        xs_sorted = np.sort(xs)
        if len(xs_sorted) > 1:
            step = float(np.median(np.diff(xs_sorted)))
            metadata["Step Size (µm)"] = step

    # Rename BC → Band Contrast, MAD → MAD (Mean Angular Deviation)
    df.rename(columns={"BC": "Band Contrast", "BS": "Band Slope"}, inplace=True)

    metadata["Source Format"] = "CTF (Oxford/HKL)"
    metadata["Rows"] = len(df)
    metadata["Columns"] = list(df.columns)

    return df, metadata


# ══════════════════════════════════════════════════════════════════════════════
#  BCF READER  (Bruker Esprit EBSD)
# ══════════════════════════════════════════════════════════════════════════════
# BCF is a ZIP-like container (SFS = AidAim Single File System).
# For EBSD data, Bruker stores it as a ZIP with:
#   - Header.xml  (scan parameters)
#   - EBSD/Area[n]/KAM.dat, Euler.dat, BC.dat, etc.  (binary float32 arrays)
# We try two strategies:
#   1. Treat as ZIP and extract binary EBSD data arrays
#   2. If not parseable as ZIP, fall back to embedded CSV/text data

BCF_EBSD_CHANNELS = {
    "phi1":            "Euler1 (phi1)",
    "phi":             "Euler2 (Phi)",
    "phi2":            "Euler3 (phi2)",
    "euler1":          "Euler1 (phi1)",
    "euler2":          "Euler2 (Phi)",
    "euler3":          "Euler3 (phi2)",
    "kam":             "KAM",
    "bc":              "Band Contrast",
    "bandcontrast":    "Band Contrast",
    "bs":              "Band Slope",
    "bandslope":       "Band Slope",
    "mad":             "MAD",
    "phase":           "Phase",
    "phaseid":         "Phase",
    "iq":              "Image Quality",
    "ci":              "Confidence Index",
    "fit":             "Fit",
    "grainid":         "Grain ID",
    "graindiameter":   "Grain Diameter",
    "grainsize":       "Grain Diameter",
    "misorientation":  "Misorientation Angle",
    "averagemisori":   "Misorientation Angle",
    "x":               "X",
    "y":               "Y",
}


def _try_zip_bcf(file_bytes: bytes) -> tuple[pd.DataFrame | None, dict]:
    """
    Attempt to open BCF as a ZIP archive and extract EBSD binary arrays.
    Returns (DataFrame or None, metadata_dict).
    """
    metadata = {}
    try:
        zf = zipfile.ZipFile(io.BytesIO(file_bytes))
    except zipfile.BadZipFile:
        return None, {}

    names = zf.namelist()

    # ── Parse XML header ──────────────────────────────────────────────────────
    xml_candidates = [n for n in names if n.lower().endswith(".xml")]
    width, height = None, None
    phase_map = {}

    for xname in xml_candidates:
        try:
            xml_content = zf.read(xname).decode("utf-8", errors="replace")
            root = ET.fromstring(xml_content)

            # Look for scan dimensions
            for tag in ["Width", "NColumns", "NCOLS", "MapWidth"]:
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    try:
                        width = int(el.text)
                    except ValueError:
                        pass
            for tag in ["Height", "NRows", "NROWS", "MapHeight"]:
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    try:
                        height = int(el.text)
                    except ValueError:
                        pass

            # Step size
            for tag in ["StepSize", "Resolution", "StepX"]:
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    try:
                        metadata["Step Size (µm)"] = float(el.text)
                    except ValueError:
                        pass

            # Phase names
            for phase_el in root.iter("Phase"):
                idx_el = phase_el.find("Index")
                name_el = phase_el.find("Name")
                if idx_el is not None and name_el is not None:
                    try:
                        phase_map[int(idx_el.text)] = name_el.text
                    except Exception:
                        pass
        except ET.ParseError:
            pass

    # ── Find and read binary data files ──────────────────────────────────────
    # Bruker stores each channel as a flat binary float32 or uint8 array
    data_files = [n for n in names if not n.lower().endswith(".xml")]
    columns_data = {}

    for fname in data_files:
        basename = fname.split("/")[-1].split("\\")[-1]
        stem = re.sub(r"\.(dat|bin|raw)$", "", basename, flags=re.IGNORECASE).lower()
        friendly = BCF_EBSD_CHANNELS.get(stem)
        if friendly is None:
            continue
        try:
            raw = zf.read(fname)
            # Try float32 first
            arr = np.frombuffer(raw, dtype=np.float32)
            if width and height and len(arr) == width * height:
                columns_data[friendly] = arr.ravel()
            elif len(arr) > 100:
                columns_data[friendly] = arr.ravel()
        except Exception:
            pass

    if not columns_data:
        return None, metadata

    # Align lengths
    min_len = min(len(v) for v in columns_data.values())
    df = pd.DataFrame({k: v[:min_len] for k, v in columns_data.items()})

    # Add X, Y grid if dimensions known and not already present
    if width and height and "X" not in df.columns and "Y" not in df.columns:
        step = metadata.get("Step Size (µm)", 1.0)
        xs = np.tile(np.arange(width) * step, height)[:min_len]
        ys = np.repeat(np.arange(height) * step, width)[:min_len]
        df.insert(0, "X", xs)
        df.insert(1, "Y", ys)

    # Map phase IDs to names
    if phase_map and "Phase" in df.columns:
        df["Phase"] = df["Phase"].map(
            lambda x: phase_map.get(int(round(x)), f"Phase {int(round(x))}")
            if pd.notna(x) else x
        )

    metadata["Source Format"] = "BCF (Bruker Esprit EBSD)"
    metadata["Rows"] = len(df)
    return df, metadata


def _try_embedded_csv_bcf(file_bytes: bytes) -> tuple[pd.DataFrame | None, dict]:
    """
    Some Bruker BCF files contain embedded CSV/text data. Scan for it.
    """
    # Look for a text block with tab/comma separated numbers
    try:
        text = file_bytes.decode("utf-8", errors="replace")
    except Exception:
        return None, {}

    # Find lines that look like data rows (mostly numeric)
    lines = text.splitlines()
    data_start = -1
    col_names = []

    for i, line in enumerate(lines):
        parts = line.strip().split("\t")
        if len(parts) >= 5:
            # Check if at least 4 of them are numeric
            numeric_count = sum(1 for p in parts if re.match(r"^-?\d+\.?\d*([eE][+-]?\d+)?$", p.strip()))
            if numeric_count >= 4:
                # Look back one line for headers
                if i > 0:
                    prev = lines[i - 1].strip().split("\t")
                    if any(re.match(r"^[A-Za-z]", p.strip()) for p in prev):
                        col_names = prev
                        data_start = i
                        break
                data_start = i
                break

    if data_start == -1:
        return None, {}

    data_text = "\n".join(lines[data_start:])
    try:
        df = pd.read_csv(
            io.StringIO(data_text),
            sep="\t",
            header=None,
            names=col_names if col_names else None,
            engine="python",
            on_bad_lines="skip",
        )
        df = _coerce_numeric(df)
        return df, {"Source Format": "BCF (embedded text data)"}
    except Exception:
        return None, {}


def read_bcf(file_bytes: bytes) -> tuple[pd.DataFrame, dict]:
    """
    Parse a BCF file (Bruker Esprit EBSD) and return (DataFrame, metadata_dict).
    Tries ZIP extraction first, then embedded-text fallback.
    """
    df, meta = _try_zip_bcf(file_bytes)
    if df is not None and len(df) > 0:
        return df, meta

    df, meta = _try_embedded_csv_bcf(file_bytes)
    if df is not None and len(df) > 0:
        return df, meta

    raise ValueError(
        "Could not extract EBSD data from this BCF file.\n\n"
        "**Why this can happen:**\n"
        "- Bruker BCF EBSD files use a proprietary binary format (SFS container) "
        "that is different from ZIP and not fully documented.\n\n"
        "**Solution — export from AztecCrystal or ESPRIT:**\n"
        "1. Open your project in **Bruker ESPRIT** (or Oxford AztecCrystal)\n"
        "2. Go to **File → Export → Text / CSV**\n"
        "3. Select columns: X, Y, Euler1, Euler2, Euler3, Band Contrast, MAD, Phase, KAM\n"
        "4. Save as `.csv` and upload that file here\n\n"
        "The app will then process it exactly the same way."
    )


# ══════════════════════════════════════════════════════════════════════════════
#  UNIFIED LOADER
# ══════════════════════════════════════════════════════════════════════════════

def load_ebsd_file(file_bytes: bytes, filename: str,
                   sep: str = ",", decimal: str = ".") -> tuple[pd.DataFrame, dict]:
    """
    Auto-detect file format and return (DataFrame, metadata).
    Supports: .csv, .txt, .ctf, .bcf
    """
    ext = filename.lower().rsplit(".", 1)[-1]

    if ext == "ctf":
        return read_ctf(file_bytes)

    elif ext == "bcf":
        return read_bcf(file_bytes)

    else:
        # CSV / TXT
        try:
            df = pd.read_csv(
                io.BytesIO(file_bytes),
                sep=sep,
                decimal=decimal,
                comment="#",
                skip_blank_lines=True,
                engine="python",
                on_bad_lines="skip",
            )
            df.columns = df.columns.str.strip()
            df = _coerce_numeric(df)
            meta = {"Source Format": "CSV/Text", "Rows": len(df)}
            return df, meta
        except Exception as e:
            raise ValueError(f"Could not read CSV file: {e}")
