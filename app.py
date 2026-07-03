"""
EBSD Analyzer — Streamlit App v3
Supports CTF (Oxford/HKL), BCF (Bruker), CSV.
Full grain segmentation from pixel-level CTF data.
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from scipy import stats
from scipy.stats import lognorm
import io, zipfile, warnings
warnings.filterwarnings("ignore")

from file_readers import load_ebsd_file
from excel_reference import parse_reference_workbook
from ctf_processing import (
    compute_kam, segment_grains,
    compute_grain_stats, compute_grain_misorientation,
    compute_gnd_from_orientations, burgers_from_lattice,
)
from ipf_plots import (
    plot_ipf_density,
    plot_ipf_3d_cube,
    compute_ipf_colormap,
    plot_ipf_legend,
)
from pf_plots import (
    plot_pole_figures,
    plot_ipf_2d_map,
    st_figure_download,
)

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(page_title="EBSD Analyzer", page_icon="🔬",
                   layout="wide", initial_sidebar_state="expanded")

# ── Publication style ─────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "font.family": "serif", "font.size": 11,
    "axes.labelsize": 12, "axes.titlesize": 13, "axes.linewidth": 1.2,
    "xtick.major.width": 1.2, "ytick.major.width": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "legend.frameon": True, "legend.fontsize": 10,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})
PAL = {"blue":"#1f77b4","orange":"#e87722","green":"#2ca02c",
       "red":"#d62728","purple":"#9467bd","gray":"#7f7f7f"}

def fig_bytes(fig, fmt="png"):
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=300, bbox_inches="tight")
    buf.seek(0); return buf.read()

# ── Column synonym detection ──────────────────────────────────────────────────
GRAIN_SYN   = ["grain diameter","grain size","ecd","equivalent circle diameter",
               "grain_diameter","grain_size","diameter","mean grain diameter"]
MISORI_SYN  = ["misorientation angle","misorientation","grain misorientation",
               "misorientation_angle","average misorientation"]
MAD_SYN     = ["mad","mean angular deviation"]
AREA_SYN    = ["area","grain area","grainarea","grain_area"]
PHASE_SYN   = ["phase","phase name","phase_name"]
KAM_SYN     = ["kam","kernel average misorientation","kernel_average_misorientation",
               "local misorientation"]
EULER_SYN   = [["phi1","euler1","euler_1","euler1 (phi1)"],
               ["phi","euler2","euler_2","euler2 (phi)","euler2 (phi)"],
               ["phi2","euler3","euler_3","euler3 (phi2)"]]
IQ_SYN      = ["iq","image quality","image_quality","fit","band contrast","band_contrast",
               "bandcontrast","bc"]
CI_SYN      = ["ci","confidence index","confidence_index","error"]

def find_col(df, synonyms):
    cl = {c.lower().strip(): c for c in df.columns}
    for s in synonyms:
        if s.lower() in cl:
            return cl[s.lower()]
    return None

def find_euler(df):
    return [find_col(df, s) for s in EULER_SYN]

# ── Outlier detection ─────────────────────────────────────────────────────────
def detect_outliers(series, method="IQR"):
    s = series.dropna()
    if method == "IQR":
        q1, q3 = s.quantile(0.25), s.quantile(0.75)
        iqr = q3 - q1
        mask = (series < q1 - 1.5*iqr) | (series > q3 + 1.5*iqr)
    elif method == "Z-score":
        z = np.abs(stats.zscore(s))
        mask = pd.Series(False, index=series.index)
        mask[s.index] = z > 3
    else:
        med = np.median(s); mad = np.median(np.abs(s - med))
        mz = 0.6745*(s - med)/(mad+1e-9)
        mask = pd.Series(False, index=series.index)
        mask[s.index] = np.abs(mz) > 3.5
    return mask

def stats_table(series, label):
    s = series.dropna()
    d = {"Count": len(s), "Mean": f"{s.mean():.4f}", "Median": f"{s.median():.4f}",
         "Std Dev": f"{s.std():.4f}", "CV (%)": f"{100*s.std()/s.mean():.2f}" if s.mean()!=0 else "—",
         "Min": f"{s.min():.4f}", "Max": f"{s.max():.4f}",
         "Skewness": f"{s.skew():.4f}", "Kurtosis": f"{s.kurtosis():.4f}",
         "P10": f"{s.quantile(0.10):.4f}", "P25 (Q1)": f"{s.quantile(0.25):.4f}",
         "P75 (Q3)": f"{s.quantile(0.75):.4f}", "P90": f"{s.quantile(0.90):.4f}"}
    return pd.DataFrame.from_dict(d, orient="index", columns=[label])

# ═══════════════════════════════════════════════════════════════════════════════
#  SIDEBAR
# ═══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.title("⚙️ Settings")
    uploaded = st.file_uploader(
        "Upload EBSD file(s)",
        type=["csv","txt","ctf","bcf"],
        accept_multiple_files=True,
        help="CTF (Oxford/HKL), BCF (Bruker), CSV/TXT",
    )
    st.caption("Separator settings apply to CSV/TXT only.")
    sep_choice = st.selectbox("Column separator", [",",";","\\t"," "], index=0)
    sep = "\t" if sep_choice == "\\t" else sep_choice
    decimal_c = st.selectbox("Decimal separator", [".","," ], index=0)

    st.divider()
    st.subheader("Reference workbook (optional)")
    ref_upload = st.file_uploader(
        "EBSD export workbook (.xlsx / .xlsm)",
        type=["xlsx", "xlsm"],
        accept_multiple_files=False,
        help="Optional Excel export from AztecCrystal / ESPRIT (Overview, Grain List, "
             "Boundary Statistics, pole & Mackenzie plots). Used for calibration and "
             "cross-checking — it does NOT replace the uploaded EBSD map.",
    )

    st.divider()
    st.subheader("CTF processing")
    grain_threshold = st.slider(
        "Grain boundary threshold (°)", 2, 20, 10,
        help="Misorientation angle used to separate grains in CTF pixel data. "
             "Typical: 10°–15° for ferritic/austenitic steels.")
    kam_threshold = st.slider(
        "KAM threshold (°)", 1, 15, 5,
        help="Maximum misorientation included in KAM calculation (suppress GBs).")

    st.divider()
    st.subheader("Analysis options")
    outlier_method = st.selectbox("Outlier detection method",
                                  ["IQR","Z-score","Modified Z-score"])
    remove_outliers = st.toggle("Remove outliers from plots", value=False)
    bins_grain  = st.slider("Grain size bins",  10, 80, 30)
    bins_misori = st.slider("Misorientation bins", 10, 90, 36)

    st.divider()
    st.subheader("Plot options")
    color_grain  = st.color_picker("Grain size color",   "#1f77b4")
    color_misori = st.color_picker("Misorientation color","#e87722")
    color_out    = st.color_picker("Outlier color",       "#d62728")
    show_fit     = st.toggle("Log-normal fit (grain size)",    value=True)
    show_mack    = st.toggle("Mackenzie curve (misorientation)",value=True)
    plot_fmt     = st.selectbox("Export format", ["PNG","SVG","PDF"])

# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════
st.title("🔬 EBSD Analyzer")
st.markdown("Upload **CTF**, **BCF**, or **CSV** — grain size, misorientation, texture, KAM, outliers, publication-ready figures.")

if not uploaded:
    st.info("""
    👈 Upload your EBSD file in the sidebar.

    **Supported formats:**
    - `.ctf` — Oxford Instruments / HKL Channel 5 *(pixel data → auto grain segmentation)*
    - `.bcf` — Bruker Esprit EBSD
    - `.csv` / `.txt` — OIM, AztecCrystal, MTEX exports

    Sample datasets in `sample_data/` folder.
    """)
    st.stop()

# ── Load files ────────────────────────────────────────────────────────────────
@st.cache_data(show_spinner=False)
def load_file(file_bytes, filename, sep, decimal):
    try:
        df, meta = load_ebsd_file(file_bytes, filename, sep=sep, decimal=decimal)
        return df, meta, None
    except Exception as e:
        return None, {}, str(e)

all_data, all_meta = {}, {}
for f in uploaded:
    raw = f.read()
    with st.spinner(f"Reading {f.name}…"):
        df, meta, err = load_file(raw, f.name, sep, decimal_c)
    if err:
        st.error(f"**{f.name}**: {err}")
    else:
        all_data[f.name] = df
        all_meta[f.name] = meta

if not all_data:
    st.stop()

# ═══════════════════════════════════════════════════════════════════════════════
#  OPTIONAL REFERENCE WORKBOOK (.xlsx / .xlsm)
# ═══════════════════════════════════════════════════════════════════════════════
@st.cache_data(show_spinner=False)
def load_reference(file_bytes):
    try:
        return parse_reference_workbook(file_bytes), None
    except Exception as e:
        return None, str(e)

ref = None
if ref_upload is not None:
    ref_bytes = ref_upload.read()
    with st.spinner(f"Reading reference workbook {ref_upload.name}…"):
        ref, ref_err = load_reference(ref_bytes)
    if ref_err:
        st.error(f"**{ref_upload.name}**: could not parse reference workbook — {ref_err}")
    elif ref:
        ov = ref.get("overview", {})
        ref_step = ov.get("step_size_um")
        msg = f"📥 Reference workbook loaded: **{ref_upload.name}**"
        extras = []
        if ref_step:
            extras.append(f"Step Size **{ref_step:.4g} µm**")
        if ov.get("pixel_count"):
            extras.append(f"Pixel Count **{int(ov['pixel_count']):,}**")
        if extras:
            msg += " — " + ", ".join(extras)
        st.success(msg)
        with st.expander("📊 Reference workbook summary (Excel)", expanded=True):
            # ── Overview / acquisition metadata ────────────────────────────────
            if ov:
                st.markdown("**Acquisition metadata (Overview)**")
                ov_rows = []
                _labels = {
                    "source_file": "Source file",
                    "step_size_um": "Step Size (µm)",
                    "pixel_count": "Pixel Count",
                    "raster": "Raster",
                    "hit_rate_pct": "Hit Rate (%)",
                    "zero_solution_count": "Zero Solution Count",
                }
                for k, lbl in _labels.items():
                    if ov.get(k) is not None:
                        ov_rows.append({"Property": lbl, "Value": ov[k]})
                if ov_rows:
                    st.dataframe(pd.DataFrame(ov_rows), use_container_width=False,
                                 hide_index=True)
                if ov.get("phases"):
                    ph_df = pd.DataFrame(ov["phases"])[["index", "name", "fraction_pct"]]
                    ph_df.columns = ["#", "Phase", "Fraction (%)"]
                    st.markdown("**Phase fractions**")
                    st.dataframe(ph_df, use_container_width=False, hide_index=True)

            # ── Grain list summary ──────────────────────────────────────────────
            grain_ref = ref.get("grain", {})
            if grain_ref.get("grain_count"):
                st.markdown(f"**Grain List** — {grain_ref['grain_count']:,} grains")
                if grain_ref.get("phase_counts"):
                    st.caption("Grains per phase: " + ", ".join(
                        f"{k}: {v}" for k, v in grain_ref["phase_counts"].items()))
                summary_stats = grain_ref.get("stats", {})
                if summary_stats:
                    label_map = {
                        "area": "Area (µm²)", "ecd": "ECD (µm)",
                        "feret": "Max Feret (µm)", "perimeter": "Perimeter (µm)",
                        "mos_mean": "Mean Orientation Spread (°)",
                        "mos_max": "Maximum Orientation Spread (°)",
                    }
                    srows = []
                    for key, lbl in label_map.items():
                        if key in summary_stats:
                            s = summary_stats[key]
                            srows.append({
                                "Attribute": lbl, "N": s["count"],
                                "Mean": round(s["mean"], 4), "Median": round(s["median"], 4),
                                "Std": round(s["std"], 4), "Min": round(s["min"], 4),
                                "Max": round(s["max"], 4),
                            })
                    if srows:
                        st.dataframe(pd.DataFrame(srows), use_container_width=True,
                                     hide_index=True)

            # ── Boundary statistics ─────────────────────────────────────────────
            bnd = ref.get("boundary", {})
            if bnd:
                st.markdown("**Boundary Statistics (LAGB 2–10° / HAGB >10°)**")
                brows = []
                for phase, d in bnd.items():
                    brows.append({
                        "Phase": phase,
                        "LAGB length (µm)": d.get("LAGB_length_um"),
                        "LAGB fraction (%)": d.get("LAGB_fraction_pct"),
                        "HAGB length (µm)": d.get("HAGB_length_um"),
                        "HAGB fraction (%)": d.get("HAGB_fraction_pct"),
                    })
                st.dataframe(pd.DataFrame(brows), use_container_width=True, hide_index=True)

            # ── Pole figure MUD peaks ───────────────────────────────────────────
            pole = ref.get("pole", [])
            if pole:
                st.markdown("**Pole-figure profiles — texture strength (m.u.d.)**")
                st.dataframe(pd.DataFrame([{
                    "Sheet": p["sheet"], "Max MUD": round(p["max_mud"], 3),
                    "Mean MUD": round(p["mean_mud"], 3),
                    "Angle at peak (°)": p["angle_at_peak_deg"],
                } for p in pole]), use_container_width=True, hide_index=True)

            # ── Mackenzie / disorientation ──────────────────────────────────────
            mack = ref.get("mackenzie", [])
            if mack:
                st.markdown("**Mackenzie / disorientation (measured neighbour pairs)**")
                st.dataframe(pd.DataFrame([{
                    "Sheet": m["sheet"],
                    "Mean disorientation (°)": round(m["mean_neighbor_disorientation_deg"], 3),
                    "Peak angle (°)": m["peak_angle_deg"],
                    "LAGB fraction (<15°)": round(m["lagb_fraction_lt15deg"], 4)
                    if m["lagb_fraction_lt15deg"] is not None else None,
                } for m in mack]), use_container_width=True, hide_index=True)

            st.caption(
                "This workbook is a **reference/cross-check** derived from an already "
                "analysed dataset. The app still computes everything from the EBSD map "
                "you uploaded; workbook values are shown for calibration and comparison.")

# preferred calibration step size from the reference workbook (if any)
ref_step_um = None
if ref and ref.get("overview"):
    ref_step_um = ref["overview"].get("step_size_um")

# ═══════════════════════════════════════════════════════════════════════════════
#  PER-FILE TABS
# ═══════════════════════════════════════════════════════════════════════════════
for fname, df_raw in all_data.items():
    st.header(f"📂 {fname}")
    meta = all_meta.get(fname, {})
    is_ctf = meta.get("Source Format","").startswith("CTF")
    is_pixel_data = is_ctf   # CTF = pixel data, needs segmentation

    # ── Metadata panel ────────────────────────────────────────────────────────
    with st.expander("📄 File info & header", expanded=False):
        st.markdown(f"**Format:** `{meta.get('Source Format','CSV')}`  |  "
                    f"**Rows:** {len(df_raw):,}  |  "
                    f"**Columns:** {len(df_raw.columns)}")
        skip = {"Source Format","Rows","Columns"}
        m2 = {k:v for k,v in meta.items() if k not in skip}
        if m2:
            st.dataframe(pd.DataFrame.from_dict(m2, orient="index", columns=["Value"]),
                         use_container_width=False)

    # ── Raw preview ───────────────────────────────────────────────────────────
    with st.expander("📋 Raw data preview", expanded=False):
        st.dataframe(df_raw.head(50), use_container_width=True)

    # ── Column mapping ────────────────────────────────────────────────────────
    num_cols = ["(none)"] + [c for c in df_raw.columns
                             if pd.api.types.is_numeric_dtype(df_raw[c])]
    cat_cols = ["(none)"] + list(df_raw.columns)

    with st.expander("🔧 Column mapping (auto-detected)", expanded=is_pixel_data):
        if is_pixel_data:
            st.info("CTF file detected — pixel-level data. The app will **segment grains automatically** "
                    "using the Euler angles and the grain boundary threshold set in the sidebar. "
                    "Select the Euler angle columns below (usually auto-detected).")

        c1, c2, c3 = st.columns(3)
        with c1:
            e1_auto = find_col(df_raw, EULER_SYN[0])
            e1_col  = st.selectbox("Euler φ₁",
                                   num_cols, index=num_cols.index(e1_auto) if e1_auto in num_cols else 0,
                                   key=f"e1_{fname}")
        with c2:
            e2_auto = find_col(df_raw, EULER_SYN[1])
            e2_col  = st.selectbox("Euler Φ",
                                   num_cols, index=num_cols.index(e2_auto) if e2_auto in num_cols else 0,
                                   key=f"e2_{fname}")
        with c3:
            e3_auto = find_col(df_raw, EULER_SYN[2])
            e3_col  = st.selectbox("Euler φ₂",
                                   num_cols, index=num_cols.index(e3_auto) if e3_auto in num_cols else 0,
                                   key=f"e3_{fname}")

        c4, c5, c6 = st.columns(3)
        with c4:
            mad_auto  = find_col(df_raw, MAD_SYN)
            mad_col   = st.selectbox("MAD (Mean Angular Dev.)",
                                     num_cols, index=num_cols.index(mad_auto) if mad_auto in num_cols else 0,
                                     key=f"mad_{fname}")
            phase_auto = find_col(df_raw, PHASE_SYN)
            phase_col  = st.selectbox("Phase",
                                      cat_cols, index=cat_cols.index(phase_auto) if phase_auto in cat_cols else 0,
                                      key=f"phase_{fname}")
        with c5:
            iq_auto = find_col(df_raw, IQ_SYN)
            iq_col  = st.selectbox("Band Contrast / IQ",
                                   num_cols, index=num_cols.index(iq_auto) if iq_auto in num_cols else 0,
                                   key=f"iq_{fname}")
            ci_auto = find_col(df_raw, CI_SYN)
            ci_col  = st.selectbox("Error / CI",
                                   num_cols, index=num_cols.index(ci_auto) if ci_auto in num_cols else 0,
                                   key=f"ci_{fname}")
        with c6:
            # Non-CTF only
            grain_auto  = find_col(df_raw, GRAIN_SYN)
            grain_col   = st.selectbox("Grain Diameter (CSV only)",
                                       num_cols, index=num_cols.index(grain_auto) if grain_auto in num_cols else 0,
                                       key=f"grain_{fname}")
            misori_auto = find_col(df_raw, MISORI_SYN)
            misori_col  = st.selectbox("Misorientation (CSV only)",
                                       num_cols, index=num_cols.index(misori_auto) if misori_auto in num_cols else 0,
                                       key=f"misori_{fname}")
            kam_auto = find_col(df_raw, KAM_SYN)
            kam_col  = st.selectbox("KAM (CSV only)",
                                    num_cols, index=num_cols.index(kam_auto) if kam_auto in num_cols else 0,
                                    key=f"kam_{fname}")

    def col_or_none(v): return None if v == "(none)" else v
    e1_col    = col_or_none(e1_col)
    e2_col    = col_or_none(e2_col)
    e3_col    = col_or_none(e3_col)
    mad_col   = col_or_none(mad_col)
    phase_col = col_or_none(phase_col)
    iq_col    = col_or_none(iq_col)
    ci_col    = col_or_none(ci_col)
    grain_col = col_or_none(grain_col)
    misori_col= col_or_none(misori_col)
    kam_col   = col_or_none(kam_col)

    # ── CI / Error filter ─────────────────────────────────────────────────────
    df = df_raw.copy()
    if ci_col:
        ci_series = df[ci_col].dropna()
        if ci_series.max() <= 1.0:
            ci_min = st.slider("Min CI (quality filter)", 0.0, 1.0, 0.0, 0.05, key=f"cif_{fname}")
        else:
            ci_min = st.slider("Max Error filter (0=indexed)", 0.0, float(ci_series.max()), float(ci_series.max()), 0.1, key=f"cif_{fname}")
        df = df[df[ci_col] <= ci_min if ci_series.max() > 1.0 else df[ci_col] >= ci_min].copy()
        st.caption(f"{len(df):,} points retained after quality filter.")

    # ════════════════════════════════════════════════════════════════════════
    #  CTF GRAIN SEGMENTATION
    # ════════════════════════════════════════════════════════════════════════
    grain_df = None   # per-grain summary DataFrame

    if is_pixel_data and e1_col and e2_col and e3_col:
        x_col_name = find_col(df, ["x","x position","x_position"]) or "X"
        y_col_name = find_col(df, ["y","y position","y_position"]) or "Y"

        @st.cache_data(show_spinner=False)
        def run_segmentation(data_bytes, phi1c, Phic, phi2c, xc, yc, thresh, kam_thresh, step):
            _df = pd.read_json(io.StringIO(data_bytes.decode()))
            gids   = segment_grains(_df, phi1c, Phic, phi2c, xc, yc, thresh)
            g_df   = compute_grain_stats(_df, gids, xc, yc, step)
            if not g_df.empty:
                g_df["Misorientation Angle"] = compute_grain_misorientation(
                    _df, g_df, phi1c, Phic, phi2c).values
            kam_arr = compute_kam(_df, phi1c, Phic, phi2c, xc, yc,
                                  kernel_order=1, threshold_deg=kam_thresh)
            return g_df, kam_arr, gids

        step_meta = meta.get("Step Size (µm)")
        if step_meta is None and ref_step_um:
            step_val = float(ref_step_um)
            st.caption(f"ℹ️ Step size not found in the EBSD file — using reference "
                       f"workbook value **{step_val:.4g} µm** for grain segmentation.")
        else:
            step_val = float(step_meta) if step_meta is not None else 0.5
            if ref_step_um and step_meta is not None and abs(step_val - ref_step_um) > 1e-4:
                st.caption(f"ℹ️ Step size — file: **{step_val:.4g} µm**, reference "
                           f"workbook: **{ref_step_um:.4g} µm** (using file value for segmentation).")

        with st.spinner("Segmenting grains and computing KAM… (may take ~30 s for large files)"):
            try:
                df_bytes = df.to_json().encode()
                grain_df, kam_arr, grain_ids = run_segmentation(
                    df_bytes, e1_col, e2_col, e3_col,
                    x_col_name, y_col_name,
                    grain_threshold, kam_threshold, step_val
                )
                df["KAM (computed)"] = kam_arr
                kam_col_use = "KAM (computed)"
                if not grain_df.empty:
                    st.success(f"Segmentation complete — **{len(grain_df):,} grains** identified "
                               f"(threshold {grain_threshold}°, step {step_val} µm).")
            except Exception as ex:
                st.warning(f"Grain segmentation error: {ex}. Pixel-level analysis only.")
                grain_df = pd.DataFrame()
                kam_col_use = mad_col
    else:
        # CSV mode: use provided columns directly
        grain_df    = None
        kam_col_use = kam_col

    # Working DataFrames
    # grain_df → per-grain stats (from CTF segmentation or None for CSV)
    # df       → pixel/point data (with computed KAM for CTF)

    def grain_series(col_name):
        """Get a series from grain_df if available, else from df."""
        if grain_df is not None and not grain_df.empty and col_name in grain_df.columns:
            return grain_df[col_name].dropna()
        if col_name in df.columns:
            return df[col_name].dropna()
        return pd.Series(dtype=float)

    def get_gs():
        if grain_df is not None and not grain_df.empty and "Grain Diameter" in grain_df.columns:
            return grain_df["Grain Diameter"].dropna()
        if grain_col:
            return df[grain_col].dropna()
        return pd.Series(dtype=float)

    def get_misori():
        if grain_df is not None and not grain_df.empty and "Misorientation Angle" in grain_df.columns:
            return grain_df["Misorientation Angle"].dropna()
        if misori_col:
            return df[misori_col].dropna()
        if mad_col:
            return df[mad_col].dropna()
        return pd.Series(dtype=float)

    # ════════════════════════════════════════════════════════════════════════
    #  ANALYSIS TABS
    # ════════════════════════════════════════════════════════════════════════
    sec = st.tabs(["📏 Grain Size", "📐 Misorientation", "🌐 Texture",
                   "⚠️ Outliers", "📊 KAM / Band Contrast",
                   "🔷 Pole Figures (IPF)", "💾 Export"])

    # ───────────────────────────────────────────────────────────────────────
    # TAB 1 — GRAIN SIZE
    # ───────────────────────────────────────────────────────────────────────
    with sec[0]:
        gs = get_gs()
        if gs.empty:
            if is_pixel_data and (not e1_col or not e2_col or not e3_col):
                st.warning("Select the three Euler angle columns above to enable grain segmentation.")
            else:
                st.warning("No grain diameter data available.")
        else:
            gs = gs[gs > 0]
            omask = detect_outliers(gs, outlier_method)
            gs_plot = gs[~omask] if remove_outliers else gs
            if remove_outliers and omask.sum():
                st.info(f"Removed {omask.sum()} outlier(s) ({100*omask.sum()/len(gs):.1f}%).")

            st.subheader("Descriptive Statistics — Grain Diameter (µm)")
            st.dataframe(stats_table(gs_plot, "Value"), use_container_width=False)

            if len(gs_plot) >= 8:
                _, p_raw = stats.shapiro(gs_plot.sample(min(5000,len(gs_plot)), random_state=42))
                _, p_log = stats.shapiro(np.log(gs_plot.clip(1e-9)).sample(min(5000,len(gs_plot)), random_state=42))
                ca, cb = st.columns(2)
                with ca: st.metric("Shapiro-Wilk p (raw)", f"{p_raw:.4f}", help="p>0.05 → normal")
                with cb: st.metric("Shapiro-Wilk p (log)", f"{p_log:.4f}", help="p>0.05 → log-normal")

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            ax = axes[0]
            cnt, _, _ = ax.hist(gs_plot, bins=bins_grain, color=color_grain,
                                edgecolor="white", lw=0.5, alpha=0.85, density=True, label="Data")
            if show_fit:
                try:
                    sh, loc, sc = lognorm.fit(gs_plot, floc=0)
                    xf = np.linspace(gs_plot.min(), gs_plot.max(), 300)
                    ax.plot(xf, lognorm.pdf(xf, sh, loc, sc), color=PAL["red"], lw=2,
                            label=f"Log-normal (σ={sh:.2f}, μ={np.log(sc):.2f})")
                except Exception: pass
            ax.set_xlabel("Grain Diameter (µm)"); ax.set_ylabel("Probability Density")
            ax.set_title("Grain Size Distribution"); ax.legend()
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

            ax2 = axes[1]
            sgs = np.sort(gs_plot)
            ax2.plot(sgs, np.arange(1,len(sgs)+1)/len(sgs)*100, color=color_grain, lw=2)
            ax2.axvline(gs_plot.median(), color=PAL["red"],    ls="--", lw=1.5,
                        label=f"Median = {gs_plot.median():.2f} µm")
            ax2.axvline(gs_plot.mean(),   color=PAL["orange"], ls=":",  lw=1.5,
                        label=f"Mean = {gs_plot.mean():.2f} µm")
            ax2.set_xlabel("Grain Diameter (µm)"); ax2.set_ylabel("Cumulative Frequency (%)")
            ax2.set_title("Cumulative Grain Size Distribution"); ax2.legend()
            ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            fig.tight_layout()
            st.pyplot(fig)
            st_figure_download(fig, f"grain_size_{fname.rsplit('.',1)[0]}",
                               fmt=plot_fmt.lower(), key=f"dl_gs_{fname}")
            st.session_state[f"fig_gs_{fname}"] = fig

            # By phase
            phase_src = grain_df if (grain_df is not None and not grain_df.empty
                                     and "Grain Diameter" in grain_df.columns) else df
            phase_col_use = "Phase" if "Phase" in phase_src.columns else phase_col
            if phase_col_use and phase_col_use in phase_src.columns and "Grain Diameter" in phase_src.columns:
                st.subheader("Grain Size by Phase")
                fig_ph, ax_ph = plt.subplots(figsize=(8, 4))
                for ph in phase_src[phase_col_use].dropna().unique():
                    sub = phase_src[phase_src[phase_col_use]==ph]["Grain Diameter"].dropna()
                    sub = sub[sub>0]
                    if len(sub) > 1:
                        ax_ph.hist(sub, bins=bins_grain, alpha=0.6, label=str(ph),
                                   density=True, edgecolor="white", lw=0.4)
                ax_ph.set_xlabel("Grain Diameter (µm)"); ax_ph.set_ylabel("Probability Density")
                ax_ph.set_title("Grain Size by Phase"); ax_ph.legend()
                fig_ph.tight_layout(); st.pyplot(fig_ph)
                st_figure_download(fig_ph, f"grain_size_by_phase_{fname.rsplit('.',1)[0]}",
                                   fmt=plot_fmt.lower(), key=f"dl_gs_ph_{fname}")
                st.session_state[f"fig_gs_ph_{fname}"] = fig_ph

    # ───────────────────────────────────────────────────────────────────────
    # TAB 2 — MISORIENTATION
    # ───────────────────────────────────────────────────────────────────────
    with sec[1]:
        mo = get_misori()
        if mo.empty:
            st.warning("No misorientation data available. "
                       "For CTF files, select Euler columns and run grain segmentation.")
        else:
            mo = mo[(mo >= 0) & (mo <= 180)]
            if mo.empty:
                st.warning("Misorientation values out of valid range (0–180°).")
            else:
                omask_mo = detect_outliers(mo, outlier_method)
                mo_plot  = mo[~omask_mo] if remove_outliers else mo

                st.subheader("Descriptive Statistics — Misorientation Angle (°)")
                st.dataframe(stats_table(mo_plot, "Value"), use_container_width=False)

                low  = (mo_plot < 15).sum() / len(mo_plot) * 100
                high = (mo_plot >= 15).sum() / len(mo_plot) * 100
                r1, r2, r3 = st.columns(3)
                with r1: st.metric("LAGB < 15°",    f"{low:.1f}%")
                with r2: st.metric("HAGB ≥ 15°",    f"{high:.1f}%")
                with r3: st.metric("Fraction 15°–65°",
                                   f"{(mo_plot[(mo_plot>=15)&(mo_plot<=65)].count()/len(mo_plot)*100):.1f}%")

                def mackenzie(a_deg):
                    a = np.deg2rad(a_deg)
                    p = np.zeros_like(a, dtype=float)
                    m1 = a <= np.deg2rad(45)
                    m2 = (a > np.deg2rad(45)) & (a <= np.deg2rad(60))
                    p[m1] = (1 - np.cos(a[m1])) * (2/np.pi) * (8/(np.sqrt(2)-1))
                    v2 = np.sqrt(2)*np.cos(a[m2]) - 2*np.cos(a[m2]) + 1
                    p[m2] = np.clip(v2, 0, None) * (2/np.pi) * (8/(np.sqrt(2)-1))
                    return p

                fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))
                ax = axes2[0]
                cnt_m, _, _ = ax.hist(mo_plot, bins=bins_misori, range=(0,65),
                                      color=color_misori, edgecolor="white", lw=0.5,
                                      alpha=0.85, density=True, label="Measured")
                if show_mack:
                    xm = np.linspace(0, 65, 500)
                    ym = mackenzie(xm)
                    ym = ym / (ym.max()+1e-9) * cnt_m.max()
                    ax.plot(xm, ym, color=PAL["blue"], lw=2, ls="--",
                            label="Mackenzie (random)")
                ax.axvline(15, color=PAL["red"], lw=1.2, ls=":",
                           label="15° LAGB/HAGB")
                ax.set_xlabel("Misorientation Angle (°)"); ax.set_ylabel("Frequency Density")
                ax.set_title("Misorientation Distribution"); ax.legend()
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

                ax2b = axes2[1]
                bins_s = np.arange(0, 66, 5)
                h_la, _ = np.histogram(mo_plot[mo_plot < 15],  bins=bins_s)
                h_ha, _ = np.histogram(mo_plot[mo_plot >= 15], bins=bins_s)
                bc2 = 0.5*(bins_s[:-1]+bins_s[1:])
                ax2b.bar(bc2, h_la, width=4.5, color=PAL["blue"],    alpha=0.85,
                         label="LAGB (<15°)", edgecolor="white")
                ax2b.bar(bc2, h_ha, width=4.5, color=color_misori, alpha=0.85,
                         bottom=h_la, label="HAGB (≥15°)", edgecolor="white")
                ax2b.set_xlabel("Misorientation Angle (°)"); ax2b.set_ylabel("Count")
                ax2b.set_title("LAGB vs. HAGB Frequency"); ax2b.legend()
                ax2b.xaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax2b.yaxis.set_minor_locator(ticker.AutoMinorLocator())
                fig2.tight_layout(); st.pyplot(fig2)
                st_figure_download(fig2, f"misorientation_{fname.rsplit('.',1)[0]}",
                                   fmt=plot_fmt.lower(), key=f"dl_mo_{fname}")
                st.session_state[f"fig_mo_{fname}"] = fig2

    # ───────────────────────────────────────────────────────────────────────
    # TAB 3 — TEXTURE
    # ───────────────────────────────────────────────────────────────────────
    with sec[2]:
        if not (e1_col and e2_col and e3_col):
            st.warning("Select the three Euler angle columns above.")
        else:
            # Use pixel data directly for texture (more statistically representative)
            phi1_s = df[e1_col].dropna()
            Phi_s  = df[e2_col].dropna()
            phi2_s = df[e3_col].dropna()
            cidx   = phi1_s.index.intersection(Phi_s.index).intersection(phi2_s.index)
            phi1_s, Phi_s, phi2_s = phi1_s[cidx], Phi_s[cidx], phi2_s[cidx]

            st.subheader("Euler Angle Summary")
            ce1, ce2, ce3 = st.columns(3)
            for col_e, (nm, ser) in zip([ce1,ce2,ce3],
                                        [("φ₁",phi1_s),("Φ",Phi_s),("φ₂",phi2_s)]):
                with col_e:
                    st.metric(f"Mean {nm}", f"{ser.mean():.2f}°")
                    st.metric(f"Std {nm}",  f"{ser.std():.2f}°")

            # Ideal orientations (BCC — ferrite/stainless)
            IDEALS = {
                "{001}<100> Cube":        (0,  0,  0),
                "{110}<001> Goss":        (0,  45, 0),
                "{111}<110> γ-fiber":     (0,  55, 45),
                "{111}<112> γ-fiber":     (30, 55, 45),
                "{112}<110> Brass (BCC)": (35, 35, 0),
                "{001}<110> Rotated Cube":(45, 0,  0),
            }
            tol = 15
            st.subheader(f"Ideal Orientation Fractions (±{tol}° tolerance)")
            rows_t = []
            for name, (p1i, Pi, p2i) in IDEALS.items():
                d1 = np.abs(phi1_s.values - p1i) % 360; d1 = np.minimum(d1, 360-d1)
                dP = np.abs(Phi_s.values  - Pi)
                d2 = np.abs(phi2_s.values - p2i) % 360; d2 = np.minimum(d2, 360-d2)
                near = ((d1<=tol)&(dP<=tol)&(d2<=tol)).sum()
                rows_t.append({"Orientation": name,
                                "Count": int(near),
                                "Fraction (%)": f"{100*near/len(phi1_s):.2f}"})
            st.dataframe(pd.DataFrame(rows_t), use_container_width=True)

            # Euler distributions
            fig3, axes3 = plt.subplots(1, 3, figsize=(13, 4))
            for ax_e, (nm, ser, xlim) in zip(axes3, [
                ("φ₁ (°)", phi1_s, (0,360)),
                ("Φ (°)",  Phi_s,  (0,90)),
                ("φ₂ (°)", phi2_s, (0,360)),
            ]):
                ax_e.hist(ser, bins=36, color=PAL["purple"],
                          edgecolor="white", lw=0.4, alpha=0.85, density=True)
                ax_e.set_xlabel(nm); ax_e.set_ylabel("Frequency Density")
                ax_e.set_title(f"{nm} Distribution"); ax_e.set_xlim(xlim)
                ax_e.xaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax_e.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            fig3.suptitle("Euler Angle Distributions", fontsize=13, y=1.02)
            fig3.tight_layout(); st.pyplot(fig3)
            st_figure_download(fig3, f"euler_distributions_{fname.rsplit('.',1)[0]}",
                               fmt=plot_fmt.lower(), key=f"dl_euler_{fname}")
            st.session_state[f"fig_euler_{fname}"] = fig3

            # ODF section Φ vs φ₂
            st.subheader("Φ vs. φ₂ Orientation Density (φ₁=0° section)")
            fig4, ax4 = plt.subplots(figsize=(6, 5))
            h2d, xe, ye = np.histogram2d(phi2_s, Phi_s, bins=90, range=[[0,90],[0,90]])
            h2d = h2d.T; h2d[h2d==0] = np.nan
            im = ax4.imshow(h2d, origin="lower", aspect="auto",
                            extent=[0,90,0,90], cmap="inferno",
                            norm=LogNorm(vmin=1))
            plt.colorbar(im, ax=ax4, label="Count (log)")
            ax4.set_xlabel("φ₂ (°)"); ax4.set_ylabel("Φ (°)")
            ax4.set_title("Orientation Density — Φ vs. φ₂")
            ax4.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            fig4.tight_layout(); st.pyplot(fig4)
            st_figure_download(fig4, f"odf_section_{fname.rsplit('.',1)[0]}",
                               fmt=plot_fmt.lower(), key=f"dl_odf_{fname}")
            st.session_state[f"fig_odf_{fname}"] = fig4

    # ───────────────────────────────────────────────────────────────────────
    # TAB 4 — OUTLIERS
    # ───────────────────────────────────────────────────────────────────────
    with sec[3]:
        st.subheader(f"Outlier Detection — {outlier_method}")
        targets = {}
        gs_s = get_gs()
        mo_s = get_misori()
        if not gs_s.empty: targets["Grain Diameter"] = gs_s
        if not mo_s.empty: targets["Misorientation"] = mo_s
        if mad_col and mad_col in df.columns:
            targets["MAD"] = df[mad_col].dropna()
        if grain_df is not None and not grain_df.empty and "Grain Area" in grain_df.columns:
            targets["Grain Area"] = grain_df["Grain Area"].dropna()

        if not targets:
            st.info("Run grain segmentation (select Euler columns) to enable outlier analysis.")
        else:
            fig_out, axes_out = plt.subplots(1, len(targets), figsize=(5*len(targets), 4.5))
            if len(targets)==1: axes_out=[axes_out]
            summary = []
            for ax_o, (label, series) in zip(axes_out, targets.items()):
                mask = detect_outliers(series, outlier_method)
                n_out = mask.sum()
                summary.append({"Column": label, "Total": len(series),
                                 "Outliers": int(n_out),
                                 "Fraction (%)": f"{100*n_out/len(series):.2f}",
                                 "Method": outlier_method})
                ax_o.boxplot(series, vert=True, patch_artist=True,
                             boxprops=dict(facecolor=PAL["blue"], alpha=0.5),
                             medianprops=dict(color=PAL["red"], lw=2),
                             flierprops=dict(marker="o", color=color_out, markersize=4, alpha=0.5))
                ov = series[mask]
                if len(ov):
                    ax_o.scatter([1]*len(ov), ov, color=color_out, s=20, zorder=5,
                                 label=f"{n_out} outliers")
                    ax_o.legend(fontsize=9)
                ax_o.set_title(label); ax_o.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            fig_out.suptitle(f"Outlier Detection ({outlier_method})", fontsize=13)
            fig_out.tight_layout(); st.pyplot(fig_out)
            st_figure_download(fig_out, f"outliers_{fname.rsplit('.',1)[0]}",
                               fmt=plot_fmt.lower(), key=f"dl_out_{fname}")
            st.session_state[f"fig_out_{fname}"] = fig_out
            st.dataframe(pd.DataFrame(summary), use_container_width=True)

    # ───────────────────────────────────────────────────────────────────────
    # TAB 5 — KAM / BAND CONTRAST
    # ───────────────────────────────────────────────────────────────────────
    with sec[4]:
        plots_k = []
        if kam_col_use and kam_col_use in df.columns:
            plots_k.append(("KAM (°)", df[kam_col_use].dropna()))
        if mad_col and mad_col in df.columns and mad_col != kam_col_use:
            plots_k.append(("MAD (°)", df[mad_col].dropna()))
        if iq_col and iq_col in df.columns:
            plots_k.append(("Band Contrast", df[iq_col].dropna()))

        if not plots_k:
            st.info("No KAM or Band Contrast data available yet. "
                    "For CTF files, select Euler columns to trigger KAM computation.")
        else:
            fig_k, axes_k = plt.subplots(1, len(plots_k), figsize=(6*len(plots_k), 4.5))
            if len(plots_k)==1: axes_k=[axes_k]
            for ax_k, (label, series) in zip(axes_k, plots_k):
                ax_k.hist(series, bins=50, color=PAL["green"],
                          edgecolor="white", lw=0.4, alpha=0.85, density=True)
                ax_k.axvline(series.mean(), color=PAL["red"], ls="--", lw=1.5,
                             label=f"Mean = {series.mean():.3f}")
                ax_k.set_xlabel(label); ax_k.set_ylabel("Frequency Density")
                ax_k.set_title(f"{label} Distribution"); ax_k.legend()
                ax_k.xaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax_k.yaxis.set_minor_locator(ticker.AutoMinorLocator())
                st.subheader(f"Statistics — {label}")
                st.dataframe(stats_table(series, "Value"), use_container_width=False)
            fig_k.tight_layout(); st.pyplot(fig_k)
            st_figure_download(fig_k, f"kam_bc_{fname.rsplit('.',1)[0]}",
                               fmt=plot_fmt.lower(), key=f"dl_kam_{fname}")
            st.session_state[f"fig_kam_{fname}"] = fig_k

            # ── KAM-derived apparent GND density ───────────────────────────
            if kam_col_use and kam_col_use in df.columns:
                st.subheader("KAM-derived apparent GND density")
                st.caption("Densidade aparente de GND estimada a partir do KAM")

                st.warning(
                    "**This is not the total dislocation density.** It is an "
                    "*apparent, lower-bound* estimate of **geometrically necessary "
                    "dislocations (GNDs)** resolved by the KAM kernel. It **excludes "
                    "statistically-stored dislocations (SSDs)** and is strongly "
                    "**method-dependent**: it changes with step size, kernel order, "
                    "angular noise, map clean-up, and grain-boundary / phase / quality "
                    "filtering. Treat it as a *KAM-derived proxy*, not a measured value.")

                # Is a real pixel map available for the distance-aware gradient?
                x_col_g = find_col(df, ["x", "x position", "x_position"])
                y_col_g = find_col(df, ["y", "y position", "y_position"])
                have_map = bool(is_pixel_data and e1_col and e2_col and e3_col
                                and x_col_g and y_col_g)

                # default step: file value, else reference workbook, else 0.5
                step_default = meta.get("Step Size (µm)")
                if step_default is None:
                    step_default = ref_step_um if ref_step_um else 0.5

                # ── Phase-aware Burgers vector ──────────────────────────────
                struct_defaults = {
                    "BCC (ferrite / martensite, Fe)": ("BCC", 0.2866),
                    "FCC (austenite / Al / Cu)":      ("FCC", 0.3595),
                    "HCP (Ti / Mg / Zn)":             ("HCP", 0.2950),
                    "Custom":                          ("BCC", 0.2866),
                }
                bc1, bc2, bc3 = st.columns(3)
                with bc1:
                    struct_key = st.selectbox(
                        "Crystal structure", list(struct_defaults.keys()),
                        index=0, key=f"struct_{fname}",
                        help="Sets how b is derived from the lattice parameter a: "
                             "BCC b=√3/2·a (½⟨111⟩), FCC b=a/√2 (½⟨110⟩), HCP b≈a. "
                             "Edit these examples to match your material — they are "
                             "not hidden assumptions.")
                    structure, a_default = struct_defaults[struct_key]
                with bc2:
                    a_nm = st.number_input(
                        "Lattice parameter a (nm)", 0.10, value=float(a_default),
                        step=0.0001, format="%.4f", key=f"a_{fname}",
                        help="Editable example. BCC Fe≈0.2866 · FCC austenite≈0.358–0.360 · "
                             "Al≈0.4050 (b≈0.286) · Ti(HCP)≈0.295.")
                with bc3:
                    b_from_a = burgers_from_lattice(a_nm, structure)  # metres
                    b_override = st.number_input(
                        "…or Burgers vector b (nm)", 0.05, value=float(b_from_a * 1e9),
                        step=0.001, format="%.4f", key=f"burg_{fname}",
                        help="Auto-filled from a and structure; edit to override directly.")
                b = float(b_override) * 1e-9   # nm → m
                b_nm = b * 1e9

                sc1, sc2, sc3 = st.columns(3)
                with sc1:
                    step_gnd = st.number_input(
                        "Step size u (µm)", 0.001, value=float(step_default),
                        step=0.01, key=f"step_{fname}",
                        help="Acquisition step. Prefilled from the EBSD file, or the "
                             "reference workbook if the file has none. Used directly "
                             "only in the KAM-mean fallback; the gradient method uses "
                             "real neighbour distances (incl. √2·u diagonals).")
                with sc2:
                    noise_deg = st.number_input(
                        "Angular noise floor θ_noise (°)", 0.0, value=0.0,
                        step=0.05, format="%.2f", key=f"noise_{fname}",
                        help="Subtracted from each pair misorientation before the "
                             "gradient. Default 0. As an order-of-magnitude only, "
                             "conventional EBSD angular noise is often ~0.2–0.5° "
                             "(instrument/step dependent — not a fixed constant).")
                with sc3:
                    noise_mode = st.selectbox(
                        "Noise subtraction", ["absolute", "rms"], index=0,
                        key=f"noisemode_{fname}",
                        help="absolute: Δθ' = max(Δθ − θ_noise, 0). "
                             "rms: Δθ' = √(max(Δθ² − θ_noise², 0)).")

                method = st.selectbox(
                    "Method",
                    ["Neighbour-gradient (distance-aware)",
                     "KAM-mean  ρ = 2θ/(α·b·L_eff)",
                     "Regression-corrected gradient"],
                    index=0, key=f"gndmethod_{fname}",
                    help="Neighbour-gradient uses g_i=Δθ_i/r_i with the true distance "
                         "to each neighbour (diagonals √2·u, higher orders). "
                         "Regression fits mean Δθ vs distance and uses the slope "
                         "dθ/du, which reduces angular-noise offset sensitivity.")

                kernel_g = 2 if "Regression" in method or "L_eff" in method else 1
                # allow user to raise kernel order for the gradient/regression
                kernel_g = st.radio(
                    "Kernel order (neighbour shells)", [1, 2], index=0,
                    horizontal=True, key=f"kernord_{fname}",
                    help="1 = 3×3 (8 neighbours: 4 axial at u, 4 diagonal at √2·u). "
                         "2 = 5×5, adds farther shells (2u, √5·u, 2√2·u) — needed "
                         "for the regression slope and distance-corrected variants.")

                # ── Optional calibration factor α (OFF by default) ──────────
                use_alpha = st.checkbox(
                    "Enable optional calibration factor α (advanced)",
                    value=False, key=f"usealpha_{fname}",
                    help="Off by default. The default simplified convention uses NO "
                         "additional calibration factor (α = 1). α is a user-defined "
                         "calibration you must justify — it is NOT a fixed scientific "
                         "standard. It sits in the denominator: ρ = 2θ/(α·b·L_eff).")
                alpha = 1.0
                alpha_note = ""
                if use_alpha:
                    ac1, ac2 = st.columns([1, 2])
                    with ac1:
                        alpha = st.number_input(
                            "α (denominator)", 0.1, value=1.0, step=0.01,
                            key=f"alpha_{fname}",
                            help="ρ = 2θ/(α·b·L_eff). α>1 lowers ρ; α<1 raises it. "
                                 "(Some older workflows used α≈1.86 — this is a "
                                 "convention choice, not a standard.)")
                    with ac2:
                        alpha_note = st.text_input(
                            "Reference / justification for α (required)",
                            value="", key=f"alpharef_{fname}",
                            help="Cite the paper/convention you are following.")
                    if not alpha_note.strip():
                        st.info("α is enabled — please record a reference/justification "
                                "above. Until then it is applied as entered.")

                if not have_map:
                    st.info("No pixel-level map (X/Y + Euler) available, so the "
                            "distance-aware gradient cannot be computed. Falling back "
                            "to the **KAM-mean** convention with L_eff ≈ step size. "
                            "Upload a CTF/CSV pixel map for the gradient/regression "
                            "methods.")

                # ── Compute ────────────────────────────────────────────────
                if have_map:
                    g_df = df[[e1_col, e2_col, e3_col, x_col_g, y_col_g]
                              + ([phase_col] if phase_col and phase_col in df.columns else [])].dropna()
                    ph_arr = g_df[phase_col].values if (phase_col and phase_col in g_df.columns) else None
                    excl_pb = st.checkbox(
                        "Exclude neighbour pairs that cross a phase boundary",
                        value=bool(ph_arr is not None), key=f"exclpb_{fname}",
                        disabled=ph_arr is None,
                        help="Available when a Phase column is mapped. Grain boundaries "
                             "are already excluded by the KAM threshold below.")
                    res = compute_gnd_from_orientations(
                        g_df, e1_col, e2_col, e3_col, x_col_g, y_col_g,
                        b_m=b, kernel_order=int(kernel_g),
                        threshold_deg=float(kam_threshold),
                        noise_deg=float(noise_deg), noise_mode=noise_mode,
                        alpha=float(alpha), phase_arr=ph_arr,
                        exclude_phase_boundaries=bool(excl_pb),
                        step_x_um=float(step_gnd), step_y_um=float(step_gnd))
                    L_eff = res["L_eff_m"]
                    rho_px = res["rho_pixel"]
                    rho_px_raw = res["rho_pixel_raw"]
                    if "Regression" in method and res["regression"]:
                        rho_primary = res["regression"]["rho"]
                        rho_series = pd.Series(rho_px).dropna()
                    elif "KAM-mean" in method:
                        kam_mean_deg = float(np.nanmean(res["kam_deg"]))
                        rho_primary = (2.0 * np.deg2rad(kam_mean_deg)) / (alpha * b * L_eff)
                        rho_series = pd.Series(rho_px).dropna()
                    else:
                        rho_series = pd.Series(rho_px).dropna()
                        rho_primary = float(rho_series.mean())
                    kam_mean_show = float(np.nanmean(res["kam_deg"]))
                else:
                    # KAM-mean fallback on the precomputed KAM column
                    kam_series = df[kam_col_use].dropna()
                    L_eff = float(step_gnd) * 1e-6
                    kam_mean_show = float(kam_series.mean())
                    rho_series = (2.0 * np.deg2rad(kam_series)) / (alpha * b * L_eff)
                    rho_primary = float(rho_series.mean())
                    res = {"frac_pixels_used": np.nan, "frac_excluded_threshold": np.nan,
                           "noise_only_rho": 0.0, "regression": None,
                           "rho_pixel_raw": rho_series.values}

                # ── Result metrics ─────────────────────────────────────────
                m1, m2, m3, m4 = st.columns(4)
                with m1: st.metric("Apparent ρ_GND (m⁻²)", f"{rho_primary:.3e}")
                with m2: st.metric("Median ρ_GND (m⁻²)", f"{float(np.nanmedian(rho_series)):.3e}"
                                   if len(rho_series) else "n/a")
                with m3: st.metric("Mean KAM θ (°)", f"{kam_mean_show:.3f}")
                with m4: st.metric("L_eff (µm)", f"{L_eff*1e6:.4g}")

                # Equation shown with the actual α placement
                if use_alpha:
                    st.latex(r"\rho_{GND}^{KAM} = \frac{2\,\theta_{KAM}}{\alpha\,b\,L_{eff}}"
                             r"\qquad g_i = \frac{\Delta\theta_i}{r_i},\;"
                             r"\rho = \frac{2}{\alpha b}\,\overline{g_i}")
                else:
                    st.latex(r"\rho_{GND}^{KAM} = \frac{2\,\theta_{KAM}}{b\,L_{eff}}"
                             r"\qquad g_i = \frac{\Delta\theta_i}{r_i},\;"
                             r"\rho = \frac{2}{b}\,\overline{g_i}")
                st.caption(
                    ("Default simplified convention: **no additional calibration "
                     "factor (α = 1)**. " if not use_alpha else
                     f"Calibration factor **α = {alpha:g}** applied in the denominator"
                     + (f" (ref: {alpha_note})" if alpha_note.strip() else "") + ". ")
                    + "θ in radians, b and L_eff in metres. The gradient method divides "
                    "each pair by its **real** distance r_i (diagonals = √2·u), not by u alone.")

                # ── Raw vs noise-corrected ─────────────────────────────────
                if noise_deg > 0:
                    rho_raw_mean = float(np.nanmean(res["rho_pixel_raw"]))
                    nn1, nn2, nn3 = st.columns(3)
                    with nn1: st.metric("ρ_GND raw (no noise corr.)", f"{rho_raw_mean:.3e}")
                    with nn2: st.metric("ρ_GND noise-corrected", f"{float(np.nanmean(rho_series)):.3e}"
                                        if len(rho_series) else "n/a")
                    with nn3: st.metric("ρ from noise floor alone", f"{res['noise_only_rho']:.3e}")
                    st.caption(f"Apparent GND produced by θ_noise = {noise_deg:.2f}° alone at "
                               f"L_eff = {L_eff*1e6:.4g} µm is **{res['noise_only_rho']:.3e} m⁻²** — "
                               "if this is comparable to your result, the signal is noise-dominated.")

                with st.expander("🔬 Method & assumptions", expanded=False):
                    st.markdown(f"""
| Quantity | Value used | Note |
|---|---|---|
| Mean KAM θ | **{kam_mean_show:.3f}°** = {np.deg2rad(kam_mean_show):.3e} rad | pairs ≤ threshold only |
| Structure / a | **{structure}**, a = {a_nm:.4f} nm | b = {'√3/2·a' if structure=='BCC' else 'a/√2' if structure=='FCC' else 'a'} |
| Burgers vector b | **{b_nm:.4f} nm** = {b:.3e} m | editable example — set per your phase |
| L_eff | **{L_eff*1e6:.4g} µm** = {L_eff:.3e} m | mean **real** neighbour distance used (not just u) |
| Step size u | **{step_gnd:.4g} µm** | acquisition step (KAM-mean fallback uses L_eff≈u) |
| Kernel / threshold | order {int(kernel_g)}, ≤ {kam_threshold}° | boundaries > {kam_threshold}° excluded |
| Calibration α | **{alpha:g}** ({'enabled' if use_alpha else 'default, none'}) | denominator: 2θ/(α·b·L_eff) |
| Noise floor | **{noise_deg:.2f}°** ({noise_mode}) | Δθ corrected per pair before gradient |
| Pixels used | **{(res['frac_pixels_used']*100):.1f}%** | fraction with ≥1 valid neighbour |
| Pairs excluded by threshold | **{(res['frac_excluded_threshold']*100):.1f}%** | treated as grain boundaries |

**Reference form (Kubin & Mortensen 2003; Calcagnotto et al. 2010):**
ρ_GND ≈ 2·θ/(b·L_eff). This is a **lower-bound / partial proxy**: only GNDs
resolved by the kernel at this step are captured; SSDs are ignored.
""")
                    if res.get("regression"):
                        rg = res["regression"]
                        st.markdown(
                            f"**Regression-corrected gradient:** slope dθ/du = "
                            f"{rg['slope_rad_per_m']:.3e} rad/m, intercept = "
                            f"{rg['intercept_rad']:.3e} rad "
                            f"(≈ angular-noise offset), giving ρ = **{rg['rho']:.3e} m⁻²**. "
                            "The intercept absorbs the noise floor, so the slope-based ρ "
                            "is less sensitive to angular noise than the raw KAM mean.")

                # ── MOS / GOS clarification (fixed) ────────────────────────
                if ref and ref.get("grain", {}).get("stats", {}).get("mos_mean"):
                    mos = ref["grain"]["stats"]["mos_mean"]
                    rel = ("lower than" if kam_mean_show < mos['mean']
                           else "higher than" if kam_mean_show > mos['mean'] else "similar to")
                    st.info(
                        f"**Reference cross-check (not a GND validation).** The workbook "
                        f"Grain List reports **Mean Orientation Spread = {mos['mean']:.3f}°** "
                        f"(median {mos['median']:.3f}°). The app's mean **KAM = "
                        f"{kam_mean_show:.3f}°** is {rel} MOS. These are *different metrics* "
                        f"and need not match: **KAM** is a *local* first/near-neighbour "
                        f"misorientation, while **MOS/GOS** is a *grain-level* spread of "
                        f"orientations about the grain's mean/reference. Neither directly "
                        f"validates a GND density.")

                # ── Sensitivity / uncertainty table ────────────────────────
                if have_map and st.checkbox(
                        "Compute sensitivity / uncertainty table (extra passes — slower)",
                        value=False, key=f"sens_{fname}"):
                    with st.spinner("Running sensitivity conditions…"):
                        conditions = [
                            ("Gradient, kernel 1, thr 2°", dict(kernel_order=1, threshold_deg=2.0, noise_deg=0.0)),
                            ("Gradient, kernel 1, thr 5°", dict(kernel_order=1, threshold_deg=5.0, noise_deg=0.0)),
                            ("Gradient, kernel 2 (dist-corr)", dict(kernel_order=2, threshold_deg=float(kam_threshold), noise_deg=0.0)),
                            ("Raw gradient (current settings)", dict(kernel_order=int(kernel_g), threshold_deg=float(kam_threshold), noise_deg=0.0)),
                        ]
                        if noise_deg > 0:
                            conditions.append((f"Noise-corrected ({noise_deg:.2f}°)",
                                               dict(kernel_order=int(kernel_g), threshold_deg=float(kam_threshold), noise_deg=float(noise_deg))))
                        if use_alpha:
                            conditions.append((f"With α={alpha:g}",
                                               dict(kernel_order=int(kernel_g), threshold_deg=float(kam_threshold), noise_deg=float(noise_deg), alpha=float(alpha))))
                        rows = []
                        for name, kw in conditions:
                            kw.setdefault("alpha", 1.0)
                            r = compute_gnd_from_orientations(
                                g_df, e1_col, e2_col, e3_col, x_col_g, y_col_g,
                                b_m=b, noise_mode=noise_mode, phase_arr=ph_arr,
                                exclude_phase_boundaries=bool(excl_pb),
                                step_x_um=float(step_gnd), step_y_um=float(step_gnd), **kw)
                            rp = pd.Series(r["rho_pixel"]).dropna()
                            if len(rp) == 0:
                                continue
                            rows.append({
                                "Condition": name,
                                "ρ mean (m⁻²)":  f"{rp.mean():.3e}",
                                "ρ median (m⁻²)": f"{rp.median():.3e}",
                                "p10 (m⁻²)": f"{rp.quantile(0.10):.3e}",
                                "p90 (m⁻²)": f"{rp.quantile(0.90):.3e}",
                                "IQR (m⁻²)": f"{(rp.quantile(0.75)-rp.quantile(0.25)):.3e}",
                                "L_eff (µm)": f"{r['L_eff_m']*1e6:.4g}",
                                "Pixels used": f"{r['frac_pixels_used']*100:.1f}%",
                                "Excl. by thr.": f"{r['frac_excluded_threshold']*100:.1f}%",
                            })
                        if rows:
                            st.dataframe(pd.DataFrame(rows), use_container_width=True)
                            st.caption("Each row is a full recomputation. Spread across "
                                       "conditions reflects the method-dependence of the "
                                       "estimate — not experimental uncertainty alone.")

                    # Log-scale histogram of per-pixel apparent GND
                    rp_all = pd.Series(rho_px).dropna()
                    rp_all = rp_all[rp_all > 0]
                    if len(rp_all) > 10:
                        fig_g, ax_g = plt.subplots(figsize=(6, 4))
                        logbins = np.logspace(np.log10(rp_all.min()),
                                              np.log10(rp_all.max()), 50)
                        ax_g.hist(rp_all, bins=logbins, color=PAL["purple"],
                                  edgecolor="white", lw=0.4, alpha=0.85)
                        ax_g.set_xscale("log")
                        ax_g.set_xlabel("Apparent ρ_GND per pixel (m⁻²)")
                        ax_g.set_ylabel("Pixel count")
                        ax_g.set_title("Per-pixel apparent GND density (log scale)")
                        ax_g.axvline(rp_all.median(), color=PAL["red"], ls="--", lw=1.5,
                                     label=f"median = {rp_all.median():.2e}")
                        ax_g.legend()
                        fig_g.tight_layout(); st.pyplot(fig_g)
                        st_figure_download(fig_g, f"gnd_hist_{fname.rsplit('.',1)[0]}",
                                           fmt=plot_fmt.lower(), key=f"dl_gndhist_{fname}")

                # ── Advanced method note ───────────────────────────────────
                with st.expander("🧭 Method 2 — curvature / Nye tensor (advanced, future)",
                                 expanded=False):
                    st.markdown(
                        "A more defensible route than a single averaged KAM is the "
                        "**lattice-curvature / Nye dislocation-density tensor** approach: "
                        "the measured orientation gradients populate components of the Nye "
                        "tensor α_ij, from which a GND density is obtained (e.g. an L1/L2 "
                        "minimisation over candidate slip systems). It is less sensitive to "
                        "arbitrary kernel choices than a scalar KAM average.\n\n"
                        "**However**, 2-D EBSD only exposes the in-plane gradients "
                        "(∂/∂x, ∂/∂y); the out-of-plane terms (∂/∂z) are unmeasured, so even "
                        "the Nye-tensor result on a 2-D map remains **incomplete and a "
                        "lower bound**. Full implementation (3-D or HR-EBSD) is out of scope "
                        "here; this app prioritises an honest, clearly-labelled "
                        "**KAM-derived apparent GND** estimate.")

    # ───────────────────────────────────────────────────────────────────────
    # TAB 6 — POLE FIGURES (IPF)
    # ───────────────────────────────────────────────────────────────────────
    with sec[5]:
        if not (e1_col and e2_col and e3_col):
            st.warning("Select the three Euler angle columns in the Column mapping panel above.")
        else:
            phi1_ipf = df[e1_col].dropna().values
            Phi_ipf  = df[e2_col].dropna().values
            phi2_ipf = df[e3_col].dropna().values
            # Align indices
            valid_ipf = (np.isfinite(phi1_ipf) & np.isfinite(Phi_ipf) & np.isfinite(phi2_ipf))
            phi1_ipf = phi1_ipf[valid_ipf]
            Phi_ipf  = Phi_ipf[valid_ipf]
            phi2_ipf = phi2_ipf[valid_ipf]

            # Phase label
            phase_lbl = "Ferrite (BCC)"
            if phase_col and phase_col in df.columns:
                top_phase = df[phase_col].value_counts().idxmax()
                phase_lbl = str(top_phase)

            # ── IPF settings ──────────────────────────────────────────────
            st.subheader("Inverse Pole Figures (IPF) — MUD Density")
            ipf_c1, ipf_c2, ipf_c3 = st.columns(3)
            with ipf_c1:
                ipf_axes = st.multiselect(
                    "Sample axes to plot",
                    options=["X", "Y", "Z"],
                    default=["X", "Y", "Z"],
                    key=f"ipf_axes_{fname}",
                )
            with ipf_c2:
                ipf_cmap = st.selectbox(
                    "Color map",
                    ["jet", "rainbow", "RdYlBu_r", "plasma", "hot"],
                    index=0,
                    key=f"ipf_cmap_{fname}",
                    help="Color map for MUD density. 'jet' matches the reference figure.",
                )
            with ipf_c3:
                ipf_smooth = st.slider(
                    "Smoothing bandwidth",
                    0.005, 0.08, 0.025, 0.005,
                    key=f"ipf_smooth_{fname}",
                    help="Gaussian smoothing bandwidth as fraction of grid size. "
                         "Larger = smoother density.",
                )

            if ipf_axes:
                with st.spinner("Computing IPF density..."):
                    ipf_title = (
                        f"Inverse Pole Figures — {phase_lbl}\n"
                        f"Directions X, Y and Z"
                    )
                    fig_ipf = plot_ipf_density(
                        phi1_ipf, Phi_ipf, phi2_ipf,
                        axes=tuple(ipf_axes),
                        title=ipf_title,
                        cmap=ipf_cmap,
                        grid_size=300,
                        sigma_frac=ipf_smooth,
                    )
                st.pyplot(fig_ipf)
                st_figure_download(fig_ipf, f"ipf_density_{fname.rsplit('.',1)[0]}",
                                   fmt=plot_fmt.lower(), key=f"dl_ipf_{fname}")
                st.session_state[f"fig_ipf_{fname}"] = fig_ipf
            else:
                st.info("Select at least one sample axis above.")

            st.divider()

            # ── 3-D IPF Cube ───────────────────────────────────────────────
            st.subheader("IPF Color Map — 3D Cube")

            x_col_3d = find_col(df, ["x", "x position", "x_position"])
            y_col_3d = find_col(df, ["y", "y position", "y_position"])

            if not x_col_3d or not y_col_3d:
                st.info(
                    "X and Y position columns not detected. "
                    "The 3-D cube requires spatial coordinates. "
                    "Make sure your file has X and Y columns."
                )
            else:
                rolling_dir = st.selectbox(
                    "Rolling / reference direction (Z of the cube)",
                    ["Z (normal direction)", "Y (rolling direction)", "X (transverse)"],
                    index=0,
                    key=f"ipf3d_roll_{fname}",
                )
                equal_aspect_3d = st.checkbox(
                    "Undistorted cube (equal edges)",
                    value=True,
                    help=(
                        "On (default): the cube is drawn with equal on-screen "
                        "edges so it is not stretched; the true physical extents "
                        "are shown on the X and Y labels. Off: edges are drawn in "
                        "physical proportion (X:Y) and a single common-length "
                        "scale bar (same µm on every axis) is added."
                    ),
                    key=f"ipf3d_equal_{fname}",
                )
                with st.spinner("Rendering 3-D IPF cube... (may take ~20 s for large files)"):
                    # Align Euler and position arrays
                    work = df[[e1_col, e2_col, e3_col, x_col_3d, y_col_3d]].dropna()
                    fig_3d = plot_ipf_3d_cube(
                        work[e1_col].values,
                        work[e2_col].values,
                        work[e3_col].values,
                        work[x_col_3d].values,
                        work[y_col_3d].values,
                        title=f"IPF Color Map — {phase_lbl} (3-D)",
                        phase_label=phase_lbl,
                        equal_aspect=equal_aspect_3d,
                    )
                st.pyplot(fig_3d)
                st_figure_download(fig_3d, f"ipf_3d_cube_{fname.rsplit('.',1)[0]}",
                                   fmt=plot_fmt.lower(), key=f"dl_ipf3d_{fname}")
                st.session_state[f"fig_ipf3d_{fname}"] = fig_3d
                st.caption(
                    "Face colors encode crystal orientation relative to each sample axis "
                    "(001 = red · 101 = green · 111 = blue). "
                    "Top face = IPF-Z, front face = IPF-Y, right face = IPF-X. "
                    "The **Undistorted cube** option (default on) draws equal on-screen "
                    "edges so the box is not stretched; the true physical extents are "
                    "printed on the X/Y labels. Turn it off to draw the box in physical "
                    "proportion (X:Y) with a common-length scale bar (same µm on every axis)."
                )

            st.divider()

            # ── Pole Figures (PF) ──────────────────────────────────────────
            st.subheader("Pole Figures (PF) — Stereographic Projections")
            st.caption(
                "Each circle shows the MUD density of crystal plane normals in "
                "the sample reference frame (X = horizontal, Y = vertical). "
                "Colormap: jet (0 = weak texture → high = strong texture)."
            )

            pf_c1, pf_c2 = st.columns(2)
            with pf_c1:
                pf_cmap = st.selectbox(
                    "PF color map",
                    ["jet", "rainbow", "RdYlBu_r", "plasma", "hot"],
                    index=0,
                    key=f"pf_cmap_{fname}",
                )
            with pf_c2:
                pf_smooth = st.slider(
                    "PF smoothing bandwidth",
                    0.005, 0.08, 0.02, 0.005,
                    key=f"pf_smooth_{fname}",
                    help="Gaussian smoothing bandwidth as fraction of grid size.",
                )

            with st.spinner("Computing Pole Figures…"):
                try:
                    fig_pf = plot_pole_figures(
                        phi1_ipf, Phi_ipf, phi2_ipf,
                        planes=((1, 0, 0), (1, 1, 0), (1, 1, 1)),
                        plane_labels=("{100}", "{110}", "{111}"),
                        title=f"Pole Figures — {phase_lbl}",
                        cmap=pf_cmap,
                        sigma_frac=pf_smooth,
                        ref_dir_x="X",
                        ref_dir_y="Y",
                    )
                    st.pyplot(fig_pf)
                    st_figure_download(
                        fig_pf,
                        f"pole_figures_{fname.rsplit('.',1)[0]}",
                        fmt=plot_fmt.lower(),
                        key=f"dl_pf_{fname}",
                    )
                    st.session_state[f"fig_pf_{fname}"] = fig_pf
                except Exception as _pf_err:
                    st.warning(f"Pole figure error: {_pf_err}")

            st.divider()

            # ── IPF 2D Color Map ───────────────────────────────────────────
            st.subheader("IPF Color Map — 2D Grain Orientation Map")
            st.caption(
                "Flat microstructure map: each pixel colored by IPF color for the "
                "selected sample axis (X/Y/Z). Reference frame uses X and Y axes."
            )

            x_col_2d = find_col(df, ["x", "x position", "x_position"])
            y_col_2d = find_col(df, ["y", "y position", "y_position"])

            if not x_col_2d or not y_col_2d:
                st.info(
                    "X and Y position columns not detected. "
                    "The 2-D map requires spatial coordinates. "
                    "Make sure your file has X and Y columns."
                )
            else:
                map_c1, map_c2 = st.columns(2)
                with map_c1:
                    map_axis = st.selectbox(
                        "IPF reference axis",
                        ["Z", "X", "Y"],
                        index=0,
                        key=f"map_axis_{fname}",
                        help="Crystal direction parallel to this sample axis defines the IPF color.",
                    )
                with map_c2:
                    map_scalebar = st.number_input(
                        "Scale bar length (µm, 0 = auto)",
                        min_value=0.0, value=0.0, step=1.0,
                        key=f"map_sb_{fname}",
                    )
                    map_scalebar = None if map_scalebar == 0 else float(map_scalebar)

                with st.spinner("Rendering 2-D IPF map…"):
                    try:
                        work2d = df[[e1_col, e2_col, e3_col,
                                     x_col_2d, y_col_2d]].dropna()
                        fig_map2d = plot_ipf_2d_map(
                            work2d[e1_col].values,
                            work2d[e2_col].values,
                            work2d[e3_col].values,
                            work2d[x_col_2d].values,
                            work2d[y_col_2d].values,
                            sample_axis=map_axis,
                            title=f"IPF Color Map // {map_axis} — {phase_lbl}",
                            phase_label=phase_lbl,
                            ref_axis_h="X",
                            ref_axis_v="Y",
                            scale_bar_um=map_scalebar,
                        )
                        st.pyplot(fig_map2d)
                        st_figure_download(
                            fig_map2d,
                            f"ipf_2d_map_{fname.rsplit('.',1)[0]}",
                            fmt=plot_fmt.lower(),
                            key=f"dl_map2d_{fname}",
                        )
                        st.session_state[f"fig_map2d_{fname}"] = fig_map2d
                    except Exception as _map_err:
                        st.warning(f"2-D IPF map error: {_map_err}")

    # ───────────────────────────────────────────────────────────────────────
    # TAB 7 — EXPORT
    # ───────────────────────────────────────────────────────────────────────
    with sec[6]:
        fmt = plot_fmt.lower()
        base = fname.rsplit('.', 1)[0]
        figs_export = {k.replace(f"_{fname}", ""): v
                       for k, v in st.session_state.items()
                       if k.endswith(f"_{fname}") and hasattr(v, "savefig")}

        # ── Build CSV exports ──────────────────────────────────────────────
        csv_exports = {}
        if grain_df is not None and not grain_df.empty:
            csv_exports["grain_summary"] = grain_df.to_csv(index=False)
        # Always export the pixel/point dataframe
        try:
            csv_exports["ebsd_data"] = df.to_csv(index=False)
        except Exception:
            pass
        # Euler angle texture summary (ideal fractions)
        if e1_col and e2_col and e3_col:
            try:
                IDEALS_EX = {
                    "{001}<100> Cube":        (0,  0,  0),
                    "{110}<001> Goss":        (0,  45, 0),
                    "{111}<110> γ-fiber":     (0,  55, 45),
                    "{111}<112> γ-fiber":     (30, 55, 45),
                    "{112}<110> Brass (BCC)": (35, 35, 0),
                    "{001}<110> Rotated Cube":(45, 0,  0),
                }
                phi1_ex = df[e1_col].dropna().values
                Phi_ex  = df[e2_col].dropna().values
                phi2_ex = df[e3_col].dropna().values
                rows_ex = []
                for iname, (p1i, Pi, p2i) in IDEALS_EX.items():
                    d1 = np.abs(phi1_ex - p1i) % 360; d1 = np.minimum(d1, 360 - d1)
                    dP = np.abs(Phi_ex  - Pi)
                    d2 = np.abs(phi2_ex - p2i) % 360; d2 = np.minimum(d2, 360 - d2)
                    near = ((d1 <= 15) & (dP <= 15) & (d2 <= 15)).sum()
                    rows_ex.append({"Orientation": iname, "Count": int(near),
                                    "Fraction (%)": round(100 * near / max(len(phi1_ex), 1), 3)})
                csv_exports["texture_summary"] = pd.DataFrame(rows_ex).to_csv(index=False)
            except Exception:
                pass

        # ── Full ZIP export ────────────────────────────────────────────────
        has_content = bool(figs_export) or bool(csv_exports)
        if not has_content:
            st.info("Generate plots in the other tabs first, then return here to export.")
        else:
            zbuf = io.BytesIO()
            with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as zf:
                # Figures
                for nm, fg in figs_export.items():
                    zf.writestr(f"figures/{nm}.{fmt}", fig_bytes(fg, fmt))
                # CSVs
                for cnm, cdata in csv_exports.items():
                    zf.writestr(f"data/{cnm}_{base}.csv", cdata)
            zbuf.seek(0)
            st.download_button(
                f"⬇️ Download complete package — figures ({fmt.upper()}) + CSV data",
                zbuf.getvalue(),
                file_name=f"EBSD_{base}_export.zip",
                mime="application/zip",
                use_container_width=True,
            )

            if figs_export:
                st.divider()
                st.subheader("Individual figure downloads")
                cols_dl = st.columns(2)
                for i, (nm, fg) in enumerate(figs_export.items()):
                    with cols_dl[i % 2]:
                        st.download_button(
                            f"⬇️ {nm}.{fmt}",
                            fig_bytes(fg, fmt),
                            file_name=f"{nm}.{fmt}",
                            key=f"dl_exp_{nm}_{fname}",
                        )

            if csv_exports:
                st.divider()
                st.subheader("Individual CSV downloads")
                for cnm, cdata in csv_exports.items():
                    st.download_button(
                        f"⬇️ {cnm}_{base}.csv",
                        cdata,
                        file_name=f"{cnm}_{base}.csv",
                        mime="text/csv",
                        key=f"dl_csv_{cnm}_{fname}",
                    )

        # Grain table preview
        if grain_df is not None and not grain_df.empty:
            st.divider()
            st.subheader("Grain summary preview (first 100 rows)")
            st.dataframe(grain_df.head(100), use_container_width=True)
