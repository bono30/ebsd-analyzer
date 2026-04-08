"""
Pole Figure (PF) and IPF 2D Color Map for EBSD Analyzer.

New functions (do NOT modify ipf_plots.py):
  - plot_pole_figures     : circular stereographic projections with MUD density
  - plot_ipf_2d_map       : 2-D grain orientation map colored by IPF (with scale bar)
  - fig_download_button   : helper — Streamlit download button for any matplotlib figure
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch
from scipy.ndimage import gaussian_filter
import io
import warnings
warnings.filterwarnings("ignore")

# Re-use from ipf_plots (already installed in same folder)
from ipf_plots import (
    _euler_to_matrix_vec,
    stereographic_projection,
    ipf_color_bcc,
    plot_ipf_legend,
    _P001, _P101, _P111,
)


# ══════════════════════════════════════════════════════════════════════════════
#  CUBIC SYMMETRY OPERATORS  (24 proper rotations)
# ══════════════════════════════════════════════════════════════════════════════
_SYM = np.array([
    [[ 1, 0, 0],[ 0, 1, 0],[ 0, 0, 1]],
    [[-1, 0, 0],[ 0,-1, 0],[ 0, 0, 1]],
    [[-1, 0, 0],[ 0, 1, 0],[ 0, 0,-1]],
    [[ 1, 0, 0],[ 0,-1, 0],[ 0, 0,-1]],
    [[ 0, 0, 1],[ 1, 0, 0],[ 0, 1, 0]],
    [[ 0, 0, 1],[-1, 0, 0],[ 0,-1, 0]],
    [[ 0, 0,-1],[-1, 0, 0],[ 0, 1, 0]],
    [[ 0, 0,-1],[ 1, 0, 0],[ 0,-1, 0]],
    [[ 0, 1, 0],[ 0, 0, 1],[ 1, 0, 0]],
    [[ 0,-1, 0],[ 0, 0, 1],[-1, 0, 0]],
    [[ 0, 1, 0],[ 0, 0,-1],[-1, 0, 0]],
    [[ 0,-1, 0],[ 0, 0,-1],[ 1, 0, 0]],
    [[ 0, 1, 0],[ 1, 0, 0],[ 0, 0,-1]],
    [[ 0,-1, 0],[-1, 0, 0],[ 0, 0,-1]],
    [[ 0, 1, 0],[-1, 0, 0],[ 0, 0, 1]],
    [[ 0,-1, 0],[ 1, 0, 0],[ 0, 0, 1]],
    [[ 1, 0, 0],[ 0, 0, 1],[ 0,-1, 0]],
    [[-1, 0, 0],[ 0, 0, 1],[ 0, 1, 0]],
    [[-1, 0, 0],[ 0, 0,-1],[ 0,-1, 0]],
    [[ 1, 0, 0],[ 0, 0,-1],[ 0, 1, 0]],
    [[ 0, 0, 1],[ 0, 1, 0],[-1, 0, 0]],
    [[ 0, 0, 1],[ 0,-1, 0],[ 1, 0, 0]],
    [[ 0, 0,-1],[ 0, 1, 0],[ 1, 0, 0]],
    [[ 0, 0,-1],[ 0,-1, 0],[-1, 0, 0]],
], dtype=np.float64)   # (24, 3, 3)


# ══════════════════════════════════════════════════════════════════════════════
#  POLE FIGURE HELPER: crystal planes → sample frame
# ══════════════════════════════════════════════════════════════════════════════

def _hkl_to_unit(hkl):
    """Normalise a (h,k,l) plane normal to unit vector."""
    v = np.asarray(hkl, dtype=float)
    return v / np.linalg.norm(v)


def crystal_to_sample_poles(phi1_arr, Phi_arr, phi2_arr, hkl=(1, 0, 0)):
    """
    For each orientation, rotate all 24 symmetry-equivalent {hkl} directions
    from crystal frame to sample frame.

    Returns
    -------
    sample_dirs : (N * n_sym_equiv, 3)  unit vectors in sample frame
    """
    R = _euler_to_matrix_vec(
        np.asarray(phi1_arr, float),
        np.asarray(Phi_arr,  float),
        np.asarray(phi2_arr, float),
    )   # (N, 3, 3)

    v0 = _hkl_to_unit(hkl)  # (3,)

    # Apply symmetry to get all equivalent crystal directions: (24, 3)
    equiv = np.einsum('sij,j->si', _SYM, v0)   # (24, 3)

    # Rotate each equivalent direction to sample frame:
    # sample_dir = R^T · v_crystal  (R maps crystal→sample, so R^T maps sample→crystal inverse)
    # Actually: R maps crystal frame to sample frame, so v_sample = R · v_crystal
    # shape: (N, 24, 3)
    # v_sample[n, s] = R[n] @ equiv[s]
    dirs = np.einsum('nij,sj->nsi', R, equiv)    # (N, 24, 3)
    dirs = dirs.reshape(-1, 3)                    # (N*24, 3)

    # Normalise
    norms = np.linalg.norm(dirs, axis=1, keepdims=True)
    dirs  = dirs / np.where(norms > 0, norms, 1.0)

    # Keep only upper hemisphere (z >= 0), flip negative-z poles
    dirs[dirs[:, 2] < 0] *= -1
    return dirs


def _pf_density_circle(x2d, y2d, grid_size=300, sigma_frac=0.02):
    """
    Compute MUD density on a unit-circle grid for a pole figure.
    Returns (mud_grid, extent).
    """
    pad = 0.05
    xmin, xmax = -1.0 - pad, 1.0 + pad
    ymin, ymax = -1.0 - pad, 1.0 + pad

    mask = ((x2d >= xmin) & (x2d <= xmax) &
            (y2d >= ymin) & (y2d <= ymax))
    xm, ym = x2d[mask], y2d[mask]
    n_total = max(len(xm), 1)

    H, _, _ = np.histogram2d(xm, ym, bins=grid_size,
                              range=[[xmin, xmax], [ymin, ymax]])
    H = H.T.astype(float)

    sigma = max(1.0, grid_size * sigma_frac)
    H_smooth = gaussian_filter(H, sigma=sigma)

    total_cells = grid_size * grid_size
    uniform_count = n_total / total_cells
    with np.errstate(invalid='ignore', divide='ignore'):
        mud = H_smooth / (uniform_count + 1e-9)
    mud = np.where(mud <= 0, np.nan, mud)

    return mud, [xmin, xmax, ymin, ymax]


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN: POLE FIGURES (PF)
# ══════════════════════════════════════════════════════════════════════════════

def plot_pole_figures(phi1_arr, Phi_arr, phi2_arr,
                      planes=((1, 0, 0), (1, 1, 0), (1, 1, 1)),
                      plane_labels=("{100}", "{110}", "{111}"),
                      title="Pole Figures",
                      cmap="jet",
                      grid_size=300,
                      sigma_frac=0.02,
                      ref_dir_x="X",
                      ref_dir_y="Y"):
    """
    Plot pole figures (equal-area / stereographic projection) for given planes.

    Each PF is a filled circle with MUD density coloring and a colorbar.
    Reference directions X and Y are labeled inside the circle.

    Parameters
    ----------
    phi1_arr, Phi_arr, phi2_arr : Euler angles (degrees)
    planes      : list of (h,k,l) tuples
    plane_labels: display labels for each plane
    title       : figure suptitle
    cmap        : colormap
    grid_size   : density grid resolution
    sigma_frac  : smoothing bandwidth
    ref_dir_x   : label for horizontal axis (default "X")
    ref_dir_y   : label for vertical axis   (default "Y")

    Returns
    -------
    fig : matplotlib Figure
    """
    phi1 = np.asarray(phi1_arr, float)
    Phi  = np.asarray(Phi_arr,  float)
    phi2 = np.asarray(phi2_arr, float)
    valid = np.isfinite(phi1) & np.isfinite(Phi) & np.isfinite(phi2)
    phi1, Phi, phi2 = phi1[valid], Phi[valid], phi2[valid]

    n_pf = len(planes)
    fig, axs = plt.subplots(1, n_pf, figsize=(4.5 * n_pf, 5.0))
    if n_pf == 1:
        axs = [axs]
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)

    for ax, hkl, label in zip(axs, planes, plane_labels):
        # 1. Compute poles in sample frame
        dirs = crystal_to_sample_poles(phi1, Phi, phi2, hkl)

        # 2. Stereographic projection
        x2d, y2d = stereographic_projection(dirs)

        # 3. MUD density
        mud, extent = _pf_density_circle(x2d, y2d, grid_size=grid_size,
                                          sigma_frac=sigma_frac)

        # 4. Circular mask
        gx = np.linspace(extent[0], extent[1], grid_size)
        gy = np.linspace(extent[2], extent[3], grid_size)
        GX, GY = np.meshgrid(gx, gy)
        in_circle = (GX**2 + GY**2) <= 1.0

        # Fill NaN holes inside circle with 0
        if mud is not None:
            mud = np.where(np.isnan(mud) & in_circle, 0.0, mud)
        mud_masked = np.where(in_circle, mud, np.nan)

        # 5. Color scale
        valid_v = mud_masked[in_circle]
        vmax = float(np.nanpercentile(valid_v, 99)) if len(valid_v) > 0 else 1.0
        vmax = max(vmax, 1.0)

        # 6. Plot density
        im = ax.imshow(mud_masked, origin="lower", aspect="equal",
                       extent=extent, cmap=cmap,
                       vmin=0, vmax=vmax, interpolation="bilinear")

        # 7. Draw outer circle
        theta = np.linspace(0, 2 * np.pi, 360)
        ax.plot(np.cos(theta), np.sin(theta), color="black", lw=1.2, zorder=5)

        # 8. Draw cross lines (reference directions)
        ax.axhline(0, color="black", lw=0.6, ls="-", alpha=0.4, zorder=4)
        ax.axvline(0, color="black", lw=0.6, ls="-", alpha=0.4, zorder=4)

        # 9. Reference direction labels
        ax.text( 1.03, 0.0,  ref_dir_x, ha="left",   va="center", fontsize=10, fontweight="bold")
        ax.text( 0.0,  1.06, ref_dir_y, ha="center", va="bottom", fontsize=10, fontweight="bold")

        # 10. Colorbar (vertical, with MUD levels matching reference figure style)
        cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.05, shrink=0.85)
        n_ticks = 8
        ticks = np.linspace(0, vmax, n_ticks)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f"{v:.2f}" for v in ticks])
        cbar.set_label("MUD", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(label, fontsize=12, fontweight="bold", pad=6)

    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  IPF 2D COLOR MAP  (flat microstructure map colored by IPF)
# ══════════════════════════════════════════════════════════════════════════════

def plot_ipf_2d_map(phi1_arr, Phi_arr, phi2_arr,
                    x_pos, y_pos,
                    sample_axis="Z",
                    title="IPF Color Map",
                    phase_label="Ferrite (BCC)",
                    ref_axis_h="X",
                    ref_axis_v="Y",
                    scale_bar_um=None):
    """
    Plot a 2-D map of the scan area colored by IPF color for the given sample axis.
    Similar to the reference image (flat microstructure map with grain colors).

    Parameters
    ----------
    phi1_arr, Phi_arr, phi2_arr : Euler angles (degrees)
    x_pos, y_pos : pixel positions in µm
    sample_axis  : 'X', 'Y', or 'Z' — which sample axis defines the IPF color
    title        : figure title
    phase_label  : phase name for legend
    ref_axis_h   : label for horizontal reference direction (default 'X')
    ref_axis_v   : label for vertical reference direction   (default 'Y')
    scale_bar_um : length of scale bar in µm (auto-computed if None)

    Returns
    -------
    fig : matplotlib Figure
    """
    phi1 = np.asarray(phi1_arr, float)
    Phi  = np.asarray(Phi_arr,  float)
    phi2 = np.asarray(phi2_arr, float)
    xp   = np.asarray(x_pos, float)
    yp   = np.asarray(y_pos, float)

    valid = (np.isfinite(phi1) & np.isfinite(Phi) & np.isfinite(phi2)
             & np.isfinite(xp) & np.isfinite(yp))
    phi1, Phi, phi2 = phi1[valid], Phi[valid], phi2[valid]
    xp, yp = xp[valid], yp[valid]

    # Reconstruct pixel grid
    ux = np.sort(np.unique(xp))
    uy = np.sort(np.unique(yp))
    sx = float(np.median(np.diff(ux))) if len(ux) > 1 else 1.0
    sy = float(np.median(np.diff(uy))) if len(uy) > 1 else 1.0
    xi = np.round((xp - xp.min()) / sx).astype(int)
    yi = np.round((yp - yp.min()) / sy).astype(int)
    ncols = xi.max() + 1
    nrows = yi.max() + 1

    # Compute IPF colors
    rgb = ipf_color_bcc(
        __import__('ipf_plots').euler_to_crystal_direction(phi1, Phi, phi2, sample_axis)
    )   # (N, 3)

    # Build image
    img = np.ones((nrows, ncols, 3), dtype=float) * 0.15   # dark background
    img[yi, xi] = rgb

    # Figure layout: map on left, legend+axes on right
    fig = plt.figure(figsize=(8, 6))
    ax_map = fig.add_axes([0.05, 0.05, 0.65, 0.88])

    # Physical extent
    x_um = ncols * sx
    y_um = nrows * sy
    ax_map.imshow(img, origin="lower",
                  extent=[0, x_um, 0, y_um],
                  aspect="equal", interpolation="nearest")

    # ── Scale bar ──────────────────────────────────────────────────────────
    if scale_bar_um is None:
        # Auto: ~15% of map width, rounded to nice number
        raw = x_um * 0.15
        mag = 10 ** np.floor(np.log10(raw + 1e-9))
        scale_bar_um = float(round(raw / mag) * mag)
        scale_bar_um = max(scale_bar_um, sx)

    bar_x0 = x_um * 0.05
    bar_y0 = y_um * 0.04
    bar_h  = y_um * 0.012
    rect = mpatches.FancyBboxPatch(
        (bar_x0, bar_y0), scale_bar_um, bar_h,
        boxstyle="square,pad=0", facecolor="white", edgecolor="white", lw=0)
    ax_map.add_patch(rect)
    ax_map.text(bar_x0 + scale_bar_um / 2, bar_y0 + bar_h * 2.2,
                f"{scale_bar_um:.4g} µm",
                ha="center", va="bottom", fontsize=9,
                color="white", fontweight="bold")

    ax_map.set_xlabel(f"{ref_axis_h} (µm)", fontsize=10)
    ax_map.set_ylabel(f"{ref_axis_v} (µm)", fontsize=10)
    ax_map.set_title(title, fontsize=12, fontweight="bold", pad=6)
    ax_map.tick_params(labelsize=9)

    # ── Reference direction axes (top right inset) ─────────────────────────
    ax_ref = fig.add_axes([0.74, 0.72, 0.22, 0.18])
    kw = dict(arrowstyle="-|>", color="black", lw=1.3,
              mutation_scale=10)
    ax_ref.annotate("", xy=(0.85, 0.5),  xytext=(0.15, 0.5),  arrowprops=kw)
    ax_ref.annotate("", xy=(0.5,  0.90), xytext=(0.5,  0.10), arrowprops=kw)
    ax_ref.text(0.90, 0.50, ref_axis_h, ha="left",   va="center",
                fontsize=11, fontweight="bold")
    ax_ref.text(0.50, 0.95, ref_axis_v, ha="center", va="bottom",
                fontsize=11, fontweight="bold")
    ax_ref.set_xlim(0, 1); ax_ref.set_ylim(0, 1)
    ax_ref.axis("off")

    # ── IPF legend triangle ────────────────────────────────────────────────
    ax_leg = fig.add_axes([0.72, 0.35, 0.26, 0.30])
    plot_ipf_legend(ax=ax_leg, title=f"{phase_label}\n001", fontsize=8)

    # ── IPF axis label ─────────────────────────────────────────────────────
    fig.text(0.84, 0.31, f"IPF // {sample_axis}",
             ha="center", va="top", fontsize=9, style="italic")

    fig.patch.set_facecolor("white")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  STREAMLIT HELPER: per-figure download button
# ══════════════════════════════════════════════════════════════════════════════

def st_figure_download(fig, filename_stem, fmt="png", label=None, key=None):
    """
    Render a Streamlit download button for a matplotlib figure.
    Call this immediately after st.pyplot(fig).

    Parameters
    ----------
    fig           : matplotlib Figure
    filename_stem : base filename without extension
    fmt           : 'png', 'svg', or 'pdf'
    label         : button label (auto-generated if None)
    key           : Streamlit widget key (auto-generated if None)
    """
    import streamlit as st
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=300, bbox_inches="tight")
    buf.seek(0)

    mime_map = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}
    mime = mime_map.get(fmt, "application/octet-stream")
    btn_label = label or f"⬇️ Download  {filename_stem}.{fmt}"
    widget_key = key or f"dl_inline_{filename_stem}_{fmt}"

    st.download_button(
        label=btn_label,
        data=buf.getvalue(),
        file_name=f"{filename_stem}.{fmt}",
        mime=mime,
        key=widget_key,
        use_container_width=False,
    )
