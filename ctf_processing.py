"""
CTF-specific processing: grain segmentation, misorientation, KAM
from pixel-level EBSD data (Oxford/HKL CTF format).

CTF files contain one row per scan pixel — NOT per grain.
All microstructural quantities (grain size, misorientation, KAM)
must be derived from the pixel data.
"""

import numpy as np
import pandas as pd
from scipy.ndimage import label as nd_label


# ── Misorientation between two orientations (Euler → rotation matrix) ────────

def euler_to_matrix(phi1_deg, Phi_deg, phi2_deg):
    """
    Bunge convention (ZXZ): φ₁, Φ, φ₂ → 3×3 rotation matrix.
    Accepts scalar or array inputs.
    """
    p1  = np.deg2rad(phi1_deg)
    P   = np.deg2rad(Phi_deg)
    p2  = np.deg2rad(phi2_deg)

    cp1, sp1 = np.cos(p1), np.sin(p1)
    cP,  sP  = np.cos(P),  np.sin(P)
    cp2, sp2 = np.cos(p2), np.sin(p2)

    # Shape: (..., 3, 3)
    R = np.stack([
        np.stack([ cp1*cp2 - sp1*sp2*cP,  sp1*cp2 + cp1*sp2*cP,  sp2*sP], axis=-1),
        np.stack([-cp1*sp2 - sp1*cp2*cP, -sp1*sp2 + cp1*cp2*cP,  cp2*sP], axis=-1),
        np.stack([ sp1*sP,               -cp1*sP,                  cP    ], axis=-1),
    ], axis=-2)
    return R


def misorientation_angle(R1, R2):
    """
    Minimum misorientation angle (degrees) between two rotation matrices.
    R1, R2: shape (..., 3, 3)
    Uses cubic symmetry operators (24 operators).
    """
    # Relative rotation
    dR = np.matmul(R1, R2.swapaxes(-1, -2))  # R1 · R2^T

    # Cubic symmetry operators (24)
    sqrt2 = np.sqrt(2) / 2
    sym_ops = np.array([
        [[1,0,0],[0,1,0],[0,0,1]],
        [[-1,0,0],[0,-1,0],[0,0,1]],
        [[-1,0,0],[0,1,0],[0,0,-1]],
        [[1,0,0],[0,-1,0],[0,0,-1]],
        [[0,0,1],[1,0,0],[0,1,0]],
        [[0,0,1],[-1,0,0],[0,-1,0]],
        [[0,0,-1],[-1,0,0],[0,1,0]],
        [[0,0,-1],[1,0,0],[0,-1,0]],
        [[0,1,0],[0,0,1],[1,0,0]],
        [[0,-1,0],[0,0,1],[-1,0,0]],
        [[0,1,0],[0,0,-1],[-1,0,0]],
        [[0,-1,0],[0,0,-1],[1,0,0]],
        [[0,1,0],[1,0,0],[0,0,-1]],
        [[0,-1,0],[-1,0,0],[0,0,-1]],
        [[0,1,0],[-1,0,0],[0,0,1]],
        [[0,-1,0],[1,0,0],[0,0,1]],
        [[1,0,0],[0,0,1],[0,-1,0]],
        [[-1,0,0],[0,0,1],[0,1,0]],
        [[-1,0,0],[0,0,-1],[0,-1,0]],
        [[1,0,0],[0,0,-1],[0,1,0]],
        [[0,0,1],[0,1,0],[-1,0,0]],
        [[0,0,1],[0,-1,0],[1,0,0]],
        [[0,0,-1],[0,1,0],[1,0,0]],
        [[0,0,-1],[0,-1,0],[-1,0,0]],
    ], dtype=np.float64)  # shape (24, 3, 3)

    # For each symmetry op, compute trace of (sym · dR)
    # dR shape: (..., 3, 3)
    # Broadcast: (24, 3, 3) · (..., 3, 3) → (..., 24, 3, 3)
    dR_exp  = dR[..., np.newaxis, :, :]           # (..., 1, 3, 3)
    sym_exp = sym_ops[np.newaxis, ...]             # (1, 24, 3, 3) if 1D
    # Need explicit einsum
    # sym_dR[..., s, i, j] = sum_k sym_ops[s,i,k] * dR[..., k, j]
    sym_dR = np.einsum('sik,...kj->...sij', sym_ops, dR)  # (..., 24, 3, 3)
    traces = np.einsum('...sii->...s', sym_dR)             # (..., 24)
    traces_clipped = np.clip(traces, -1.0, 3.0)
    angles = np.arccos(np.clip((traces_clipped - 1.0) / 2.0, -1.0, 1.0))  # (..., 24)
    min_angle = np.min(angles, axis=-1)  # (...,)
    return np.rad2deg(min_angle)


# ── KAM from pixel map ────────────────────────────────────────────────────────

def compute_kam(df: pd.DataFrame,
                phi1_col: str, Phi_col: str, phi2_col: str,
                x_col: str = "X", y_col: str = "Y",
                kernel_order: int = 1,
                threshold_deg: float = 5.0) -> np.ndarray:
    """
    Compute Kernel Average Misorientation (KAM) for each pixel.

    Parameters
    ----------
    df           : CTF DataFrame (one row per pixel)
    phi1_col     : column name for Euler φ₁
    Phi_col      : column name for Euler Φ
    phi2_col     : column name for Euler φ₂
    x_col, y_col : position columns
    kernel_order : 1 = nearest neighbours only (3×3 kernel)
    threshold_deg: ignore pairs with misorientation > threshold (default 5°)

    Returns
    -------
    kam_values : np.ndarray, length = len(df), in degrees
    """
    if x_col not in df.columns or y_col not in df.columns:
        return np.full(len(df), np.nan)

    xs = df[x_col].values
    ys = df[y_col].values

    # Determine step size
    ux = np.sort(np.unique(xs))
    uy = np.sort(np.unique(ys))
    if len(ux) < 2 or len(uy) < 2:
        return np.full(len(df), np.nan)
    step_x = float(np.median(np.diff(ux)))
    step_y = float(np.median(np.diff(uy)))

    # Build index map: (row, col) → dataframe index
    x_idx = np.round((xs - xs.min()) / step_x).astype(int)
    y_idx = np.round((ys - ys.min()) / step_y).astype(int)
    nrows = y_idx.max() + 1
    ncols = x_idx.max() + 1

    grid = np.full((nrows, ncols), -1, dtype=int)
    for df_i, (xi, yi) in enumerate(zip(x_idx, y_idx)):
        grid[yi, xi] = df_i

    # Precompute rotation matrices (vectorized)
    phi1 = df[phi1_col].values
    Phi  = df[Phi_col].values
    phi2 = df[phi2_col].values
    R_all = euler_to_matrix(phi1, Phi, phi2)  # shape (N, 3, 3)

    kam_values = np.zeros(len(df))
    count_arr  = np.zeros(len(df), dtype=int)

    # 8-neighbourhood (kernel_order=1)
    offsets = [(dr, dc)
               for dr in range(-kernel_order, kernel_order + 1)
               for dc in range(-kernel_order, kernel_order + 1)
               if not (dr == 0 and dc == 0)]

    for dr, dc in offsets:
        # Shift grid
        r0 = max(0, -dr);  r1 = nrows - max(0, dr)
        c0 = max(0, -dc);  c1 = ncols - max(0, dc)
        r0n = max(0, dr);  r1n = nrows - max(0, -dr)
        c0n = max(0, dc);  c1n = ncols - max(0, -dc)

        center_idx  = grid[r0:r1, c0:c1].ravel()
        neigh_idx   = grid[r0n:r1n, c0n:c1n].ravel()

        valid = (center_idx >= 0) & (neigh_idx >= 0)
        ci = center_idx[valid]
        ni = neigh_idx[valid]

        if len(ci) == 0:
            continue

        angles = misorientation_angle(R_all[ci], R_all[ni])

        # Apply threshold
        below = angles <= threshold_deg
        ci_b = ci[below]
        ang_b = angles[below]

        np.add.at(kam_values, ci_b, ang_b)
        np.add.at(count_arr,  ci_b, 1)

    with np.errstate(invalid='ignore'):
        kam_out = np.where(count_arr > 0, kam_values / count_arr, 0.0)
    return kam_out


# ── Burgers vector from lattice parameter ─────────────────────────────────────

def burgers_from_lattice(a_nm: float, structure: str = "BCC") -> float:
    """
    Shortest lattice Burgers vector length (in metres) from the cubic lattice
    parameter ``a_nm`` (nm).

        BCC: b = (√3 / 2) · a   (½<111> slip)
        FCC: b = a / √2         (½<110> slip)

    Any other/unknown ``structure`` falls back to the FCC expression and the
    caller should treat the value as an editable example, not a fixed truth.
    """
    a_m = float(a_nm) * 1e-9
    s = (structure or "").strip().upper()
    if s.startswith("BCC"):
        return (np.sqrt(3.0) / 2.0) * a_m
    if s.startswith("FCC"):
        return a_m / np.sqrt(2.0)
    if s.startswith("HCP"):
        # For HCP the a-parameter itself is the basal <a> Burgers vector.
        return a_m
    return a_m / np.sqrt(2.0)


# ── KAM-derived apparent GND density (distance-aware) ─────────────────────────

def compute_gnd_from_orientations(df: pd.DataFrame,
                                   phi1_col: str, Phi_col: str, phi2_col: str,
                                   x_col: str, y_col: str,
                                   b_m: float,
                                   kernel_order: int = 1,
                                   threshold_deg: float = 5.0,
                                   noise_deg: float = 0.0,
                                   noise_mode: str = "absolute",
                                   alpha: float = 1.0,
                                   phase_arr=None,
                                   exclude_phase_boundaries: bool = False,
                                   step_x_um: float = None,
                                   step_y_um: float = None) -> dict:
    """
    Estimate the **KAM-derived apparent GND density** from pixel-level
    orientation data using the real centre→neighbour distances (including
    diagonals √2·u and higher-order offsets), instead of dividing a single
    averaged KAM by the axial step size only.

    Per valid neighbour pair the local orientation gradient is

        g_i = Δθ_i / r_i            (Δθ in radians, r_i in metres)

    and the per-pixel apparent GND density is

        ρ_GND ≈ (2 / (α·b)) · mean_i(g_i)

    which reduces to ρ = 2·θ / (α·b·L_eff) for the simple KAM-mean variant,
    where L_eff is the mean included neighbour distance (NOT necessarily u).

    Noise correction (when ``noise_deg`` > 0) is applied per pair before the
    gradient:  absolute → Δθ' = max(Δθ − θ_noise, 0);
               rms      → Δθ' = sqrt(max(Δθ² − θ_noise², 0)).

    Returns a dict of arrays/scalars (see keys assembled at the end). All ρ
    values are in m⁻². This is a lower-bound / partial proxy: it captures only
    the GNDs resolved by the kernel at this step and ignores SSDs.
    """
    N = len(df)
    nan_ret = {
        "rho_pixel": np.full(N, np.nan),
        "rho_pixel_raw": np.full(N, np.nan),
        "kam_deg": np.full(N, np.nan),
        "n_used": np.zeros(N, dtype=int),
        "L_eff_m": np.nan,
        "frac_excluded_threshold": np.nan,
        "frac_pixels_used": 0.0,
        "regression": None,
        "noise_only_rho": np.nan,
        "step_x_um": np.nan, "step_y_um": np.nan,
    }
    if x_col not in df.columns or y_col not in df.columns:
        return nan_ret

    xs = df[x_col].values
    ys = df[y_col].values
    ux = np.sort(np.unique(xs))
    uy = np.sort(np.unique(ys))
    if len(ux) < 2 or len(uy) < 2:
        return nan_ret

    step_x = float(step_x_um) if step_x_um else float(np.median(np.diff(ux)))
    step_y = float(step_y_um) if step_y_um else float(np.median(np.diff(uy)))
    sxm = step_x * 1e-6
    sym = step_y * 1e-6

    x_idx = np.round((xs - xs.min()) / step_x).astype(int)
    y_idx = np.round((ys - ys.min()) / step_y).astype(int)
    nrows = y_idx.max() + 1
    ncols = x_idx.max() + 1

    grid = np.full((nrows, ncols), -1, dtype=int)
    grid[y_idx, x_idx] = np.arange(N)

    R_all = euler_to_matrix(df[phi1_col].values, df[Phi_col].values, df[phi2_col].values)

    ph = None
    if phase_arr is not None and exclude_phase_boundaries:
        ph = np.asarray(phase_arr)

    noise_rad = np.deg2rad(max(0.0, float(noise_deg)))

    grad_sum      = np.zeros(N)   # noise-corrected Σ g_i (rad/m)
    grad_sum_raw  = np.zeros(N)   # raw Σ g_i (rad/m)
    grad_cnt      = np.zeros(N, dtype=int)
    kam_sum       = np.zeros(N)   # Σ Δθ (deg), below threshold
    r_sum         = np.zeros(N)   # Σ r_i (m), for L_eff
    total_pairs   = np.zeros(N, dtype=int)
    excl_thresh   = np.zeros(N, dtype=int)

    # Regression accumulators keyed by rounded distance (m)
    reg_ang = {}   # dist_m -> Σ Δθ (rad, raw, below threshold)
    reg_cnt = {}

    offsets = [(dr, dc)
               for dr in range(-kernel_order, kernel_order + 1)
               for dc in range(-kernel_order, kernel_order + 1)
               if not (dr == 0 and dc == 0)]

    for dr, dc in offsets:
        r_i = np.sqrt((dr * sym) ** 2 + (dc * sxm) ** 2)  # metres
        if r_i <= 0:
            continue
        r0 = max(0, -dr);  r1 = nrows - max(0, dr)
        c0 = max(0, -dc);  c1 = ncols - max(0, dc)
        r0n = max(0, dr);  r1n = nrows - max(0, -dr)
        c0n = max(0, dc);  c1n = ncols - max(0, -dc)

        center_idx = grid[r0:r1, c0:c1].ravel()
        neigh_idx  = grid[r0n:r1n, c0n:c1n].ravel()
        valid = (center_idx >= 0) & (neigh_idx >= 0)
        ci = center_idx[valid]; ni = neigh_idx[valid]
        if len(ci) == 0:
            continue
        if ph is not None:
            same_phase = ph[ci] == ph[ni]
            ci = ci[same_phase]; ni = ni[same_phase]
            if len(ci) == 0:
                continue

        angles = misorientation_angle(R_all[ci], R_all[ni])  # degrees
        np.add.at(total_pairs, ci, 1)

        below = angles <= threshold_deg
        np.add.at(excl_thresh, ci[~below], 1)
        ci_b = ci[below]; ang_b = angles[below]
        if len(ci_b) == 0:
            continue

        ang_rad = np.deg2rad(ang_b)
        # Noise correction
        if noise_rad > 0:
            if noise_mode == "rms":
                ang_corr = np.sqrt(np.clip(ang_rad**2 - noise_rad**2, 0.0, None))
            else:
                ang_corr = np.clip(ang_rad - noise_rad, 0.0, None)
        else:
            ang_corr = ang_rad

        np.add.at(grad_sum,     ci_b, ang_corr / r_i)
        np.add.at(grad_sum_raw, ci_b, ang_rad  / r_i)
        np.add.at(grad_cnt,     ci_b, 1)
        np.add.at(kam_sum,      ci_b, ang_b)
        np.add.at(r_sum,        ci_b, r_i)

        key = round(r_i, 15)
        reg_ang[key] = reg_ang.get(key, 0.0) + float(ang_rad.sum())
        reg_cnt[key] = reg_cnt.get(key, 0)   + int(len(ang_rad))

    with np.errstate(invalid='ignore', divide='ignore'):
        mean_grad     = np.where(grad_cnt > 0, grad_sum / grad_cnt, np.nan)
        mean_grad_raw = np.where(grad_cnt > 0, grad_sum_raw / grad_cnt, np.nan)
        kam_deg       = np.where(grad_cnt > 0, kam_sum / grad_cnt, np.nan)

    rho_pixel     = (2.0 / (alpha * b_m)) * mean_grad
    rho_pixel_raw = (2.0 / (alpha * b_m)) * mean_grad_raw

    L_eff = float(np.nansum(r_sum) / max(1, int(grad_cnt.sum())))
    frac_excl = float(excl_thresh.sum()) / max(1, int(total_pairs.sum()))
    frac_used = float((grad_cnt > 0).sum()) / max(1, N)

    # Apparent GND from the noise floor alone at this L_eff
    noise_only_rho = float((2.0 / (alpha * b_m)) * (noise_rad / L_eff)) if (noise_rad > 0 and L_eff > 0) else 0.0

    # ── Regression method: mean Δθ(r) vs r → slope dθ/du ──────────────────────
    regression = None
    if len(reg_cnt) >= 2:
        dists = np.array(sorted(reg_ang.keys()))
        mean_dtheta = np.array([reg_ang[d] / reg_cnt[d] for d in dists])  # rad
        A = np.vstack([dists, np.ones_like(dists)]).T
        slope, intercept = np.linalg.lstsq(A, mean_dtheta, rcond=None)[0]
        rho_reg = float((2.0 / (alpha * b_m)) * slope) if slope > 0 else 0.0
        regression = {
            "distances_m": dists,
            "mean_dtheta_rad": mean_dtheta,
            "slope_rad_per_m": float(slope),
            "intercept_rad": float(intercept),
            "rho": rho_reg,
        }

    return {
        "rho_pixel": rho_pixel,
        "rho_pixel_raw": rho_pixel_raw,
        "kam_deg": kam_deg,
        "n_used": grad_cnt,
        "L_eff_m": L_eff,
        "frac_excluded_threshold": frac_excl,
        "frac_pixels_used": frac_used,
        "regression": regression,
        "noise_only_rho": noise_only_rho,
        "step_x_um": step_x, "step_y_um": step_y,
    }


# ── Grain segmentation from pixel map ────────────────────────────────────────

def segment_grains(df: pd.DataFrame,
                   phi1_col: str, Phi_col: str, phi2_col: str,
                   x_col: str = "X", y_col: str = "Y",
                   threshold_deg: float = 10.0) -> np.ndarray:
    """
    Assign a grain ID to each pixel using flood-fill segmentation.
    Two adjacent pixels belong to the same grain if their misorientation < threshold_deg.

    Returns grain_ids array (length = len(df)).
    """
    if x_col not in df.columns or y_col not in df.columns:
        return np.arange(len(df))

    xs = df[x_col].values
    ys = df[y_col].values

    ux = np.sort(np.unique(xs))
    uy = np.sort(np.unique(ys))
    if len(ux) < 2 or len(uy) < 2:
        return np.arange(len(df))

    step_x = float(np.median(np.diff(ux)))
    step_y = float(np.median(np.diff(uy)))

    x_idx = np.round((xs - xs.min()) / step_x).astype(int)
    y_idx = np.round((ys - ys.min()) / step_y).astype(int)
    nrows = y_idx.max() + 1
    ncols = x_idx.max() + 1

    grid = np.full((nrows, ncols), -1, dtype=int)
    for df_i, (xi, yi) in enumerate(zip(x_idx, y_idx)):
        grid[yi, xi] = df_i

    phi1 = df[phi1_col].values
    Phi  = df[Phi_col].values
    phi2 = df[phi2_col].values
    R_all = euler_to_matrix(phi1, Phi, phi2)

    # Build boundary map: pixel is a boundary if it differs > threshold from right or bottom neighbour
    is_boundary = np.zeros((nrows, ncols), dtype=bool)

    for dr, dc in [(0, 1), (1, 0)]:
        r0 = 0; r1 = nrows - dr
        c0 = 0; c1 = ncols - dc
        ci_grid = grid[r0:r1, c0:c1]
        ni_grid = grid[r0+dr:r1+dr, c0+dc:c1+dc]
        valid = (ci_grid >= 0) & (ni_grid >= 0)
        rows_v, cols_v = np.where(valid)
        ci = ci_grid[rows_v, cols_v]
        ni = ni_grid[rows_v, cols_v]
        if len(ci) == 0:
            continue
        angles = misorientation_angle(R_all[ci], R_all[ni])
        is_gb = angles >= threshold_deg
        is_boundary[rows_v[is_gb], cols_v[is_gb]] = True
        is_boundary[rows_v[is_gb] + dr, cols_v[is_gb] + dc] = True

    # Label connected regions (interior pixels)
    interior = ~is_boundary & (grid >= 0)
    labeled, n_grains = nd_label(interior)

    # Assign grain IDs back to df indices
    grain_ids = np.zeros(len(df), dtype=int)
    for row in range(nrows):
        for col in range(ncols):
            df_i = grid[row, col]
            if df_i >= 0:
                grain_ids[df_i] = labeled[row, col]

    return grain_ids


def compute_grain_stats(df: pd.DataFrame,
                        grain_ids: np.ndarray,
                        x_col: str = "X", y_col: str = "Y",
                        step_size_um: float = None) -> pd.DataFrame:
    """
    Given pixel-level data and grain IDs, compute per-grain statistics:
    - Grain area (µm²)
    - Grain diameter (ECD, µm)
    - Mean Euler angles per grain
    - Mean MAD per grain
    - Mean Band Contrast per grain
    - Number of pixels

    Returns a DataFrame with one row per grain (grain_id > 0).
    """
    df2 = df.copy()
    df2["__grain_id__"] = grain_ids

    # Estimate step size if not provided
    if step_size_um is None:
        xs = df[x_col].dropna().values if x_col in df.columns else None
        if xs is not None and len(xs) > 1:
            ux = np.sort(np.unique(xs))
            if len(ux) > 1:
                step_size_um = float(np.median(np.diff(ux)))
            else:
                step_size_um = 1.0
        else:
            step_size_um = 1.0

    pixel_area = step_size_um ** 2  # µm² per pixel

    rows = []
    for gid, grp in df2[df2["__grain_id__"] > 0].groupby("__grain_id__"):
        n_pix = len(grp)
        area  = n_pix * pixel_area
        diam  = 2.0 * np.sqrt(area / np.pi)  # ECD

        row = {
            "Grain ID":       gid,
            "Pixel Count":    n_pix,
            "Grain Area":     round(area, 4),
            "Grain Diameter": round(diam, 4),
        }
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                row[f"Mean {col}"] = round(grp[col].mean(), 4)
        rows.append(row)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def compute_grain_misorientation(df: pd.DataFrame,
                                 grain_stats: pd.DataFrame,
                                 phi1_col: str, Phi_col: str, phi2_col: str) -> pd.Series:
    """
    Compute average grain-to-grain misorientation for each grain
    relative to the mean orientation of all grains.
    Returns a Series indexed like grain_stats.
    """
    mean_p1_col = f"Mean {phi1_col}"
    mean_P_col  = f"Mean {Phi_col}"
    mean_p2_col = f"Mean {phi2_col}"

    if not all(c in grain_stats.columns for c in [mean_p1_col, mean_P_col, mean_p2_col]):
        return pd.Series(np.nan, index=grain_stats.index)

    R_grains = euler_to_matrix(
        grain_stats[mean_p1_col].values,
        grain_stats[mean_P_col].values,
        grain_stats[mean_p2_col].values,
    )  # (N_grains, 3, 3)

    # Mean orientation across all grains (simple average of rotation matrices, then re-orthogonalise)
    R_mean = R_grains.mean(axis=0)
    U, _, Vt = np.linalg.svd(R_mean)
    R_mean_ortho = U @ Vt  # nearest rotation matrix

    R_mean_rep = np.tile(R_mean_ortho[np.newaxis], (len(grain_stats), 1, 1))
    angles = misorientation_angle(R_grains, R_mean_rep)
    return pd.Series(angles, index=grain_stats.index, name="Misorientation Angle")
