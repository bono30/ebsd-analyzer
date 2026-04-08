"""
Inverse Pole Figure (IPF) module for EBSD Analyzer.

Provides:
  - euler_to_crystal_direction : rotate sample axis into crystal frame
  - ipf_color_bcc              : map crystal direction to RGB (BCC triangle)
  - ipf_triangle_density       : plot IPF with MUD density (Multiples of Uniform Distribution)
  - ipf_color_map_3d           : 3-D cube with IPF-colored faces (X, Y, Z)
  - ipf_legend_triangle        : small legend showing the BCC color key

References:
  Bunge, H.J. (1982) Texture Analysis in Materials Science.
  Nolze & Hielscher (2016) J. Appl. Cryst. 49, 1786–1802.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings("ignore")


# ══════════════════════════════════════════════════════════════════════════════
#  EULER → ROTATION MATRIX  (Bunge ZXZ, same convention as ctf_processing)
# ══════════════════════════════════════════════════════════════════════════════
def _euler_to_matrix_vec(phi1, Phi, phi2):
    """Vectorised Euler (Bunge ZXZ) → rotation matrices. Shape: (N,3,3)."""
    p1 = np.deg2rad(phi1); P = np.deg2rad(Phi); p2 = np.deg2rad(phi2)
    cp1, sp1 = np.cos(p1), np.sin(p1)
    cP,  sP  = np.cos(P),  np.sin(P)
    cp2, sp2 = np.cos(p2), np.sin(p2)
    R = np.stack([
        np.stack([ cp1*cp2 - sp1*sp2*cP,  sp1*cp2 + cp1*sp2*cP,  sp2*sP], axis=-1),
        np.stack([-cp1*sp2 - sp1*cp2*cP, -sp1*sp2 + cp1*cp2*cP,  cp2*sP], axis=-1),
        np.stack([ sp1*sP,               -cp1*sP,                  cP   ], axis=-1),
    ], axis=-2)
    return R   # (N,3,3)


# ══════════════════════════════════════════════════════════════════════════════
#  SAMPLE AXIS → CRYSTAL DIRECTION  (applying rotation matrices)
# ══════════════════════════════════════════════════════════════════════════════
_SAMPLE_AXES = {
    "X": np.array([1., 0., 0.]),
    "Y": np.array([0., 1., 0.]),
    "Z": np.array([0., 0., 1.]),
}

def euler_to_crystal_direction(phi1_arr, Phi_arr, phi2_arr, sample_axis="Z"):
    """
    For each pixel/grain, rotate the given sample axis into the crystal frame.
    Returns unit vectors in crystal coordinates, shape (N, 3).

    g_crystal = R · v_sample
    R is the orientation matrix (crystal → sample), so we use R^T = R^-1
    to go sample → crystal.
    """
    R = _euler_to_matrix_vec(np.asarray(phi1_arr, float),
                              np.asarray(Phi_arr,  float),
                              np.asarray(phi2_arr, float))  # (N,3,3)
    v = _SAMPLE_AXES[sample_axis.upper()]                   # (3,)
    # crystal direction = R^T · v  →  einsum 'nij,j->ni'
    cryst = np.einsum('nji,j->ni', R, v)                    # (N,3) — R^T via swapped indices
    # Normalise
    norms = np.linalg.norm(cryst, axis=1, keepdims=True)
    cryst = cryst / np.where(norms > 0, norms, 1.0)
    # Map to upper hemisphere (all z ≥ 0 by convention)
    cryst[cryst[:, 2] < 0] *= -1
    return cryst


# ══════════════════════════════════════════════════════════════════════════════
#  STEREOGRAPHIC PROJECTION  (upper hemisphere → 2-D)
# ══════════════════════════════════════════════════════════════════════════════
def stereographic_projection(xyz):
    """
    Projects unit vectors from upper hemisphere onto the equatorial plane.
    Returns (x2d, y2d) arrays.
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    z = np.clip(z, -1 + 1e-9, 1.0)
    denom = 1.0 + z
    return x / denom, y / denom


# ══════════════════════════════════════════════════════════════════════════════
#  BCC FUNDAMENTAL ZONE  (stereographic triangle 001-101-111)
# ══════════════════════════════════════════════════════════════════════════════
# Corners in crystal coordinates (unit vectors):
_C001 = np.array([0., 0., 1.])
_C101 = np.array([1., 0., 1.]) / np.sqrt(2)
_C111 = np.array([1., 1., 1.]) / np.sqrt(3)

# Stereographic projections of the three corners:
def _sp(v):
    return v[0] / (1 + v[2]), v[1] / (1 + v[2])

_P001 = _sp(_C001)   # (0, 0)
_P101 = _sp(_C101)   # (sqrt2/(1+sqrt2/sqrt2), 0) ≈ (0.414, 0)
_P111 = _sp(_C111)   # both coords equal

# Pre-compute triangle boundary for masking
def _in_fundamental_triangle(x2d, y2d):
    """
    Return boolean mask: True if the stereographic point lies inside the
    BCC fundamental triangle (001)-(101)-(111).
    """
    # Three edges defined by great circles; approximate using barycentric coords
    # Vertices:
    v0 = np.array(_P001)
    v1 = np.array(_P101)
    v2 = np.array(_P111)

    def sign(p, a, b):
        return (p[0]-b[0])*(a[1]-b[1]) - (a[0]-b[0])*(p[1]-b[1])

    pts = np.stack([x2d, y2d], axis=1)  # (N,2)
    d1 = (pts[:,0]-v1[0])*(v0[1]-v1[1]) - (v0[0]-v1[0])*(pts[:,1]-v1[1])
    d2 = (pts[:,0]-v2[0])*(v1[1]-v2[1]) - (v1[0]-v2[0])*(pts[:,1]-v2[1])
    d3 = (pts[:,0]-v0[0])*(v2[1]-v0[1]) - (v2[0]-v0[0])*(pts[:,1]-v0[1])

    has_neg = (d1 < 0) | (d2 < 0) | (d3 < 0)
    has_pos = (d1 > 0) | (d2 > 0) | (d3 > 0)
    return ~(has_neg & has_pos)


# ══════════════════════════════════════════════════════════════════════════════
#  IPF COLOR  (RGB mapping inside BCC triangle)
# ══════════════════════════════════════════════════════════════════════════════
def ipf_color_bcc(xyz):
    """
    Map an array of crystal unit vectors to IPF RGB colors (BCC).
    Standard color key:
        001 → red   [1, 0, 0]
        101 → green [0, 1, 0]
        111 → blue  [0, 0, 1]

    Uses barycentric interpolation inside the fundamental triangle.

    Parameters
    ----------
    xyz : (N, 3) unit vectors in crystal frame (upper hemisphere)

    Returns
    -------
    rgb : (N, 3) float array in [0, 1]
    """
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]

    # Reduce to fundamental triangle by applying cubic symmetry:
    # Sort components so that |x| ≤ |y| ≤ |z| and then into the sector z ≥ y ≥ x ≥ 0
    coords = np.sort(np.abs(xyz), axis=1)   # ascending: (a, b, c) with a≤b≤c
    a, b, c = coords[:, 0], coords[:, 1], coords[:, 2]

    # Normalise so that c = 1 (max component)
    with np.errstate(divide='ignore', invalid='ignore'):
        scale = np.where(c > 0, 1.0 / c, 1.0)
    a, b, c = a*scale, b*scale, np.ones_like(a)

    # Barycentric weights for triangle 001-101-111
    # 001: (a=0, b=0)  101: (a=0, b=1)  111: (a=1, b=1)
    # w001 = 1 - b
    # w101 = b - a
    # w111 = a
    w001 = np.clip(1.0 - b, 0, 1)
    w101 = np.clip(b - a,   0, 1)
    w111 = np.clip(a,        0, 1)
    total = w001 + w101 + w111 + 1e-9
    w001, w101, w111 = w001/total, w101/total, w111/total

    # Colors: 001=red, 101=green, 111=blue
    R = w001 * 1.0 + w101 * 0.0 + w111 * 0.0
    G = w001 * 0.0 + w101 * 1.0 + w111 * 0.0
    B = w001 * 0.0 + w101 * 0.0 + w111 * 1.0

    # Boost saturation slightly (as in most EBSD software)
    rgb = np.stack([R, G, B], axis=1)
    mn  = rgb.min(axis=1, keepdims=True)
    rgb = rgb - mn
    mx  = rgb.max(axis=1, keepdims=True)
    rgb = np.where(mx > 0, rgb / mx, rgb)

    return np.clip(rgb, 0, 1)


# ══════════════════════════════════════════════════════════════════════════════
#  MUD (Multiples of Uniform Distribution) DENSITY
# ══════════════════════════════════════════════════════════════════════════════
def _compute_mud(x2d, y2d, grid_size=256, sigma_frac=0.03):
    """
    Estimate MUD (Multiples of Uniform Distribution) on a regular grid.
    MUD = local density / random (uniform) density.
    MUD = 1.0 means random texture; MUD > 1 means preferred orientation.
    """
    # Grid bounds from triangle vertices with small padding
    x_all = [_P001[0], _P101[0], _P111[0]]
    y_all = [_P001[1], _P101[1], _P111[1]]
    pad = 0.015
    xmin, xmax = min(x_all) - pad, max(x_all) + pad
    ymin, ymax = min(y_all) - pad, max(y_all) + pad

    # Keep only points inside the bounding box
    mask = ((x2d >= xmin) & (x2d <= xmax) &
            (y2d >= ymin) & (y2d <= ymax))
    xm, ym = x2d[mask], y2d[mask]
    n_total = len(xm)
    if n_total == 0:
        return np.ones((grid_size, grid_size)), [xmin, xmax, ymin, ymax]

    # Build histogram
    H, _, _ = np.histogram2d(xm, ym, bins=grid_size,
                              range=[[xmin, xmax], [ymin, ymax]])
    H = H.T.astype(float)  # (grid_size, grid_size), rows=y, cols=x

    # Smooth with Gaussian
    sigma = max(1.0, grid_size * sigma_frac)
    H_smooth = gaussian_filter(H, sigma=sigma)

    # MUD = observed density / expected uniform density
    # Expected count per cell = n_total / (number of cells inside triangle)
    # Approximate: use total cells as denominator for simplicity
    total_cells = grid_size * grid_size
    uniform_count_per_cell = n_total / total_cells
    with np.errstate(invalid='ignore', divide='ignore'):
        mud = H_smooth / (uniform_count_per_cell + 1e-9)

    # Values ≤ 0 → NaN (will be masked anyway by triangle mask)
    mud = np.where(mud <= 0, np.nan, mud)
    return mud, [xmin, xmax, ymin, ymax]


# ══════════════════════════════════════════════════════════════════════════════
#  TRIANGLE PATCH  (draw the fundamental triangle as a filled polygon)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_triangle_frame(ax, linewidth=1.2, color="black"):
    """Draw the boundary of the BCC fundamental triangle on ax."""
    pts = [_P001, _P101, _P111, _P001]
    xs  = [p[0] for p in pts]
    ys  = [p[1] for p in pts]
    ax.plot(xs, ys, color=color, lw=linewidth, zorder=10)


def _annotate_corners(ax, fontsize=11):
    """Label the three corners of the fundamental triangle."""
    ax.text(_P001[0], _P001[1] - 0.018, "(001)",
            ha="center", va="top",  fontsize=fontsize, fontweight="bold")
    ax.text(_P101[0], _P101[1] - 0.018, "(101)",
            ha="center", va="top",  fontsize=fontsize, fontweight="bold")
    ax.text(_P111[0], _P111[1] + 0.018, "(111)",
            ha="center", va="bottom", fontsize=fontsize, fontweight="bold")


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN PUBLIC FUNCTION: IPF DENSITY PLOT (X, Y, Z)
# ══════════════════════════════════════════════════════════════════════════════
def plot_ipf_density(phi1_arr, Phi_arr, phi2_arr,
                     axes=("X", "Y", "Z"),
                     title="Inverse Pole Figures — Ferrite (BCC)",
                     cmap="jet",
                     grid_size=300,
                     sigma_frac=0.025):
    """
    Plot IPF density (MUD) triangles for the requested sample axes.

    Parameters
    ----------
    phi1_arr, Phi_arr, phi2_arr : array-like, Euler angles in degrees
    axes     : tuple of sample axes to plot, e.g. ("X","Y","Z")
    title    : figure title
    cmap     : matplotlib colormap name
    grid_size: density grid resolution
    sigma_frac: smoothing bandwidth

    Returns
    -------
    fig : matplotlib Figure
    """
    n_axes = len(axes)
    fig, axs = plt.subplots(1, n_axes, figsize=(5.5 * n_axes, 5.5))
    if n_axes == 1:
        axs = [axs]
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.02)

    phi1 = np.asarray(phi1_arr, float)
    Phi  = np.asarray(Phi_arr,  float)
    phi2 = np.asarray(phi2_arr, float)

    # Remove NaN rows
    valid = np.isfinite(phi1) & np.isfinite(Phi) & np.isfinite(phi2)
    phi1, Phi, phi2 = phi1[valid], Phi[valid], phi2[valid]

    for ax, sample_axis in zip(axs, axes):
        # 1. Crystal directions
        cryst = euler_to_crystal_direction(phi1, Phi, phi2, sample_axis)

        # 2. Stereographic projection
        x2d, y2d = stereographic_projection(cryst)

        # 3. MUD grid
        mud, extent = _compute_mud(x2d, y2d, grid_size=grid_size,
                                   sigma_frac=sigma_frac)

        # 4. Build triangle mask
        gx = np.linspace(extent[0], extent[1], grid_size)
        gy = np.linspace(extent[2], extent[3], grid_size)
        GX, GY = np.meshgrid(gx, gy)
        in_tri = _in_fundamental_triangle(GX.ravel(), GY.ravel()).reshape(grid_size, grid_size)

        # Replace NaN holes INSIDE the triangle with 0 (empty cells = 0 MUD)
        if mud is not None:
            mud = np.where(np.isnan(mud) & in_tri, 0.0, mud)

        # Mask outside triangle
        mud_masked = np.where(in_tri, mud, np.nan)

        # 5. Plot
        valid_vals = mud_masked[in_tri]
        vmax = float(np.nanpercentile(valid_vals, 99)) if len(valid_vals) > 0 else 1.0
        vmax = max(vmax, 1.0)
        im = ax.imshow(mud_masked, origin="lower", aspect="equal",
                       extent=extent, cmap=cmap,
                       vmin=0, vmax=vmax, interpolation="bilinear")

        # 6. Triangle border and corner labels
        _draw_triangle_frame(ax)
        _annotate_corners(ax, fontsize=10)

        # 7. Colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.85)
        # Show MUD levels like the reference figure
        levels = np.linspace(0, vmax, 10)
        cbar.set_ticks(levels)
        cbar.set_ticklabels([f"{v:.2f}" for v in levels])
        cbar.set_label("MUD", fontsize=9)

        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2] - 0.02, extent[3] + 0.03)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(sample_axis, fontsize=13, fontweight="bold", loc="left", pad=8)

    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  IPF COLOR MAP PER PIXEL
# ══════════════════════════════════════════════════════════════════════════════
def compute_ipf_colormap(phi1_arr, Phi_arr, phi2_arr, sample_axis="Z"):
    """
    Return per-pixel/grain IPF RGB color for the given sample axis.

    Parameters
    ----------
    phi1_arr, Phi_arr, phi2_arr : Euler angles in degrees (length N)
    sample_axis : 'X', 'Y', or 'Z'

    Returns
    -------
    rgb : (N, 3) float array in [0, 1]
    """
    phi1 = np.asarray(phi1_arr, float)
    Phi  = np.asarray(Phi_arr,  float)
    phi2 = np.asarray(phi2_arr, float)
    cryst = euler_to_crystal_direction(phi1, Phi, phi2, sample_axis)
    return ipf_color_bcc(cryst)


# ══════════════════════════════════════════════════════════════════════════════
#  IPF COLOR LEGEND TRIANGLE
# ══════════════════════════════════════════════════════════════════════════════
def plot_ipf_legend(ax=None, title="", fontsize=9):
    """
    Draw the BCC IPF color key (001=red, 101=green, 111=blue) on ax.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(2.5, 2.5))
    else:
        fig = ax.get_figure()

    # Sample a dense grid of stereographic points inside the triangle
    gx = np.linspace(_P001[0] - 0.01, _P111[0] + 0.02, 300)
    gy = np.linspace(_P001[1] - 0.01, _P111[1] + 0.02, 300)
    GX, GY = np.meshgrid(gx, gy)
    pts2d = np.stack([GX.ravel(), GY.ravel()], axis=1)

    # Inverse stereographic to get unit vectors
    r2 = pts2d[:, 0]**2 + pts2d[:, 1]**2
    z  = (1 - r2) / (1 + r2)
    scale = 1.0 / (1 + r2 + 1e-9)
    x  = 2 * pts2d[:, 0] * scale * (1 + r2) / 2   # simplified
    x  = pts2d[:, 0] * (1 + z)
    y  = pts2d[:, 1] * (1 + z)
    xyz = np.stack([x, y, z], axis=1)
    nrm = np.linalg.norm(xyz, axis=1, keepdims=True)
    xyz = xyz / np.where(nrm > 0, nrm, 1.0)

    # IPF colors
    rgb = ipf_color_bcc(xyz)

    # Mask to triangle
    in_tri = _in_fundamental_triangle(pts2d[:, 0], pts2d[:, 1])
    img_r = np.where(in_tri, rgb[:, 0], np.nan)
    img_g = np.where(in_tri, rgb[:, 1], np.nan)
    img_b = np.where(in_tri, rgb[:, 2], np.nan)

    # Combine to RGBA
    alpha = np.where(in_tri, 1.0, 0.0)
    img = np.stack([
        img_r.reshape(300, 300),
        img_g.reshape(300, 300),
        img_b.reshape(300, 300),
        alpha.reshape(300, 300),
    ], axis=-1)
    img = np.nan_to_num(img, nan=0.0)

    extent = [gx.min(), gx.max(), gy.min(), gy.max()]
    ax.imshow(img, origin="lower", extent=extent, aspect="equal",
              interpolation="bilinear")
    _draw_triangle_frame(ax, linewidth=1.0)
    _annotate_corners(ax, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize, fontweight="bold")
    ax.axis("off")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
#  3-D IPF COLOR MAP CUBE
# ══════════════════════════════════════════════════════════════════════════════
def plot_ipf_3d_cube(phi1_arr, Phi_arr, phi2_arr,
                     x_pos, y_pos,
                     step_x=None, step_y=None,
                     title="IPF Color Map — 3D",
                     phase_label="Ferrite"):
    """
    Generate a 3-D cube figure with three IPF-colored faces (X, Y, Z planes),
    similar to the reference image showing grain microstructure in 3-D.

    The three visible faces of the cube correspond to:
      - Top  face (XY plane) → IPF colors for Z axis
      - Front face (XZ plane) → IPF colors for Y axis
      - Right face (YZ plane) → IPF colors for X axis

    Parameters
    ----------
    phi1_arr, Phi_arr, phi2_arr : Euler angles (degrees), length N
    x_pos, y_pos : pixel X and Y positions (µm)
    step_x, step_y : grid step sizes (µm); estimated if None
    title : figure title
    phase_label : label for the legend triangle

    Returns
    -------
    fig : matplotlib Figure
    """
    phi1 = np.asarray(phi1_arr, float)
    Phi  = np.asarray(Phi_arr,  float)
    phi2 = np.asarray(phi2_arr, float)
    xp   = np.asarray(x_pos, float)
    yp   = np.asarray(y_pos, float)

    # ── Rebuild pixel grid ────────────────────────────────────────────────────
    valid = np.isfinite(phi1) & np.isfinite(Phi) & np.isfinite(phi2) \
          & np.isfinite(xp)  & np.isfinite(yp)
    phi1, Phi, phi2 = phi1[valid], Phi[valid], phi2[valid]
    xp, yp = xp[valid], yp[valid]

    ux = np.sort(np.unique(xp))
    uy = np.sort(np.unique(yp))
    sx = float(np.median(np.diff(ux))) if step_x is None and len(ux) > 1 else (step_x or 1.0)
    sy = float(np.median(np.diff(uy))) if step_y is None and len(uy) > 1 else (step_y or 1.0)

    xi = np.round((xp - xp.min()) / sx).astype(int)
    yi = np.round((yp - yp.min()) / sy).astype(int)
    ncols = xi.max() + 1
    nrows = yi.max() + 1

    # ── Compute IPF colors for all three axes ─────────────────────────────────
    rgb_z = compute_ipf_colormap(phi1, Phi, phi2, "Z")   # top face
    rgb_y = compute_ipf_colormap(phi1, Phi, phi2, "Y")   # front face
    rgb_x = compute_ipf_colormap(phi1, Phi, phi2, "X")   # right face

    def make_image(rgb, row_idx, col_idx, nr, nc):
        """Place N RGB values on a (nr, nc, 3) image grid."""
        img = np.ones((nr, nc, 3), dtype=float) * 0.9
        img[row_idx, col_idx] = rgb
        return img

    img_z = make_image(rgb_z, yi, xi, nrows, ncols)   # top  (rows=Y, cols=X)
    img_y = make_image(rgb_y, yi, xi, nrows, ncols)   # front (rows=Y, cols=X)
    img_x = make_image(rgb_x, yi, xi, nrows, ncols)   # right (rows=Y, cols=X)

    # ── Build 3-D figure ──────────────────────────────────────────────────────
    fig = plt.figure(figsize=(9, 8))
    ax3 = fig.add_axes([0.0, 0.05, 0.78, 0.90], projection="3d")

    W = ncols - 1  # width  (X direction)
    H = nrows - 1  # height (Y direction)
    D = max(W, H)  # depth  (Z direction, use same scale for aesthetics)

    # ── Top face (Z plane): IPF-Z ─────────────────────────────────────────────
    # Surface at z = D; x varies 0→W, y varies 0→H
    xs_t = np.linspace(0, W, ncols)
    ys_t = np.linspace(0, H, nrows)
    XS, YS = np.meshgrid(xs_t, ys_t)
    ZS = np.full_like(XS, D)

    # imshow → face colors need (nr-1, nc-1, 4) or use plot_surface with facecolors
    # We use plot_surface with facecolors
    fc_top = img_z                          # (nrows, ncols, 3)
    # plot_surface facecolors expects (nrows-1, ncols-1, 4) — use subsampled RGBA
    def _to_face_rgba(img):
        """Down-sample image to face-centre colors: (nr-1, nc-1, 4)."""
        r = (img[:-1,:-1,:] + img[1:,:-1,:] + img[:-1,1:,:] + img[1:,1:,:]) / 4
        a = np.ones((*r.shape[:2], 1))
        return np.concatenate([r, a], axis=-1)

    ax3.plot_surface(XS, YS, ZS,
                     facecolors=_to_face_rgba(fc_top),
                     rstride=1, cstride=1,
                     linewidth=0, antialiased=False, shade=False)

    # ── Front face (Y=0 plane): IPF-Y ─────────────────────────────────────────
    xs_f = np.linspace(0, W, ncols)
    zs_f = np.linspace(0, D, nrows)
    XF, ZF = np.meshgrid(xs_f, zs_f)
    YF = np.zeros_like(XF)

    # Front face: columns=X, rows=Z (bottom of map = z=0, top = z=D)
    # Flip Y image vertically so "bottom" of map is at z=0
    fc_front = img_y[::-1, :, :]
    ax3.plot_surface(XF, YF, ZF,
                     facecolors=_to_face_rgba(fc_front),
                     rstride=1, cstride=1,
                     linewidth=0, antialiased=False, shade=False)

    # ── Right face (X=W plane): IPF-X ─────────────────────────────────────────
    ys_r = np.linspace(0, H, nrows)
    zs_r = np.linspace(0, D, nrows)
    YR, ZR = np.meshgrid(ys_r, zs_r)
    XR = np.full_like(YR, W)

    fc_right = img_x[::-1, :, :]
    ax3.plot_surface(XR, YR, ZR,
                     facecolors=_to_face_rgba(fc_right),
                     rstride=1, cstride=1,
                     linewidth=0, antialiased=False, shade=False)

    # ── Axes styling ──────────────────────────────────────────────────────────
    ax3.set_xlim(0, W); ax3.set_ylim(0, H); ax3.set_zlim(0, D)
    ax3.set_xlabel(f"X  ({W*sx:.1f} µm)", fontsize=9, labelpad=4)
    ax3.set_ylabel(f"Y  ({H*sy:.1f} µm)", fontsize=9, labelpad=4)
    ax3.set_zlabel(f"Z", fontsize=9, labelpad=4)
    ax3.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax3.view_init(elev=28, azim=-50)

    # ── Add grey border edges ─────────────────────────────────────────────────
    edge_kw = dict(color="grey", lw=0.8, alpha=0.7)
    # Bottom rectangle
    for xs, ys, zs in [
        ([0,W,W,0,0],[0,0,H,H,0],[0,0,0,0,0]),
        ([0,W,W,0,0],[0,0,H,H,0],[D,D,D,D,D]),
        ([0,0],[0,0],[0,D]), ([W,W],[0,0],[0,D]),
        ([W,W],[H,H],[0,D]), ([0,0],[H,H],[0,D]),
    ]:
        ax3.plot(xs, ys, zs, **edge_kw)

    # ── Legend triangle (inset axes) ─────────────────────────────────────────
    ax_leg = fig.add_axes([0.78, 0.55, 0.20, 0.22])
    plot_ipf_legend(ax=ax_leg, title=f"{phase_label}\n001", fontsize=8)

    # ── Sample reference frame arrow ─────────────────────────────────────────
    ax_ref = fig.add_axes([0.78, 0.30, 0.20, 0.18])
    ax_ref.annotate("", xy=(0.7, 0.5), xytext=(0.0, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))
    ax_ref.annotate("", xy=(0.5, 0.95), xytext=(0.5, 0.2),
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))
    ax_ref.annotate("", xy=(0.15, 0.15), xytext=(0.5, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))
    ax_ref.text(0.72, 0.46, "X", fontsize=9, fontweight="bold", ha="left")
    ax_ref.text(0.50, 0.98, "Y", fontsize=9, fontweight="bold", ha="center", va="bottom")
    ax_ref.text(0.10, 0.10, "Z", fontsize=9, fontweight="bold", ha="right")
    ax_ref.set_xlim(0, 1); ax_ref.set_ylim(0, 1); ax_ref.axis("off")
    ax_ref.set_title("Z = Rolling direction", fontsize=7, pad=2)

    fig.patch.set_facecolor("white")
    return fig

