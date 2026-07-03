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
def _in_fundamental_triangle(x2d, y2d, expand=0.0):
    """
    Return boolean mask: True if the stereographic point lies inside the
    BCC fundamental triangle (001)-(101)-(111).

    ``expand`` slightly grows the triangle outward from its centroid (in the
    same units as the stereographic coords). A tiny positive value lets the
    filled color reach fully under the drawn boundary line so no white sliver
    is left at the edges.
    """
    # Three edges defined by great circles; approximate using barycentric coords
    # Vertices (optionally pushed outward from the centroid):
    v0 = np.array(_P001, float)
    v1 = np.array(_P101, float)
    v2 = np.array(_P111, float)
    if expand:
        cen = (v0 + v1 + v2) / 3.0
        v0 = cen + (v0 - cen) * (1.0 + expand)
        v1 = cen + (v1 - cen) * (1.0 + expand)
        v2 = cen + (v2 - cen) * (1.0 + expand)

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

    # Sample a dense grid covering ALL three corners (101 is the widest in x,
    # 111 is the tallest in y). Padding a little past each corner guarantees the
    # bilinear-interpolated fill reaches the triangle edges with no white sliver.
    pad = 0.02
    x_corners = [_P001[0], _P101[0], _P111[0]]
    y_corners = [_P001[1], _P101[1], _P111[1]]
    gx = np.linspace(min(x_corners) - pad, max(x_corners) + pad, 400)
    gy = np.linspace(min(y_corners) - pad, max(y_corners) + pad, 400)
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

    # Mask to triangle (expand a hair so color fills under the border line)
    in_tri = _in_fundamental_triangle(pts2d[:, 0], pts2d[:, 1], expand=0.015)
    img_r = np.where(in_tri, rgb[:, 0], np.nan)
    img_g = np.where(in_tri, rgb[:, 1], np.nan)
    img_b = np.where(in_tri, rgb[:, 2], np.nan)

    # Combine to RGBA
    alpha = np.where(in_tri, 1.0, 0.0)
    img = np.stack([
        img_r.reshape(400, 400),
        img_g.reshape(400, 400),
        img_b.reshape(400, 400),
        alpha.reshape(400, 400),
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
                     phase_label="Ferrite",
                     equal_aspect=True):
    """
    3-D isometric cube — clean IPF faces, no internal diagonal lines,
    axis labels drawn OUTSIDE the cube, and a sample-axes key that matches
    the isometric view directions.

    Parameters
    ----------
    equal_aspect : bool, default True
        If True the cube is drawn as a visually undistorted cube (all three
        edges the same length on screen) regardless of the map's physical
        aspect ratio; the true physical extents are reported in the axis
        labels. This is the scientifically clearest default for a map whose
        X and Y spans differ a lot. If False the cube edges are drawn in
        physical proportion (X:Y), and a single common-length scale bar (same
        µm value on every axis) is added.
    """
    from scipy.ndimage import map_coordinates

    phi1 = np.asarray(phi1_arr, float)
    Phi  = np.asarray(Phi_arr,  float)
    phi2 = np.asarray(phi2_arr, float)
    xp   = np.asarray(x_pos,   float)
    yp   = np.asarray(y_pos,   float)

    # ── Grid ──────────────────────────────────────────────────────────────────
    valid = (np.isfinite(phi1) & np.isfinite(Phi) & np.isfinite(phi2)
             & np.isfinite(xp) & np.isfinite(yp))
    phi1, Phi, phi2 = phi1[valid], Phi[valid], phi2[valid]
    xp, yp = xp[valid], yp[valid]

    ux = np.sort(np.unique(xp))
    uy = np.sort(np.unique(yp))
    sx = float(np.median(np.diff(ux))) if step_x is None and len(ux)>1 else (step_x or 1.)
    sy = float(np.median(np.diff(uy))) if step_y is None and len(uy)>1 else (step_y or 1.)

    xi = np.round((xp - xp.min()) / sx).astype(int)
    yi = np.round((yp - yp.min()) / sy).astype(int)
    ncols = xi.max() + 1
    nrows = yi.max() + 1

    MAX_PX = 300
    ds = max(1, max(nrows, ncols) // MAX_PX)
    if ds > 1:
        keep = (xi % ds == 0) & (yi % ds == 0)
        phi1, Phi, phi2 = phi1[keep], Phi[keep], phi2[keep]
        xi = xi[keep] // ds; yi = yi[keep] // ds
        ncols = xi.max() + 1; nrows = yi.max() + 1

    W, H = ncols, nrows   # W = X cols, H = Y rows

    # ── IPF images ────────────────────────────────────────────────────────────
    def make_img(rgb, ri, ci, nr, nc):
        img = np.full((nr, nc, 3), 0.82, dtype=np.float32)
        img[ri, ci] = np.clip(rgb, 0, 1).astype(np.float32)
        return img

    img_z = make_img(compute_ipf_colormap(phi1,Phi,phi2,"Z"), yi,xi,H,W)
    img_y = make_img(compute_ipf_colormap(phi1,Phi,phi2,"Y"), yi,xi,H,W)
    img_x = make_img(compute_ipf_colormap(phi1,Phi,phi2,"X"), yi,xi,H,W)

    # ── Isometric geometry ────────────────────────────────────────────────────
    # Edge lengths (canvas pixels) along each cube axis.
    SCALE = 4
    if equal_aspect:
        # Visually undistorted cube: every edge the same on-screen length.
        EDGE = SCALE * max(W, H)
        LX = LY = LZ = float(EDGE)
    else:
        # Physical proportion: X:Y match the map; Z reuses the Y edge length.
        LX = float(W * SCALE)   # X edge (front width)
        LY = float(H * SCALE)   # Y edge (depth)
        LZ = float(H * SCALE)   # Z edge (height)

    a30 = np.radians(30)
    # Full-edge isometric vectors (canvas pixels, y=down):
    #   +X → right-down,  +Y → left-down (depth),  +Z → straight up.
    EX = np.array([ np.cos(a30),  np.sin(a30)]) * LX
    EY = np.array([-np.cos(a30),  np.sin(a30)]) * LY
    EZ = np.array([ 0.0,         -1.0         ]) * LZ

    def P(fx, fy, fz):
        """Fractional cube coordinate (each in [0,1]) → canvas pixel."""
        return fx*EX + fy*EY + fz*EZ

    # Raw corners (fractional along each edge)
    raw_corners = {
        "BLF": P(0,1,0), "BRF": P(1,1,0), "BLB": P(0,0,0), "BRB": P(1,0,0),
        "TLF": P(0,1,1), "TRF": P(1,1,1), "TLB": P(0,0,1), "TRB": P(1,0,1),
    }
    all_pts = np.array(list(raw_corners.values()))
    PAD = 24
    shift = np.array([PAD - all_pts[:,0].min(), PAD - all_pts[:,1].min()])
    C = {k: v + shift for k, v in raw_corners.items()}

    iso_max = all_pts + shift
    canvas_W = int(np.ceil(iso_max[:,0].max())) + PAD
    canvas_H = int(np.ceil(iso_max[:,1].max())) + PAD

    BLF,BRF,BLB,BRB = C["BLF"],C["BRF"],C["BLB"],C["BRB"]
    TLF,TRF,TLB,TRB = C["TLF"],C["TRF"],C["TLB"],C["TRB"]

    # ── Warp helper ───────────────────────────────────────────────────────────
    def warp_face(src, bl, br, tl, shade=1.0):
        """
        Warp src (H_s×W_s×3 float32) into parallelogram defined by
        bottom-left=bl, bottom-right=br, top-left=tl (canvas pixel coords).
        right = br-bl  maps u in [0,1] → src cols
        up    = tl-bl  maps v in [0,1] → src rows (v=1 → row 0 = top of image)
        """
        H_s, W_s = src.shape[:2]
        bl = np.asarray(bl, float)
        right = np.asarray(br, float) - bl
        up    = np.asarray(tl, float) - bl

        A  = np.stack([right, up], axis=1)
        Ai = np.linalg.inv(A)

        gc = np.tile(np.arange(canvas_W, dtype=np.float32), (canvas_H, 1))
        gr = np.tile(np.arange(canvas_H, dtype=np.float32)[:,None], (1, canvas_W))

        dc = gc - bl[0]; dr = gr - bl[1]
        u  = Ai[0,0]*dc + Ai[0,1]*dr
        v  = Ai[1,0]*dc + Ai[1,1]*dr

        # Mask: strictly inside the parallelogram
        EPS = 3e-3
        mask = (u >= -EPS) & (u <= 1+EPS) & (v >= -EPS) & (v <= 1+EPS)
        u = np.clip(u, 0, 1); v = np.clip(v, 0, 1)

        src_col = u * (W_s - 1)
        src_row = (1.0 - v) * (H_s - 1)   # v=1 → row 0 (top of image = top of face)

        face = np.zeros((canvas_H, canvas_W, 3), dtype=np.float32)
        for ch in range(3):
            vals = map_coordinates(src[:,:,ch], [src_row, src_col],
                                   order=1, mode='nearest')
            face[:,:,ch] = np.clip(vals * shade, 0, 1)
        return face, mask

    # ── Render canvas ─────────────────────────────────────────────────────────
    canvas = np.ones((canvas_H, canvas_W, 3), dtype=np.float32)

    # Top face   (IPF-Z): bottom=back edge, top=front edge
    f, m = warp_face(img_z[::-1], TLB, TRB, TLF, shade=1.00)
    canvas[m] = f[m]
    # Front face (IPF-Y): bottom=bottom-front edge, top=top-front edge
    f, m = warp_face(img_y[::-1], BLF, BRF, TLF, shade=0.78)
    canvas[m] = f[m]
    # Right face (IPF-X): bottom=bottom-right edge, top=top-right edge
    # For right face: "columns" go along Y-depth (BRF→BRB), rows along Z
    f, m = warp_face(img_x[::-1], BRF, BRB, TRF, shade=0.60)
    canvas[m] = f[m]

    # ── Scale bars (drawn onto canvas array) ──────────────────────────────────
    def draw_scalebar_on_canvas(canvas, p_start, p_end, label,
                                bar_thickness=3, color=(1.,1.,1.)):
        """Draw a scale bar line between two canvas pixel positions."""
        from PIL import Image as PILImg, ImageDraw, ImageFont
        # Draw line
        img_pil = PILImg.fromarray((canvas * 255).astype(np.uint8))
        draw = ImageDraw.Draw(img_pil)
        x0,y0 = int(round(p_start[0])), int(round(p_start[1]))
        x1,y1 = int(round(p_end[0])),   int(round(p_end[1]))
        col8 = tuple(int(c*255) for c in color)
        draw.line([(x0,y0),(x1,y1)], fill=col8, width=bar_thickness)
        # End ticks (perpendicular, 4px each side)
        dx, dy = x1-x0, y1-y0
        L = max(1, np.sqrt(dx**2+dy**2))
        nx, ny = -dy/L*4, dx/L*4
        for px, py in [(x0,y0),(x1,y1)]:
            draw.line([(int(px-nx),int(py-ny)),(int(px+nx),int(py+ny))],
                      fill=col8, width=bar_thickness)
        return np.array(img_pil).astype(np.float32)/255.0

    # Compute a nice round scale bar length in µm (~25% of the given span)
    def nice_bar(total_um, frac=0.25):
        raw = total_um * frac
        mag = 10**np.floor(np.log10(raw))
        for m in [1,2,5,10]:
            if m*mag >= raw*0.5:
                return m*mag
        return mag*10

    X_um = W * sx * ds
    Y_um = H * sy * ds

    # Physical pixels-per-µm along each drawn edge (differ when equal_aspect).
    ppu_x = LX / X_um
    ppu_y = LY / Y_um
    ppu_z = LZ / Y_um   # Z is a projection direction, not a physical depth

    # Scale bars are only meaningful when the cube is in physical proportion.
    # We use ONE common µm value on all three axes ("same scale convention"),
    # which — because EBSD grids are square (sx≈sy) — renders as the same
    # visible length on every face. In equal-aspect (normalized) mode we omit
    # the bars and report the true physical extents in the axis labels instead.
    scalebar_um = None
    if not equal_aspect:
        scalebar_um = nice_bar(max(X_um, Y_um))

        # Front face bar (along +X), near the bottom-left of the front face
        sb_fx_start = BLF + (BRF-BLF)*0.07 + (TLF-BLF)*0.08
        bar_x_vec   = (BRF - BLF) / np.linalg.norm(BRF - BLF)
        sb_fx_end   = sb_fx_start + bar_x_vec * (scalebar_um * ppu_x)

        # Right face bar (along +Y depth)
        sb_ry_start = BRF + (BRB-BRF)*0.07 + (TRF-BRF)*0.08
        bar_y_vec   = (BRB - BRF) / np.linalg.norm(BRB - BRF)
        sb_ry_end   = sb_ry_start + bar_y_vec * (scalebar_um * ppu_y)

        # Top face bar (along +X, near the front edge)
        sb_tx_start = TLF + (TRF-TLF)*0.07 + (TLB-TLF)*0.05
        bar_tx_vec  = (TRF - TLF) / np.linalg.norm(TRF - TLF)
        sb_tx_end   = sb_tx_start + bar_tx_vec * (scalebar_um * ppu_x)

        # Draw the three bars (thicker so they read clearly on the map)
        for p0, p1 in [(sb_fx_start,sb_fx_end),
                       (sb_ry_start,sb_ry_end),
                       (sb_tx_start,sb_tx_end)]:
            canvas = draw_scalebar_on_canvas(canvas, p0, p1, "", bar_thickness=4)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_axes([0.05, 0.03, 0.66, 0.93])
    ax.imshow(canvas, origin="upper", interpolation="nearest", aspect="equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── Cube edges — ONLY the 12 outer edges, NO internal diagonals ───────────
    ekw = dict(color="black", lw=1.0, zorder=5)
    outer_edges = [
        # Bottom: only the 3 visible edges
        (BLF,BRF),               # front bottom
        (BRF,BRB),               # right bottom
        (BLF,BLB),               # left bottom (visible as left edge of front face)
        # Top: all 4 visible
        (TLF,TRF),               # front top
        (TRF,TRB),               # right top
        (TRB,TLB),               # back top
        (TLB,TLF),               # left top
        # Vertical pillars: only 3 visible
        (BLF,TLF),               # front-left pillar
        (BRF,TRF),               # front-right pillar
        (BRB,TRB),               # back-right pillar
    ]
    for a,b in outer_edges:
        ax.plot([a[0],b[0]],[a[1],b[1]],**ekw)

    # ── Scale bar labels (only in physical mode; placed beside each bar) ──────
    if not equal_aspect and scalebar_um is not None:
        sblkw = dict(fontsize=8, ha="center", va="center", zorder=12,
                     color="black", fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.15", fc="white",
                               alpha=0.85, lw=0))
        for p0, p1 in [(sb_fx_start, sb_fx_end),
                       (sb_ry_start, sb_ry_end),
                       (sb_tx_start, sb_tx_end)]:
            mid  = (p0 + p1) / 2
            perp = np.array([-(p1-p0)[1], (p1-p0)[0]])
            perp /= (np.linalg.norm(perp) + 1e-9)
            ax.text(mid[0]+perp[0]*11, mid[1]+perp[1]*11,
                    f"{scalebar_um:.0f} µm", **sblkw)

    # ── Axis labels — OUTSIDE the cube, with annotation arrows ───────────────
    akw = dict(fontsize=10, fontweight="bold", ha="center", va="center",
               zorder=12,
               bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#aaa",
                         alpha=0.95, lw=0.7))
    arro = dict(arrowstyle="-|>", color="black", lw=1.0)

    # X label: below the bottom-front edge midpoint
    mX = (BLF + BRF) / 2
    ax.annotate(f"X\n({X_um:.0f} µm)", xy=(mX[0], mX[1]),
                xytext=(mX[0], mX[1]+30),
                arrowprops=arro, fontsize=9, fontweight="bold",
                ha="center", va="top", zorder=12, color="darkred",
                bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#ccc",alpha=0.95,lw=0.5))

    # Y label: to the left of the bottom-left edge midpoint
    mY = (BLF + BLB) / 2
    ax.annotate(f"Y\n({Y_um:.0f} µm)", xy=(mY[0], mY[1]),
                xytext=(mY[0]-38, mY[1]),
                arrowprops=arro, fontsize=9, fontweight="bold",
                ha="right", va="center", zorder=12, color="#006600",
                bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#ccc",alpha=0.95,lw=0.5))

    # Z label: to the left of the left-vertical edge midpoint
    mZ = (BLF + TLF) / 2
    ax.annotate("Z", xy=(mZ[0], mZ[1]),
                xytext=(mZ[0]-32, mZ[1]),
                arrowprops=arro, fontsize=10, fontweight="bold",
                ha="right", va="center", zorder=12, color="navy",
                bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#ccc",alpha=0.95,lw=0.5))

    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax.set_xlim(-70, canvas_W + 2)
    ax.set_ylim(canvas_H + 55, -2)   # extra room at bottom for labels/note

    # Short note clarifying the display aspect (kept OUTSIDE the cube, at the
    # very bottom of the map axes).
    if equal_aspect:
        aspect_note = ("Cube drawn with equal edges (display-normalized); "
                       "true extents are given on the X/Y labels.")
    else:
        aspect_note = ("Cube drawn in physical proportion (X:Y). "
                       f"Scale bar = {scalebar_um:.0f} µm on every axis.")
    ax.text(canvas_W/2, canvas_H + 48, aspect_note,
            ha="center", va="bottom", fontsize=7, color="#444", zorder=12)

    # ── Legend triangle (IPF color key) ───────────────────────────────────────
    ax_leg = fig.add_axes([0.72, 0.60, 0.25, 0.23])
    plot_ipf_legend(ax=ax_leg, title=str(phase_label), fontsize=8)

    # ── Sample-axes key — matches the isometric view directions exactly ───────
    # Reuse the same edge vectors that build the cube so the arrows point the
    # same way. Canvas is y-down, so flip the y-component for this y-up axes.
    ax_ref = fig.add_axes([0.72, 0.35, 0.25, 0.20])
    dX = np.array([ EX[0], -EX[1]]); dX /= np.linalg.norm(dX)
    dY = np.array([ EY[0], -EY[1]]); dY /= np.linalg.norm(dY)
    dZ = np.array([ EZ[0], -EZ[1]]); dZ /= np.linalg.norm(dZ)
    origin = np.array([0.42, 0.40])
    L = 0.42
    for dvec, lbl, col in [(dX,"X","darkred"), (dY,"Y","#006600"), (dZ,"Z","navy")]:
        head = origin + dvec * L
        ax_ref.annotate("", xy=head, xytext=origin,
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
        lp = origin + dvec * (L + 0.14)
        ax_ref.text(lp[0], lp[1], lbl, fontsize=10, fontweight="bold",
                    color=col, ha="center", va="center")
    ax_ref.set_xlim(-0.25, 1.15); ax_ref.set_ylim(-0.25, 1.15)
    ax_ref.set_aspect("equal")
    ax_ref.axis("off")
    ax_ref.set_title("Sample axes (view)", fontsize=7, pad=2)

    # ── Face-identity key — which IPF projection each face shows (OUTSIDE) ────
    ax_key = fig.add_axes([0.72, 0.14, 0.25, 0.16])
    ax_key.axis("off")
    ax_key.set_title("Faces", fontsize=7, pad=2)
    face_key = [
        ("Top face",   "IPF ∥ Z", 1.00),
        ("Front face", "IPF ∥ Y", 0.78),
        ("Right face", "IPF ∥ X", 0.60),
    ]
    for i, (face, proj, shade) in enumerate(face_key):
        yy = 0.80 - i * 0.32
        ax_key.add_patch(plt.Rectangle((0.02, yy-0.09), 0.16, 0.18,
                                       facecolor=(shade*0.6, shade*0.6, shade*0.6),
                                       edgecolor="black", lw=0.5))
        ax_key.text(0.24, yy, f"{face} = {proj}", fontsize=7.5,
                    va="center", ha="left")
    ax_key.set_xlim(0, 1); ax_key.set_ylim(0, 1)

    return fig
