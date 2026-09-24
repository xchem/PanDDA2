"""Local unmask: build a small density box from the sparse representation WITHOUT
allocating the full unit cell.

This is the foundation of the ``PANDDA_LOCAL_AUTOBUILD`` path. PanDDA stores per-event/model maps sparsely
(``reference_frame.mask.indicies`` = a 3-tuple ``(U,V,W)`` of native-grid indices,
``sparse.data[i]`` the value at ``(U[i],V[i],W[i])`` on the native ``(nu,nv,nw)``
P1 grid). ``reference_frame.unmask`` densifies the WHOLE cell -- catastrophic when
the cell has a ~190 A axis and it runs per (model x event x conformer) task.

``cut_local_grid_from_sparse`` instead resamples an orthonormal n^3 box about the
event centroid directly from the sparse points in the box's native footprint. The
returned grid lives in a LOCAL frame: box corner -> (0,0,0). Translate structures
by ``-box_origin`` to score/fit against it, and add ``box_origin`` back to map a
fitted pose into the native frame. Scoring (RSCC, CNN, masks) is
translation-invariant, so the local frame changes nothing but the memory.
"""

from __future__ import annotations

import numpy as np
import gemmi


def _frac_matrix(cell: gemmi.UnitCell) -> np.ndarray:
    """3x3 fractionalisation matrix F such that frac = F @ cartesian. Taken from
    gemmi so it matches its cell geometry exactly (handles non-orthogonal cells)."""
    cols = [cell.fractionalize(gemmi.Position(1, 0, 0)),
            cell.fractionalize(gemmi.Position(0, 1, 0)),
            cell.fractionalize(gemmi.Position(0, 0, 1))]
    return np.array([[c.x, c.y, c.z] for c in cols]).T


def _trilinear_window(window: np.ndarray, coords: np.ndarray, shape: tuple) -> np.ndarray:
    """Trilinear sample ``window`` (a periodic native-lattice block) at fractional
    window-index ``coords`` (N,3); periodic wrap via modulo on the window shape's
    parent lattice is handled by the caller folding indices in. Out-of-window
    contributions are treated as 0 (unmasked native points are 0 in PanDDA)."""
    n0, n1, n2 = window.shape
    i0 = np.floor(coords).astype(np.int64)
    f = coords - i0
    out = np.zeros(coords.shape[0], dtype=np.float64)
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                ii = i0[:, 0] + di
                jj = i0[:, 1] + dj
                kk = i0[:, 2] + dk
                inb = (ii >= 0) & (ii < n0) & (jj >= 0) & (jj < n1) & \
                      (kk >= 0) & (kk < n2)
                wt = (np.where(di, f[:, 0], 1 - f[:, 0]) *
                      np.where(dj, f[:, 1], 1 - f[:, 1]) *
                      np.where(dk, f[:, 2], 1 - f[:, 2]))
                idx = np.where(inb)[0]
                out[idx] += wt[idx] * window[ii[idx], jj[idx], kk[idx]]
    return out


def box_origin_for(centroid, n: int, spacing: float) -> np.ndarray:
    """Native Cartesian position of box voxel (0,0,0), snapped to the box
    lattice so that every cut about the same centroid shares one local frame."""
    centroid = np.asarray(centroid, dtype=np.float64)
    half = (n / 2.0) * spacing
    return np.round((centroid - half) / spacing) * spacing


def _box_cart(box_origin: np.ndarray, n: int, spacing: float) -> np.ndarray:
    """(n^3, 3) native Cartesian positions of the box voxels."""
    ax = np.arange(n, dtype=np.float64) * spacing
    grid_ijk = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), axis=-1)
    return (grid_ijk + box_origin[None, None, None, :]).reshape(-1, 3)


def _as_local_grid(values: np.ndarray, n: int, spacing: float) -> gemmi.FloatGrid:
    local = gemmi.FloatGrid(n, n, n)
    local.set_unit_cell(gemmi.UnitCell(n * spacing, n * spacing, n * spacing,
                                       90.0, 90.0, 90.0))
    local.spacegroup = gemmi.SpaceGroup("P 1")
    np.array(local, copy=False)[:, :, :] = values
    return local


def _trilinear_periodic(arr: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Trilinear sample of a full-cell dense array at fractional grid ``coords``
    (N,3), wrapping periodically -- the array spans the whole cell, so an index
    outside it is the same density one cell over."""
    n0, n1, n2 = arr.shape
    i0 = np.floor(coords).astype(np.int64)
    f = coords - i0
    out = np.zeros(coords.shape[0], dtype=np.float64)
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                wt = (np.where(di, f[:, 0], 1 - f[:, 0]) *
                      np.where(dj, f[:, 1], 1 - f[:, 1]) *
                      np.where(dk, f[:, 2], 1 - f[:, 2]))
                out += wt * arr[(i0[:, 0] + di) % n0,
                                (i0[:, 1] + dj) % n1,
                                (i0[:, 2] + dk) % n2]
    return out


def cut_local_grid_from_dense(dense_array, unit_cell, box_origin,
                              n: int, spacing: float) -> gemmi.FloatGrid:
    """Cut the same local box from a DENSE full-cell array that has its own
    sampling, independent of the reference frame's.

    Needed for the raw xmap: it is sampled at ``sample_rate=3`` while the
    reference frame uses ``resolution/0.4999``, so the two grids have different
    shapes and the frame's mask indices do not address this array at all. Only
    the unit cell is shared -- which is enough, because the box is defined in
    Cartesian space.
    """
    arr = np.asarray(dense_array, dtype=np.float32)
    cell = gemmi.UnitCell(*unit_cell)
    box_origin = np.asarray(box_origin, dtype=np.float64)

    frac = _box_cart(box_origin, n, spacing) @ _frac_matrix(cell).T
    gi = frac * np.array(arr.shape, dtype=np.float64)
    vals = _trilinear_periodic(arr, gi).reshape(n, n, n).astype(np.float32)
    return _as_local_grid(vals, n, spacing)


def cut_local_grid_from_sparse(reference_frame, sparse_data, centroid,
                               n: int, spacing: float):
    """Resample an orthonormal n^3 box (A spacing) about ``centroid`` from the
    sparse density, without densifying the full cell.

    Returns ``(local_grid, box_origin)`` -- a P1 gemmi FloatGrid with cubic cell
    ``n*spacing`` holding the density in box-local frame, and ``box_origin`` (the
    native Cartesian position of box voxel (0,0,0)). Native pos of voxel (i,j,k)
    = box_origin + (i,j,k)*spacing; the grid stores that density at box-frame
    (i,j,k)*spacing, so sample/score with structures translated by -box_origin.
    """
    cell = gemmi.UnitCell(*reference_frame.unit_cell)
    nu, nv, nw = reference_frame.spacing
    U, V, W = reference_frame.mask.indicies
    data = np.asarray(sparse_data, dtype=np.float32)

    box_origin = box_origin_for(centroid, n, spacing)

    # Native fractional -> native grid coords for every box voxel.
    frac = _box_cart(box_origin, n, spacing) @ _frac_matrix(cell).T
    gi = frac * np.array([nu, nv, nw])

    # Native-index footprint of the box (with margin), then build the window.
    lo = np.floor(gi.min(0)).astype(np.int64) - 2
    hi = np.ceil(gi.max(0)).astype(np.int64) + 3
    shp = (int(hi[0] - lo[0]), int(hi[1] - lo[1]), int(hi[2] - lo[2]))
    window = np.zeros(shp, dtype=np.float32)
    # Scatter sparse points whose (periodically folded) index lands in the window.
    uu = (U.astype(np.int64) - lo[0]) % nu
    vv = (V.astype(np.int64) - lo[1]) % nv
    ww = (W.astype(np.int64) - lo[2]) % nw
    sel = (uu < shp[0]) & (vv < shp[1]) & (ww < shp[2])
    window[uu[sel], vv[sel], ww[sel]] = data[sel]

    # Sample the box from the window (coords relative to window origin `lo`).
    vals = _trilinear_window(window, gi - lo[None, :], shp).reshape(n, n, n).astype(np.float32)
    return _as_local_grid(vals, n, spacing), box_origin
