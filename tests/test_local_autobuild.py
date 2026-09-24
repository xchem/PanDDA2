"""Validate cut_local_grid_from_sparse against the full-cell cut it replaces.

The local cut must reproduce, in the box region, exactly what densifying the whole
cell and sampling there would give -- on a non-orthogonal (monoclinic) cell, which
is the hard case. Tested two ways: exactly for a linear field (trilinear is exact,
so any discrepancy is a geometry/index bug), and against gemmi's own full-grid
interpolation for a smooth field.
"""

import numpy as np
import gemmi
import pytest

from pandda_gemmi.autobuild.local_grid import (
    cut_local_grid_from_sparse, cut_local_grid_from_dense, box_origin_for,
    _frac_matrix,
)


class _MockMask:
    def __init__(self, indicies):
        self.indicies = indicies


class _MockFrame:
    """Stand-in for DFrame exposing only what the local cut needs."""
    def __init__(self, unit_cell, spacing, indicies):
        self.unit_cell = unit_cell
        self.spacing = spacing
        self.mask = _MockMask(indicies)


def _orth_matrix(cell):
    cols = [cell.orthogonalize(gemmi.Fractional(1, 0, 0)),
            cell.orthogonalize(gemmi.Fractional(0, 1, 0)),
            cell.orthogonalize(gemmi.Fractional(0, 0, 1))]
    return np.array([[c.x, c.y, c.z] for c in cols]).T


def _synthetic(cell_params, nu, nv, nw, fn):
    """A native gemmi FloatGrid filled with fn(cartesian), plus a MockFrame whose
    sparse 'mask' is every grid point (so the local cut sees the full density)."""
    cell = gemmi.UnitCell(*cell_params)
    grid = gemmi.FloatGrid(nu, nv, nw)
    grid.set_unit_cell(cell)
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    M = _orth_matrix(cell)
    iu = np.arange(nu) / nu
    iv = np.arange(nv) / nv
    iw = np.arange(nw) / nw
    fu, fv, fw = np.meshgrid(iu, iv, iw, indexing="ij")
    frac = np.stack([fu, fv, fw], axis=-1)
    cart = frac @ M.T
    arr = fn(cart).astype(np.float32)
    np.array(grid, copy=False)[:, :, :] = arr
    idx = np.nonzero(np.ones((nu, nv, nw), dtype=np.int8))
    sparse = arr[idx]
    frame = _MockFrame(cell_params, (nu, nv, nw), idx)
    return grid, frame, sparse, cell


MONO = (50.0, 55.0, 60.0, 90.0, 95.0, 90.0)


def test_local_cut_linear_exact():
    """Linear field -> trilinear is exact, so the local box must equal the field
    sampled at the box voxels to machine-ish precision (catches geometry/index
    errors on a non-orthogonal cell)."""
    def fn(c):
        return 1.0 + 0.30 * c[..., 0] - 0.20 * c[..., 1] + 0.15 * c[..., 2]
    grid, frame, sparse, cell = _synthetic(MONO, 80, 88, 96, fn)
    centroid = np.array(cell.orthogonalize(gemmi.Fractional(0.5, 0.5, 0.5)).tolist())

    n, sp = 16, 0.5
    local, box_origin = cut_local_grid_from_sparse(frame, sparse, centroid, n, sp)
    got = np.array(local, copy=False)

    ax = np.arange(n) * sp
    gi = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), axis=-1) + box_origin
    expected = fn(gi)
    assert np.allclose(got, expected, atol=2e-3), \
        f"max dev {np.abs(got - expected).max():.4f}"


def test_local_cut_matches_full_unmask():
    """Smooth field -> the local cut must match gemmi's own full-grid
    interpolation at the same box positions (i.e. local cut == full cut)."""
    def fn(c):
        return (np.sin(0.25 * c[..., 0]) + np.cos(0.20 * c[..., 1])
                + 0.5 * np.sin(0.15 * c[..., 2]))
    grid, frame, sparse, cell = _synthetic(MONO, 80, 88, 96, fn)
    centroid = np.array(cell.orthogonalize(gemmi.Fractional(0.45, 0.5, 0.55)).tolist())

    n, sp = 24, 0.5
    local, box_origin = cut_local_grid_from_sparse(frame, sparse, centroid, n, sp)
    got = np.array(local, copy=False)

    # reference: sample the FULL native grid at the same native box positions
    transform = gemmi.Transform()
    transform.mat.fromlist((np.eye(3) * sp).tolist())
    transform.vec.fromlist([float(x) for x in box_origin])
    ref = np.zeros((n, n, n), dtype=np.float32)
    grid.interpolate_values(ref, transform)

    assert np.corrcoef(got.ravel(), ref.ravel())[0, 1] > 0.999
    assert np.sqrt(((got - ref) ** 2).mean()) < 0.02


def test_frac_matrix_inverts_orth():
    """Sanity: F (frac<-cart) is the inverse of the orthogonalisation matrix."""
    cell = gemmi.UnitCell(*MONO)
    F = _frac_matrix(cell)
    O = _orth_matrix(cell)
    assert np.allclose(F @ O, np.eye(3), atol=1e-6)


def test_dense_cut_matches_gemmi_on_a_DIFFERENT_grid():
    """The raw xmap is sampled independently of the reference frame (sample_rate=3
    vs resolution/0.4999), so it has a different shape and the frame's mask
    indices do not address it at all. cut_local_grid_from_dense must cut from
    that array on its OWN geometry and match gemmi's interpolation of it.

    Regression: cutting the raw xmap with cut_local_grid_from_sparse instead
    (i.e. via mask indices belonging to another grid) silently read the wrong
    voxels, corrupted the build CNN's xmap channel, and moved real poses by up
    to 22 A -- while every synthetic test still passed, because they all used a
    single grid for everything.
    """
    def fn(c):
        return (np.sin(0.21 * c[..., 0]) + np.cos(0.17 * c[..., 1])
                + 0.5 * np.sin(0.13 * c[..., 2]))

    # Reference frame at one sampling ...
    _, frame, _, cell = _synthetic(MONO, 80, 88, 96, fn)
    # ... and the "raw xmap" on a coarser, differently-shaped grid, same cell.
    raw_grid, _, _, _ = _synthetic(MONO, 50, 55, 60, fn)
    raw_dense = np.array(raw_grid, copy=False)
    assert raw_dense.shape != tuple(frame.spacing), "shapes must differ for this test"

    centroid = np.array(cell.orthogonalize(gemmi.Fractional(0.45, 0.5, 0.55)).tolist())
    n, sp = 24, 0.5
    box_origin = box_origin_for(centroid, n, sp)

    got = np.array(cut_local_grid_from_dense(raw_dense, MONO, box_origin, n, sp), copy=False)

    # Reference: gemmi interpolating the raw grid itself at the box positions.
    transform = gemmi.Transform()
    transform.mat.fromlist((np.eye(3) * sp).tolist())
    transform.vec.fromlist([float(x) for x in box_origin])
    ref = np.zeros((n, n, n), dtype=np.float32)
    raw_grid.interpolate_values(ref, transform)

    assert np.corrcoef(got.ravel(), ref.ravel())[0, 1] > 0.999
    assert np.sqrt(((got - ref) ** 2).mean()) < 0.02


def test_box_origin_is_shared_across_cuts():
    """Every cut about one centroid must land on the same local frame, or the
    channels handed to the scorer are mutually offset."""
    def fn(c):
        return c[..., 0] * 0.1
    _, frame, sparse, cell = _synthetic(MONO, 80, 88, 96, fn)
    centroid = np.array(cell.orthogonalize(gemmi.Fractional(0.5, 0.5, 0.5)).tolist())
    n, sp = 16, 0.5
    _, origin_sparse = cut_local_grid_from_sparse(frame, sparse, centroid, n, sp)
    assert np.allclose(origin_sparse, box_origin_for(centroid, n, sp))
