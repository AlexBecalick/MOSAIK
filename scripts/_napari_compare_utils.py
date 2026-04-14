#!/usr/bin/env python3
"""Utility helpers for the Xenium/MERSCOPE Napari comparison viewer."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

log = logging.getLogger("napari_compare_utils")


def first_existing_col(df_like, candidates: Iterable[str]) -> str | None:
    """Return the first existing column name from candidates."""
    cols = set(map(str, list(df_like.columns)))
    for col in candidates:
        if col in cols:
            return col
    return None


def to_pandas(df_like) -> pd.DataFrame:
    """Convert a pandas/dask-like dataframe to an in-memory pandas DataFrame."""
    if isinstance(df_like, pd.DataFrame):
        return df_like.copy()
    if hasattr(df_like, "compute"):
        return df_like.compute()
    return pd.DataFrame(df_like).copy()


def assignment_mask(series_like) -> pd.Series:
    """Infer assigned/unassigned transcript status from a column."""
    series = pd.Series(series_like)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float) > 0

    as_str = series.astype("string")
    bad_values = {"", "0", "-1", "nan", "None", "<NA>"}
    return as_str.notna() & ~as_str.isin(bad_values)


def get_scale0_dataarray(image_elem):
    """Return the underlying DataArray from a SpatialData image element."""
    if hasattr(image_elem, "keys") and "scale0" in image_elem:
        node = image_elem["scale0"]
        if hasattr(node, "ds"):
            if "image" in node.ds:
                return node.ds["image"]
            if len(node.ds.data_vars) > 0:
                return next(iter(node.ds.data_vars.values()))

    if hasattr(image_elem, "ds"):
        if "image" in image_elem.ds:
            return image_elem.ds["image"]
        if len(image_elem.ds.data_vars) > 0:
            return next(iter(image_elem.ds.data_vars.values()))

    return image_elem


def ensure_cyx(image_da):
    """Normalize a DataArray to (c, y, x) dimensions."""
    dims = tuple(str(d) for d in image_da.dims)

    if all(d in dims for d in ("c", "y", "x")):
        return image_da.transpose("c", "y", "x")
    if all(d in dims for d in ("y", "x", "c")):
        return image_da.transpose("c", "y", "x")
    if all(d in dims for d in ("y", "x")):
        return image_da.expand_dims(c=["c0"]).transpose("c", "y", "x")

    raise ValueError(f"Unsupported image dims for channel extraction: {dims}")


def channel_labels(image_cyx) -> list[str]:
    """Get channel labels from a (c, y, x) image DataArray."""
    if "c" in image_cyx.coords:
        return [str(c) for c in image_cyx.coords["c"].values]
    return [f"c{i}" for i in range(int(image_cyx.sizes.get("c", 1)))]


def _coords_origin_step(coords) -> tuple[float, float]:
    """Infer origin and step for a monotonic coordinate array."""
    if coords is None:
        return 0.0, 1.0

    arr = np.asarray(coords, dtype=float)
    if arr.size == 0:
        return 0.0, 1.0
    if arr.size == 1:
        return float(arr[0]), 1.0

    diffs = np.diff(arr)
    step = float(np.median(diffs))
    if not np.allclose(diffs, step, rtol=1e-3, atol=1e-6):
        step = float(diffs[0])

    return float(arr[0]), float(step)


def build_napari_affine_from_px_to_um(
    x_transform: tuple[float, float, float],
    y_transform: tuple[float, float, float],
    x_coords=None,
    y_coords=None,
) -> np.ndarray:
    """Build a 3x3 affine for napari (row/col -> y/x) from px->um transforms."""
    a, b, c = map(float, x_transform)
    d, e, f = map(float, y_transform)

    x_origin, x_step = _coords_origin_step(x_coords)
    y_origin, y_step = _coords_origin_step(y_coords)

    # Napari uses (row=y_idx, col=x_idx) order.
    # Convert idx -> pixel coords -> micron coords in one matrix.
    return np.array(
        [
            [e * y_step, d * x_step, d * x_origin + e * y_origin + f],  # y_um
            [b * y_step, a * x_step, a * x_origin + b * y_origin + c],  # x_um
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def resolve_dataset_mask_affine(
    dataset_name: str,
    merscope_transform_path: str | Path | None = None,
    xenium_spec_path: str | Path | None = None,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Resolve pixel->micron affine transform for a dataset."""
    ds = str(dataset_name).upper()

    if ds == "MERSCOPE":
        path = Path(merscope_transform_path) if merscope_transform_path else None
        if path is not None and path.exists():
            matrix = np.loadtxt(path)
            inv = np.linalg.inv(matrix)
            return (
                (float(inv[0, 0]), float(inv[0, 1]), float(inv[0, 2])),
                (float(inv[1, 0]), float(inv[1, 1]), float(inv[1, 2])),
            )

        log.warning(
            "[MERSCOPE] Transform file missing; using fallback 0.108 um/px isotropic scale."
        )
        return (0.108, 0.0, 0.0), (0.0, 0.108, 0.0)

    if ds == "XENIUM":
        mpp = None
        path = Path(xenium_spec_path) if xenium_spec_path else None
        if path is not None and path.exists():
            try:
                spec = json.loads(path.read_text())
                if "pixel_size" in spec:
                    mpp = float(spec["pixel_size"])
            except Exception as exc:
                log.warning("[XENIUM] Failed to parse spec file %s (%s)", path, exc)

        if mpp is None:
            mpp = 0.2125
            log.warning(
                "[XENIUM] Spec file missing/unreadable; using fallback pixel_size=%s um.",
                mpp,
            )

        return (float(mpp), 0.0, 0.0), (0.0, float(mpp), 0.0)

    raise ValueError(f"Unknown dataset: {dataset_name}")


def load_points_dataframe(
    points_obj,
    x_col: str,
    y_col: str,
    assignment_col: str | None = None,
    max_points: int | None = None,
    random_state: int = 42,
) -> pd.DataFrame:
    """Load selected points columns into pandas with optional unbiased sampling."""
    cols = [x_col, y_col] + ([assignment_col] if assignment_col is not None else [])
    work = points_obj[cols]

    if max_points is None:
        return to_pandas(work)

    if hasattr(work, "npartitions") and hasattr(work, "compute"):
        total = int(work.map_partitions(len, meta=("n", "i8")).sum().compute())
        if total <= max_points:
            pdf = work.compute()
        else:
            frac = float(max_points) / float(total)
            pdf = work.sample(frac=frac, random_state=random_state).compute()
        if len(pdf) > max_points:
            pdf = pdf.sample(n=max_points, random_state=random_state)
        return pdf

    pdf = to_pandas(work)
    if len(pdf) > max_points:
        pdf = pdf.sample(n=max_points, random_state=random_state)
    return pdf


def geometry_to_napari_polygons(
    geometries,
    max_shapes: int | None = None,
) -> list[np.ndarray]:
    """Convert shapely Polygon/MultiPolygon geometries to napari polygon arrays."""
    out: list[np.ndarray] = []
    n_added = 0

    for geom in geometries:
        if max_shapes is not None and n_added >= max_shapes:
            break
        if geom is None or geom.is_empty:
            continue

        if geom.geom_type == "Polygon":
            arr = np.asarray(geom.exterior.coords, dtype=np.float32)
            if arr.shape[0] >= 3:
                out.append(arr[:, [1, 0]])  # y, x order for napari
                n_added += 1
            continue

        if geom.geom_type == "MultiPolygon":
            for part in geom.geoms:
                if max_shapes is not None and n_added >= max_shapes:
                    break
                arr = np.asarray(part.exterior.coords, dtype=np.float32)
                if arr.shape[0] >= 3:
                    out.append(arr[:, [1, 0]])  # y, x order
                    n_added += 1

    return out
