"""
Shared helpers for ProSeg sweep, napari annotation, and scoring scripts.

All data-loading, coordinate-transform, and config utilities live here so
that ``proseg_grid_sweep.py``, ``napari_annotate_cells.py``, and
``score_proseg_sweep.py`` share a single implementation.
"""

import json
import logging
import pickle
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src" / "resegmentation"))

log = logging.getLogger("sweep_utils")


# ---------------------------------------------------------------------------
# Coordinate transforms
# ---------------------------------------------------------------------------

def global_to_pixel(x_global, y_global, M):
    """Transform global micron coordinates to mosaic pixel coordinates."""
    ones = np.ones_like(x_global)
    pts = np.vstack([x_global, y_global, ones])
    pix = M @ pts
    return pix[0], pix[1]


def micron_to_local_px(geom, M_transform, x0_off, y0_off):
    """Transform a Shapely geometry from global microns to local ROI pixels."""
    from shapely.ops import transform as shp_transform

    def _tx(x, y):
        arr = np.vstack([x, y, np.ones_like(x)])
        out = M_transform @ arr
        return out[0] - x0_off, out[1] - y0_off

    return shp_transform(_tx, geom)


def draw_boundaries(ax, gdf_like, color="cyan", linewidth=0.5):
    """Plot polygon boundaries from a GeoDataFrame onto a matplotlib Axes."""
    geoms = gdf_like.geometry if hasattr(gdf_like, "geometry") else gdf_like
    for geom in geoms:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "Polygon":
            x, y = geom.exterior.xy
            ax.plot(x, y, color=color, linewidth=linewidth)
        elif geom.geom_type == "MultiPolygon":
            for part in geom.geoms:
                x, y = part.exterior.xy
                ax.plot(x, y, color=color, linewidth=linewidth)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(yaml_path: str) -> dict:
    """Load and validate a sweep YAML config file."""
    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    for k in ("paths", "roi"):
        if k not in config:
            raise ValueError(f"Missing required config key: {k}")
    return config


def checkpoint_name(roi_cfg: dict) -> str:
    """Build the checkpoint base-name string from an ROI config dict."""
    return (
        f"roi_x{roi_cfg['x_min']:.0f}_{roi_cfg['x_max']:.0f}"
        f"_y{roi_cfg['y_min']:.0f}_{roi_cfg['y_max']:.0f}"
        f"_z{roi_cfg['z_start']}_{roi_cfg['z_end']}"
    )


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_checkpoint(config: dict) -> dict:
    """Load checkpoint arrays and metadata (no zarr / SpatialData needed)."""
    ckpt_dir = Path(config["paths"]["checkpoint_dir"])
    roi_cfg = config["roi"]
    base = ckpt_dir / checkpoint_name(roi_cfg)

    arr = np.load(f"{base}_arrays.npz")
    filtered_masks = arr["filtered_masks"]
    masks_arr = arr["masks"]

    cells_local = gpd.read_parquet(f"{base}_cells_local.parquet")

    with open(f"{base}_meta.pkl", "rb") as f:
        meta = pickle.load(f)

    return {
        "filtered_masks": filtered_masks,
        "masks": masks_arr,
        "cells_local": cells_local,
        "meta": meta,
        "ckpt_base": base,
    }


def build_roi_image(config: dict, meta: dict) -> Tuple[np.ndarray, np.ndarray]:
    """Build the z-range max-projection ROI image.

    Returns ``(crop_raw, img8)`` where *crop_raw* is the raw projection
    ``(y, x, c)`` and *img8* is the contrast-stretched uint8 version.
    """
    import dask.array as da
    import spatialdata as sd

    zarr_path = Path(config["paths"]["zarr_path"])
    sdata = sd.read_zarr(zarr_path)

    x0, x1 = meta["x0"], meta["x1"]
    y0, y1 = meta["y0"], meta["y1"]
    z_start, z_end = meta["z_start"], meta["z_end"]
    channel_names_to_use = meta.get(
        "channel_names_to_use", config["roi"].get("channels")
    )

    pat = re.compile(r"^(?P<prefix>.+)_z(?P<z>\d+)$")
    plane_keys = []
    for key in sorted(sdata.images.keys()):
        m = pat.match(key)
        if m:
            plane_keys.append((int(m.group("z")), key))

    selected_keys = [k for z, k in plane_keys if z_start <= z <= z_end]

    lazy_planes = []
    for key in selected_keys:
        img_obj = sdata.images[key]
        img_xr = (
            img_obj["scale0"].ds["image"]
            if hasattr(img_obj, "__contains__") and "scale0" in img_obj
            else img_obj
        )
        roi_xr = img_xr.isel(y=slice(y0, y1), x=slice(x0, x1))
        c_all = (
            [str(c) for c in img_xr.coords["c"].values]
            if "c" in img_xr.coords
            else None
        )
        if channel_names_to_use and c_all:
            keep = [c for c in channel_names_to_use if c in c_all]
            if keep:
                roi_xr = roi_xr.sel(c=keep)
        lazy_planes.append(roi_xr.transpose("y", "x", "c").data)

    crop_raw = da.max(da.stack(lazy_planes, axis=0), axis=0).compute()

    img = crop_raw.astype(np.float32)
    p2, p98 = np.percentile(img, (2, 98))
    img = np.clip((img - p2) / (p98 - p2 + 1e-8), 0, 1)
    img8 = (img * 255).astype(np.uint8)
    if img8.ndim == 2:
        img8 = np.stack([img8] * 3, axis=-1)
    elif img8.shape[-1] == 1:
        img8 = np.repeat(img8, 3, axis=-1)
    elif img8.shape[-1] == 2:
        img8 = np.concatenate([img8, np.zeros_like(img8[..., :1])], axis=-1)
    elif img8.shape[-1] > 3:
        img8 = img8[..., :3]

    return crop_raw, img8


def prepare_data(config: dict) -> dict:
    """Full data preparation: checkpoint + image + transcripts + transforms.

    This is the heavy-weight loader used by the sweep runner and scoring
    scripts.  For lighter use (napari annotation) call ``load_checkpoint``
    and ``build_roi_image`` directly.
    """
    import spatialdata as sd

    paths = config["paths"]
    roi_cfg = config["roi"]

    zarr_path = Path(paths["zarr_path"])
    M = np.loadtxt(zarr_path / "micron_to_mosaic_pixel_transform.csv")
    Minv = np.linalg.inv(M)
    log.info("M (micron->pixel):\n%s", M)

    ckpt = load_checkpoint(config)
    filtered_masks = ckpt["filtered_masks"]
    masks_arr = ckpt["masks"]
    cells_local = ckpt["cells_local"]
    meta = ckpt["meta"]

    x0, x1 = meta["x0"], meta["x1"]
    y0, y1 = meta["y0"], meta["y1"]
    poly_scale = meta["poly_scale"]
    points_key = meta["points_key"]

    log.info(
        "Checkpoint: masks %s, %d cells, pixel bbox (%d,%d,%d,%d)",
        filtered_masks.shape, len(cells_local), x0, x1, y0, y1,
    )

    _, img8 = build_roi_image(config, meta)
    log.info("img8 shape: %s", img8.shape)

    # --- Transcripts ---
    log.info("Preparing transcripts")
    sdata = sd.read_zarr(zarr_path)
    ROI = roi_cfg
    pts_full = sdata.points[points_key][["x", "y", "z", "gene"]].compute()
    pts_roi = pts_full[
        (pts_full["x"] >= ROI["x_min"])
        & (pts_full["x"] <= ROI["x_max"])
        & (pts_full["y"] >= ROI["y_min"])
        & (pts_full["y"] <= ROI["y_max"])
    ].copy()

    xpix, ypix = global_to_pixel(
        pts_roi["x"].to_numpy(), pts_roi["y"].to_numpy(), M
    )
    pts_roi["x_local_px"] = xpix - x0
    pts_roi["y_local_px"] = ypix - y0

    pts_gdf = gpd.GeoDataFrame(
        pts_roi,
        geometry=gpd.points_from_xy(pts_roi["x_local_px"], pts_roi["y_local_px"]),
    )

    if len(cells_local) > 0:
        join = gpd.sjoin(
            pts_gdf,
            cells_local[["cell_id", "geometry"]],
            how="left",
            predicate="within",
        )
        pts_gdf["cell_id"] = (
            join.groupby(join.index)["cell_id"].first().reindex(pts_gdf.index)
        )
    else:
        pts_gdf["cell_id"] = np.nan

    assigned = pts_gdf["cell_id"].notna().sum()
    log.info(
        "Transcripts: %d total, %d assigned (%.1f%%)",
        len(pts_gdf), assigned, 100 * assigned / max(len(pts_gdf), 1),
    )

    z_vals = pd.to_numeric(pts_gdf["z"], errors="coerce").dropna().unique()
    voxel_layers = int(len(np.unique(z_vals))) if len(z_vals) > 0 else 1
    log.info("voxel_layers: %d", voxel_layers)

    trans = pts_gdf[["x", "y", "z", "gene", "cell_id"]].copy()
    trans["x_micron"] = trans["x"].astype(float)
    trans["y_micron"] = trans["y"].astype(float)
    trans["feature_name"] = trans["gene"].astype(str)
    trans["cell_id"] = trans["cell_id"].fillna("0").astype(str)
    trans["z_micron"] = (
        pd.to_numeric(trans["z"], errors="coerce").fillna(0).astype(float)
    )

    cellpose_masks = filtered_masks if filtered_masks is not None else masks_arr
    sf = float(poly_scale) if poly_scale > 1.0 else 1.0
    S = np.array(
        [[sf, 0.0, float(x0)], [0.0, sf, float(y0)], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    T = Minv @ S
    cellpose_x_transform = (float(T[0, 0]), float(T[0, 1]), float(T[0, 2]))
    cellpose_y_transform = (float(T[1, 0]), float(T[1, 1]), float(T[1, 2]))
    log.info("Cellpose x_transform: %s", cellpose_x_transform)
    log.info("Cellpose y_transform: %s", cellpose_y_transform)

    h, w = img8.shape[:2]
    zoom_w, zoom_h = min(2200, w), min(2200, h)
    cx, cy = w // 2, h // 2
    zx0 = max(0, cx - zoom_w // 2)
    zx1 = min(w, cx + zoom_w // 2)
    zy0 = max(0, cy - zoom_h // 2)
    zy1 = min(h, cy + zoom_h // 2)

    pts_zoom = pts_gdf[
        (pts_gdf["x_local_px"] >= zx0)
        & (pts_gdf["x_local_px"] <= zx1)
        & (pts_gdf["y_local_px"] >= zy0)
        & (pts_gdf["y_local_px"] <= zy1)
    ]

    return {
        "M": M,
        "Minv": Minv,
        "x0": x0, "y0": y0, "x1": x1, "y1": y1,
        "img8": img8,
        "filtered_masks": filtered_masks,
        "cells_local": cells_local,
        "pts_gdf": pts_gdf,
        "pts_zoom": pts_zoom,
        "trans": trans,
        "cellpose_masks": cellpose_masks,
        "cellpose_x_transform": cellpose_x_transform,
        "cellpose_y_transform": cellpose_y_transform,
        "voxel_layers": voxel_layers,
        "zoom_bounds": (zx0, zx1, zy0, zy1),
    }


# ---------------------------------------------------------------------------
# ProSeg output loading
# ---------------------------------------------------------------------------

def load_refined_shapes(
    path: Path,
    M_transform=None,
    x0_off: float = 0,
    y0_off: float = 0,
) -> Tuple[gpd.GeoDataFrame, str]:
    """Load ProSeg cell boundaries and optionally transform to local pixels."""
    import spatialdata as sd

    sdata_out = sd.read_zarr(path)
    if len(sdata_out.shapes) == 0:
        raise ValueError(f"No shapes in {path}")
    shape_key = list(sdata_out.shapes.keys())[0]
    shapes = sdata_out.shapes[shape_key]

    try:
        gdf = shapes[["geometry"]].copy()
    except Exception:
        gdf = gpd.GeoDataFrame({"geometry": shapes.geometry})

    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()

    if M_transform is not None:
        gdf["geometry"] = gdf["geometry"].apply(
            lambda g: micron_to_local_px(g, M_transform, x0_off, y0_off)
        )
    return gdf, shape_key
