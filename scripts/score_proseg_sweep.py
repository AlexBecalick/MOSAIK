#!/usr/bin/env python3
"""
Score ProSeg sweep outputs against human-drawn ground-truth masks.

Loads human annotations (GeoJSON from ``napari_annotate_cells.py``), matches
each human polygon to a Cellpose nucleus via centroid proximity, then
compares the corresponding ProSeg polygons across all sweep conditions.

Metrics per matched pair: IoU, Dice, symmetric mean boundary distance,
Hausdorff distance.  The primary ranking loss is ``1 - mean_IoU``.

Usage:
    python scripts/score_proseg_sweep.py \\
        --config scripts/proseg_sweep_config.yaml \\
        --human-masks notebooks/_roi_checkpoints/roi_x5000_6200_y3500_4700_z0_6_human_masks.geojson
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon, shape

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _sweep_utils import (
    load_config,
    load_checkpoint,
    load_refined_shapes,
    build_roi_image,
    draw_boundaries,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("score_proseg")


# ---------------------------------------------------------------------------
# Human mask loading
# ---------------------------------------------------------------------------

def load_human_masks(geojson_path: str) -> List[Polygon]:
    """Load human-drawn polygons from GeoJSON in local ROI pixel coords."""
    with open(geojson_path) as f:
        data = json.load(f)
    polys = []
    for feat in data.get("features", []):
        geom = shape(feat["geometry"])
        if not geom.is_valid:
            geom = geom.buffer(0)
        if not geom.is_empty:
            polys.append(geom)
    log.info("Loaded %d human polygons from %s", len(polys), geojson_path)
    return polys


# ---------------------------------------------------------------------------
# Nucleus centroid matching
# ---------------------------------------------------------------------------

def compute_nucleus_centroids(filtered_masks: np.ndarray) -> Dict[int, Tuple[float, float]]:
    """Compute centroid (row, col) for each labelled region in the mask array."""
    from skimage.measure import regionprops

    centroids = {}
    for prop in regionprops(filtered_masks):
        centroids[prop.label] = prop.centroid  # (row, col) = (y, x)
    log.info("Found %d nucleus centroids", len(centroids))
    return centroids


def match_human_to_nuclei(
    human_polys: List[Polygon],
    centroids: Dict[int, Tuple[float, float]],
) -> Dict[int, int]:
    """Map each human polygon index to the closest nucleus label.

    Strategy: for each human polygon, find the nucleus whose centroid falls
    inside it.  If none (or multiple), use the nearest centroid to the
    polygon's own centroid.

    Returns {human_idx: nucleus_label}.
    """
    from shapely.geometry import Point

    label_arr = np.array(list(centroids.keys()))
    # centroids are (row, col) = (y, x), but GeoJSON polygons use (x, y)
    centroid_xy = np.array(
        [[centroids[l][1], centroids[l][0]] for l in label_arr]
    )

    mapping: Dict[int, int] = {}
    used_labels = set()

    for hi, hpoly in enumerate(human_polys):
        inside = []
        for li, label in enumerate(label_arr):
            pt = Point(centroid_xy[li])
            if hpoly.contains(pt):
                inside.append((label, li))

        if len(inside) == 1:
            chosen = inside[0][0]
        else:
            hc = np.array([hpoly.centroid.x, hpoly.centroid.y])
            dists = np.linalg.norm(centroid_xy - hc, axis=1)
            # Prefer centroids not already used
            order = np.argsort(dists)
            chosen = None
            for oi in order:
                if label_arr[oi] not in used_labels:
                    chosen = label_arr[oi]
                    break
            if chosen is None:
                chosen = label_arr[order[0]]

        mapping[hi] = int(chosen)
        used_labels.add(int(chosen))

    log.info("Matched %d human polygons to nucleus labels", len(mapping))
    return mapping


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def _sample_boundary(geom: Polygon, n_points: int = 200) -> np.ndarray:
    """Sample evenly-spaced points along a polygon exterior."""
    boundary = geom.exterior
    distances = np.linspace(0, boundary.length, n_points, endpoint=False)
    pts = np.array([boundary.interpolate(d).coords[0] for d in distances])
    return pts


def _symmetric_boundary_distance(a: Polygon, b: Polygon, n_points: int = 200) -> float:
    """Mean of min-distances from boundary A to B and B to A."""
    pts_a = _sample_boundary(a, n_points)
    pts_b = _sample_boundary(b, n_points)

    from scipy.spatial import cKDTree

    tree_b = cKDTree(pts_b)
    d_a2b, _ = tree_b.query(pts_a)
    tree_a = cKDTree(pts_a)
    d_b2a, _ = tree_a.query(pts_b)

    return float((d_a2b.mean() + d_b2a.mean()) / 2.0)


def compute_pair_metrics(human: Polygon, proseg: Polygon) -> dict:
    """Compute IoU, Dice, symmetric boundary distance, Hausdorff for one pair."""
    intersection = human.intersection(proseg).area
    union = human.union(proseg).area
    iou = intersection / union if union > 0 else 0.0
    dice = 2 * intersection / (human.area + proseg.area) if (human.area + proseg.area) > 0 else 0.0
    hausdorff = human.hausdorff_distance(proseg)
    sym_bd = _symmetric_boundary_distance(human, proseg)

    return {
        "iou": iou,
        "dice": dice,
        "boundary_dist": sym_bd,
        "hausdorff": hausdorff,
    }


# ---------------------------------------------------------------------------
# Per-condition scoring
# ---------------------------------------------------------------------------

def _find_proseg_polygon_for_label(
    gdf: gpd.GeoDataFrame,
    label: int,
    filtered_masks: np.ndarray,
    x0: int,
    y0: int,
) -> Optional[Polygon]:
    """Find the ProSeg polygon that best corresponds to a given nucleus label.

    Since ProSeg outputs don't carry the original Cellpose label, we match by
    finding the ProSeg polygon whose centroid is closest to the nucleus
    centroid in local-pixel space.
    """
    from skimage.measure import regionprops
    from shapely.geometry import Point

    props = {p.label: p.centroid for p in regionprops(filtered_masks)}
    if label not in props:
        return None
    nuc_y, nuc_x = props[label]  # local pixel coords

    if gdf.empty:
        return None

    nuc_pt = np.array([nuc_x, nuc_y])
    best_dist = float("inf")
    best_geom = None

    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        c = np.array([geom.centroid.x, geom.centroid.y])
        d = np.linalg.norm(c - nuc_pt)
        if d < best_dist:
            best_dist = d
            best_geom = geom

    return best_geom


def score_condition(
    zarr_path: Path,
    human_polys: List[Polygon],
    matching: Dict[int, int],
    filtered_masks: np.ndarray,
    M: np.ndarray,
    x0: int,
    y0: int,
) -> Tuple[dict, List[dict]]:
    """Score a single ProSeg condition against human masks.

    Returns (aggregate_dict, list_of_pair_dicts).
    """
    gdf, _ = load_refined_shapes(zarr_path, M_transform=M, x0_off=x0, y0_off=y0)

    pair_results = []
    for hi, nuc_label in matching.items():
        proseg_geom = _find_proseg_polygon_for_label(
            gdf, nuc_label, filtered_masks, x0, y0
        )
        if proseg_geom is None:
            continue

        if isinstance(proseg_geom, MultiPolygon):
            proseg_geom = max(proseg_geom.geoms, key=lambda g: g.area)

        human_geom = human_polys[hi]
        metrics = compute_pair_metrics(human_geom, proseg_geom)
        metrics["human_idx"] = hi
        metrics["nuc_label"] = nuc_label
        pair_results.append(metrics)

    if pair_results:
        agg = {
            "mean_iou": np.mean([r["iou"] for r in pair_results]),
            "mean_dice": np.mean([r["dice"] for r in pair_results]),
            "mean_boundary_dist": np.mean([r["boundary_dist"] for r in pair_results]),
            "mean_hausdorff": np.mean([r["hausdorff"] for r in pair_results]),
            "n_matched": len(pair_results),
        }
        agg["loss"] = 1.0 - agg["mean_iou"]
    else:
        agg = {
            "mean_iou": np.nan,
            "mean_dice": np.nan,
            "mean_boundary_dist": np.nan,
            "mean_hausdorff": np.nan,
            "n_matched": 0,
            "loss": np.nan,
        }

    return agg, pair_results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_ranking_bar(scores_df: pd.DataFrame, plot_dir: Path):
    """Bar chart of all conditions ranked by loss (ascending = best first)."""
    df = scores_df.dropna(subset=["loss"]).sort_values("loss")
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(max(12, len(df) * 0.3), 6))
    x = np.arange(len(df))
    bars = ax.bar(x, df["loss"], color="steelblue", edgecolor="none")

    # Highlight top-5
    for i in range(min(5, len(bars))):
        bars[i].set_color("seagreen")

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=90, fontsize=5)
    ax.set_ylabel("Loss (1 - mean IoU)")
    ax.set_title("ProSeg Conditions Ranked by Loss")
    fig.tight_layout()

    fname = plot_dir / "ranking_bar.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    log.info("Saved %s", fname)


def plot_best_worst_overlay(
    scores_df: pd.DataFrame,
    human_polys: List[Polygon],
    matching: Dict[int, int],
    filtered_masks: np.ndarray,
    data: dict,
    plot_dir: Path,
):
    """Side-by-side overlays of best and worst conditions vs human masks."""
    df = scores_df.dropna(subset=["loss"]).sort_values("loss")
    if len(df) < 2:
        return

    best_row = df.iloc[0]
    worst_row = df.iloc[-1]

    img8 = data["img8"]
    M = data["M"]
    x0, y0 = data["x0"], data["y0"]

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    for ax, row, title in [
        (axes[0], best_row, f"BEST: {best_row['label']}\nloss={best_row['loss']:.4f}"),
        (axes[1], worst_row, f"WORST: {worst_row['label']}\nloss={worst_row['loss']:.4f}"),
    ]:
        ax.imshow(img8[..., 0], cmap="gray")

        # Human masks in yellow
        for hi in matching:
            geom = human_polys[hi]
            if geom.geom_type == "Polygon":
                xc, yc = geom.exterior.xy
                ax.plot(xc, yc, color="yellow", linewidth=1.0, label="human" if hi == 0 else "")

        # ProSeg masks in cyan
        zp = Path(row["zarr_path"])
        try:
            gdf, _ = load_refined_shapes(zp, M_transform=M, x0_off=x0, y0_off=y0)
            draw_boundaries(ax, gdf, color="cyan", linewidth=0.6)
        except Exception as e:
            log.warning("Could not load %s: %s", zp, e)

        ax.set_title(title, fontsize=9)
        ax.axis("off")

    fig.tight_layout()
    fname = plot_dir / "best_worst_overlay.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    log.info("Saved %s", fname)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Score ProSeg sweep outputs against human ground-truth masks."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument("--human-masks", required=True, help="Path to human masks GeoJSON")
    args = parser.parse_args()

    config = load_config(args.config)

    # -- Load data --
    zarr_path = Path(config["paths"]["zarr_path"])
    M = np.loadtxt(zarr_path / "micron_to_mosaic_pixel_transform.csv")

    ckpt = load_checkpoint(config)
    filtered_masks = ckpt["filtered_masks"]
    meta = ckpt["meta"]
    x0, y0 = meta["x0"], meta["y0"]

    _, img8 = build_roi_image(config, meta)

    data = {"M": M, "x0": x0, "y0": y0, "img8": img8}

    # -- Load human masks & match to nuclei --
    human_polys = load_human_masks(args.human_masks)
    if not human_polys:
        log.error("No human polygons found. Annotate cells first.")
        sys.exit(1)

    centroids = compute_nucleus_centroids(filtered_masks)
    matching = match_human_to_nuclei(human_polys, centroids)

    # -- Load sweep results --
    output_dir = Path(config["paths"]["output_dir"])
    results_csv = output_dir / "results.csv"
    if not results_csv.exists():
        log.error("No results.csv found at %s. Run the sweep first.", results_csv)
        sys.exit(1)

    results_df = pd.read_csv(results_csv)
    ok_runs = results_df[results_df["status"] == "ok"].copy()
    log.info("Scoring %d successful conditions", len(ok_runs))

    # -- Score each condition --
    sweep_keys = list(config.get("sweep", {}).keys())
    score_rows = []

    for idx, (_, row) in enumerate(ok_runs.iterrows(), 1):
        label = row["label"]
        zp = Path(row["zarr_path"])

        if not zp.exists():
            log.warning("[%d/%d] SKIP (missing zarr): %s", idx, len(ok_runs), label)
            continue

        log.info("[%d/%d] Scoring: %s", idx, len(ok_runs), label)
        try:
            agg, _ = score_condition(
                zp, human_polys, matching, filtered_masks, M, x0, y0
            )
        except Exception as e:
            log.error("  FAILED: %s", e)
            agg = {
                "mean_iou": np.nan, "mean_dice": np.nan,
                "mean_boundary_dist": np.nan, "mean_hausdorff": np.nan,
                "n_matched": 0, "loss": np.nan,
            }

        rec = {"label": label, "zarr_path": str(zp)}
        rec.update(agg)
        for k in sweep_keys:
            rec[k] = row.get(k)
        score_rows.append(rec)

    scores_df = pd.DataFrame(score_rows)

    # -- Save CSV --
    scores_csv = output_dir / "scores.csv"
    scores_df.to_csv(scores_csv, index=False)
    log.info("Scores saved to %s", scores_csv)

    # -- Print top-5 / bottom-5 --
    ranked = scores_df.dropna(subset=["loss"]).sort_values("loss")
    if len(ranked) > 0:
        print("\n=== TOP 5 (lowest loss) ===")
        print(ranked.head(5)[["label", "loss", "mean_iou", "mean_dice", "mean_boundary_dist", "n_matched"]].to_string(index=False))
        print("\n=== BOTTOM 5 (highest loss) ===")
        print(ranked.tail(5)[["label", "loss", "mean_iou", "mean_dice", "mean_boundary_dist", "n_matched"]].to_string(index=False))

    # -- Generate plots --
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_ranking_bar(scores_df, plot_dir)
    plot_best_worst_overlay(
        scores_df, human_polys, matching, filtered_masks, data, plot_dir
    )

    log.info("Done.")


if __name__ == "__main__":
    main()
