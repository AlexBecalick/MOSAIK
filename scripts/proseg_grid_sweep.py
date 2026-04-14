#!/usr/bin/env python3
"""
ProSeg parameter grid sweep for MERSCOPE resegmentation tuning.

Reads a YAML config describing sweep axes and fixed parameters, runs every
combination of the sweep grid through ProSeg, then generates image-overlay
grid plots for each pair of sweep axes.

Usage:
    python scripts/proseg_grid_sweep.py scripts/proseg_sweep_config.yaml
    python scripts/proseg_grid_sweep.py scripts/proseg_sweep_config.yaml --plot-only
"""

import argparse
import itertools
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _sweep_utils import (
    REPO,
    draw_boundaries,
    load_config,
    load_refined_shapes,
    prepare_data,
)

sys.path.insert(0, str(REPO / "src" / "resegmentation"))
from proseg_wrapper import run_proseg_refinement  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("proseg_sweep")


# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------

def build_param_grid(config: dict) -> Tuple[List[dict], Dict[str, np.ndarray]]:
    """Return (list-of-param-dicts, {axis_name: values_array})."""
    sweep = config["sweep"]
    fixed = config["fixed_params"]

    axes: Dict[str, np.ndarray] = {}
    for name, spec in sweep.items():
        vals = np.arange(spec["start"], spec["stop"] + spec["step"] / 2, spec["step"])
        axes[name] = np.round(vals, 10)

    axis_names = list(axes.keys())
    axis_values = [axes[n] for n in axis_names]

    param_grid: List[dict] = []
    for combo in itertools.product(*axis_values):
        params = dict(zip(axis_names, [float(v) for v in combo]))
        params["samples"] = int(fixed["samples"])
        params["voxel_size"] = float(fixed["voxel_size"])
        params["burnin_voxel_size"] = float(fixed["burnin_voxel_size"])
        params["num_threads"] = int(fixed.get("num_threads", 64))
        parts = [f"{k[:4]}_{v:g}" for k, v in zip(axis_names, combo)]
        params["label"] = "__".join(parts)
        param_grid.append(params)

    return param_grid, axes


def _label_to_dir(label: str) -> str:
    return label.replace(" ", "_")


# ---------------------------------------------------------------------------
# Sweep runner
# ---------------------------------------------------------------------------

def run_sweep(
    config: dict,
    param_grid: List[dict],
    data: dict,
) -> pd.DataFrame:
    """Run ProSeg for every parameter combination. Returns summary DataFrame."""
    output_dir = Path(config["paths"]["output_dir"])
    proseg_binary = config["paths"]["proseg_binary"]
    output_dir.mkdir(parents=True, exist_ok=True)
    results_csv = output_dir / "results.csv"

    existing_results = {}
    if results_csv.exists():
        prev = pd.read_csv(results_csv)
        for _, row in prev.iterrows():
            existing_results[row["label"]] = row.to_dict()

    trans = data["trans"]
    voxel_layers = data["voxel_layers"]
    cellpose_masks = data["cellpose_masks"]
    cx = data["cellpose_x_transform"]
    cy = data["cellpose_y_transform"]

    results: List[dict] = []
    total = len(param_grid)
    cumulative_seconds = 0.0

    for idx, params in enumerate(param_grid, 1):
        label = params["label"]
        run_dir = output_dir / "outputs" / _label_to_dir(label)
        zarr_path = run_dir / "proseg_output.zarr"

        if label in existing_results and zarr_path.exists():
            log.info("[%d/%d] SKIP (exists): %s", idx, total, label)
            results.append(existing_results[label])
            continue

        run_dir.mkdir(parents=True, exist_ok=True)
        log.info("[%d/%d] Running: %s", idx, total, label)
        t0 = time.time()

        try:
            run_proseg_refinement(
                transcripts_df=trans,
                output_path=zarr_path,
                proseg_binary=proseg_binary,
                x_col="x_micron",
                y_col="y_micron",
                z_col="z_micron",
                gene_col="feature_name",
                cell_id_col="cell_id",
                samples=params["samples"],
                burnin_voxel_size=params["burnin_voxel_size"],
                voxel_size=params["voxel_size"],
                voxel_layers=voxel_layers,
                nuclear_reassignment_prob=params["nuclear_reassignment_prob"],
                diffusion_probability=params["diffusion_probability"],
                cell_compactness=params.get("cell_compactness"),
                cellpose_masks=cellpose_masks,
                cellpose_x_transform=cx,
                cellpose_y_transform=cy,
                diffusion_sigma_far=params.get("diffusion_sigma_far"),
                num_threads=params["num_threads"],
                overwrite=True,
                logger=None,
            )
            elapsed = time.time() - t0
            cumulative_seconds += elapsed

            n_cells, total_tx, mean_tx = _extract_metrics(zarr_path)
            status = "ok"
        except Exception as e:
            elapsed = time.time() - t0
            cumulative_seconds += elapsed
            log.error("FAILED %s: %s", label, e)
            n_cells, total_tx, mean_tx = np.nan, np.nan, np.nan
            status = f"error: {e}"

        rec = {
            "label": label,
            "status": status,
            "n_cells": n_cells,
            "total_transcripts": total_tx,
            "mean_tx_per_cell": mean_tx,
            "elapsed_s": round(elapsed, 1),
            "zarr_path": str(zarr_path),
        }
        for k in config["sweep"]:
            rec[k] = params.get(k)
        results.append(rec)

        pd.DataFrame(results).to_csv(results_csv, index=False)

        remaining = total - idx
        avg = cumulative_seconds / idx
        eta_h = remaining * avg / 3600
        log.info(
            "  done %s  cells=%s  elapsed=%.0fs  avg=%.0fs  ETA=%.1fh",
            label,
            n_cells if not np.isnan(n_cells) else "?",
            elapsed,
            avg,
            eta_h,
        )

    df = pd.DataFrame(results)
    df.to_csv(results_csv, index=False)
    log.info("Results saved to %s", results_csv)
    return df


def _extract_metrics(zarr_path: Path) -> Tuple[float, float, float]:
    """Pull cell count + transcript stats from a ProSeg zarr output."""
    import spatialdata as sd

    try:
        sdata_out = sd.read_zarr(zarr_path)
        shape_key = list(sdata_out.shapes.keys())[0]
        shapes = sdata_out.shapes[shape_key]
        n_cells = len(shapes)
    except Exception:
        n_cells = np.nan

    try:
        from anndata import read_zarr as read_anndata_zarr

        table_path = zarr_path / "tables" / "table"
        if table_path.exists():
            adata = read_anndata_zarr(table_path)
            total_tx = float(adata.X.sum())
            mean_tx = float(total_tx / max(adata.n_obs, 1))
        else:
            total_tx, mean_tx = np.nan, np.nan
    except Exception:
        total_tx, mean_tx = np.nan, np.nan

    return n_cells, total_tx, mean_tx


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def generate_plots(
    config: dict,
    results_df: pd.DataFrame,
    axes: Dict[str, np.ndarray],
    data: dict,
):
    """Generate all comparison image-overlay grid plots."""
    output_dir = Path(config["paths"]["output_dir"])
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    ok = results_df[results_df["status"] == "ok"].copy()
    if ok.empty:
        log.warning("No successful runs to plot.")
        return

    axis_names = list(axes.keys())

    ref_values = {}
    for name, vals in axes.items():
        ref_values[name] = vals[len(vals) // 2]
    log.info("Reference values for grid slices: %s", ref_values)

    pairs = list(itertools.combinations(axis_names, 2))
    for ax_a, ax_b in pairs:
        _plot_image_grid_pair(ok, axes, ax_a, ax_b, ref_values, data, plot_dir)

    log.info("All plots saved to %s", plot_dir)


def _plot_image_grid_pair(
    df: pd.DataFrame,
    axes: Dict[str, np.ndarray],
    ax_a: str,
    ax_b: str,
    ref_values: Dict[str, float],
    data: dict,
    plot_dir: Path,
):
    axis_names = list(axes.keys())
    other_axes = [n for n in axis_names if n not in (ax_a, ax_b)]

    subset = df.copy()
    for oa in other_axes:
        closest = axes[oa][np.argmin(np.abs(axes[oa] - ref_values[oa]))]
        subset = subset[np.isclose(subset[oa], closest)]

    if subset.empty:
        log.warning("No runs for grid plot %s vs %s", ax_a, ax_b)
        return

    vals_a = np.sort(subset[ax_a].unique())
    vals_b = np.sort(subset[ax_b].unique())
    ncols = len(vals_a)
    nrows = len(vals_b)

    M = data["M"]
    x0, y0 = data["x0"], data["y0"]
    img8 = data["img8"]
    cells_local = data["cells_local"]
    pts_zoom = data["pts_zoom"]
    zx0, zx1, zy0, zy1 = data["zoom_bounds"]

    fig, axarr = plt.subplots(
        nrows, ncols,
        figsize=(4 * ncols, 4 * nrows),
        constrained_layout=True,
    )
    if nrows == 1 and ncols == 1:
        axarr = np.array([[axarr]])
    elif nrows == 1:
        axarr = axarr[np.newaxis, :]
    elif ncols == 1:
        axarr = axarr[:, np.newaxis]

    fixed_str = ", ".join(f"{oa}={ref_values[oa]:g}" for oa in other_axes)
    fig.suptitle(f"{ax_a} vs {ax_b}  (fixed: {fixed_str})", fontsize=14)

    for ri, vb in enumerate(vals_b):
        for ci, va in enumerate(vals_a):
            ax = axarr[ri, ci]
            row = subset[np.isclose(subset[ax_a], va) & np.isclose(subset[ax_b], vb)]
            if row.empty or pd.isna(row.iloc[0]["zarr_path"]):
                ax.axis("off")
                ax.set_title(f"{ax_a}={va:g}\n{ax_b}={vb:g}\n(missing)", fontsize=7)
                continue

            zp = Path(row.iloc[0]["zarr_path"])
            n_cells = row.iloc[0]["n_cells"]
            try:
                gdf, _ = load_refined_shapes(zp, M_transform=M, x0_off=x0, y0_off=y0)
            except Exception as e:
                ax.axis("off")
                ax.set_title(f"error: {e}", fontsize=6)
                continue

            ax.imshow(img8[..., 0], cmap="gray")

            if len(cells_local) > 0:
                draw_boundaries(ax, cells_local, color="red", linewidth=0.3)
            if len(gdf) > 0:
                draw_boundaries(ax, gdf, color="cyan", linewidth=0.4)

            if pts_zoom is not None and len(pts_zoom) > 0:
                ax.scatter(
                    pts_zoom["x_local_px"],
                    pts_zoom["y_local_px"],
                    s=1, c="lime", alpha=0.25, linewidths=0, rasterized=True,
                )

            ax.set_xlim(zx0, zx1)
            ax.set_ylim(zy1, zy0)
            ax.set_aspect("equal")
            ax.axis("off")
            n_str = f"{n_cells:.0f}" if not np.isnan(n_cells) else "?"
            ax.set_title(f"{ax_a}={va:g}  {ax_b}={vb:g}\nn={n_str}", fontsize=7)

    fname = plot_dir / f"grid__{ax_a}__vs__{ax_b}.png"
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved %s", fname)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ProSeg parameter grid sweep for MERSCOPE resegmentation tuning."
    )
    parser.add_argument("config", help="Path to YAML config file")
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Skip the sweep, only regenerate plots from existing results.csv",
    )
    parser.add_argument(
        "--estimate",
        action="store_true",
        help="Run only the first combination, print ETA, then exit.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    param_grid, axes = build_param_grid(config)

    total = len(param_grid)
    log.info("Sweep axes: %s", {k: len(v) for k, v in axes.items()})
    log.info("Total combinations: %d", total)
    for name, vals in axes.items():
        log.info("  %s: %s", name, [f"{v:g}" for v in vals])

    if args.estimate:
        log.info("--estimate: running first combination only")
        data = prepare_data(config)
        param_grid = param_grid[:1]
        run_sweep(config, param_grid, data)
        return

    if args.plot_only:
        results_csv = Path(config["paths"]["output_dir"]) / "results.csv"
        if not results_csv.exists():
            log.error("No results.csv found at %s", results_csv)
            sys.exit(1)
        results_df = pd.read_csv(results_csv)
        log.info("Loaded %d results from %s", len(results_df), results_csv)
        data = prepare_data(config)
        generate_plots(config, results_df, axes, data)
        return

    data = prepare_data(config)
    results_df = run_sweep(config, param_grid, data)
    log.info("Sweep complete. Generating plots...")
    generate_plots(config, results_df, axes, data)
    log.info("Done.")


if __name__ == "__main__":
    main()
