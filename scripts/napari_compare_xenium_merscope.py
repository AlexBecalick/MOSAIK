#!/usr/bin/env python3
"""
Interactive Napari viewer for comparing MERSCOPE and Xenium SpatialData outputs.

Features:
- Single viewer with dataset switcher (MERSCOPE/XENIUM).
- All image channels as separate image layers in micron coordinates.
- All shape keys as separate polygon layers with deterministic colors.
- Assigned vs unassigned transcript layers from points['assignment']-style columns.
"""

from __future__ import annotations

import argparse
import colorsys
import gc
import hashlib
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import spatialdata as sd

try:
    import napari
except ImportError:
    print(
        "ERROR: napari is not installed in this environment.\n"
        "Install it with: pip install 'napari[all]' or conda install napari",
        file=sys.stderr,
    )
    sys.exit(1)

try:
    from qtpy.QtWidgets import (
        QComboBox,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QVBoxLayout,
        QWidget,
    )
except ImportError:
    print(
        "ERROR: qtpy is required for the Napari dock widget.\n"
        "Install it with: pip install qtpy",
        file=sys.stderr,
    )
    sys.exit(1)

try:
    import psutil
except Exception:
    psutil = None

try:
    from _napari_compare_utils import (
        assignment_mask,
        build_napari_affine_from_px_to_um,
        channel_labels,
        ensure_cyx,
        first_existing_col,
        geometry_to_napari_polygons,
        get_scale0_dataarray,
        load_points_dataframe,
        resolve_dataset_mask_affine,
    )
except ImportError:
    from scripts._napari_compare_utils import (
        assignment_mask,
        build_napari_affine_from_px_to_um,
        channel_labels,
        ensure_cyx,
        first_existing_col,
        geometry_to_napari_polygons,
        get_scale0_dataarray,
        load_points_dataframe,
        resolve_dataset_mask_affine,
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("napari_compare")


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    zarr_path: Path
    merscope_transform_path: Path | None = None
    xenium_spec_path: Path | None = None


def memory_snapshot_gb() -> dict[str, float]:
    """Return current process/system memory snapshot in GB."""
    if psutil is None:
        return {"rss_gb": float("nan"), "sys_used_gb": float("nan")}

    proc = psutil.Process()
    rss = proc.memory_info().rss / (1024**3)
    used = psutil.virtual_memory().used / (1024**3)
    return {"rss_gb": float(rss), "sys_used_gb": float(used)}


def stable_layer_color(key: str, alpha: float = 1.0) -> tuple[float, float, float, float]:
    """Generate deterministic RGBA color for a layer key."""
    digest = hashlib.md5(key.encode("utf-8")).hexdigest()
    hue = int(digest[:8], 16) / 0xFFFFFFFF
    sat = 0.75
    val = 0.95
    r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
    return (float(r), float(g), float(b), float(alpha))


def image_colormap_for_channel(channel_name: str, fallback_index: int = 0) -> str:
    """Select a napari colormap based on channel name."""
    name = str(channel_name).lower()
    if "dapi" in name:
        return "blue"
    if ("polyt" in name) or ("18s" in name) or ("rna" in name):
        return "green"

    fallback = ["gray", "magenta", "cyan", "yellow", "red", "orange"]
    return fallback[fallback_index % len(fallback)]


class DatasetSwitcherWidget(QWidget):
    """Simple dock widget with dataset selector + reload button."""

    def __init__(self, datasets: list[str], load_callback, initial_dataset: str = "MERSCOPE"):
        super().__init__()
        self._load_callback = load_callback

        self._dataset_combo = QComboBox()
        self._dataset_combo.addItems(datasets)
        if initial_dataset in datasets:
            self._dataset_combo.setCurrentText(initial_dataset)
        self._dataset_combo.currentTextChanged.connect(self._on_dataset_changed)

        self._reload_button = QPushButton("Reload Dataset")
        self._reload_button.clicked.connect(self._on_reload_clicked)

        self._status_label = QLabel("Ready")
        self._status_label.setWordWrap(True)

        row = QHBoxLayout()
        row.addWidget(QLabel("Dataset:"))
        row.addWidget(self._dataset_combo, stretch=1)
        row.addWidget(self._reload_button)

        root = QVBoxLayout()
        root.addLayout(row)
        root.addWidget(self._status_label)
        self.setLayout(root)

    @property
    def current_dataset(self) -> str:
        return str(self._dataset_combo.currentText())

    def set_status(self, text: str):
        self._status_label.setText(str(text))

    def _on_dataset_changed(self, text: str):
        if not text:
            return
        self._load_callback(str(text), False)

    def _on_reload_clicked(self):
        self._load_callback(self.current_dataset, True)


class ComparisonViewerController:
    """Coordinate loading/clearing napari layers for each dataset."""

    def __init__(self, viewer: napari.Viewer, datasets: dict[str, DatasetConfig], args):
        self.viewer = viewer
        self.datasets = datasets
        self.args = args
        self.active_dataset: str | None = None
        self._status_callback = None

    def set_status_callback(self, fn):
        self._status_callback = fn

    def _set_status(self, text: str):
        if self._status_callback is not None:
            self._status_callback(text)

    def _clear_layers(self):
        for layer in list(self.viewer.layers):
            self.viewer.layers.remove(layer)

    def _resolve_optional_transform_paths(self, ds: str, cfg: DatasetConfig) -> tuple[Path | None, Path | None]:
        """Resolve optional transform/spec paths from explicit args or common defaults."""
        merscope_transform_path = cfg.merscope_transform_path
        xenium_spec_path = cfg.xenium_spec_path

        if ds == "MERSCOPE" and merscope_transform_path is None:
            candidate = cfg.zarr_path / "micron_to_mosaic_pixel_transform.csv"
            if candidate.exists():
                merscope_transform_path = candidate

        if ds == "XENIUM" and xenium_spec_path is None:
            candidates = [
                cfg.zarr_path / "experiment.xenium",
                cfg.zarr_path.parent / "experiment.xenium",
            ]
            for candidate in candidates:
                if candidate.exists():
                    xenium_spec_path = candidate
                    break

        return merscope_transform_path, xenium_spec_path

    def load_dataset(self, dataset_name: str, force: bool = False):
        ds = str(dataset_name).upper()
        if ds not in self.datasets:
            self._set_status(f"Unknown dataset: {dataset_name}")
            return

        if (self.active_dataset == ds) and (not force):
            return

        mem_before = memory_snapshot_gb()
        t0 = time.time()
        cfg = self.datasets[ds]

        try:
            self._set_status(f"Loading {ds}...")
            log.info("[%s] Loading dataset from %s", ds, cfg.zarr_path)

            self._clear_layers()
            gc.collect()

            sdata = sd.read_zarr(str(cfg.zarr_path))
            ms_tf_path, xe_spec_path = self._resolve_optional_transform_paths(ds, cfg)
            x_transform, y_transform = resolve_dataset_mask_affine(
                ds,
                merscope_transform_path=ms_tf_path,
                xenium_spec_path=xe_spec_path,
            )
            log.info("[%s] Using image px->um transform x=%s y=%s", ds, x_transform, y_transform)

            img_stats = self._add_image_layers(ds, sdata, x_transform, y_transform)
            shape_stats = self._add_shape_layers(ds, sdata)
            tx_stats = self._add_transcript_layers(ds, sdata)

            self.active_dataset = ds
            elapsed = time.time() - t0
            mem_after = memory_snapshot_gb()

            image_summary = f"images={img_stats['layers']}"
            if img_stats.get("skipped", False):
                image_summary = "images=skipped"
            elif img_stats.get("failed_keys", 0) > 0:
                image_summary = (
                    f"images={img_stats['layers']} "
                    f"(failed_keys={img_stats['failed_keys']})"
                )

            summary = (
                f"{ds} loaded in {elapsed:.1f}s | "
                f"{image_summary} | "
                f"shape_layers={shape_stats['layers']} polygons={shape_stats['polygons']:,} | "
                f"tx_total={tx_stats['total']:,} assigned={tx_stats['assigned']:,} "
                f"unassigned={tx_stats['unassigned']:,} | "
                f"RSS {mem_before['rss_gb']:.1f}->{mem_after['rss_gb']:.1f} GB"
            )
            self._set_status(summary)
            log.info(summary)

        except Exception as exc:
            self._set_status(f"Failed to load {ds}: {exc}")
            log.exception("[%s] Failed to load dataset", ds)

    def _add_image_layers(self, ds: str, sdata, x_transform, y_transform) -> dict[str, int]:
        if getattr(self.args, "skip_images", False):
            log.info("[%s] Image loading skipped (--skip-images).", ds)
            return {"layers": 0, "failed_keys": 0, "skipped": True}

        visible = not self.args.hide_images
        total_layers = 0
        failed_keys = 0

        try:
            image_keys = list(sdata.images.keys())
        except Exception as exc:
            log.warning("[%s] Could not enumerate images; skipping image loading (%s)", ds, exc)
            return {"layers": 0, "failed_keys": 0, "skipped": True}

        if len(image_keys) == 0:
            log.info("[%s] No images found in SpatialData; continuing without image layers.", ds)
            return {"layers": 0, "failed_keys": 0, "skipped": True}

        for image_key in image_keys:
            try:
                da = get_scale0_dataarray(sdata.images[image_key])
                image_cyx = ensure_cyx(da)
                labels = channel_labels(image_cyx)

                x_coords = (
                    np.asarray(image_cyx.coords["x"].values)
                    if "x" in image_cyx.coords
                    else None
                )
                y_coords = (
                    np.asarray(image_cyx.coords["y"].values)
                    if "y" in image_cyx.coords
                    else None
                )
                affine = build_napari_affine_from_px_to_um(
                    x_transform=x_transform,
                    y_transform=y_transform,
                    x_coords=x_coords,
                    y_coords=y_coords,
                )

                for chan_idx, chan_name in enumerate(labels):
                    layer_name = f"{ds} | image | {image_key} | {chan_name}"
                    ch_data = image_cyx.isel(c=chan_idx).data
                    cmap = image_colormap_for_channel(chan_name, chan_idx)

                    self.viewer.add_image(
                        ch_data,
                        name=layer_name,
                        affine=affine,
                        colormap=cmap,
                        blending="additive",
                        opacity=self.args.image_opacity,
                        visible=visible,
                    )
                    total_layers += 1
            except Exception as exc:
                failed_keys += 1
                log.warning(
                    "[%s] Skipping image key '%s' due to load error (%s)",
                    ds,
                    image_key,
                    exc,
                )
                continue

        return {"layers": total_layers, "failed_keys": failed_keys, "skipped": False}

    def _add_shape_layers(self, ds: str, sdata) -> dict[str, int]:
        visible = not self.args.hide_shapes
        total_layers = 0
        total_polygons = 0

        for shape_key in sorted(sdata.shapes.keys()):
            gdf = sdata.shapes[shape_key]
            polygons = geometry_to_napari_polygons(
                gdf.geometry,
                max_shapes=self.args.max_shapes_per_layer,
            )
            if len(polygons) == 0:
                log.info("[%s] shapes[%s] has no drawable polygons", ds, shape_key)
                continue

            edge_color = stable_layer_color(shape_key, alpha=self.args.shape_opacity)
            layer_name = f"{ds} | shapes | {shape_key}"

            self.viewer.add_shapes(
                polygons,
                shape_type="polygon",
                name=layer_name,
                edge_color=edge_color,
                edge_width=self.args.shape_edge_width,
                face_color="transparent",
                visible=visible,
            )
            total_layers += 1
            total_polygons += len(polygons)

            limit_note = (
                f", capped at {self.args.max_shapes_per_layer:,}"
                if self.args.max_shapes_per_layer is not None
                else ""
            )
            log.info(
                "[%s] Added shape layer %s with %s polygons%s",
                ds,
                shape_key,
                f"{len(polygons):,}",
                limit_note,
            )

        return {"layers": total_layers, "polygons": total_polygons}

    def _add_transcript_layers(self, ds: str, sdata) -> dict[str, int]:
        if len(sdata.points) == 0:
            log.warning("[%s] No points available in SpatialData.", ds)
            return {"total": 0, "assigned": 0, "unassigned": 0}

        points_key = list(sdata.points.keys())[0]
        points_obj = sdata.points[points_key]

        x_col = first_existing_col(points_obj, ["x", "x_micron", "global_x", "x_location", "observed_x"])
        y_col = first_existing_col(points_obj, ["y", "y_micron", "global_y", "y_location", "observed_y"])
        assignment_col = first_existing_col(points_obj, ["assignment", "cell", "cell_id"])

        if x_col is None or y_col is None:
            raise KeyError(f"Could not resolve x/y columns in points[{points_key}]")

        points_pdf = load_points_dataframe(
            points_obj,
            x_col=x_col,
            y_col=y_col,
            assignment_col=assignment_col,
            max_points=self.args.max_transcripts,
            random_state=self.args.random_state,
        )

        if assignment_col is not None and assignment_col in points_pdf.columns:
            assigned_mask = assignment_mask(points_pdf[assignment_col]).to_numpy(dtype=bool, copy=False)
        else:
            assigned_mask = np.ones(len(points_pdf), dtype=bool)

        x_vals = points_pdf[x_col].to_numpy(dtype=np.float32, copy=False)
        y_vals = points_pdf[y_col].to_numpy(dtype=np.float32, copy=False)

        assigned_coords = np.column_stack([y_vals[assigned_mask], x_vals[assigned_mask]])
        unassigned_coords = np.column_stack([y_vals[~assigned_mask], x_vals[~assigned_mask]])

        visible = not self.args.hide_transcripts

        if len(unassigned_coords) > 0:
            self.viewer.add_points(
                unassigned_coords,
                name=f"{ds} | transcripts | unassigned",
                face_color=self.args.unassigned_color,
                edge_color=self.args.unassigned_color,
                size=self.args.point_size,
                opacity=self.args.point_opacity,
                visible=visible,
            )

        if len(assigned_coords) > 0:
            self.viewer.add_points(
                assigned_coords,
                name=f"{ds} | transcripts | assigned",
                face_color=self.args.assigned_color,
                edge_color=self.args.assigned_color,
                size=self.args.point_size,
                opacity=self.args.point_opacity,
                visible=visible,
            )

        total = int(len(points_pdf))
        assigned = int(assigned_mask.sum())
        unassigned = total - assigned

        limit_note = ""
        if self.args.max_transcripts is not None:
            limit_note = f" (sampled/capped to {self.args.max_transcripts:,})"

        log.info(
            "[%s] Added transcripts from points[%s]: total=%s assigned=%s unassigned=%s%s",
            ds,
            points_key,
            f"{total:,}",
            f"{assigned:,}",
            f"{unassigned:,}",
            limit_note,
        )

        return {"total": total, "assigned": assigned, "unassigned": unassigned}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Napari comparison viewer for MERSCOPE and Xenium SpatialData outputs."
    )
    parser.add_argument("--merscope-zarr", required=True, type=Path, help="Path to MERSCOPE .zarr output")
    parser.add_argument("--xenium-zarr", required=True, type=Path, help="Path to Xenium .zarr output")

    parser.add_argument(
        "--merscope-transform-path",
        default=None,
        type=Path,
        help="Optional MERSCOPE micron_to_mosaic_pixel_transform.csv path",
    )
    parser.add_argument(
        "--xenium-spec-path",
        default=None,
        type=Path,
        help="Optional Xenium experiment.xenium path (for pixel_size)",
    )

    parser.add_argument(
        "--initial-dataset",
        default="MERSCOPE",
        choices=["MERSCOPE", "XENIUM"],
        help="Dataset loaded first when the viewer opens",
    )

    parser.add_argument(
        "--max-transcripts",
        default=None,
        type=int,
        help="Optional cap/sampling limit for transcripts per dataset load (default: full)",
    )
    parser.add_argument(
        "--max-shapes-per-layer",
        default=None,
        type=int,
        help="Optional cap for polygons per shape layer (default: full)",
    )
    parser.add_argument("--random-state", default=42, type=int, help="Random seed used for sampling")

    parser.add_argument("--point-size", default=2.0, type=float, help="Napari point size for transcript layers")
    parser.add_argument("--point-opacity", default=0.50, type=float, help="Napari point opacity")
    parser.add_argument("--shape-edge-width", default=0.75, type=float, help="Shape layer edge width")
    parser.add_argument("--shape-opacity", default=0.95, type=float, help="Shape layer edge alpha")
    parser.add_argument("--image-opacity", default=1.0, type=float, help="Image layer opacity")
    parser.add_argument("--assigned-color", default="yellow", help="Color for assigned transcript points")
    parser.add_argument("--unassigned-color", default="#d62728", help="Color for unassigned transcript points")

    parser.add_argument("--hide-images", action="store_true", help="Start with image layers hidden")
    parser.add_argument("--hide-shapes", action="store_true", help="Start with shape layers hidden")
    parser.add_argument("--hide-transcripts", action="store_true", help="Start with transcript layers hidden")
    parser.add_argument(
        "--skip-images",
        action="store_true",
        help="Do not attempt to load image layers (useful for image-less zarr transfers).",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if not args.merscope_zarr.exists():
        raise FileNotFoundError(f"MERSCOPE zarr path not found: {args.merscope_zarr}")
    if not args.xenium_zarr.exists():
        raise FileNotFoundError(f"Xenium zarr path not found: {args.xenium_zarr}")

    datasets = {
        "MERSCOPE": DatasetConfig(
            name="MERSCOPE",
            zarr_path=args.merscope_zarr,
            merscope_transform_path=args.merscope_transform_path,
            xenium_spec_path=args.xenium_spec_path,
        ),
        "XENIUM": DatasetConfig(
            name="XENIUM",
            zarr_path=args.xenium_zarr,
            merscope_transform_path=args.merscope_transform_path,
            xenium_spec_path=args.xenium_spec_path,
        ),
    }

    viewer = napari.Viewer(title="MOSAIK Xenium vs MERSCOPE Comparison")
    controller = ComparisonViewerController(viewer=viewer, datasets=datasets, args=args)

    switcher = DatasetSwitcherWidget(
        datasets=["MERSCOPE", "XENIUM"],
        load_callback=controller.load_dataset,
        initial_dataset=args.initial_dataset,
    )
    controller.set_status_callback(switcher.set_status)
    viewer.window.add_dock_widget(switcher, area="right", name="Dataset Switcher")

    controller.load_dataset(args.initial_dataset, force=True)
    napari.run()


if __name__ == "__main__":
    main()
