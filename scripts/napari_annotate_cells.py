#!/usr/bin/env python3
"""
Napari-based tool for manually annotating cell boundaries on an ROI crop.

Opens the ROI image (DAPI as magenta, PolyT as green) with existing Cellpose
masks shown for reference, and provides a Shapes layer where the user can
draw polygons around perceived cell boundaries.  Polygons are saved to a
GeoJSON file on close and reloaded on the next session.

Usage:
    python scripts/napari_annotate_cells.py scripts/proseg_sweep_config.yaml
    python scripts/napari_annotate_cells.py scripts/proseg_sweep_config.yaml --output my_masks.geojson
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("napari_annotate")

try:
    import napari
except ImportError:
    print(
        "ERROR: napari is not installed in this environment.\n"
        "Install it with:  pip install 'napari[all]'  or  conda install napari\n"
        "Then re-run this script.",
        file=sys.stderr,
    )
    sys.exit(1)

# Ensure scripts/ is on sys.path for _sweep_utils
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _sweep_utils import (
    build_roi_image,
    checkpoint_name,
    load_checkpoint,
    load_config,
)


# ---------------------------------------------------------------------------
# GeoJSON I/O for the Shapes layer
# ---------------------------------------------------------------------------

def _shapes_to_geojson(shapes_data: list, shape_types: list) -> dict:
    """Convert napari Shapes layer data to a GeoJSON FeatureCollection.

    Each polygon is stored with its 1-based integer ID.
    Coordinates are in local ROI pixel space (row=y, col=x).
    """
    features = []
    for idx, (coords, stype) in enumerate(zip(shapes_data, shape_types), 1):
        if stype != "polygon":
            continue
        ring = [[float(c[1]), float(c[0])] for c in coords]
        if ring and ring[0] != ring[-1]:
            ring.append(ring[0])
        feature = {
            "type": "Feature",
            "properties": {"id": idx},
            "geometry": {"type": "Polygon", "coordinates": [ring]},
        }
        features.append(feature)
    return {"type": "FeatureCollection", "features": features}


def _geojson_to_shapes(geojson: dict) -> list:
    """Convert GeoJSON FeatureCollection back to napari polygon arrays.

    Returns a list of (N, 2) arrays in napari row-col order.
    """
    polygons = []
    for feat in geojson.get("features", []):
        geom = feat.get("geometry", {})
        if geom.get("type") != "Polygon":
            continue
        ring = geom["coordinates"][0]
        arr = np.array([[pt[1], pt[0]] for pt in ring], dtype=np.float64)
        polygons.append(arr)
    return polygons


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Annotate cell boundaries on an ROI crop using napari."
    )
    parser.add_argument("config", help="Path to YAML config file")
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Path for the output GeoJSON file.  Defaults to "
            "<checkpoint_dir>/<ckpt_name>_human_masks.geojson"
        ),
    )
    args = parser.parse_args()

    config = load_config(args.config)

    # -- Determine output path --
    ckpt = load_checkpoint(config)
    meta = ckpt["meta"]
    ckpt_base = ckpt["ckpt_base"]

    if args.output:
        geojson_path = Path(args.output)
    else:
        geojson_path = Path(f"{ckpt_base}_human_masks.geojson")
    log.info("GeoJSON path: %s", geojson_path)

    # -- Build ROI image --
    log.info("Building ROI image (this may take a moment)...")
    crop_raw, img8 = build_roi_image(config, meta)
    log.info("img8 shape: %s, crop_raw shape: %s", img8.shape, crop_raw.shape)

    # Separate channels: assume index 0 = DAPI, index 1 = PolyT
    ch_names = meta.get("channel_names_to_use", ["ch0", "ch1"])
    if crop_raw.shape[-1] >= 2:
        dapi = crop_raw[..., 0].astype(np.float32)
        polyt = crop_raw[..., 1].astype(np.float32)
    else:
        dapi = crop_raw[..., 0].astype(np.float32)
        polyt = np.zeros_like(dapi)

    filtered_masks = ckpt["filtered_masks"]

    # -- Launch napari --
    viewer = napari.Viewer(title="Cell Boundary Annotation")

    dapi_name = ch_names[0] if len(ch_names) > 0 else "DAPI"
    polyt_name = ch_names[1] if len(ch_names) > 1 else "PolyT"

    viewer.add_image(
        dapi,
        name=dapi_name,
        colormap="magenta",
        blending="additive",
        contrast_limits=[np.percentile(dapi, 2), np.percentile(dapi, 99.5)],
    )
    viewer.add_image(
        polyt,
        name=polyt_name,
        colormap="green",
        blending="additive",
        contrast_limits=[np.percentile(polyt, 2), np.percentile(polyt, 99.5)],
    )

    viewer.add_labels(
        filtered_masks,
        name="cellpose_masks (reference)",
        opacity=0.25,
    )

    # -- Shapes layer: load existing or create new --
    existing_polygons = []
    if geojson_path.exists():
        log.info("Loading existing annotations from %s", geojson_path)
        with open(geojson_path) as f:
            existing_polygons = _geojson_to_shapes(json.load(f))
        log.info("Loaded %d polygons", len(existing_polygons))

    if existing_polygons:
        shapes_layer = viewer.add_shapes(
            existing_polygons,
            shape_type="polygon",
            name="human_cells",
            edge_color="yellow",
            edge_width=2,
            face_color="transparent",
        )
    else:
        shapes_layer = viewer.add_shapes(
            name="human_cells",
            shape_type="polygon",
            edge_color="yellow",
            edge_width=2,
            face_color="transparent",
        )

    shapes_layer.mode = "add_polygon"

    # -- Auto-save on close --
    def _save_on_close(_event=None):
        data = shapes_layer.data
        types = shapes_layer.shape_type
        if len(data) == 0:
            log.info("No polygons drawn; nothing to save.")
            return
        geojson = _shapes_to_geojson(data, types)
        geojson_path.parent.mkdir(parents=True, exist_ok=True)
        with open(geojson_path, "w") as f:
            json.dump(geojson, f, indent=2)
        log.info("Saved %d polygons to %s", len(geojson["features"]), geojson_path)

    viewer.window._qt_window.destroyed.connect(_save_on_close)

    log.info(
        "Napari is open.  Draw polygons on the 'human_cells' layer.\n"
        "  - Use the polygon tool (top toolbar) to click vertices.\n"
        "  - Close the window when done; annotations auto-save to:\n"
        "    %s",
        geojson_path,
    )

    napari.run()

    # Belt-and-suspenders: also try saving after run() returns
    try:
        _save_on_close()
    except Exception:
        pass


if __name__ == "__main__":
    main()
