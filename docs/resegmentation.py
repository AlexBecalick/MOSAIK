#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import gc
import numpy as np
import pandas as pd
import geopandas as gpd
import anndata as ad
from PIL import Image
from collections import Counter, OrderedDict
from shapely.geometry import Polygon
from shapely.affinity import scale
from skimage.measure import label, regionprops, find_contours
from skimage.transform import resize
import tensorflow as tf
from dask.distributed import Client
from cellpose import io, models
import spatialdata as sd
from spatialdata.models import Labels2DModel, ShapesModel


class Resegmentation:
    """Pipeline for spatial omics segmentation and data integration.

    This class encapsulates preprocessing, segmentation, mask filtering,
    shape/label creation, transcript assignment, and writing results
    back to a SpatialData Zarr file.
    """

    def __init__(self, zarr_dir: str, zarr_name: str, fov: int):
        self.zarr_dir = zarr_dir
        self.zarr_name = zarr_name
        self.fov = fov
        self.sdata = None
        self.image = None
        self.masks = None
        self.flows = None
        self.styles = None
        self.gdf_polygons = None
        self.shapes_model = None
        self.labels_model = None
        self.gdf_points = None
        self.vdata = None

    # ---------------- GPU UTILITIES ---------------- #
    @staticmethod
    def check_gpu():
        """Print Python environment, TensorFlow version, and GPU availability."""
        print("Python path:", sys.executable)
        print("TensorFlow version:", tf.__version__)
        print("GPU available:", tf.config.list_physical_devices("GPU"))

    @staticmethod
    def clear_gpu_memory():
        """Clear TensorFlow GPU memory and trigger garbage collection."""
        try:
            if tf.config.list_physical_devices("GPU"):
                tf.keras.backend.clear_session()
            gc.collect()
        except Exception as e:
            print(f"Error clearing GPU memory: {e}")

    # ---------------- IMAGE PREPROCESS ---------------- #
    def preprocess_image(
        self,
        channel_names: list[str],
        channels_to_use: list[str] = ["Y", "U"],
        output_path: str | None = None,
        thumbnail_size: tuple[int, int] = (768, 768),
    ):
        """Load and preprocess Zarr image for segmentation."""
        zarr_path = os.path.join(self.zarr_dir, self.zarr_name)
        self.sdata = sd.read_zarr(zarr_path)

        img_xr = self.sdata.images[f"{self.fov}_image"]
        img = img_xr.transpose("y", "x", "c").data.compute()

        channel_indices = [channel_names.index(ch) for ch in channels_to_use]
        img = img[..., channel_indices]

        p2, p98 = np.percentile(img, (2, 98))
        img_stretched = np.clip((img - p2) / (p98 - p2), 0, 1).astype(np.float32)
        img_8bit = (img_stretched * 255).astype(np.uint8)

        if img_8bit.ndim == 3 and img_8bit.shape[-1] == 2:
            zero_channel = np.zeros_like(img_8bit[..., :1])
            img_8bit = np.concatenate([img_8bit, zero_channel], axis=-1)
        elif img_8bit.ndim == 2:
            img_8bit = np.stack([img_8bit] * 3, axis=-1)

        seg_pil = Image.fromarray(img_8bit, mode="RGB")
        seg_pil.thumbnail(thumbnail_size, Image.Resampling.LANCZOS)

        if output_path is not None:
            seg_pil.save(output_path)

        self.image = seg_pil
        return self.sdata

    # ---------------- SEGMENTATION ---------------- #
    def run_cellpose(
        self,
        img_path: str,
        model_type: str = "cyto3",
        gpu: bool = True,
        channels: list[int] = [0, 0],
        diameter: float | None = None,
        flow_threshold: float = 1,
        cellprob_threshold: float = -3,
    ):
        """Run Cellpose-SAM segmentation on the preprocessed image."""
        print(img_path)
        img = io.imread(img_path)
        model = models.CellposeModel(model_type=model_type, gpu=gpu)

        self.masks, self.flows, self.styles = model.eval(
            img,
            diameter=diameter,
            channels=channels,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )

        self.clear_gpu_memory()
        return self.masks, self.flows, self.styles

    # ---------------- MASKS & SHAPES ---------------- #
    @staticmethod
    def filter_cell_by_regionprops(seg_masks, max_eccentricity=0.95):
        """Filter segmented cell masks by size and eccentricity."""
        labeled = label(seg_masks)
        regions = regionprops(labeled)
        if not regions:
            return np.zeros_like(seg_masks, dtype=np.int32)

        areas = [r.area for r in regions]
        min_area = np.median(areas)
        cleaned_mask = np.zeros_like(seg_masks, dtype=np.int32)
        current_label = 1

        for region in regions:
            if region.area < min_area:
                continue
            if region.eccentricity > max_eccentricity:
                continue
            cleaned_mask[labeled == region.label] = current_label
            current_label += 1

        return cleaned_mask

    @staticmethod
    def masks_to_polygons(seg_masks):
        """Convert segmentation masks into scaled polygons."""
        polygons = []
        for region in regionprops(seg_masks):
            mask = (seg_masks == region.label).astype(int)
            contours = find_contours(np.array(mask))
            if contours:
                contour = max(contours, key=lambda x: len(x))
                if not ((contour[0] == contour[-1]).all()):
                    contour = np.vstack([contour, contour[0]])
                poly = Polygon(contour[:, [1, 0]]).buffer(0)
                if poly.is_valid and not poly.is_empty:
                    polygons.append(poly)

        h, w = seg_masks.shape
        if h < 4256 and w < 4256:
            scale_factor = 4256 / h
            polygons = [
                scale(poly, xfact=scale_factor, yfact=scale_factor, origin=(0, 0))
                for poly in polygons
            ]
        return polygons

    def process_masks_to_shapes(self, max_eccentricity: float = 0.95):
        """Filter masks, convert to polygons, and create shapes model."""
        self.masks = self.filter_cell_by_regionprops(self.masks, max_eccentricity)
        polygons = self.masks_to_polygons(self.masks)

        if not polygons:
            raise RuntimeError("No polygons extracted from mask — check segmentation.")

        self.gdf_polygons = gpd.GeoDataFrame(
            {"cell_id": [f"cell_{i+1}" for i in range(len(polygons))], "geometry": polygons}
        )
        shapes_df = gpd.GeoDataFrame(
            {"geometry": polygons, "region": [f"cell_{i+1}" for i in range(len(polygons))]}
        )
        shapes_df.set_geometry("geometry", inplace=True)
        self.shapes_model = ShapesModel.parse(shapes_df)
        return self.masks, self.gdf_polygons, self.shapes_model

    # ---------------- LABELS ---------------- #
    def process_labels(self, target_shape: tuple = (4256, 4256)):
        """Upscale masks and create Labels2DModel."""
        self.masks = resize(
            self.masks,
            target_shape,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        ).astype(np.int32)

        self.labels_model = Labels2DModel.parse(data=np.squeeze(self.masks), dims=("y", "x"))
        self.labels_model.name = f"{self.fov}_labels"
        return self.masks, self.labels_model

    # ---------------- POINTS & TABLES ---------------- #
    def process_points_and_tables(self):
        """Assign transcript points to segmented cells and build AnnData table."""
        points = self.sdata.points[f"{self.fov}_points"]
        points_df = points.compute()

        self.gdf_points = gpd.GeoDataFrame(
            points_df,
            geometry=gpd.points_from_xy(points_df["x"], points_df["y"]),
            crs=self.gdf_polygons.crs,
        )

        joined = gpd.sjoin(
            self.gdf_points,
            self.gdf_polygons[["cell_id", "geometry"]],
            how="left",
            predicate="within",
        )
        self.gdf_points["cell_ID"] = joined["cell_id"].values

        all_genes = self.gdf_points["target"].unique().tolist()
        cell_gene_counts = {}

        for cell in self.gdf_points["cell_ID"].dropna().unique():
            cell_data = self.gdf_points[self.gdf_points["cell_ID"] == cell]
            gene_count = Counter(cell_data["target"])
            ordered_counts = OrderedDict((g, gene_count.get(g, 0)) for g in all_genes)
            cell_gene_counts[str(cell)] = ordered_counts

        df = pd.DataFrame.from_dict(cell_gene_counts, orient="index", columns=all_genes)
        df = df.reindex(columns=all_genes, fill_value=0)

        self.vdata = ad.AnnData(X=df.values)
        self.vdata.obs_names = df.index
        self.vdata.var_names = df.columns
        self.vdata.obs["cells"] = self.vdata.obs_names
        self.vdata.var["genes"] = self.vdata.var_names
        return self.gdf_points, self.vdata

    # ---------------- SAVE ---------------- #
    def update_spatialdata(self):
        """Save processed SpatialData into a new Zarr file."""
        self.process_masks_to_shapes()
        self.process_labels()
        self.process_points_and_tables()
        
        new_sdata = sd.SpatialData()
        dst_path = os.path.join(self.zarr_dir, f"{self.fov}.zarr")

        new_sdata.images[f"{self.fov}_image"] = self.sdata.images[f"{self.fov}_image"]
        new_sdata.labels[f"{self.fov}_labels"] = self.labels_model
        new_sdata.shapes[f"{self.fov}_shapes"] = self.shapes_model
        new_sdata.points[f"{self.fov}_points"] = self.sdata.points[f"{self.fov}_points"]
        # new_sdata.tables = self.vdata
        # new_sdata.tables = {'vdata': self.vdata}

        try:
            client = Client.current()
            client.close()
        except ValueError:
            pass

        new_sdata.write(dst_path, overwrite=True)

        for var in (self.masks, self.flows, self.styles):
            try:
                del var
            except Exception:
                pass
        gc.collect()

        return new_sdata
